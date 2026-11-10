#include "mesh_adaptivity/amr_mesh_adapter.h"

#include <algorithm>
#include <cmath>
#include <set>

#include "node_base.h"

namespace mpm {
namespace amr {
namespace {

bool set_error(std::string* error, const std::string& message) {
  if (error != nullptr) *error = message;
  return false;
}

double rounded_coordinate(double value) {
  return std::round(value * 1.0e10) / 1.0e10;
}

bool near(double lhs, double rhs) {
  return std::abs(lhs - rhs) <= 1.0e-10;
}

}  // namespace

bool AmrMeshAdapter::populate(const AmrexAmrBackend& backend, Mesh<3>* mesh,
                              const std::string& node_type,
                              const std::shared_ptr<Element<3>>& element,
                              std::string* error) {
  if (mesh == nullptr) return set_error(error, "MPM mesh output is null");
  if (element == nullptr) return set_error(error, "MPM element is null");
  if (mesh->nnodes() != 0 || mesh->ncells() != 0)
    return set_error(error, "AMR adapter requires an empty MPM mesh");

  using Vector3d = Eigen::Matrix<double, 3, 1>;
  std::vector<Vector3d> node_coordinates;
  std::map<NodePosition, Index> node_ids;
  std::vector<std::vector<Index>> cell_connectivity;
  cell_keys_.clear();
  cell_ids_.clear();
  hanging_node_constraints_.clear();

  auto node_id = [&](double x, double y, double z) -> Index {
    const NodePosition position{rounded_coordinate(x), rounded_coordinate(y),
                                rounded_coordinate(z)};
    const auto found = node_ids.find(position);
    if (found != node_ids.end()) return found->second;
    const Index id = static_cast<Index>(node_coordinates.size());
    node_ids.emplace(position, id);
    node_coordinates.emplace_back(Vector3d(x, y, z));
    return id;
  };

  for (int level = 0; level < backend.number_of_levels(); ++level) {
    const auto domain_min = backend.domain_minimum();
    const auto spacing = backend.cell_size(level);
    for (const auto& patch : backend.patches(level)) {
      for (int k = patch.lo[2]; k <= patch.hi[2]; ++k) {
        for (int j = patch.lo[1]; j <= patch.hi[1]; ++j) {
          for (int i = patch.lo[0]; i <= patch.hi[0]; ++i) {
            const std::array<int, 3> index{i, j, k};
            if (!backend.is_finest_valid_cell(level, index)) continue;

            const double x0 = domain_min[0] + i * spacing[0];
            const double x1 = x0 + spacing[0];
            const double y0 = domain_min[1] + j * spacing[1];
            const double y1 = y0 + spacing[1];
            const double z0 = domain_min[2] + k * spacing[2];
            const double z1 = z0 + spacing[2];

            const std::array<Index, 8> nodes = {
                node_id(x0, y0, z0), node_id(x1, y0, z0),
                node_id(x1, y1, z0), node_id(x0, y1, z0),
                node_id(x0, y0, z1), node_id(x1, y0, z1),
                node_id(x1, y1, z1), node_id(x0, y1, z1)};

            const Index cell_id = static_cast<Index>(cell_keys_.size());
            cell_keys_.push_back(AmrCellKey{level, index});
            cell_ids_.emplace(std::array<int, 4>{level, i, j, k}, cell_id);
            cell_connectivity.emplace_back(nodes.begin(), nodes.end());
          }
        }
      }
    }
  }

  if (node_coordinates.empty() || cell_connectivity.empty())
    return set_error(error, "AMR hierarchy produced no finest valid cells");

  std::vector<std::set<int>> incident_levels(node_coordinates.size());
  for (Index cell_id = 0; cell_id < cell_connectivity.size(); ++cell_id) {
    for (const Index node : cell_connectivity[cell_id])
      incident_levels[node].insert(cell_keys_[cell_id].level);
  }

  const auto domain_min = backend.domain_minimum();
  for (Index node = 0; node < node_coordinates.size(); ++node) {
    const auto& levels = incident_levels[node];
    if (levels.empty()) continue;
    const int coarsest_level = *levels.begin();
    const int finest_level = *levels.rbegin();
    for (int level = coarsest_level + 1; level < finest_level; ++level) {
      if (levels.find(level) == levels.end())
        return set_error(error,
                         "AMR transition violates the 2:1 level-balance rule");
    }
    if (finest_level == 0) continue;

    const int parent_level = finest_level - 1;
    const auto parent_spacing = backend.cell_size(parent_level);
    std::array<std::set<int>, 3> candidate_indices;
    for (int axis = 0; axis < 3; ++axis) {
      const double relative =
          (node_coordinates[node][axis] - domain_min[axis]) /
          parent_spacing[axis];
      const int base = static_cast<int>(std::floor(relative));
      for (int offset = -1; offset <= 1; ++offset)
        candidate_indices[axis].insert(base + offset);
    }

    std::optional<AmrCellKey> parent_key;
    std::optional<Index> parent_cell_id;
    std::array<double, 3> parent_min{};
    std::array<double, 3> parent_max{};
    for (const int k : candidate_indices[2]) {
      for (const int j : candidate_indices[1]) {
        for (const int i : candidate_indices[0]) {
          const std::array<int, 3> index{i, j, k};
          if (!backend.is_finest_valid_cell(parent_level, index)) continue;

          std::array<double, 3> lower{};
          std::array<double, 3> upper{};
          bool contains = true;
          int boundary_axes = 0;
          for (int axis = 0; axis < 3; ++axis) {
            lower[axis] =
                domain_min[axis] + index[axis] * parent_spacing[axis];
            upper[axis] = lower[axis] + parent_spacing[axis];
            const double coordinate = node_coordinates[node][axis];
            if (coordinate < lower[axis] - 1.0e-10 ||
                coordinate > upper[axis] + 1.0e-10) {
              contains = false;
              break;
            }
            if (near(coordinate, lower[axis]) ||
                near(coordinate, upper[axis]))
              ++boundary_axes;
          }

          // A hanging node lies on a coarse face or edge, but is not a vertex.
          if (!contains || boundary_axes < 1 || boundary_axes > 2) continue;
          const auto found =
              cell_ids_.find({parent_level, i, j, k});
          if (found == cell_ids_.end()) continue;
          if (!parent_cell_id.has_value() || found->second > *parent_cell_id) {
            parent_key = AmrCellKey{parent_level, index};
            parent_cell_id = found->second;
            parent_min = lower;
            parent_max = upper;
          }
        }
      }
    }

    if (!parent_key.has_value()) continue;

    const std::array<std::array<double, 3>, 8> parent_vertices = {{
        {{parent_min[0], parent_min[1], parent_min[2]}},
        {{parent_max[0], parent_min[1], parent_min[2]}},
        {{parent_max[0], parent_max[1], parent_min[2]}},
        {{parent_min[0], parent_max[1], parent_min[2]}},
        {{parent_min[0], parent_min[1], parent_max[2]}},
        {{parent_max[0], parent_min[1], parent_max[2]}},
        {{parent_max[0], parent_max[1], parent_max[2]}},
        {{parent_min[0], parent_max[1], parent_max[2]}},
    }};

    std::array<Index, 8> support_nodes{};
    for (unsigned i = 0; i < parent_vertices.size(); ++i) {
      const NodePosition position{
          rounded_coordinate(parent_vertices[i][0]),
          rounded_coordinate(parent_vertices[i][1]),
          rounded_coordinate(parent_vertices[i][2])};
      const auto found = node_ids.find(position);
      if (found == node_ids.end())
        return set_error(
            error, "A coarse parent-cell vertex is absent from the MPM mesh");
      support_nodes[i] = found->second;
    }

    const double xi =
        2.0 * (node_coordinates[node][0] - parent_min[0]) /
            (parent_max[0] - parent_min[0]) -
        1.0;
    const double eta =
        2.0 * (node_coordinates[node][1] - parent_min[1]) /
            (parent_max[1] - parent_min[1]) -
        1.0;
    const double zeta =
        2.0 * (node_coordinates[node][2] - parent_min[2]) /
            (parent_max[2] - parent_min[2]) -
        1.0;
    std::array<double, 8> weights = {
        0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta),
        0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta),
        0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta),
        0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta),
        0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta),
        0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta),
        0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta),
        0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta)};
    double weight_sum = 0.0;
    for (const double weight : weights) weight_sum += weight;
    if (!near(weight_sum, 1.0))
      return set_error(error, "Failed to construct a partition-of-unity "
                              "hanging-node interpolation");

    hanging_node_constraints_.push_back(
        HangingNodeConstraint{node, *parent_key, support_nodes, weights});
  }

  const std::set<Index> hanging_nodes = [&]() {
    std::set<Index> result;
    for (const auto& constraint : hanging_node_constraints_)
      result.insert(constraint.hanging_node);
    return result;
  }();
  for (const auto& constraint : hanging_node_constraints_)
    for (const Index support : constraint.support_nodes)
      if (hanging_nodes.find(support) != hanging_nodes.end())
        return set_error(error,
                         "Chained hanging-node constraints are unsupported");

  if (!mesh->create_nodes(0, node_type, node_coordinates, false))
    return set_error(error, "Failed to create MPM nodes from AMR hierarchy");
  if (!mesh->create_cells(0, element, cell_connectivity, false))
    return set_error(error, "Failed to create MPM cells from AMR hierarchy");
  for (Index cell_id = 0; cell_id < cell_keys_.size(); ++cell_id)
    mesh->assign_amr_cell_key(cell_id, cell_keys_[cell_id]);
  mesh->assign_hanging_node_constraints(hanging_node_constraints_);
  return true;
}

std::optional<Index> AmrMeshAdapter::cell_id(
    int level, const std::array<int, 3>& index) const {
  const auto found =
      cell_ids_.find(std::array<int, 4>{level, index[0], index[1], index[2]});
  if (found == cell_ids_.end()) return std::nullopt;
  return found->second;
}

}  // namespace amr
}  // namespace mpm
