#include <cmath>
#include <iostream>

#include "data_structures/factory.h"
#include "mesh_adaptivity/amr_mesh_adapter.h"
#include "mesh_adaptivity/amrex_runtime.h"

int main(int argc, char** argv) {
  mpm::amr::AmrexRuntime runtime(argc, argv);

  mpm::amr::AmrConfig config;
  config.domain_minimum = {0.0, 0.0, 0.0};
  config.domain_maximum = {8000.0, 8000.0, 8000.0};
  config.coarse_cells = {10, 10, 10};
  config.refinement_ratio = {2, 2, 2};
  config.max_grid_size = 4;
  config.levels = {
      {1, {{{0.0, 0.0, 1600.0}, {4000.0, 4000.0, 6400.0}}}},
      {2, {{{0.0, 0.0, 3200.0}, {1600.0, 1600.0, 4800.0}}}}};

  auto backend = std::make_shared<mpm::amr::AmrexAmrBackend>();
  std::string error;
  if (!backend->initialise(config, &error)) {
    std::cerr << error << '\n';
    return 1;
  }

  auto element = Factory<mpm::Element<3>>::instance()->create("ED3H8");
  mpm::Mesh<3> mesh(0, false);
  mesh.assign_amr_backend(backend);
  mpm::amr::AmrMeshAdapter adapter;
  if (!adapter.populate(*backend, &mesh, "N3D", element, &error)) {
    std::cerr << error << '\n';
    return 1;
  }

  if (mesh.ncells() == 0 || mesh.nnodes() == 0) return 1;
  if (mesh.ncells() != adapter.cell_keys().size()) return 1;
  if (mesh.hanging_node_constraints().empty()) return 1;
  if (mesh.hanging_node_constraints().size() !=
      adapter.hanging_node_constraints().size())
    return 1;
  for (const auto& constraint : mesh.hanging_node_constraints()) {
    double weight_sum = 0.0;
    for (const double weight : constraint.weights) weight_sum += weight;
    if (std::abs(weight_sum - 1.0) > 1.0e-10) return 1;
  }
  for (const auto& key : adapter.cell_keys()) {
    const auto id = adapter.cell_id(key.level, key.index);
    if (!id.has_value() || mesh.cell(id.value()) == nullptr) return 1;
  }

  const auto& constraint = mesh.hanging_node_constraints().front();
  const auto hanging = mesh.node(constraint.hanging_node);
  hanging->update_mass(false, 0, 8.0);
  hanging->update_momentum(false, 0, mpm::Mesh<3>::VectorDim(8.0, 16.0, 24.0));
  mesh.condense_hanging_node_mass_momentum(0);
  if (std::abs(hanging->mass(0)) > 1.0e-14) return 1;

  double condensed_mass = 0.0;
  mpm::Mesh<3>::VectorDim condensed_momentum =
      mpm::Mesh<3>::VectorDim::Zero();
  for (const auto support_id : constraint.support_nodes) {
    const auto support = mesh.node(support_id);
    condensed_mass += support->mass(0);
    condensed_momentum += support->momentum(0);
  }
  if (std::abs(condensed_mass - 8.0) > 1.0e-10) return 1;
  if ((condensed_momentum -
       mpm::Mesh<3>::VectorDim(8.0, 16.0, 24.0))
          .norm() >
      1.0e-10)
    return 1;

  hanging->update_external_force(
      false, 0, mpm::Mesh<3>::VectorDim(8.0, 16.0, 24.0));
  hanging->update_internal_force(
      false, 0, mpm::Mesh<3>::VectorDim(4.0, 8.0, 12.0));
  mesh.project_hanging_node_forces(0);
  mpm::Mesh<3>::VectorDim external_force =
      mpm::Mesh<3>::VectorDim::Zero();
  mpm::Mesh<3>::VectorDim internal_force =
      mpm::Mesh<3>::VectorDim::Zero();
  mpm::Mesh<3>::VectorDim expected_velocity =
      mpm::Mesh<3>::VectorDim::Zero();
  mpm::Mesh<3>::VectorDim expected_acceleration =
      mpm::Mesh<3>::VectorDim::Zero();
  for (unsigned i = 0; i < constraint.support_nodes.size(); ++i) {
    const auto support = mesh.node(constraint.support_nodes[i]);
    external_force += support->external_force(0);
    internal_force += support->internal_force(0);
    const mpm::Mesh<3>::VectorDim velocity(
        static_cast<double>(i + 1), static_cast<double>(2 * i + 1),
        static_cast<double>(3 * i + 1));
    const mpm::Mesh<3>::VectorDim acceleration(
        static_cast<double>(2 * i + 1), static_cast<double>(3 * i + 1),
        static_cast<double>(4 * i + 1));
    support->assign_velocity_acceleration(0, velocity, acceleration);
    expected_velocity += constraint.weights[i] * velocity;
    expected_acceleration += constraint.weights[i] * acceleration;
  }
  if ((external_force -
       mpm::Mesh<3>::VectorDim(8.0, 16.0, 24.0))
          .norm() >
      1.0e-10)
    return 1;
  if ((internal_force -
       mpm::Mesh<3>::VectorDim(4.0, 8.0, 12.0))
          .norm() >
      1.0e-10)
    return 1;

  mesh.reconstruct_hanging_node_kinematics(0);
  if ((hanging->velocity(0) - expected_velocity).norm() > 1.0e-10) return 1;
  if ((hanging->acceleration(0) - expected_acceleration).norm() > 1.0e-10)
    return 1;

  const mpm::Mesh<3>::VectorDim coordinate(500.0, 500.0, 3500.0);
  auto particle =
      Factory<mpm::ParticleBase<3>, mpm::Index,
              const mpm::Mesh<3>::VectorDim&>::instance()
          ->create("P3D", 0, coordinate);
  if (!mesh.add_particle(particle, false)) return 1;
  if (!mesh.locate_particles_mesh().empty()) return 1;

  const auto finest = backend->locate_finest({500.0, 500.0, 3500.0});
  const auto expected_cell =
      adapter.cell_id(finest->level, finest->index);
  if (!expected_cell.has_value() || mesh.particle(0) == nullptr ||
      mesh.particle(0)->cell_id() != expected_cell.value())
    return 1;

  auto unbalanced = config;
  unbalanced.levels[0].boxes[0] = {{0.0, 0.0, 0.0},
                                    {1600.0, 1600.0, 1600.0}};
  unbalanced.levels[1].boxes[0] = {{0.0, 0.0, 0.0},
                                    {1600.0, 1600.0, 1600.0}};
  auto unbalanced_backend = std::make_shared<mpm::amr::AmrexAmrBackend>();
  if (!unbalanced_backend->initialise(unbalanced, &error)) return 1;
  mpm::Mesh<3> unbalanced_mesh(1, false);
  mpm::amr::AmrMeshAdapter unbalanced_adapter;
  if (unbalanced_adapter.populate(*unbalanced_backend, &unbalanced_mesh,
                                  "N3D", element, &error))
    return 1;
  if (error.find("2:1") == std::string::npos) return 1;

  std::cout << "AMR mesh adapter test passed: " << mesh.ncells()
            << " cells, " << mesh.nnodes() << " nodes\n";
  return 0;
}
