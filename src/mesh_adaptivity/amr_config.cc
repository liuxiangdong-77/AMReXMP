#include "mesh_adaptivity/amr_config.h"

#include <algorithm>
#include <cmath>
#include <set>

namespace mpm {
namespace amr {
namespace {

bool set_error(std::string* error, const std::string& message) {
  if (error != nullptr) *error = message;
  return false;
}

bool parse_triplet(const AmrJson& value, const std::string& name,
                   std::array<double, 3>* result, std::string* error) {
  if (!value.is_array() || value.size() != 3) {
    return set_error(error, name + " must contain three numbers");
  }
  for (int i = 0; i < 3; ++i) {
    if (!value[i].is_number()) {
      return set_error(error, name + " entries must be numbers");
    }
    (*result)[i] = value[i].template get<double>();
    if (!std::isfinite((*result)[i])) {
      return set_error(error, name + " entries must be finite");
    }
  }
  return true;
}

bool parse_int_triplet(const AmrJson& value, const std::string& name,
                       std::array<int, 3>* result, std::string* error) {
  if (!value.is_array() || value.size() != 3) {
    return set_error(error, name + " must contain three integers");
  }
  for (int i = 0; i < 3; ++i) {
    if (!value[i].is_number_integer()) {
      return set_error(error, name + " entries must be integers");
    }
    (*result)[i] = value[i].template get<int>();
    if ((*result)[i] <= 0) {
      return set_error(error, name + " entries must be positive");
    }
  }
  return true;
}

}  // namespace

bool parse_amr_config(const AmrJson& json, AmrConfig* config,
                      std::string* error) {
  if (!json.is_object()) return set_error(error, "mesh.amr must be an object");
  if (config == nullptr) return set_error(error, "AMR config output is null");

  AmrConfig parsed;
  if (!parse_triplet(json.at("domain_min"), "domain_min",
                     &parsed.domain_minimum, error) ||
      !parse_triplet(json.at("domain_max"), "domain_max",
                     &parsed.domain_maximum, error) ||
      !parse_int_triplet(json.at("coarse_cells"), "coarse_cells",
                         &parsed.coarse_cells, error)) {
    return false;
  }

  if (json.contains("refinement_ratio") &&
      !parse_int_triplet(json.at("refinement_ratio"), "refinement_ratio",
                         &parsed.refinement_ratio, error)) {
    return false;
  }

  if (json.contains("max_grid_size")) {
    parsed.max_grid_size = json.at("max_grid_size").template get<int>();
    if (parsed.max_grid_size <= 0)
      return set_error(error, "max_grid_size must be positive");
  }

  for (int axis = 0; axis < 3; ++axis) {
    if (parsed.domain_maximum[axis] <= parsed.domain_minimum[axis])
      return set_error(error, "domain_max must be greater than domain_min");
  }

  const auto levels = json.at("levels");
  if (!levels.is_array() || levels.empty())
    return set_error(error, "levels must be a non-empty array");

  std::set<int> level_ids;
  for (const auto& level_json : levels) {
    if (!level_json.is_object())
      return set_error(error, "each AMR level must be an object");

    AmrLevelRegion region;
    region.level = level_json.at("level").template get<int>();
    if (region.level <= 0)
      return set_error(error, "AMR refinement levels start at 1");
    if (!level_ids.insert(region.level).second)
      return set_error(error, "AMR level ids must be unique");

    const auto boxes = level_json.at("boxes");
    if (!boxes.is_array() || boxes.empty())
      return set_error(error, "AMR level boxes must be non-empty");

    for (const auto& box_json : boxes) {
      AmrPhysicalBox box;
      if (!parse_triplet(box_json.at("min"), "box.min", &box.minimum, error) ||
          !parse_triplet(box_json.at("max"), "box.max", &box.maximum, error)) {
        return false;
      }
      for (int axis = 0; axis < 3; ++axis) {
        if (box.maximum[axis] <= box.minimum[axis])
          return set_error(error, "box.max must be greater than box.min");
      }
      region.boxes.emplace_back(box);
    }
    parsed.levels.emplace_back(std::move(region));
  }

  std::sort(parsed.levels.begin(), parsed.levels.end(),
            [](const AmrLevelRegion& lhs, const AmrLevelRegion& rhs) {
              return lhs.level < rhs.level;
            });

  for (std::size_t i = 0; i < parsed.levels.size(); ++i) {
    if (parsed.levels[i].level != static_cast<int>(i + 1))
      return set_error(error, "AMR levels must be consecutive from 1");
  }

  *config = std::move(parsed);
  return true;
}

}  // namespace amr
}  // namespace mpm
