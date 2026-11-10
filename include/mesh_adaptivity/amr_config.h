#ifndef MPM_MESH_ADAPTIVITY_AMR_CONFIG_H_
#define MPM_MESH_ADAPTIVITY_AMR_CONFIG_H_

#include <array>
#include <string>
#include <vector>

#include "json.hpp"

namespace mpm {
namespace amr {

using AmrJson = nlohmann::json;

struct AmrPhysicalBox {
  std::array<double, 3> minimum{};
  std::array<double, 3> maximum{};
};

struct AmrLevelRegion {
  int level{0};
  std::vector<AmrPhysicalBox> boxes;
};

struct AmrConfig {
  std::array<double, 3> domain_minimum{};
  std::array<double, 3> domain_maximum{};
  std::array<int, 3> coarse_cells{};
  std::array<int, 3> refinement_ratio{2, 2, 2};
  int max_grid_size{32};
  std::vector<AmrLevelRegion> levels;
};

//! Parse the mesh.amr JSON object.
bool parse_amr_config(const AmrJson& json, AmrConfig* config,
                      std::string* error = nullptr);

}  // namespace amr
}  // namespace mpm

#endif  // MPM_MESH_ADAPTIVITY_AMR_CONFIG_H_
