#ifndef MPM_MESH_ADAPTIVITY_AMR_HIERARCHY_H_
#define MPM_MESH_ADAPTIVITY_AMR_HIERARCHY_H_

#include <array>

#include "data_types.h"

namespace mpm {
namespace amr {

struct AmrCellKey {
  int level{0};
  std::array<int, 3> index{};
};

struct AmrPatchInfo {
  int level{0};
  int patch{0};
  int rank{0};
  std::array<int, 3> lo{};
  std::array<int, 3> hi{};
};

//! Direct hanging-node constraint: u_h = sum_i weight_i * u_i.
struct HangingNodeConstraint {
  Index hanging_node{0};
  AmrCellKey parent_cell{};
  std::array<Index, 8> support_nodes{};
  std::array<double, 8> weights{};
};

}  // namespace amr
}  // namespace mpm

#endif  // MPM_MESH_ADAPTIVITY_AMR_HIERARCHY_H_
