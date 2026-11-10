#ifndef MPM_MESH_ADAPTIVITY_AMR_MESH_ADAPTER_H_
#define MPM_MESH_ADAPTIVITY_AMR_MESH_ADAPTER_H_

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "amr_hierarchy.h"
#include "amrex_amr_backend.h"
#include "element.h"
#include "mesh.h"

namespace mpm {
namespace amr {

//! Populate a three-dimensional MPM Mesh from the finest valid AMR cells.
class AmrMeshAdapter {
 public:
  bool populate(const AmrexAmrBackend& backend, Mesh<3>* mesh,
                const std::string& node_type,
                const std::shared_ptr<Element<3>>& element,
                std::string* error = nullptr);

  const std::vector<AmrCellKey>& cell_keys() const { return cell_keys_; }
  const std::vector<HangingNodeConstraint>& hanging_node_constraints() const {
    return hanging_node_constraints_;
  }
  std::optional<Index> cell_id(int level,
                               const std::array<int, 3>& index) const;

 private:
  struct NodePosition {
    double x{0.0};
    double y{0.0};
    double z{0.0};

    bool operator<(const NodePosition& other) const {
      if (x != other.x) return x < other.x;
      if (y != other.y) return y < other.y;
      return z < other.z;
    }
  };

  std::vector<AmrCellKey> cell_keys_;
  std::map<std::array<int, 4>, Index> cell_ids_;
  std::vector<HangingNodeConstraint> hanging_node_constraints_;
};

}  // namespace amr
}  // namespace mpm

#endif  // MPM_MESH_ADAPTIVITY_AMR_MESH_ADAPTER_H_
