#ifndef MPM_MESH_ADAPTIVITY_AMREX_AMR_BACKEND_H_
#define MPM_MESH_ADAPTIVITY_AMREX_AMR_BACKEND_H_

#include <array>
#include <optional>
#include <string>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>

#include "amr_config.h"
#include "amr_hierarchy.h"

namespace mpm {
namespace amr {

//! Static block-structured AMR backend for three-dimensional Hex8 MPM.
class AmrexAmrBackend {
 public:
  bool initialise(const AmrConfig& config, std::string* error = nullptr);

  int number_of_levels() const;
  const amrex::BoxArray& box_array(int level) const;
  const amrex::DistributionMapping& distribution_mapping(int level) const;
  amrex::Box domain_box(int level) const;
  std::array<double, 3> domain_minimum() const;
  std::array<double, 3> cell_size(int level) const;

  std::optional<AmrCellKey> locate_finest(
      const std::array<double, 3>& point) const;
  bool is_finest_valid_cell(
      int level, const std::array<int, 3>& index) const;

  std::vector<AmrPatchInfo> patches(int level) const;

 private:
  std::optional<amrex::IntVect> level_index(
      const std::array<double, 3>& point, int level) const;

  AmrConfig config_;
  std::vector<amrex::Box> level_domains_;
  std::vector<amrex::BoxArray> level_boxes_;
  std::vector<amrex::DistributionMapping> level_distributions_;
  std::vector<std::array<double, 3>> level_cell_sizes_;
  std::vector<std::array<int, 3>> level_cell_counts_;
};

}  // namespace amr
}  // namespace mpm

#endif  // MPM_MESH_ADAPTIVITY_AMREX_AMR_BACKEND_H_
