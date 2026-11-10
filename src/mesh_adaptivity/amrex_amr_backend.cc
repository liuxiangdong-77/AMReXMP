#include "mesh_adaptivity/amrex_amr_backend.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_Box.H>
#include <AMReX_BoxList.H>
#include <AMReX_IntVect.H>

namespace mpm {
namespace amr {
namespace {

constexpr double tolerance = 1.0e-10;

bool set_error(std::string* error, const std::string& message) {
  if (error != nullptr) *error = message;
  return false;
}

amrex::IntVect to_int_vect(const std::array<int, 3>& values) {
  return amrex::IntVect(AMREX_D_DECL(values[0], values[1], values[2]));
}

}  // namespace

bool AmrexAmrBackend::initialise(const AmrConfig& config, std::string* error) {
  config_ = config;
  level_domains_.clear();
  level_boxes_.clear();
  level_distributions_.clear();
  level_cell_sizes_.clear();
  level_cell_counts_.clear();

  for (int axis = 0; axis < 3; ++axis) {
    if (config_.domain_maximum[axis] <= config_.domain_minimum[axis])
      return set_error(error, "AMR domain bounds are invalid");
    if (config_.coarse_cells[axis] <= 0 ||
        config_.refinement_ratio[axis] <= 0) {
      return set_error(error, "AMR cell counts and ratios must be positive");
    }
  }

  const int levels = static_cast<int>(config_.levels.size()) + 1;
  level_domains_.reserve(levels);
  level_boxes_.reserve(levels);
  level_distributions_.reserve(levels);
  level_cell_sizes_.reserve(levels);
  level_cell_counts_.reserve(levels);

  std::array<double, 3> cell_size{};
  std::array<int, 3> cell_counts = config_.coarse_cells;
  for (int axis = 0; axis < 3; ++axis) {
    cell_size[axis] = (config_.domain_maximum[axis] -
                       config_.domain_minimum[axis]) /
                      static_cast<double>(config_.coarse_cells[axis]);
  }

  amrex::BoxArray coarse_boxes(
      amrex::Box(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                 to_int_vect({config_.coarse_cells[0] - 1,
                              config_.coarse_cells[1] - 1,
                              config_.coarse_cells[2] - 1})));
  coarse_boxes.maxSize(config_.max_grid_size);
  level_domains_.emplace_back(coarse_boxes.minimalBox());
  level_boxes_.emplace_back(coarse_boxes);
  level_distributions_.emplace_back(level_boxes_.back());
  level_cell_sizes_.emplace_back(cell_size);
  level_cell_counts_.emplace_back(cell_counts);

  for (const auto& level_region : config_.levels) {
    std::array<double, 3> child_cell_size = cell_size;
    std::array<int, 3> child_cell_counts = cell_counts;
    for (int axis = 0; axis < 3; ++axis) {
      child_cell_size[axis] /=
          static_cast<double>(config_.refinement_ratio[axis]);
      child_cell_counts[axis] *= config_.refinement_ratio[axis];
    }

    std::vector<amrex::Box> raw_boxes;
    for (const auto& physical_box : level_region.boxes) {
      std::array<int, 3> lo{};
      std::array<int, 3> hi{};
      for (int axis = 0; axis < 3; ++axis) {
        const double lower =
            (physical_box.minimum[axis] - config_.domain_minimum[axis]) /
            child_cell_size[axis];
        const double upper =
            (physical_box.maximum[axis] - config_.domain_minimum[axis]) /
            child_cell_size[axis];
        if (lower < -tolerance ||
            upper > child_cell_counts[axis] + tolerance)
          return set_error(error, "AMR refinement box is outside the domain");

        lo[axis] = std::max(
            0, static_cast<int>(std::floor(lower + tolerance)));
        hi[axis] = std::min(
            child_cell_counts[axis] - 1,
            static_cast<int>(std::ceil(upper - tolerance)) - 1);
        if (lo[axis] > hi[axis])
          return set_error(error, "AMR refinement box is thinner than one cell");
      }
      raw_boxes.emplace_back(to_int_vect(lo), to_int_vect(hi));
    }

    for (const auto& raw_box : raw_boxes) {
      for (int axis = 0; axis < 3; ++axis) {
        const int ratio = config_.refinement_ratio[axis];
        if (raw_box.smallEnd(axis) % ratio != 0 ||
            (raw_box.bigEnd(axis) + 1) % ratio != 0)
          return set_error(
              error,
              "AMR refinement boxes must align with parent-cell boundaries");
      }
    }

    for (int i = 0; i < raw_boxes.size(); ++i) {
      for (int j = i + 1; j < raw_boxes.size(); ++j) {
        if (raw_boxes[i].intersects(raw_boxes[j]))
          return set_error(error, "AMR boxes on one level overlap");
      }
    }

    const auto& parent = level_boxes_.back();
    for (const auto& raw_box : raw_boxes) {
      const amrex::Box coarse_box =
          amrex::coarsen(raw_box, to_int_vect(config_.refinement_ratio));
      if (!parent.contains(coarse_box))
        return set_error(error, "AMR level is not properly nested");
    }

    amrex::BoxList box_list;
    for (const auto& raw_box : raw_boxes) box_list.push_back(raw_box);
    amrex::BoxArray boxes(box_list);
    boxes.maxSize(config_.max_grid_size);
    level_domains_.emplace_back(boxes.minimalBox());
    level_boxes_.emplace_back(boxes);
    level_distributions_.emplace_back(level_boxes_.back());

    level_cell_sizes_.emplace_back(child_cell_size);
    level_cell_counts_.emplace_back(child_cell_counts);
    cell_size = child_cell_size;
    cell_counts = child_cell_counts;
  }

  return true;
}

int AmrexAmrBackend::number_of_levels() const {
  return static_cast<int>(level_boxes_.size());
}

const amrex::BoxArray& AmrexAmrBackend::box_array(int level) const {
  return level_boxes_.at(level);
}

const amrex::DistributionMapping& AmrexAmrBackend::distribution_mapping(
    int level) const {
  return level_distributions_.at(level);
}

amrex::Box AmrexAmrBackend::domain_box(int level) const {
  return level_domains_.at(level);
}

std::array<double, 3> AmrexAmrBackend::domain_minimum() const {
  return config_.domain_minimum;
}

std::array<double, 3> AmrexAmrBackend::cell_size(int level) const {
  return level_cell_sizes_.at(level);
}

std::optional<amrex::IntVect> AmrexAmrBackend::level_index(
    const std::array<double, 3>& point, int level) const {
  if (level < 0 || level >= number_of_levels()) return std::nullopt;
  const auto& dx = level_cell_sizes_[level];
  const auto& counts = level_cell_counts_[level];
  std::array<int, 3> index{};
  for (int axis = 0; axis < 3; ++axis) {
    const double relative =
        (point[axis] - config_.domain_minimum[axis]) / dx[axis];
    if (relative < -tolerance || relative > counts[axis] + tolerance)
      return std::nullopt;
    index[axis] = static_cast<int>(std::floor(relative + tolerance));
    index[axis] = std::min(counts[axis] - 1, std::max(0, index[axis]));
  }
  return to_int_vect(index);
}

std::optional<AmrCellKey> AmrexAmrBackend::locate_finest(
    const std::array<double, 3>& point) const {
  for (int level = number_of_levels() - 1; level >= 0; --level) {
    const auto index = this->level_index(point, level);
    if (index.has_value() && level_boxes_[level].contains(index.value())) {
      return AmrCellKey{level,
                        {AMREX_D_DECL((*index)[0], (*index)[1],
                                      (*index)[2])}};
    }
  }
  return std::nullopt;
}

bool AmrexAmrBackend::is_finest_valid_cell(
    int level, const std::array<int, 3>& index) const {
  if (level < 0 || level >= number_of_levels()) return false;
  const amrex::IntVect iv = to_int_vect(index);
  if (!level_boxes_[level].contains(iv)) return false;

  std::array<double, 3> center{};
  const auto& dx = level_cell_sizes_[level];
  for (int axis = 0; axis < 3; ++axis) {
    center[axis] = config_.domain_minimum[axis] +
                   (static_cast<double>(index[axis]) + 0.5) * dx[axis];
  }
  const auto finest = this->locate_finest(center);
  return finest.has_value() && finest->level == level;
}

std::vector<AmrPatchInfo> AmrexAmrBackend::patches(int level) const {
  std::vector<AmrPatchInfo> result;
  if (level < 0 || level >= number_of_levels()) return result;
  const auto& boxes = level_boxes_[level];
  const auto& distribution = level_distributions_[level];
  result.reserve(boxes.size());
  for (int i = 0; i < boxes.size(); ++i) {
    const auto& box = boxes[i];
    result.push_back(
        AmrPatchInfo{level,
                     i,
                     distribution[i],
                     {AMREX_D_DECL(box.smallEnd(0), box.smallEnd(1),
                                   box.smallEnd(2))},
                     {AMREX_D_DECL(box.bigEnd(0), box.bigEnd(1),
                                   box.bigEnd(2))}});
  }
  return result;
}

}  // namespace amr
}  // namespace mpm
