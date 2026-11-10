#include <cmath>
#include <iostream>

#include "mesh_adaptivity/amrex_amr_backend.h"
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

  mpm::amr::AmrexAmrBackend backend;
  std::string error;
  if (!backend.initialise(config, &error)) {
    std::cerr << error << '\n';
    return 1;
  }

  if (backend.number_of_levels() != 3) {
    std::cerr << "Unexpected level count: " << backend.number_of_levels()
              << '\n';
    return 1;
  }
  if (!backend.box_array(0).ok() || !backend.box_array(1).ok() ||
      !backend.box_array(2).ok()) {
    std::cerr << "Invalid BoxArray\n";
    return 1;
  }

  const auto level1 =
      backend.locate_finest({3000.0, 3000.0, 5500.0});
  if (!level1.has_value() || level1->level != 1) {
    std::cerr << "Expected level 1, got "
              << (level1.has_value() ? level1->level : -1) << '\n';
    return 1;
  }

  const auto level2 =
      backend.locate_finest({500.0, 500.0, 3500.0});
  if (!level2.has_value() || level2->level != 2) {
    std::cerr << "Expected level 2, got "
              << (level2.has_value() ? level2->level : -1) << '\n';
    return 1;
  }

  if (backend.locate_finest({9000.0, 500.0, 500.0}).has_value()) {
    std::cerr << "Out-of-domain point was located\n";
    return 1;
  }

  if (backend.is_finest_valid_cell(0, {0, 0, 4})) {
    std::cerr << "Covered coarse cell was treated as finest valid\n";
    return 1;
  }
  if (!backend.is_finest_valid_cell(2, {2, 2, 17})) {
    std::cerr << "Expected finest cell was rejected\n";
    return 1;
  }

  for (int level = 0; level < backend.number_of_levels(); ++level) {
    const auto patches = backend.patches(level);
    if (patches.empty()) {
      std::cerr << "No patches on level " << level << '\n';
      return 1;
    }
    for (const auto& patch : patches) {
      if (patch.rank < 0) {
        std::cerr << "Invalid patch rank\n";
        return 1;
      }
    }
  }

  auto misaligned = config;
  misaligned.levels[0].boxes[0].minimum[2] = 2000.0;
  mpm::amr::AmrexAmrBackend invalid_backend;
  if (invalid_backend.initialise(misaligned, &error)) {
    std::cerr << "A parent-cell-misaligned refinement box was accepted\n";
    return 1;
  }

  std::cout << "AMReX AMR backend test passed\n";
  return 0;
}
