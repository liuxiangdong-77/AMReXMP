#include <fstream>
#include <iostream>

#include "mesh_adaptivity/amrex_amr_backend.h"
#include "mesh_adaptivity/amrex_runtime.h"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " AMR_CONFIG.json\n";
    return 2;
  }

  mpm::amr::AmrexRuntime runtime(argc, argv);
  std::ifstream input(argv[1]);
  if (!input) {
    std::cerr << "Cannot open " << argv[1] << '\n';
    return 2;
  }

  mpm::amr::AmrJson json;
  input >> json;
  mpm::amr::AmrConfig config;
  std::string error;
  if (!mpm::amr::parse_amr_config(json, &config, &error)) {
    std::cerr << "Invalid AMR config: " << error << '\n';
    return 1;
  }

  mpm::amr::AmrexAmrBackend backend;
  if (!backend.initialise(config, &error)) {
    std::cerr << "Cannot initialise AMR backend: " << error << '\n';
    return 1;
  }

  std::cout << "AMReX AMR levels: " << backend.number_of_levels() << '\n';
  for (int level = 0; level < backend.number_of_levels(); ++level) {
    const auto patches = backend.patches(level);
    std::cout << "Level " << level << ": " << patches.size() << " patches\n";
    for (const auto& patch : patches) {
      std::cout << "  patch " << patch.patch << " rank " << patch.rank
                << " lo=(" << patch.lo[0] << ',' << patch.lo[1] << ','
                << patch.lo[2] << ") hi=(" << patch.hi[0] << ','
                << patch.hi[1] << ',' << patch.hi[2] << ")\n";
    }
  }
  return 0;
}
