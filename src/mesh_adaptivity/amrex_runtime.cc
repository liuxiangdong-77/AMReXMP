#include "mesh_adaptivity/amrex_runtime.h"

#include <stdexcept>

#ifdef MPM_WITH_AMREX
#include <AMReX.H>
#endif

namespace mpm {
namespace amr {

AmrexRuntime::AmrexRuntime(int& argc, char**& argv)
    : instance_(nullptr) {
#ifdef MPM_WITH_AMREX
#ifdef AMREX_USE_MPI
  instance_ = amrex::Initialize(argc, argv, false, MPI_COMM_WORLD);
#else
  instance_ = amrex::Initialize(argc, argv, false);
#endif
  amrex::SetVerbose(0);
#else
  (void)argc;
  (void)argv;
  throw std::runtime_error(
      "AMReX runtime requested, but MPM was built with MPM_WITH_AMREX=OFF");
#endif
}

AmrexRuntime::~AmrexRuntime() {
#ifdef MPM_WITH_AMREX
  if (instance_ != nullptr) amrex::Finalize(instance_);
#endif
}

bool AmrexRuntime::active() const noexcept {
#ifdef MPM_WITH_AMREX
  return instance_ != nullptr;
#else
  return false;
#endif
}

std::string AmrexRuntime::version() {
#ifdef MPM_WITH_AMREX
  return amrex::Version();
#else
  return "AMReX backend disabled";
#endif
}

}  // namespace amr
}  // namespace mpm
