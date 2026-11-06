#include "mesh_adaptivity/amrex_runtime.h"

#include <stdexcept>

#ifdef MPM_WITH_AMREX
#include <AMReX.H>
#endif

namespace mpm {
namespace amr {

class AmrexRuntime::Impl {
 public:
#ifdef MPM_WITH_AMREX
  amrex::AMReX* instance{nullptr};
#endif
};

AmrexRuntime::AmrexRuntime(int& argc, char**& argv)
    : impl_(std::make_unique<Impl>()) {
#ifdef MPM_WITH_AMREX
#ifdef AMREX_USE_MPI
  impl_->instance = amrex::Initialize(argc, argv, false, MPI_COMM_WORLD);
#else
  impl_->instance = amrex::Initialize(argc, argv, false);
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
  if (impl_ != nullptr && impl_->instance != nullptr)
    amrex::Finalize(impl_->instance);
#endif
}

bool AmrexRuntime::active() const noexcept {
#ifdef MPM_WITH_AMREX
  return impl_ != nullptr && impl_->instance != nullptr;
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
