#ifndef MPM_MESH_ADAPTIVITY_AMREX_RUNTIME_H_
#define MPM_MESH_ADAPTIVITY_AMREX_RUNTIME_H_

#include <memory>
#include <string>

namespace mpm {
namespace amr {

//! Owns the AMReX runtime for the duration of an MPM analysis.
//! The class keeps AMReX headers out of the public MPM mesh interface.
class AmrexRuntime {
 public:
  AmrexRuntime(int& argc, char**& argv);
  ~AmrexRuntime();

  AmrexRuntime(const AmrexRuntime&) = delete;
  AmrexRuntime& operator=(const AmrexRuntime&) = delete;

  //! Return whether AMReX was compiled in and initialized successfully.
  bool active() const noexcept;

  //! Return the AMReX version used by this build.
  static std::string version();

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace amr
}  // namespace mpm

#endif  // MPM_MESH_ADAPTIVITY_AMREX_RUNTIME_H_
