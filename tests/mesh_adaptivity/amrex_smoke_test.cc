#include <iostream>

#include <AMReX.H>
#include <AMReX_Box.H>

#include "mesh_adaptivity/amrex_runtime.h"

int main(int argc, char** argv) {
  mpm::amr::AmrexRuntime runtime(argc, argv);

  amrex::Box box(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                 amrex::IntVect(AMREX_D_DECL(3, 3, 3)));

  std::cout << "AMReX version: " << mpm::amr::AmrexRuntime::version() << '\n';
  std::cout << "AMReX smoke box: " << box << '\n';
  std::cout << "AMReX runtime active: " << std::boolalpha << runtime.active()
            << '\n';
  return (runtime.active() && box.ok()) ? 0 : 1;
}
