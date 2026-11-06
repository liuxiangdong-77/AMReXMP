#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: $0 AMREX_SOURCE_DIR [INSTALL_PREFIX] [BUILD_DIR]" >&2
  exit 2
fi

source_dir=$(cd "$1" && pwd)
install_prefix=${2:-"$HOME/opt/amrex-26.11"}
build_dir=${3:-"$HOME/opt/build-amrex-26.11"}

cc=${AMREX_CC:-gcc-13}
cxx=${AMREX_CXX:-g++-13}
release_flags=${AMREX_RELEASE_FLAGS:-"-O1 -DNDEBUG"}
build_jobs=${AMREX_BUILD_JOBS:-4}
linear_solvers=${AMREX_LINEAR_SOLVERS:-OFF}

if [[ ! -f "$source_dir/CMakeLists.txt" ]]; then
  echo "Not an AMReX source directory: $source_dir" >&2
  exit 2
fi

cmake -S "$source_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$cc" \
  -DCMAKE_CXX_COMPILER="$cxx" \
  -DCMAKE_CXX_FLAGS_RELEASE="$release_flags" \
  -DCMAKE_INSTALL_PREFIX="$install_prefix" \
  -DAMReX_SPACEDIM=3 \
  -DAMReX_MPI=ON \
  -DAMReX_OMP=ON \
  -DAMReX_PARTICLES=ON \
  -DAMReX_LINEAR_SOLVERS="$linear_solvers"

cmake --build "$build_dir" --parallel "$build_jobs"
cmake --install "$build_dir"

echo
echo "AMReX installed in: $install_prefix"
echo "Configure MPM with:"
echo "  -DMPM_WITH_AMREX=ON"
echo "  -DCMAKE_PREFIX_PATH=$install_prefix"
