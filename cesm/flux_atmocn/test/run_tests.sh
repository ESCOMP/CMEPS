#!/usr/bin/env bash
#
# Build and run the cesm/flux_atmocn unit tests.
#
# Assumes a Fortran compiler, CMake (3.18+) and pFUnit 4 are already
# available in the environment (e.g. via modules or a local install).
#
# Usage: run_tests.sh [build_dir] [extra cmake args...]
#   build_dir  default: ./build-flux
#   Any further arguments are passed to the CMake configure step,
#   e.g. -DCESM_SHARE_TAG=<tag/sha>.
#
# Environment variables honored by CMake:
#   FC                  Fortran compiler to use (e.g. ifx, gfortran). The
#                       vendor names GNU and Intel (any case) are also
#                       accepted and mapped to gfortran and ifx.
#   CMAKE_PREFIX_PATH   or PFUNIT_DIR, if pFUnit is not found automatically
#
# All unit tests are run.

set -euo pipefail

build_dir=${1:-build-flux}
shift || true
test_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Map vendor names to compiler executables.
if [[ -n "${FC:-}" ]]; then
   case ${FC,,} in
      gnu)   export FC=gfortran ;;
      intel) export FC=ifx ;;
   esac
fi

cmake -S "${test_dir}" -B "${build_dir}" "$@"
cmake --build "${build_dir}" -j 8
ctest --test-dir "${build_dir}" --output-on-failure
