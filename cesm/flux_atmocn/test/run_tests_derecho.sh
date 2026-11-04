#!/usr/bin/env bash
#
# Build and run the cesm/flux_atmocn unit tests on Derecho.
#
# Usage: run_tests_derecho.sh [intel|gnu] [build_dir]
#   compiler   intel (default, ifx) or gnu (gfortran)
#   build_dir  default: ./build-flux-<compiler>
#
# All unit tests are run.

set -euo pipefail

compiler=${1:-intel}
build_dir=${2:-build-flux-${compiler}}
test_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if ! type module >/dev/null 2>&1; then
   source /etc/profile
fi

case ${compiler} in
   intel)
      module load intel/2025.2.1 cmake pfunit/4.14.0
      fc=ifx
      ;;
   gnu)
      module load gcc/14.3.0 cmake pfunit/4.14.0
      fc=gfortran
      ;;
   *)
      echo "Unknown compiler '${compiler}'; expected intel or gnu" >&2
      exit 1
      ;;
esac

cmake -S "${test_dir}" -B "${build_dir}" -DCMAKE_Fortran_COMPILER=${fc}
cmake --build "${build_dir}" -j 8
ctest --test-dir "${build_dir}" --output-on-failure
