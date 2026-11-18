cmake_minimum_required(VERSION 3.18)
project(flux_atmocn_tests LANGUAGES Fortran)

# pFUnit unit tests for the CESM atm/ocn bulk flux schemes in cesm/flux_atmocn.
# Builds without CIME, ESMF, PIO or MPI. See README.md for how to run.

find_package(PFUNIT 4 REQUIRED)
enable_testing()

#-----------------------------------------------------------------------------
# CESM_share sources
#
# Only a few self-contained modules are compiled from CESM_share. The sources are
# fetched but never added as a subdirectory, because CESM_share's own CMakeLists.txt
# requires CIME, PIO and ESMF (SOURCE_SUBDIR points at a directory with no
# CMakeLists.txt, so FetchContent_MakeAvailable only downloads).
#
# To test against a local CESM_share checkout instead of fetching:
#   -DFETCHCONTENT_SOURCE_DIR_CESM_SHARE=/path/to/CESM_share
#-----------------------------------------------------------------------------
include(FetchContent)
# Defaults to the head of the CESM_share main branch; set -DCESM_SHARE_TAG=<tag/sha>
# to test against a specific version instead.
set(CESM_SHARE_TAG "main" CACHE STRING "CESM_share git branch/tag/sha to test against")
FetchContent_Declare(cesm_share
  GIT_REPOSITORY https://github.com/ESCOMP/CESM_share.git
  GIT_TAG        ${CESM_SHARE_TAG}
  GIT_SHALLOW    TRUE
  SOURCE_SUBDIR  do_not_add_subdirectory)
FetchContent_MakeAvailable(cesm_share)
set(SHARE_SRC ${cesm_share_SOURCE_DIR}/src)

#-----------------------------------------------------------------------------
# Library of the code under test, plus stubs and shared test helpers
#-----------------------------------------------------------------------------
set(FLUX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_library(flux_under_test STATIC
  ${SHARE_SRC}/shr_kind_mod.F90
  ${SHARE_SRC}/shr_const_mod.F90
  ${SHARE_SRC}/shr_wv_sat_mod.F90
  stubs/shr_sys_mod_stub.F90
  ${FLUX_DIR}/shr_flux_mod.F90
  ${FLUX_DIR}/flux_atmocn_Large.F90
  ${FLUX_DIR}/flux_atmocn_COARE_mod.F90
  ${FLUX_DIR}/flux_atmocn_UA_mod.F90
  ${FLUX_DIR}/flux_atmocn_Diurnal_mod.F90
  helpers/flux_test_fixtures.F90)
set_target_properties(flux_under_test PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(flux_under_test PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_link_libraries(flux_under_test PUBLIC PFUNIT::funit)
target_compile_options(flux_under_test PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -g>
  $<$<Fortran_COMPILER_ID:IntelLLVM>:-check all -fpe0 -traceback -init=snan -g>)

#-----------------------------------------------------------------------------
# Tests: each .pf file in a subdirectory of test/ (e.g. unit_tests/) becomes its own
# ctest executable, and all of them run on every `ctest` invocation.
#
# The .pf paths must be relative: pFUnit writes the generated .F90 next to an
# absolute-path .pf file (i.e. into the source tree), but under
# CMAKE_CURRENT_BINARY_DIR for a relative one.
#-----------------------------------------------------------------------------
file(GLOB pf_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
     ${CMAKE_CURRENT_SOURCE_DIR}/*/*.pf)
foreach(pf ${pf_files})
  get_filename_component(name ${pf} NAME_WE)
  add_pfunit_ctest(${name} TEST_SOURCES ${pf} LINK_LIBRARIES flux_under_test)
endforeach()
