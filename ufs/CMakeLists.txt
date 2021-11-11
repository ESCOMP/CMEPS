project(CMEPS_share Fortran)
include(ExternalProject)

add_library(cmeps_share flux_atmocn_mod.F90 glc_elevclass_mod.F90 perf_mod.F90 ufs_const_mod.F90)

target_include_directories (cmeps_share PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ESMF_F90COMPILEPATHS} ${PIO_Fortran_INCLUDE_DIRS})
