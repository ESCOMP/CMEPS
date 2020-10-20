project(CMEPS_share Fortran)
include(ExternalProject)

add_library(cmeps_share  shr_abort_mod.F90  shr_flux_mod.F90  shr_log_mod.F90  shr_mpi_mod.F90  shr_sys_mod.F90
	   glc_elevclass_mod.F90  perf_mod.F90  shr_const_mod.F90  shr_kind_mod.F90  shr_mem_mod.F90  shr_pio_mod.F90)

target_include_directories (cmeps_share PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ESMF_F90COMPILEPATHS} ${PIO_Fortran_INCLUDE_DIRS})
