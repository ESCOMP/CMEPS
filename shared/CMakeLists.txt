project(shared Fortran)
include(ExternalProject)

add_library(shared shr_infnan_mod.F90 shr_kind_mod.F90 shr_orb_mod.F90 shr_sys_mod.F90 shr_sys_mod.F90 shr_sys_mod.F90 shr_sys_mod.F90 shr_log_mod.F90 shr_strconvert_mod.F90 shr_abort_mod.F90 shr_const_mod.F90 nuopc_shr_methods.F90)

target_include_directories (shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ESMF_F90COMPILEPATHS} ${PIO_Fortran_INCLUDE_DIRS})

target_compile_definitions(shared PUBLIC HAVE_IEEE_ARITHMETIC)
