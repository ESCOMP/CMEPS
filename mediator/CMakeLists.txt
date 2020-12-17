project(cmeps Fortran)

set(SRCFILES esmFldsExchange_cesm_mod.F90 med_fraction_mod.F90
              med_methods_mod.F90 med_phases_prep_ice_mod.F90
              med_phases_restart_mod.F90 esmFldsExchange_hafs_mod.F90
              med_internalstate_mod.F90 med_phases_aofluxes_mod.F90
              med_phases_prep_lnd_mod.F90 med_time_mod.F90
              esmFldsExchange_nems_mod.F90 med_io_mod.F90
              med_phases_history_mod.F90 med_phases_prep_ocn_mod.F90
              med_utils_mod.F90 esmFlds.F90 med_kind_mod.F90
              med_phases_ocnalb_mod.F90 med_phases_prep_rof_mod.F90
              med_constants_mod.F90 med_map_mod.F90
              med_phases_prep_atm_mod.F90 med_phases_prep_wav_mod.F90
              med.F90 med_merge_mod.F90 med_phases_prep_glc_mod.F90
              med_phases_profile_mod.F90 med_diag_mod.F90
              med_phases_post_ocn_mod.F90
              med_phases_post_atm_mod.F90 med_phases_post_ice_mod.F90
              med_phases_post_lnd_mod.F90 med_phases_post_glc_mod.F90
              med_phases_post_rof_mod.F90 med_phases_post_wav_mod.F90)

foreach(FILE ${SRCFILES})
  if(EXISTS "${CASEROOT}/SourceMods/src.cmeps/${FILE}")
    list(REMOVE_ITEM SRCFILES ${FILE})
    list(APPEND SRCFILES "${CASEROOT}/SourceMods/src.cmeps/${FILE}")
    message("Using ${FILE} from ${CASEROOT}/SourceMods/src.cmeps")
  endif()
endforeach()
add_library(cmeps ${SRCFILES})

if(BLD_STANDALONE)
  add_dependencies(cmeps cmeps_share)
endif()

target_include_directories (cmeps PUBLIC ${ESMF_F90COMPILEPATHS})
target_include_directories (cmeps PUBLIC "${CMAKE_BINARY_DIR}/nems/util")
target_include_directories (cmeps PUBLIC ${PIO_Fortran_INCLUDE_DIR})

install(TARGETS cmeps
        LIBRARY DESTINATION lib)
