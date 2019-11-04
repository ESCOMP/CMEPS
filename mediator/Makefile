ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)

ifndef PIO_INCLUDE_DIR
$(error PIO_INCLUDE_DIR not set)
endif

ifdef INTERNAL_PIO_INIT
CPPDEFS := -DINTERNAL_PIO_INIT
endif

LIBRARY := libcmeps.a

OBJ = $(patsubst %.F90, %.o, $(wildcard *.F90))

all default: $(LIBRARY)

$(LIBRARY): $(OBJ)
	$(AR) $(ARFLAGS) $@ $?

%.o: %.F90 
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(CPPDEFS) -I${PIO_INCLUDE_DIR} -I../nems/util $*.F90

clean:
	$(RM) -f $(LIBRARY) *.i90 *.o *.mod

med_kind_mod.o : 
med_constants_mod.o : med_kind_mod.o
esmFlds.o : med_kind_mod.o
esmFldsExchange.o : med_kind_mod.o med_methods_mod.o esmFlds.o med_internalstate_mod.o med_utils_mod.o   
med.o : med_kind_mod.o med_phases_profile_mod.o med_utils_mod.o med_phases_prep_rof_mod.o med_phases_aofluxes_mod.o \
   med_phases_prep_ice_mod.o med_fraction_mod.o med_map_mod.o med_constants_mod.o med_phases_prep_wav_mod.o \
   med_phases_prep_lnd_mod.o med_phases_history_mod.o med_phases_ocnalb_mod.o med_phases_restart_mod.o \
   med_time_mod.o med_internalstate_mod.o med_phases_prep_atm_mod.o esmFldsExchange.o \
   med_phases_prep_glc_mod.o esmFlds.o med_io_mod.o med_methods_mod.o med_phases_prep_ocn_mod.o   
med_fraction_mod.o : med_kind_mod.o med_utils_mod.o med_internalstate_mod.o med_constants_mod.o med_map_mod.o med_methods_mod.o esmFlds.o   
med_internalstate_mod.o : med_kind_mod.o esmFlds.o   
med_io_mod.o : med_kind_mod.o med_methods_mod.o med_constants_mod.o med_internalstate_mod.o med_utils_mod.o   
med_map_mod.o : med_kind_mod.o med_internalstate_mod.o med_constants_mod.o med_methods_mod.o esmFlds.o med_utils_mod.o   
med_merge_mod.o : med_kind_mod.o med_constants_mod.o med_internalstate_mod.o esmFlds.o med_methods_mod.o med_utils_mod.o   
med_methods_mod.o : med_kind_mod.o med_utils_mod.o med_constants_mod.o   
med_phases_aofluxes_mod.o : med_kind_mod.o med_utils_mod.o med_map_mod.o med_constants_mod.o med_internalstate_mod.o esmFlds.o med_methods_mod.o   
med_phases_history_mod.o : med_kind_mod.o med_utils_mod.o med_time_mod.o med_internalstate_mod.o med_constants_mod.o med_map_mod.o med_methods_mod.o med_io_mod.o esmFlds.o   
med_phases_ocnalb_mod.o : med_kind_mod.o med_utils_mod.o med_map_mod.o med_constants_mod.o med_internalstate_mod.o esmFlds.o med_methods_mod.o   
med_phases_prep_atm_mod.o : med_kind_mod.o esmFlds.o med_methods_mod.o med_merge_mod.o med_map_mod.o med_constants_mod.o med_phases_ocnalb_mod.o med_internalstate_mod.o med_utils_mod.o   
med_phases_prep_glc_mod.o : med_kind_mod.o med_utils_mod.o med_internalstate_mod.o med_map_mod.o med_constants_mod.o med_methods_mod.o esmFlds.o   
med_phases_prep_ice_mod.o : med_kind_mod.o med_utils_mod.o med_methods_mod.o med_merge_mod.o esmFlds.o med_internalstate_mod.o med_constants_mod.o med_map_mod.o   
med_phases_prep_lnd_mod.o : med_kind_mod.o med_internalstate_mod.o med_map_mod.o med_constants_mod.o med_merge_mod.o med_methods_mod.o esmFlds.o med_utils_mod.o   
med_phases_prep_ocn_mod.o : med_kind_mod.o med_internalstate_mod.o med_constants_mod.o med_map_mod.o med_merge_mod.o med_methods_mod.o esmFlds.o med_utils_mod.o   
med_phases_prep_rof_mod.o : med_kind_mod.o med_internalstate_mod.o med_map_mod.o med_constants_mod.o med_merge_mod.o med_methods_mod.o esmFlds.o med_utils_mod.o   
med_phases_prep_wav_mod.o : med_kind_mod.o med_utils_mod.o med_internalstate_mod.o med_constants_mod.o med_map_mod.o med_methods_mod.o med_merge_mod.o esmFlds.o   
med_phases_profile_mod.o : med_kind_mod.o med_utils_mod.o med_constants_mod.o med_internalstate_mod.o   
med_phases_restart_mod.o : med_kind_mod.o med_utils_mod.o med_constants_mod.o med_internalstate_mod.o esmFlds.o med_io_mod.o   
med_time_mod.o : med_kind_mod.o med_utils_mod.o med_constants_mod.o   
med_utils_mod.o : med_kind_mod.o med_utils_mod.F90    
