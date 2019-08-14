include ../../../conf/configure.nems

# remove openmp flag from ESMF_F90COMPILEOPTS if NEMS is not compiled with openmp flag
ifeq ($(OPENMP),N)
ESMF_F90COMPILEOPTS:=$(filter-out $(FFLAGS_OPENMP),$(ESMF_F90COMPILEOPTS))
ESMF_F90COMPILEOPTS:=$(filter-out -threads,$(ESMF_F90COMPILEOPTS))
endif

MEDDIR = ../../mediator
LIBRARY = libcmeps.a

OBJ1= \
esmFlds.o \
esmFldsExchange.o \
med.o \
med_constants_mod.o \
med_fraction_mod.o \
med_io_mod.o \
med_internalstate_mod.o \
med_map_mod.o \
med_merge_mod.o \
med_phases_aofluxes_mod.o \
med_phases_history_mod.o \
med_phases_ocnalb_mod.o \
med_phases_prep_atm_mod.o \
med_phases_prep_glc_mod.o \
med_phases_prep_ice_mod.o \
med_phases_prep_lnd_mod.o \
med_phases_prep_ocn_mod.o \
med_phases_prep_rof_mod.o \
med_phases_prep_wav_mod.o \
med_phases_profile_mod.o \
med_phases_restart_mod.o \
shr_nuopc_time_mod.o \
shr_nuopc_methods_mod.o \
shr_nuopc_utils_mod.o
 
all default: $(LIBRARY)

$(LIBRARY): $(OBJ1)
	$(AR) $(ARFLAGS) $@ $?

%.o: $(MEDDIR)/%.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) -I${PIO_INCDIR} -I../../util $(MEDDIR)/$*.F90

clean:
	$(RM) -f $(LIBRARY) *.i90 *.o *.mod

esmFlds.o: med_constants_mod.o
esmFldsExchange.o: med_constants_mod.o shr_nuopc_utils_mod.o shr_nuopc_methods_mod.o \
   med_internalstate_mod.o esmFlds.o
med_constants_mod.o:
med_phases_history_mod.o: med_constants_mod.o esmFlds.o shr_nuopc_utils_mod.o \
   shr_nuopc_methods_mod.o shr_nuopc_time_mod.o med_map_mod.o med_internalstate_mod.o \
   med_io_mod.o 
med_phases_restart_mod.o: med_constants_mod.o shr_nuopc_utils_mod.o esmFlds.o \
   med_internalstate_mod.o med_io_mod.o
med_phases_ocnalb_mod.o: med_constants_mod.o shr_nuopc_utils_mod.o shr_nuopc_methods_mod.o \
   med_internalstate_mod.o esmFlds.o med_map_mod.o
med_phases_prep_lnd_mod.o: esmFlds.o shr_nuopc_utils_mod.o shr_nuopc_methods_mod.o \
   med_constants_mod.o med_merge_mod.o med_map_mod.o med_internalstate_mod.o
med_internalstate_mod.o: esmFlds.o med_constants_mod.o
med_phases_prep_rof_mod.o: esmFlds.o med_constants_mod.o shr_nuopc_utils_mod.o \
   shr_nuopc_methods_mod.o med_internalstate_mod.o med_merge_mod.o med_map_mod.o
med.o: med_constants_mod.o shr_nuopc_utils_mod.o shr_nuopc_methods_mod.o shr_nuopc_time_mod.o \
   med_phases_history_mod.o med_phases_restart_mod.o med_phases_prep_atm_mod.o med_phases_prep_ice_mod.o \
   med_phases_prep_lnd_mod.o med_phases_prep_wav_mod.o med_phases_prep_glc_mod.o med_phases_prep_rof_mod.o \
   med_phases_prep_ocn_mod.o med_phases_ocnalb_mod.o med_phases_aofluxes_mod.o med_fraction_mod.o \
   med_phases_profile_mod.o med_internalstate_mod.o esmFlds.o esmFldsExchange.o med_map_mod.o \
   med_io_mod.o
med_phases_prep_glc_mod.o: esmFlds.o med_constants_mod.o shr_nuopc_utils_mod.o \
   shr_nuopc_methods_mod.o med_merge_mod.o med_map_mod.o med_internalstate_mod.o
med_phases_prep_ice_mod.o: esmFlds.o shr_nuopc_utils_mod.o shr_nuopc_methods_mod.o \
   med_constants_mod.o med_merge_mod.o med_map_mod.o med_internalstate_mod.o
med_merge_mod.o: med_constants_mod.o shr_nuopc_utils_mod.o shr_nuopc_methods_mod.o \
   med_internalstate_mod.o esmFlds.o
med_phases_prep_wav_mod.o: esmFlds.o med_constants_mod.o shr_nuopc_utils_mod.o \
   shr_nuopc_methods_mod.o med_merge_mod.o med_map_mod.o med_internalstate_mod.o
med_map_mod.o: esmFlds.o med_constants_mod.o shr_nuopc_utils_mod.o shr_nuopc_methods_mod.o \
   med_internalstate_mod.o
med_io_mod.o: med_constants_mod.o med_internalstate_mod.o shr_nuopc_methods_mod.o \
   shr_nuopc_utils_mod.o esmFlds.o
med_phases_profile_mod.o: med_constants_mod.o shr_nuopc_utils_mod.o med_internalstate_mod.o
med_phases_aofluxes_mod.o: med_constants_mod.o med_internalstate_mod.o shr_nuopc_utils_mod.o \
   shr_nuopc_methods_mod.o med_map_mod.o esmFlds.o
med_phases_prep_ocn_mod.o: med_internalstate_mod.o med_constants_mod.o shr_nuopc_utils_mod.o \
   shr_nuopc_methods_mod.o med_map_mod.o esmFlds.o med_merge_mod.o
med_phases_prep_atm_mod.o: esmFlds.o med_constants_mod.o shr_nuopc_utils_mod.o \
   shr_nuopc_methods_mod.o med_merge_mod.o med_map_mod.o med_internalstate_mod.o med_phases_ocnalb_mod.o
med_fraction_mod.o: esmFlds.o med_constants_mod.o shr_nuopc_utils_mod.o shr_nuopc_methods_mod.o \
   med_map_mod.o med_internalstate_mod.o
shr_nuopc_methods_mod.o: med_constants_mod.o shr_nuopc_utils_mod.o
shr_nuopc_time_mod.o: med_constants_mod.o shr_nuopc_utils_mod.o
shr_nuopc_utils_mod.o: med_constants_mod.o
