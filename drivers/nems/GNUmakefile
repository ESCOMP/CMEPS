ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)

ifneq ($(origin PIO_INC), environment)
$(error PIO_INC should point to PIO include directory.)
endif

MEDDIR = ../../mediator
LIBRARY = libcmeps.a
LIBRARY_UTIL = ../../util/libcmeps_util.a

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

all default: install

install: $(LIBRARY)
	rm -f cmeps.mk.install
	@echo "# ESMF self-describing build dependency makefile fragment" > cmeps.mk.install
	@echo "# src location: $(PWD)" >> cmeps.mk.install
	@echo  >> cmeps.mk.install
	@echo "ESMF_DEP_FRONT     = MED" >> cmeps.mk.install
	@echo "ESMF_DEP_INCPATH   = $(INSTALLDIR)" >> cmeps.mk.install
	@echo "ESMF_DEP_CMPL_OBJS = " >> cmeps.mk.install
	@echo "ESMF_DEP_LINK_OBJS = $(INSTALLDIR)/libcmeps.a $(INSTALLDIR)/libcmeps_util.a $(PIO_LIB)" >> cmeps.mk.install
	mkdir -p $(INSTALLDIR)
	cp -f $(LIBRARY_UTIL) $(INSTALLDIR)
	cp -f $(LIBRARY) $(INSTALLDIR)
	cp -f cmeps.mk.install $(INSTALLDIR)/cmeps.mk

$(LIBRARY): $(OBJ1)
	$(AR) $(ARFLAGS) $@ $?

$(LIBRARY_UTIL):
	cd ../../util ;\
	exec $(MAKE) 

%.o: $(MEDDIR)/%.F90 $(LIBRARY_UTIL)
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEPATHS) -I${PIO_INC} -I../../util $(MEDDIR)/$*.F90

clean:
	$(RM) -f $(LIBRARY) *.i90 *.o *.mod cmeps.mk.install
	cd ../../util ;\
	exec $(MAKE) clean

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
