ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)

ifndef PIO_INC
$(error PIO_INC should point to PIO include directory.)
endif

LIBRARY  = libcmeps_util.a

OBJ1= \
perf_mod.o \
ufs_const_mod.o \
flux_atmocn_mod.o \
glc_elevclass_mod.o

all default: $(LIBRARY)

$(LIBRARY): $(OBJ1)
	$(AR) $(ARFLAGS) $@ $(OBJ1)

%.o: %.F90.in
	perl genf90.pl $< > $(@:.o=.F90)
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) -I$(PIO_INC) -DFORTRANUNDERSCORE -DCPRINTEL $(@:.o=.F90)
%.o: %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) -I$(PIO_INC) -I. $*.F90

clean:
	$(RM) -f $(LIBRARY) *.f90 *.i90 *.o *.mod
	$(RM) -f ufs_const_mod.*

flux_atmocn_mod.o: ufs_const_mod.o
