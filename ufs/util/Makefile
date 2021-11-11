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
shr_abort_mod.o \
shr_const_mod.o \
shr_flux_mod.o \
shr_kind_mod.o \
shr_log_mod.o \
shr_mem_mod.o \
shr_mpi_mod.o \
shr_pio_mod.o \
shr_sys_mod.o \
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
	$(RM) -f med_constants_mod.*

perf_mod.o: shr_kind_mod.o 
shr_abort_mod.o: shr_kind_mod.o shr_mpi_mod.o shr_log_mod.o
shr_const_mod.o: shr_kind_mod.o shr_const_mod.F90
shr_flux_mod.o: shr_kind_mod.o shr_const_mod.o shr_sys_mod.o shr_log_mod.o
shr_log_mod.o: shr_kind_mod.o
shr_mem_mod.o: shr_kind_mod.o shr_log_mod.o
shr_mpi_mod.o: shr_kind_mod.o shr_log_mod.o
shr_pio_mod.o: shr_kind_mod.o shr_log_mod.o shr_mpi_mod.o shr_sys_mod.o
shr_sys_mod.o: shr_kind_mod.o shr_log_mod.o shr_abort_mod.o
glc_elevclass_mod.o: shr_kind_mod.o
