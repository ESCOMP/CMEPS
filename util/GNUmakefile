-include $(MKFILE) 

LIBRARY  = libcmeps_util.a

OBJ1= \
med_constants_mod.o

OBJ2= \
perf_mod.o \
perf_utils.o \
shr_abort_mod.o \
shr_cal_mod.o \
shr_const_mod.o \
shr_file_mod.o \
shr_flux_mod.o \
shr_kind_mod.o \
shr_log_mod.o \
shr_mem_mod.o \
shr_mpi_mod.o \
shr_pio_mod.o \
shr_string_mod.o \
shr_strconvert_mod.o \
shr_sys_mod.o \
shr_timer_mod.o \
water_isotopes.o

OBJ3 = \
shr_infnan_mod.o \
shr_assert_mod.o
 
all default: $(LIBRARY)

$(LIBRARY): $(OBJ1) $(OBJ2) $(OBJ3)
	$(AR) $(ARFLAGS) $@ $(OBJ2) $(OBJ3)

%.o: %.F90.in
	perl genf90.pl $< > $(@:.o=.F90)
	$(FC) $(FFLAGS) $(FREE) $(ESMF_INC) ${PIOLIB_INC} -c -DFORTRANUNDERSCORE -DCPRINTEL $(@:.o=.F90)

%.o: %.F90
	$(FC) $(FFLAGS) $(FREE) $(ESMF_INC) ${PIOLIB_INC} -I. -c $*.F90

clean:
	$(RM) -f $(LIBRARY) *.f90 *.i90 *.o *.mod 
	$(RM) -f shr_infnan_mod.F90 shr_assert_mod.F90
	$(RM) -f med_constants_mod.*

med_constants_mod.o: shr_kind_mod.o shr_const_mod.o
	cp -f ../mediator/med_constants_mod.F90 .
	$(FC) $(FFLAGS) $(FREE) $(ESMF_INC) ${PIOLIB_INC} -I. -c $*.F90
perf_mod.o: perf_utils.o shr_sys_mod.o shr_kind_mod.o shr_mpi_mod.o gptl.inc
perf_utils.o: gptl.inc
shr_abort_mod.o: shr_kind_mod.o shr_mpi_mod.o shr_log_mod.o
shr_assert.h: shr_assert_mod.o
shr_assert_mod.o: dtypes.h shr_kind_mod.o shr_sys_mod.o shr_log_mod.o shr_infnan_mod.o \
   shr_strconvert_mod.o
shr_cal_mod.o: shr_kind_mod.o shr_const_mod.o shr_sys_mod.o shr_string_mod.o \
   shr_log_mod.o
shr_const_mod.o: shr_kind_mod.o
shr_file_mod: shr_kind_mod.o shr_sys_mod.o shr_log_mod.o
shr_flux_mod.o: shr_kind_mod.o shr_const_mod.o shr_sys_mod.o shr_log_mod.o water_isotopes.o
shr_infnan_mod.o: dtypes.h shr_kind_mod.o
shr_log_mod.o: shr_kind_mod.o shr_strconvert_mod.o
shr_mem_mod.o: shr_kind_mod.o shr_log_mod.o 
shr_mpi_mod.o: shr_kind_mod.o shr_log_mod.o
shr_pio_mod.o: shr_kind_mod.o shr_file_mod.o shr_log_mod.o shr_mpi_mod.o shr_sys_mod.o shr_string_mod.o
shr_string_mod.o: shr_assert.h shr_kind_mod.o shr_sys_mod.o shr_timer_mod.o \
   shr_log_mod.o
shr_strconvert_mod.o: shr_kind_mod.o shr_infnan_mod.o
shr_sys_mod.o: shr_kind_mod.o shr_log_mod.o shr_abort_mod.o
shr_timer_mod.o: shr_kind_mod.o shr_log_mod.o shr_sys_mod.o
water_isotopes.o: shr_kind_mod.o shr_const_mod.o
