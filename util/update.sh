#!/bin/bash

# list of files that are needed to create libcmeps_util.a
lst="dtypes.h
genf90.pl
gptl.inc
perf_mod.F90
perf_utils.F90
shr_abort_mod.F90
shr_assert.h
shr_assert_mod.F90.in
shr_cal_mod.F90
shr_const_mod.F90
shr_file_mod.F90
shr_flux_mod.F90
shr_infnan_mod.F90.in
shr_kind_mod.F90
shr_log_mod.F90
shr_mem_mod.F90
shr_mpi_mod.F90
shr_pio_mod.F90
shr_strconvert_mod.F90
shr_string_mod.F90
shr_sys_mod.F90
shr_timer_mod.F90
water_isotopes.F90"

# clone cime
git clone https://github.com/ESMCI/cime .cime
cd .cime
git checkout cmeps_v0.4.2
cd -

# copy files
for i in $lst
do
  f=`find .cime/. -not -path '*/*/*/pio*' -not -path '*/tools*' -name "$i*"`
  if [ -z "$f" ]; then
    f=`find .cime/. -name "$i*"`
  fi  
  cp -f $f .
done

# remove cime clone
rm -rf .cime
