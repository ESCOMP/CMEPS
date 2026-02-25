module wtracers_mod

   !-----------------------------------------------------------------------------
   ! This module wraps shr_wtracers_mod from the CESM_share repository to avoid direct
   ! dependencies on this share code from CMEPS.
   !
   ! It is acceptable for CESM-specific CMEPS code (e.g., in the cesm directory) to call
   ! directly into shr_wtracers_mod, but CMEPS code shared between CESM and other modeling
   ! systems (e.g., UFS) should go through this module.
   !
   ! See also the version of wtracers_mod in the ufs directory for when we do not have
   ! access to the CESM_share library.
   !-----------------------------------------------------------------------------

   use shr_wtracers_mod, only : wtracers_present             => shr_wtracers_present
   use shr_wtracers_mod, only : wtracers_get_num_tracers     => shr_wtracers_get_num_tracers
   use shr_wtracers_mod, only : wtracers_is_wtracer_field    => shr_wtracers_is_wtracer_field
   use shr_wtracers_mod, only : wtracers_get_bulk_fieldname  => shr_wtracers_get_bulk_fieldname
   use shr_wtracers_mod, only : wtracers_check_tracer_ratios => shr_wtracers_check_tracer_ratios
   use shr_wtracers_mod, only : WTRACERS_SUFFIX

   implicit none
   private

   public :: wtracers_present             ! return true if there are water tracers in this simulation
   public :: wtracers_get_num_tracers     ! get number of water tracers in this simulation
   public :: wtracers_is_wtracer_field    ! return true if the given field name is a water tracer field
   public :: wtracers_get_bulk_fieldname  ! return the name of the equivalent bulk field corresponding to a water tracer field
   public :: wtracers_check_tracer_ratios ! check tracer ratios against expectations

   public :: WTRACERS_SUFFIX ! suffix for water tracer field names

end module wtracers_mod
