module wtracers_mod

   !-----------------------------------------------------------------------------
   ! This module wraps shr_wtracers_mod from the CESM_share repository to avoid direct
   ! dependencies on this share code from CMEPS.
   !
   ! See also the version of wtracers_mod in the ufs directory for when we do not have
   ! access to the CESM_share library.
   !-----------------------------------------------------------------------------

   use shr_wtracers_mod, only : wtracers_is_wtracer_field => shr_wtracers_is_wtracer_field

   implicit none
   private

   public :: wtracers_is_wtracer_field  ! return true if the given field name is a water tracer field

end module wtracers_mod
