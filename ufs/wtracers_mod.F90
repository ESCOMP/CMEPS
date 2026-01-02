module wtracers_mod

   !-----------------------------------------------------------------------------
   ! This module provides stub implementations for the shr_wtracers_mod code for when we
   ! do not have access to the CESM_share library.
   !
   ! See also the version of wtracers_mod in the cesm directory for when we have access to
   ! the CESM_share library.
   !-----------------------------------------------------------------------------

   implicit none
   private

   public :: wtracers_get_num_tracers   ! get number of water tracers in this simulation
   public :: wtracers_is_wtracer_field  ! return true if the given field name is a water tracer field

contains

   !-----------------------------------------------------------------------
   function wtracers_get_num_tracers()
      !
      ! !DESCRIPTION:
      ! Get number of water tracers in this simulation
      !
      ! In this stub implementation, we always return 0, since water tracers are not
      ! implemented here.
      !
      ! !ARGUMENTS
      integer :: wtracers_get_num_tracers  ! function result
      !-----------------------------------------------------------------------
      wtracers_get_num_tracers = 0
   end function wtracers_get_num_tracers

   !-----------------------------------------------------------------------
   function wtracers_is_wtracer_field(fieldname)
      !
      ! !DESCRIPTION:
      ! Return true if the given field name is a water tracer field
      !
      ! In this stub implementation, we always return false, since water tracers are not
      ! implemented here.
      !
      ! !ARGUMENTS
      character(len=*), intent(in) :: fieldname
      logical :: wtracers_is_wtracer_field
      !-----------------------------------------------------------------------
      wtracers_is_wtracer_field = .false.
   end function wtracers_is_wtracer_field

end module wtracers_mod
