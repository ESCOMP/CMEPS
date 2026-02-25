module wtracers_mod

   !-----------------------------------------------------------------------------
   ! This module provides stub implementations for the shr_wtracers_mod code for when we
   ! do not have access to the CESM_share library.
   !
   ! See also the version of wtracers_mod in the cesm directory for when we have access to
   ! the CESM_share library.
   !-----------------------------------------------------------------------------

   use shr_kind_mod      , only : r8=>SHR_KIND_R8

   implicit none
   private

   public :: wtracers_present             ! return true if there are water tracers in this simulation
   public :: wtracers_get_num_tracers     ! get number of water tracers in this simulation
   public :: wtracers_is_wtracer_field    ! return true if the given field name is a water tracer field
   public :: wtracers_get_bulk_fieldname  ! return the name of the equivalent bulk field corresponding to a water tracer field
   public :: wtracers_check_tracer_ratios ! check tracer ratios against expectations

   interface wtracers_check_tracer_ratios
      module procedure wtracers_check_tracer_ratios_1d
      module procedure wtracers_check_tracer_ratios_2d
   end interface wtracers_check_tracer_ratios

   ! Suffix for water tracer field names
   character(len=*), parameter, public :: WTRACERS_SUFFIX = "_wtracers"

contains

   !-----------------------------------------------------------------------
   function wtracers_present()
      !
      ! !DESCRIPTION:
      ! Return true if there are water tracers in this simulation
      !
      ! In this stub implementation, we always return false, since water tracers are not
      ! implemented here.
      !
      ! !ARGUMENTS
      logical :: wtracers_present  ! function result
      !-----------------------------------------------------------------------
      wtracers_present = .false.
   end function wtracers_present

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

   !-----------------------------------------------------------------------
   subroutine wtracers_get_bulk_fieldname(fieldname, is_wtracer_field, bulk_fieldname)
      !
      ! !DESCRIPTION:
      ! Return the name of the equivalent bulk field corresponding to a water tracer field
      !
      ! In this stub implementation, we always return false for is_wtracer_field, and set
      ! bulk_fieldname equal to fieldname, since water tracers are not implemented here.
      !
      ! !ARGUMENTS
      character(len=*), intent(in)  :: fieldname
      logical         , intent(out) :: is_wtracer_field
      character(len=*), intent(out) :: bulk_fieldname
      !-----------------------------------------------------------------------
      is_wtracer_field = .false.
      bulk_fieldname = fieldname
   end subroutine wtracers_get_bulk_fieldname

   !-----------------------------------------------------------------------
   subroutine wtracers_check_tracer_ratios_1d(tracers, bulk, name)
      !
      ! !DESCRIPTION:
      ! Check tracer ratios (tracer/bulk) against expectations
      !
      ! In this stub implementation, we simply return without doing anything
      !
      ! !ARGUMENTS
      real(r8), intent(in) :: tracers(:,:)  ! dimensioned [tracerNum, gridcell]
      real(r8), intent(in) :: bulk(:)
      character(len=*), intent(in) :: name  ! for diagnostic output
      !-----------------------------------------------------------------------

      ! Do nothing
   end subroutine wtracers_check_tracer_ratios_1d

   !-----------------------------------------------------------------------
   subroutine wtracers_check_tracer_ratios_2d(tracers, bulk, name)
      !
      ! !DESCRIPTION:
      ! Check tracer ratios (tracer/bulk) against expectations for 2-d bulk arrays
      !
      ! In this stub implementation, we simply return without doing anything
      !
      ! !ARGUMENTS
      real(r8), intent(in) :: tracers(:,:,:)  ! dimensioned [tracerNum, ungriddedDim, gridcell]
      real(r8), intent(in) :: bulk(:,:)       ! dimensioned [ungriddedDim, gridcell]
      character(len=*), intent(in) :: name  ! for diagnostic output
      !-----------------------------------------------------------------------

      ! Do nothing
   end subroutine wtracers_check_tracer_ratios_2d

end module wtracers_mod
