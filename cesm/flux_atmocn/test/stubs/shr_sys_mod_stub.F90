module shr_sys_mod

  !-----------------------------------------------------------------------------
  ! Unit-test stub for CESM_share's shr_sys_mod.
  !
  ! The real module depends on ESMF and MPI. This stub only provides shr_sys_abort,
  ! which raises a pFUnit exception instead of stopping the program, so that abort
  ! paths can be tested with @assertExceptionRaised.
  !-----------------------------------------------------------------------------

  use funit, only : throw

  implicit none
  private

  public :: shr_sys_abort

contains

  subroutine shr_sys_abort(string, rc, line, file)
    ! Same interface as the real shr_sys_abort (shr_abort_mod::shr_abort_abort)
    character(len=*), intent(in), optional :: string  ! error message string
    integer         , intent(in), optional :: rc      ! error code
    integer         , intent(in), optional :: line
    character(len=*), intent(in), optional :: file

    if (present(string)) then
       call throw(trim(string))
    else
       call throw('shr_sys_abort called')
    end if

  end subroutine shr_sys_abort

end module shr_sys_mod
