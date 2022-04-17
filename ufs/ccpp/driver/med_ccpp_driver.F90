module med_ccpp_driver

  use ccpp_types,          only: ccpp_t
  use ccpp_static_api_med, only: ccpp_physics_init
  use ccpp_static_api_med, only: ccpp_physics_run
  use ccpp_static_api_med, only: ccpp_physics_finalize

  use MED_data,        only: physics, cdata 

  implicit none

  private ! default private

  public :: med_ccpp_driver_init
  public :: med_ccpp_driver_run
  public :: med_ccpp_driver_finalize

!===============================================================================
contains
!===============================================================================

  subroutine med_ccpp_driver_init(ccpp_suite)
    implicit none

    !--- input arguments --------------------------------
    character(len=*), intent(in) :: ccpp_suite

    !--- local variables --------------------------------
    integer :: ierr

    ! initialize CCPP physics (run all _init routines)
    call ccpp_physics_init(cdata, suite_name=trim(ccpp_suite), ierr=ierr)
    if (ierr /= 0) then
       write(0,'(a)') "An error occurred in ccpp_physics_init"
       write(0,'(a)') trim(cdata%errmsg)
       return
    end if

  end subroutine med_ccpp_driver_init

  !=============================================================================
  subroutine med_ccpp_driver_run(ccpp_suite, group)
    implicit none

    !--- input arguments --------------------------------
    character(len=*), intent(in) :: ccpp_suite
    character(len=*), optional, intent(in) :: group

    !--- local variables --------------------------------
    integer :: ierr

    ! run CCPP physics (run all _run routines)
    if (present(group)) then
       call ccpp_physics_run(cdata, suite_name=trim(ccpp_suite), group_name=trim(group), ierr=ierr)
    else
       call ccpp_physics_run(cdata, suite_name=trim(ccpp_suite), ierr=ierr)
    end if
    if (ierr /= 0) then
       write(0,'(a)') "An error occurred in ccpp_physics_run"
       write(0,'(a)') trim(cdata%errmsg)
       return
    end if

  end subroutine med_ccpp_driver_run

  !=============================================================================
  subroutine med_ccpp_driver_finalize(ccpp_suite)
    implicit none

    !--- input arguments --------------------------------
    character(len=*), intent(in) :: ccpp_suite

    !--- local variables --------------------------------
    integer :: ierr

    ! finalize CCPP physics (run all _finalize routines)
    call ccpp_physics_finalize(cdata, suite_name=trim(ccpp_suite), ierr=ierr)
    if (ierr /= 0) then
       write(0,'(a)') "An error occurred in ccpp_physics_finalize"
       write(0,'(a)') trim(cdata%errmsg)
       return
    end if

  end subroutine med_ccpp_driver_finalize

end module med_ccpp_driver
