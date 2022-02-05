module med_ccpp_driver

  use ccpp_api,        only: ccpp_t
  use ccpp_static_api, only: ccpp_physics_init
  use ccpp_static_api, only: ccpp_physics_run
  use ccpp_static_api, only: ccpp_physics_finalize

  use med_typedefs   , only: physics, cdata 

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

    ! init
    print*, "call ccpp_physics_init for suite "//trim(ccpp_suite)
    call ccpp_physics_init(cdata, suite_name=trim(ccpp_suite), ierr=ierr)
    if (ierr /= 0) then
       write(0,'(a)') "An error occurred in ccpp_physics_init"
       write(0,'(a)') trim(cdata%errmsg)
       return
    end if

  end subroutine med_ccpp_driver_init

  !=============================================================================
  subroutine med_ccpp_driver_run(ccpp_suite_name, group)
    implicit none

    !--- input arguments --------------------------------
    character(len=*), intent(in) :: ccpp_suite_name
    character(len=*), optional, intent(in) :: group

    !--- local variables --------------------------------
    integer :: ierr

  end subroutine med_ccpp_driver_run

  !=============================================================================
  subroutine med_ccpp_driver_finalize(ccpp_suite_name)
    implicit none

    !--- input arguments --------------------------------
    character(len=*), intent(in) :: ccpp_suite_name

    !--- local variables --------------------------------
    integer :: ierr

  end subroutine med_ccpp_driver_finalize

end module med_ccpp_driver
