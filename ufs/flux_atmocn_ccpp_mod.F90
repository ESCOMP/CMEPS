module flux_atmocn_ccpp_mod

  use ccpp_api,        only: ccpp_t
  use ccpp_static_api, only: ccpp_physics_init
  use ccpp_static_api, only: ccpp_physics_run
  use ccpp_static_api, only: ccpp_physics_finalize

  implicit none

  private ! default private

  public :: flux_atmOcn_init
  public :: flux_atmOcn_run
  public :: flux_atmOcn_finalize

!===============================================================================
contains
!===============================================================================

  subroutine flux_atmOcn_init(ccpp_suite_name)
    implicit none

    !--- input arguments --------------------------------
    character(len=*), intent(in) :: ccpp_suite_name

    !--- local variables --------------------------------
    integer :: ierr

  end subroutine flux_atmOcn_init

  !=============================================================================
  subroutine flux_atmOcn_run(ccpp_suite_name, group)
    implicit none

    !--- input arguments --------------------------------
    character(len=*), intent(in) :: ccpp_suite_name
    character(len=*), optional, intent(in) :: group

    !--- local variables --------------------------------
    integer :: ierr

  end subroutine flux_atmOcn_run

  !=============================================================================
  subroutine flux_atmOcn_finalize(ccpp_suite_name)
    implicit none

    !--- input arguments --------------------------------
    character(len=*), intent(in) :: ccpp_suite_name

    !--- local variables --------------------------------
    integer :: ierr

  end subroutine flux_atmOcn_finalize

end module flux_atmocn_ccpp_mod
