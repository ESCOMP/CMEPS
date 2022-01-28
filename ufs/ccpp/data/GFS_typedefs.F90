module GFS_typedefs
  use machine, only: kind_phys

  implicit none

  !--- parameter constants used for default initializations
  real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
  real(kind=kind_phys), parameter :: clear_val = zero

  !--- data containers
  type GFS_statein_type
    real (kind=kind_phys), pointer :: prsl(:) => null() !< model layer mean pressure Pa
    real (kind=kind_phys), pointer :: tgrs(:) => null() !< model layer mean temperature in k
  contains
    procedure :: create  => statein_create  !<   allocate array data
  end type GFS_statein_type

!------------------------------------------------------------------------------------
! combined type of all of the above except GFS_control_type and GFS_interstitial_type
!------------------------------------------------------------------------------------
!! \section arg_table_GFS_data_type
!! \htmlinclude GFS_data_type.html
!!
  type GFS_data_type
     type(GFS_statein_type) :: statein
  end type GFS_data_type

  contains

  subroutine statein_create(statein, im)
    class(GFS_statein_type) :: statein
    integer, intent(in)     :: im

    allocate(statein%prsl(im))
    statein%prsl = clear_val
    allocate(statein%tgrs(im))
    statein%tgrs = clear_val

  end subroutine statein_create

end module GFS_typedefs
