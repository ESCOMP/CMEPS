module glc_elevclass_mod

  !---------------------------------------------------------------------
  ! This module contains the interfaces needed by mediator code - but
  ! is not used by the UFS system
  !---------------------------------------------------------------------

  use ufs_kind_mod , only : r8=>shr_kind_r8

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: glc_get_num_elevation_classes ! get the number of elevation classes
  public :: glc_get_elevation_classes     ! get elevation class of each grid cell on the glc grid.
  public :: glc_mean_elevation_virtual    ! get the mean elevation of a virtual elevation class
  public :: glc_get_fractional_icecov     ! get the fractional ice cover for each glc elevation class

  interface glc_get_elevation_classes
     module procedure glc_get_elevation_classes_with_bareland
     module procedure glc_get_elevation_classes_without_bareland
  end interface glc_get_elevation_classes

contains

  !-----------------------------------------------------------------------
  function glc_get_num_elevation_classes() result(num_elevation_classes)
    integer :: num_elevation_classes  ! function result
    num_elevation_classes = 0
  end function glc_get_num_elevation_classes

  !-----------------------------------------------------------------------
  subroutine glc_get_elevation_classes_without_bareland(glc_topo, glc_elevclass, logunit)
    real(r8), intent(in)  :: glc_topo(:)      ! topographic height
    integer , intent(out) :: glc_elevclass(:) ! elevation class
    integer , intent(in)  :: logunit
  end subroutine glc_get_elevation_classes_without_bareland

  !-----------------------------------------------------------------------
  subroutine glc_get_elevation_classes_with_bareland(glc_ice_covered, glc_topo, glc_elevclass, logunit)
    real(r8), intent(in)  :: glc_ice_covered(:) ! ice-covered (1) vs. ice-free (0)
    real(r8), intent(in)  :: glc_topo(:)        ! ice topographic height
    integer , intent(out) :: glc_elevclass(:)   ! elevation class
    integer , intent(in)  :: logunit
  end subroutine glc_get_elevation_classes_with_bareland

  !-----------------------------------------------------------------------
  function glc_mean_elevation_virtual(elevation_class, logunit) result(mean_elevation)
    real(r8) :: mean_elevation  ! function result
    integer, intent(in) :: elevation_class
    integer, optional, intent(in) :: logunit
    mean_elevation = 0.0_r8
  end function glc_mean_elevation_virtual

  !-----------------------------------------------------------------------
  subroutine glc_get_fractional_icecov(nec, glc_topo, glc_icefrac, glc_icefrac_ec, logunit)
    integer , intent(in)  :: nec              ! number of elevation classes 
    real(r8), intent(in)  :: glc_topo(:)      ! topographic height
    real(r8), intent(in)  :: glc_icefrac(:)
    real(r8), intent(out) :: glc_icefrac_ec(:,:)
    integer , intent(in)  :: logunit
  end subroutine glc_get_fractional_icecov

end module glc_elevclass_mod
