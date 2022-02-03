!> \file med_type_defs.F90
!!  Contains type definitions for CMEPS-related and physics-related variables

module med_type_defs 

!> \section arg_table_med_type_defs
!! \htmlinclude med_type_defs.html
!!

  use GFS_typedefs, only: GFS_statein_type
  use GFS_typedefs, only: GFS_init_type
  use GFS_typedefs, only: GFS_interstitial_type
  use GFS_typedefs, only: GFS_control_type
  use GFS_typedefs, only: GFS_coupling_type
  use ccpp_api,     only: ccpp_t

  implicit none

!! \section arg_table_physics_type
!! \htmlinclude physics_type.html
!!
  type physics_type
    type(GFS_init_type)         :: init
    type(GFS_statein_type)      :: statein
    type(GFS_interstitial_type) :: interstitial
    type(GFS_control_type)      :: model
    type(GFS_coupling_type)     :: coupling 
  end type physics_type

  type(physics_type), target :: physics
  type(ccpp_t),       target :: cdata

contains

end module med_type_defs
