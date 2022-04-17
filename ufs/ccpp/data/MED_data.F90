!> \file MED_data.F90
!!  Contains type definitions for CMEPS-related and physics-related variables

module MED_data 

!> \section arg_table_MED_data
!! \htmlinclude MED_data.html
!!

  use MED_typedefs, only: MED_statein_type
  use MED_typedefs, only: MED_init_type
  use MED_typedefs, only: MED_interstitial_type
  use MED_typedefs, only: MED_control_type
  use MED_typedefs, only: MED_coupling_type
  use MED_typedefs, only: MED_grid_type
  use MED_typedefs, only: MED_sfcprop_type  
  use MED_typedefs, only: MED_diag_type
  use ccpp_types,   only: ccpp_t

  implicit none

  public physics

!! \section arg_table_physics_type
!! \htmlinclude physics_type.html
!!
  type physics_type
    type(MED_init_type)         :: init
    type(MED_statein_type)      :: statein
    type(MED_interstitial_type) :: interstitial
    type(MED_control_type)      :: model
    type(MED_coupling_type)     :: coupling 
    type(MED_grid_type)         :: grid
    type(MED_sfcprop_type)      :: sfcprop
    type(MED_diag_type)         :: diag
  end type physics_type

  type(physics_type), save, target :: physics
  type(ccpp_t),       save, target :: cdata

contains

end module MED_data
