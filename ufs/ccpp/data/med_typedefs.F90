!> \file med_type_defs.F90
!!  Contains type definitions for CMEPS-related and physics-related variables

module med_type_defs 

  use GFS_typedefs, only: GFS_statein_type
  use machine,      only: kind_phys
  use ccpp_api,     only: ccpp_t

  implicit none

  type physics_type
    ype(GFS_statein_type) :: statein
  end type physics_type

  type(physics_type), target :: physics
  type(ccpp_t),       target :: cdata

contains

end module med_type_defs
