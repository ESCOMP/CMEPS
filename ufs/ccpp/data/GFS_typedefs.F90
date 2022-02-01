module GFS_typedefs
  use machine,  only: kind_phys
  use physcons, only: con_hvap, con_cp, con_rd, con_eps
  use physcons, only: con_epsm1, con_fvirt 

  implicit none

  !--- parameter constants used for default initializations
  real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
  real(kind=kind_phys), parameter :: clear_val = zero
  real(kind=kind_phys), parameter :: huge = 9.9692099683868690E36

  !--- data containers

!! \section arg_table_GFS_init_type
!! \htmlinclude GFS_init_type.html
!!
  type GFS_init_type
    integer, pointer              :: im                     !< horizontal loop extent
  end type GFS_init_type

!! \section arg_table_GFS_statein_type
!! \htmlinclude GFS_statein_type.html
!!
  type GFS_statein_type
    real(kind=kind_phys), pointer :: pgr(:)        => null() !< surface pressure (Pa)
    real(kind=kind_phys), pointer :: ugrs(:)       => null() !< u component of layer wind (m/s)
    real(kind=kind_phys), pointer :: vgrs(:)       => null() !< v component of layer wind (m/s)
    real(kind=kind_phys), pointer :: tgrs(:)       => null() !< model layer mean temperature (K)
    real(kind=kind_phys), pointer :: qgrs(:)       => null() !< layer mean tracer concentration (kg/kg)
    real(kind=kind_phys), pointer :: prsl(:)       => null() !< model layer mean pressure (Pa)
    contains
      procedure :: create => statein_create      !< allocate array data
  end type GFS_statein_type

!! \section arg_table_GFS_interstitial_type
!! \htmlinclude GFS_interstitial_type.html
!!
  type GFS_interstitial_type
    real(kind=kind_phys), pointer :: tsfc_water(:) => null() !< surface skin temperature over water (K)
    real(kind=kind_phys), pointer :: cd_water(:)   => null() !< surface exchange coeff for momentum over water
    real(kind=kind_phys), pointer :: cdq_water(:)  => null() !< surface exchange coeff heat surface exchange coeff heat & moisture over ocean moisture over water
    real(kind=kind_phys), pointer :: ffmm_water(:) => null() !< Monin-Obukhov similarity function for momentum over water
    real(kind=kind_phys), pointer :: fm10_water(:) => null() !< Monin-Obukhov similarity parameter for momentum at 10m over water
    real(kind=kind_phys), pointer :: prslki(:)     => null() !< Exner function ratio bt midlayer and interface at 1st layer
    real(kind=kind_phys), pointer :: wet(:)        => null() !< flag indicating presence of some ocean or lake surface area fraction
    real(kind=kind_phys), pointer :: use_flake(:)  => null() !< flag indicating lake points using flake model
    real(kind=kind_phys), pointer :: wind(:)       => null() !< wind speed at lowest model level (m/s)
    logical,              pointer :: flag_iter(:)  => null() !< flag for iteration
    real(kind=kind_phys), pointer :: qss_water(:)  => null() !< surface air saturation specific humidity over water (kg/kg)
    real(kind=kind_phys), pointer :: cmm_water(:)  => null() !< momentum exchange coefficient over water (m/s)
    real(kind=kind_phys), pointer :: chh_water(:)  => null() !< thermal exchange coefficient over water (kg/m2s)
    real(kind=kind_phys), pointer :: gflx_water(:) => null() !< soil heat flux over water (W/m2)
    real(kind=kind_phys), pointer :: evap_water(:) => null() !< kinematic surface upward latent heat flux over water (m/s)
    real(kind=kind_phys), pointer :: hflx_water(:) => null() !< kinematic surface upward sensible heat flux over water (Km/s)
    real(kind=kind_phys), pointer :: ep1d_water(:) => null() !< surface upward potential latent heat flux over water (W/m2)
    contains
      procedure :: create => interstitial_create !< allocate array data
  end type GFS_interstitial_type

!! \section arg_table_GFS_control_type
!! \htmlinclude GFS_control_type.html
!!
  type GFS_control_type
    !--- tuning parameters for physical parameterizations
    logical                       :: lseaspray               !< flag for sea spray parameterization    
    !--- coupling parameters
    logical                       :: use_med_flux            !< flag for using atmosphere-ocean fluxes form mediator
    contains
      procedure :: init  => control_initialize
  end type GFS_control_type

!! \section arg_table_GFS_coupling_type
!! \htmlinclude GFS_coupling_type.html
!!
  type GFS_coupling_type
    real(kind=kind_phys), pointer :: dtsfcino_cpl(:) => null() !< sfc latent heat flux over ocean
    real(kind=kind_phys), pointer :: dqsfcino_cpl(:) => null() !< sfc sensible heat flux over ocean
    contains
      procedure :: create  => coupling_create !< allocate array data
  end type GFS_coupling_type

  contains

  subroutine statein_create(statein, im)
    implicit none
    class(GFS_statein_type) :: statein
    integer, intent(in)     :: im

    allocate(statein%pgr(im))
    statein%pgr = clear_val
    allocate(statein%ugrs(im))
    statein%ugrs = clear_val
    allocate(statein%vgrs(im))
    statein%vgrs = clear_val
    allocate(statein%tgrs(im))
    statein%tgrs = clear_val
    allocate(statein%qgrs(im))
    statein%qgrs = clear_val
    allocate(statein%prsl(im))
    statein%prsl = clear_val

  end subroutine statein_create

  subroutine interstitial_create(interstitial, im)
    implicit none
    class(GFS_interstitial_type) :: interstitial
    integer, intent(in)          :: im

    allocate(interstitial%tsfc_water(im))
    interstitial%tsfc_water = huge
    allocate(interstitial%cd_water(im))
    interstitial%cd_water = huge
    allocate(interstitial%cdq_water(im))
    interstitial%cdq_water = huge
    allocate(interstitial%ffmm_water(im))
    interstitial%ffmm_water = huge
    allocate(interstitial%fm10_water(im))
    interstitial%fm10_water = huge
    allocate(interstitial%prslki(im))
    interstitial%prslki = clear_val
    allocate(interstitial%wet(im))
    interstitial%wet = .false.
    allocate(interstitial%use_flake(im))
    interstitial%use_flake = .false.
    allocate(interstitial%wind(im))
    interstitial%wind = huge
    allocate(interstitial%flag_iter(im))
    interstitial%flag_iter = .true.
    allocate(interstitial%qss_water(im))
    interstitial%qss_water = huge
    allocate(interstitial%cmm_water(im))
    interstitial%cmm_water = huge
    allocate(interstitial%chh_water(im))
    interstitial%chh_water = huge
    allocate(interstitial%gflx_water(im))
    interstitial%gflx_water = clear_val
    allocate(interstitial%evap_water(im))
    interstitial%evap_water = huge
    allocate(interstitial%hflx_water(im))
    interstitial%hflx_water = huge
    allocate(interstitial%ep1d_water(im))
    interstitial%ep1d_water = huge

  end subroutine interstitial_create

  subroutine control_initialize(model)
    implicit none
    class(GFS_control_type) :: model

    logical :: lseaspray = .false.
    logical :: use_med_flux = .false.

  end subroutine control_initialize

  subroutine coupling_create(coupling, im)
    implicit none
    class(GFS_coupling_type) :: coupling
    integer, intent(in)      :: im

    allocate(coupling%dtsfcino_cpl(im))
    coupling%dtsfcino_cpl = clear_val
    allocate(coupling%dqsfcino_cpl(im))
    coupling%dqsfcino_cpl = clear_val

  end subroutine coupling_create
end module GFS_typedefs
