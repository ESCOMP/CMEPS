module MED_typedefs

!> \section arg_table_MED_typedefs
!! \htmlinclude MED_typedefs.html
!!
  use machine,  only: kind_phys
  use physcons, only: con_hvap, con_cp, con_rd, con_eps
  use physcons, only: con_epsm1, con_fvirt, con_g 

  implicit none

  !--- parameter constants used for default initializations
  real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
  real(kind=kind_phys), parameter :: clear_val = zero
  real(kind=kind_phys), parameter :: huge = 9.9692099683868690E36

  !--- data containers

!! \section arg_table_MED_init_type
!! \htmlinclude MED_init_type.html
!!
  type MED_init_type
    integer                       :: im                     !< horizontal loop extent
  end type MED_init_type

!! \section arg_table_MED_statein_type
!! \htmlinclude MED_statein_type.html
!!
  type MED_statein_type
    real(kind=kind_phys), pointer :: pgr(:)        => null() !< surface pressure (Pa)
    real(kind=kind_phys), pointer :: ugrs(:)       => null() !< u component of layer wind (m/s)
    real(kind=kind_phys), pointer :: vgrs(:)       => null() !< v component of layer wind (m/s)
    real(kind=kind_phys), pointer :: tgrs(:)       => null() !< model layer mean temperature (K)
    real(kind=kind_phys), pointer :: qgrs(:)       => null() !< layer mean tracer concentration (kg/kg)
    real(kind=kind_phys), pointer :: prsl(:)       => null() !< model layer mean pressure (Pa)
    real(kind=kind_phys), pointer :: zlvl(:)       => null() !< layer 1 height above ground (m)
    real(kind=kind_phys), pointer :: prsik(:)      => null() !< dimensionless Exner function at lowest model interface
    real(kind=kind_phys), pointer :: prslk(:)      => null() !< dimensionless Exner function at lowest model layer
    real(kind=kind_phys), pointer :: u10m(:)       => null() !< 10 meter u wind speed
    real(kind=kind_phys), pointer :: v10m(:)       => null() !< 10 meter v wind speed
    contains
      procedure :: create => statein_create      !< allocate array data
  end type MED_statein_type

!! \section arg_table_MED_interstitial_type
!! \htmlinclude MED_interstitial_type.html
!!
  type MED_interstitial_type
    ! water
    real(kind=kind_phys), pointer :: tsfc_water(:)   => null() !< surface skin temperature over water (K)
    real(kind=kind_phys), pointer :: cd_water(:)     => null() !< surface exchange coeff for momentum over water
    real(kind=kind_phys), pointer :: cdq_water(:)    => null() !< surface exchange coeff heat surface exchange coeff heat & moisture over ocean moisture over water
    real(kind=kind_phys), pointer :: ffmm_water(:)   => null() !< Monin-Obukhov similarity function for momentum over water
    real(kind=kind_phys), pointer :: fm10_water(:)   => null() !< Monin-Obukhov similarity parameter for momentum at 10m over water
    real(kind=kind_phys), pointer :: prslki(:)       => null() !< Exner function ratio bt midlayer and interface at 1st layer
    logical,              pointer :: wet(:)          => null() !< flag indicating presence of some ocean or lake surface area fraction
    logical,              pointer :: use_flake(:)    => null() !< flag indicating lake points using flake model
    real(kind=kind_phys), pointer :: wind(:)         => null() !< wind speed at lowest model level (m/s)
    logical,              pointer :: flag_iter(:)    => null() !< flag for iteration
    real(kind=kind_phys), pointer :: qss_water(:)    => null() !< surface air saturation specific humidity over water (kg/kg)
    real(kind=kind_phys), pointer :: cmm_water(:)    => null() !< momentum exchange coefficient over water (m/s)
    real(kind=kind_phys), pointer :: chh_water(:)    => null() !< thermal exchange coefficient over water (kg/m2s)
    real(kind=kind_phys), pointer :: gflx_water(:)   => null() !< soil heat flux over water (W/m2)
    real(kind=kind_phys), pointer :: evap_water(:)   => null() !< kinematic surface upward latent heat flux over water (m/s)
    real(kind=kind_phys), pointer :: hflx_water(:)   => null() !< kinematic surface upward sensible heat flux over water (Km/s)
    real(kind=kind_phys), pointer :: ep1d_water(:)   => null() !< surface upward potential latent heat flux over water (W/m2)
    real(kind=kind_phys), pointer :: tsurf_water(:)  => null() !< surface skin temperature after iteration over water (K)
    real(kind=kind_phys), pointer :: uustar_water(:) => null() !< surface friction velocity over water (m/s)
    real(kind=kind_phys), pointer :: rb_water(:)     => null() !< bulk Richardson number at the surface over water 
    real(kind=kind_phys), pointer :: stress_water(:) => null() !< surface wind stress over water
    real(kind=kind_phys), pointer :: ffhh_water(:)   => null() !< Monin-Obukhov similarity function for heat over water
    real(kind=kind_phys), pointer :: fh2_water(:)    => null() !< Monin-Obukhov similarity parameter for heat at 2m over water
    real(kind=kind_phys), pointer :: ztmax_water(:)  => null() !< bounded surface roughness length for heat over water (m)

    ! land, not used to calculate aofluxes
    real(kind=kind_phys), pointer :: zvfun(:)        => null() !< function of surface roughness length and green vegetation fraction
    real(kind=kind_phys), pointer :: sigmaf(:)       => null() !< areal fractional cover of green vegetation bounded on the bottom
    logical,              pointer :: dry(:)          => null() !< flag indicating presence of some land surface area fraction
    real(kind=kind_phys), pointer :: tsfcl(:)        => null() !< surface skin temperature over land (K)
    real(kind=kind_phys), pointer :: tsurf_land(:)   => null() !< surface skin temperature after iteration over land (K) 
    real(kind=kind_phys), pointer :: uustar_land(:)  => null() !< surface friction velocity over land (m/s)
    real(kind=kind_phys), pointer :: cd_land(:)      => null() !< surface exchange coeff for momentum over land 
    real(kind=kind_phys), pointer :: cdq_land(:)     => null() !< surface exchange coeff heat surface exchange coeff heat & moisture over ocean moisture over land 
    real(kind=kind_phys), pointer :: rb_land(:)      => null() !< bulk Richardson number at the surface over land
    real(kind=kind_phys), pointer :: stress_land(:)  => null() !< surface wind stress over land 
    real(kind=kind_phys), pointer :: ffmm_land(:)    => null() !< Monin-Obukhov similarity function for momentum over land 
    real(kind=kind_phys), pointer :: ffhh_land(:)    => null() !< Monin-Obukhov similarity function for heat over land
    real(kind=kind_phys), pointer :: fm10_land(:)    => null() !< Monin-Obukhov similarity parameter for momentum at 10m over land
    real(kind=kind_phys), pointer :: fh2_land(:)     => null() !< Monin-Obukhov similarity parameter for heat at 2m over land
    real(kind=kind_phys), pointer :: ztmax_land(:)   => null() !< bounded surface roughness length for heat over land (m)

    ! ice, not used to calculate aofluxes
    logical,              pointer :: icy(:)          => null() !< flag indicating presence of some sea ice surface area fraction
    real(kind=kind_phys), pointer :: tisfc(:)        => null() !< surface skin temperature over ice (K)
    real(kind=kind_phys), pointer :: tsurf_ice(:)    => null() !< surface skin temperature after iteration over ice (K) 
    real(kind=kind_phys), pointer :: uustar_ice(:)   => null() !< surface friction velocity over ice (m/s)
    real(kind=kind_phys), pointer :: cd_ice(:)       => null() !< surface exchange coeff for momentum over ice 
    real(kind=kind_phys), pointer :: cdq_ice(:)      => null() !< surface exchange coeff heat surface exchange coeff heat & moisture over ocean moisture over ice 
    real(kind=kind_phys), pointer :: rb_ice(:)       => null() !< bulk Richardson number at the surface over ice
    real(kind=kind_phys), pointer :: stress_ice(:)   => null() !< surface wind stress over ice
    real(kind=kind_phys), pointer :: ffmm_ice(:)     => null() !< Monin-Obukhov similarity function for momentum over ice
    real(kind=kind_phys), pointer :: ffhh_ice(:)     => null() !< Monin-Obukhov similarity function for heat over ice
    real(kind=kind_phys), pointer :: fm10_ice(:)     => null() !< Monin-Obukhov similarity parameter for momentum at 10m over ice
    real(kind=kind_phys), pointer :: fh2_ice(:)      => null() !< Monin-Obukhov similarity parameter for heat at 2m over ice
    real(kind=kind_phys), pointer :: ztmax_ice(:)    => null() !< bounded surface roughness length for heat over ice (m)

    ! others
    real(kind=kind_phys), pointer :: z01d(:)         => null() !< perturbation of momentum roughness length
    real(kind=kind_phys), pointer :: zt1d(:)         => null() !< perturbation of heat to momentum roughness length ratio
    logical,              pointer :: flag_guess(:)   => null() !< flag for guess run
    contains
      procedure :: create => interstitial_create !< allocate array data
      procedure :: phys_reset  => interstitial_phys_reset !<   reset array data for physics
  end type MED_interstitial_type

!! \section arg_table_MED_control_type
!! \htmlinclude MED_control_type.html
!!
  type MED_control_type
    !--- tuning parameters for physical parameterizations
    logical                       :: lseaspray               !< flag for sea spray parameterization    
    !--- coupling parameters
    logical                       :: use_med_flux            !< flag for using atmosphere-ocean fluxes form mediator
    !--- land/surface model parameters, not used to calculate aofluxes
    integer                       :: ivegsrc                 !< land use dataset choice 0 => USGS, 1 => IGBP, 2 => UMD  
    integer                       :: lsm                     !< flag for land surface model
    integer                       :: lsm_noahmp              !< flag for NOAH MP land surface model
    !--- tuning parameters for physical parameterizations
    logical                       :: redrag                  !< flag for reduced drag coeff. over sea
    !--- surface layer z0 scheme
    integer                       :: sfc_z0_type             !< surface roughness options over water
    !--- potential temperature definition in surface layer physics
    logical                       :: thsfc_loc               !< flag for reference pressure in theta calculation
    !--- near surface temperature model
    integer                       :: nstf_name(5)            !< NSSTM flag: off/uncoupled/coupled=0/1/2
    contains
      procedure :: init  => control_initialize
  end type MED_control_type

!! \section arg_table_MED_coupling_type
!! \htmlinclude MED_coupling_type.html
!!
  type MED_coupling_type
    real(kind=kind_phys), pointer :: dtsfcino_cpl(:) => null() !< sfc latent heat flux over ocean
    real(kind=kind_phys), pointer :: dqsfcino_cpl(:) => null() !< sfc sensible heat flux over ocean
    contains
      procedure :: create  => coupling_create !< allocate array data
  end type MED_coupling_type

!! \section arg_table_MED_grid_type
!! \htmlinclude MED_grid_type.html
!!
  type MED_grid_type
    real(kind=kind_phys), pointer :: area(:)         => null() !< area of the grid cell
    contains
      procedure :: create  => grid_create !< allocate array data
  end type MED_grid_type

!! \section arg_table_MED_sfcprop_type
!! \htmlinclude MED_sfcprop_type.html
!!
  type MED_sfcprop_type
    ! water
    real(kind=kind_phys), pointer :: zorlw(:)        => null()  !< surface roughness length over water (cm)

    ! land, not used to calculate aofluxes
    integer,              pointer :: vtype(:)        => null()  !< vegetation type
    real(kind=kind_phys), pointer :: shdmax(:)       => null()  !< max fractional coverage of green vegetation
    real(kind=kind_phys), pointer :: zorll(:)        => null()  !< surface roughness length over land (cm)

    ! ice, not used to calculate aofluxes
    real(kind=kind_phys), pointer :: zorli(:)        => null()  !< surface roughness length over ice (cm)

    ! wave
    real(kind=kind_phys), pointer :: zorlwav(:)      => null()  !< surface roughness length from wave model (cm)

    ! other
    real(kind=kind_phys), pointer :: zorl(:)         => null()  !< surface roughness length (cm)

    contains
      procedure :: create  => sfcprop_create !< allocate array data
  end type MED_sfcprop_type

  public MED_init_type
  public MED_statein_type
  public MED_coupling_type
  public MED_control_type
  public MED_interstitial_type
  public MED_grid_type

  contains

  subroutine statein_create(statein, im)
    implicit none
    class(MED_statein_type) :: statein
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
    allocate(statein%zlvl(im))
    statein%zlvl = clear_val
    allocate(statein%prsik(im))
    statein%prsik = clear_val
    allocate(statein%prslk(im))
    statein%prslk = clear_val
    allocate(statein%u10m(im))
    statein%u10m = clear_val
    allocate(statein%v10m(im))
    statein%v10m = clear_val

  end subroutine statein_create

  subroutine interstitial_create(interstitial, im)
    implicit none
    class(MED_interstitial_type) :: interstitial
    integer, intent(in)          :: im

    ! water
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
    allocate(interstitial%tsurf_water(im))
    interstitial%tsurf_water = huge
    allocate(interstitial%uustar_water(im))
    interstitial%uustar_water = huge
    allocate(interstitial%rb_water(im))
    interstitial%rb_water = huge
    allocate(interstitial%stress_water(im))
    interstitial%stress_water = huge
    allocate(interstitial%ffhh_water(im))
    interstitial%ffhh_water = huge
    allocate(interstitial%fh2_water(im))
    interstitial%fh2_water = huge
    allocate(interstitial%ztmax_water(im))
    interstitial%ztmax_water = clear_val

    ! land
    allocate(interstitial%zvfun(im))
    interstitial%zvfun = clear_val
    allocate(interstitial%sigmaf(im))
    interstitial%sigmaf = clear_val
    allocate(interstitial%dry(im))
    interstitial%dry = .false.
    allocate(interstitial%tsfcl(im))
    interstitial%tsfcl = clear_val
    allocate(interstitial%tsurf_land(im))
    interstitial%tsurf_land = huge
    allocate(interstitial%uustar_land(im))
    interstitial%uustar_land = huge
    allocate(interstitial%cd_land(im))
    interstitial%cd_land = huge
    allocate(interstitial%cdq_land(im))
    interstitial%cdq_land = huge
    allocate(interstitial%rb_land(im))
    interstitial%rb_land = huge
    allocate(interstitial%stress_land(im))
    interstitial%stress_land = huge
    allocate(interstitial%ffmm_land(im))
    interstitial%ffmm_land = huge
    allocate(interstitial%ffhh_land(im))
    interstitial%ffhh_land = huge
    allocate(interstitial%fm10_land(im))
    interstitial%fm10_land = huge
    allocate(interstitial%fh2_land(im))
    interstitial%fh2_land = huge
    allocate(interstitial%ztmax_land(im))
    interstitial%ztmax_land = clear_val

    ! ice
    allocate(interstitial%icy(im))
    interstitial%icy = .false.
    allocate(interstitial%tisfc(im))
    interstitial%tisfc = clear_val
    allocate(interstitial%tsurf_ice(im))
    interstitial%tsurf_ice = huge
    allocate(interstitial%uustar_ice(im))
    interstitial%uustar_ice = huge
    allocate(interstitial%cd_ice(im))
    interstitial%cd_ice = huge
    allocate(interstitial%cdq_ice(im))
    interstitial%cdq_ice = huge
    allocate(interstitial%rb_ice(im))
    interstitial%rb_ice = huge
    allocate(interstitial%stress_ice(im))
    interstitial%stress_ice = huge
    allocate(interstitial%ffmm_ice(im))
    interstitial%ffmm_ice = huge
    allocate(interstitial%ffhh_ice(im))
    interstitial%ffhh_ice = huge
    allocate(interstitial%fm10_ice(im))
    interstitial%fm10_ice = huge
    allocate(interstitial%fh2_ice(im))
    interstitial%fh2_ice = huge
    allocate(interstitial%ztmax_ice(im))
    interstitial%ztmax_ice = clear_val

    ! others
    allocate(interstitial%z01d(im))
    interstitial%z01d = clear_val
    allocate(interstitial%zt1d(im))
    interstitial%zt1d = clear_val
    allocate(interstitial%flag_guess(im))
    interstitial%flag_guess = .false.

  end subroutine interstitial_create

  subroutine interstitial_phys_reset(interstitial)
    implicit none
    class(MED_interstitial_type) :: interstitial

    interstitial%cd_ice = huge
    interstitial%cd_land = huge
    interstitial%cd_water = huge
    interstitial%cdq_ice = huge
    interstitial%cdq_land = huge
    interstitial%cdq_water = huge
    interstitial%chh_water = huge
    interstitial%cmm_water = huge
    interstitial%dry = .false.
    interstitial%ep1d_water = huge
    interstitial%evap_water = huge
    interstitial%ffhh_ice = huge
    interstitial%ffhh_land = huge
    interstitial%ffhh_water = huge
    interstitial%ffmm_ice = huge
    interstitial%ffmm_land = huge
    interstitial%ffmm_water = huge
    interstitial%fh2_ice = huge
    interstitial%fh2_land = huge
    interstitial%fh2_water = huge
    interstitial%flag_guess = .false.
    interstitial%flag_iter = .true.
    interstitial%fm10_ice = huge
    interstitial%fm10_land = huge
    interstitial%fm10_water = huge
    interstitial%gflx_water = clear_val
    interstitial%hflx_water = huge
    interstitial%icy = .false.
    interstitial%prslki = clear_val
    interstitial%qss_water = huge
    interstitial%rb_ice = huge
    interstitial%rb_land = huge
    interstitial%rb_water = huge
    interstitial%sigmaf = clear_val
    interstitial%stress_ice = huge
    interstitial%stress_land = huge
    interstitial%stress_water = huge
    interstitial%tisfc = clear_val
    interstitial%tsfc_water = huge
    interstitial%tsfcl = clear_val
    interstitial%tsurf_ice = huge
    interstitial%tsurf_land = huge
    interstitial%tsurf_water = huge
    interstitial%use_flake = .false.
    interstitial%uustar_ice = huge
    interstitial%uustar_land = huge
    interstitial%uustar_water = huge
    interstitial%wet = .false.
    interstitial%wind = huge
    interstitial%z01d = clear_val
    interstitial%zt1d = clear_val
    interstitial%ztmax_ice = clear_val
    interstitial%ztmax_land = clear_val
    interstitial%ztmax_water = clear_val
    interstitial%zvfun = clear_val

  end subroutine interstitial_phys_reset

  subroutine control_initialize(model)
    implicit none
    class(MED_control_type) :: model

    model%lseaspray = .false.
    model%use_med_flux = .false.
    model%ivegsrc = 2
    model%redrag = .false.
    model%sfc_z0_type = 0
    model%thsfc_loc = .true.
    model%lsm = 1
    model%lsm_noahmp = 2
    model%nstf_name = (/0,0,1,0,5/)

  end subroutine control_initialize

  subroutine coupling_create(coupling, im)
    implicit none
    class(MED_coupling_type) :: coupling
    integer, intent(in)      :: im

    allocate(coupling%dtsfcino_cpl(im))
    coupling%dtsfcino_cpl = clear_val
    allocate(coupling%dqsfcino_cpl(im))
    coupling%dqsfcino_cpl = clear_val

  end subroutine coupling_create

  subroutine grid_create(grid, im)
    implicit none
    class(MED_grid_type) :: grid
    integer, intent(in)  :: im

    allocate(grid%area(im))
    grid%area = clear_val

  end subroutine grid_create

  subroutine sfcprop_create(sfcprop, im)
    implicit none
    class(MED_sfcprop_type) :: sfcprop
    integer, intent(in)  :: im

    allocate(sfcprop%vtype(im))
    sfcprop%vtype = zero
    allocate(sfcprop%shdmax(im))
    sfcprop%shdmax = clear_val
    allocate(sfcprop%zorl(im))
    sfcprop%zorl = clear_val
    allocate(sfcprop%zorlw(im))
    sfcprop%zorlw = clear_val
    allocate(sfcprop%zorll(im))
    sfcprop%zorll = clear_val
    allocate(sfcprop%zorli(im))
    sfcprop%zorli = clear_val
    allocate(sfcprop%zorlwav(im))
    sfcprop%zorlwav = clear_val

  end subroutine sfcprop_create

end module MED_typedefs
