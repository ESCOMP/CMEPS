module MED_typedefs

!> \section arg_table_MED_typedefs
!! \htmlinclude MED_typedefs.html
!!
  use machine,  only: kind_phys
  use physcons, only: con_hvap, con_cp, con_rd, con_eps, con_rocp
  use physcons, only: con_epsm1, con_fvirt, con_g 
  use physcons, only: con_tice, karman

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
    real(kind=kind_phys), pointer :: u10m(:)       => null() !< 10 meter u wind speed (m/s)
    real(kind=kind_phys), pointer :: v10m(:)       => null() !< 10 meter v wind speed (m/s)
    real(kind=kind_phys), pointer :: stc(:,:)      => null() !< soil temperature (K)
    contains
      procedure :: create => statein_create      !< allocate array data
  end type MED_statein_type

!! \section arg_table_MED_stateout_type
!! \htmlinclude MED_stateout_type.html
!!
  type MED_stateout_type
    real(kind=kind_phys), pointer :: gu0(:)   => null()  !< updated zonal wind
    real(kind=kind_phys), pointer :: gv0(:)   => null()  !< updated meridional wind
    real(kind=kind_phys), pointer :: gt0(:)   => null()  !< updated temperature
    real(kind=kind_phys), pointer :: gq0(:)   => null()  !< updated tracers
    contains
      procedure :: create  => stateout_create  !<   allocate array data
  end type MED_stateout_type

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
    integer,              pointer :: use_lake_model(:)=>null() !< 0 for points that don't use a lake model, lkm for points that do
    real (kind=kind_phys),pointer :: lake_t2m (:)    => null() !< 2 meter temperature from CLM Lake model 
    real (kind=kind_phys),pointer :: lake_q2m (:)    => null() !< 2 meter humidity from CLM Lake model
    real(kind=kind_phys), pointer :: wind(:)         => null() !< wind speed at lowest model level (m/s)
    logical,              pointer :: flag_iter(:)    => null() !< flag for iteration
    logical,              pointer :: flag_lakefreeze(:) => null() !< flag for lake freeze
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
    logical,              pointer :: lake(:)         => null() !< flag indicating presence of some lake surface area fraction
    real(kind=kind_phys), pointer :: tprcp_water(:)  => null() !< total precipitation amount in each time step over water

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
    real(kind=kind_phys), pointer :: frland(:)       => null() !< land area fraction used in microphysics schemes
    real(kind=kind_phys), pointer :: tprcp_land(:)   => null() !< total precipitation amount in each time step over land
    real(kind=kind_phys), pointer :: qss_land(:)     => null() !< surface air saturation specific humidity over land (kg/kg)
    real(kind=kind_phys), pointer :: evap_land(:)    => null() !< kinematic surface upward latent heat flux over land (m/s)
    real(kind=kind_phys), pointer :: hflx_land(:)    => null() !< kinematic surface upward sensible heat flux over land (Km/s)
    real(kind=kind_phys), pointer :: hflxq(:)        => null() !< kinematic surface upward sensible heat flux reduced by surface roughness and vegetation
    real(kind=kind_phys), pointer :: chh_land(:)     => null() !< thermal exchange coefficient over land (kg/m2s)
    real(kind=kind_phys), pointer :: cmm_land(:)     => null() !< momentum exchange coefficient over land (m/s)
    real(kind=kind_phys), pointer :: gflx_land(:)    => null() !< soil heat flux over land (W/m2)
    real(kind=kind_phys), pointer :: ep1d_land(:)    => null() !< surface upward potential latent heat flux over land (W/m2)

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
    logical,              pointer :: flag_cice(:)    => null() !< flag for cice
    real(kind=kind_phys), pointer :: tprcp_ice(:)    => null() !< total precipitation amount in each time step over ice
    integer,              pointer :: islmsk(:)       => null() !< sea/land/ice mask (=0/1/2)
    integer,              pointer :: islmsk_cice(:)  => null() !< sea/land/ice mask cice (=0/1/2)
    real(kind=kind_phys), pointer :: ep1d_ice(:)     => null() !< surface upward potential latent heat flux over ice (W/m2)
    real(kind=kind_phys), pointer :: gflx_ice(:)     => null() !< soil heat flux over ice
    real(kind=kind_phys), pointer :: qss_ice(:)      => null() !< surface air saturation specific humidity over ice (kg/kg)
    real(kind=kind_phys), pointer :: evap_ice(:)     => null() !< kinematic surface upward latent heat flux over ice (m/s)
    real(kind=kind_phys), pointer :: hflx_ice(:)     => null() !< kinematic surface upward sensible heat flux over ice (Km/s)
    real(kind=kind_phys), pointer :: chh_ice(:)      => null() !< thermal exchange coefficient over ice (kg/m2s)
    real(kind=kind_phys), pointer :: cmm_ice(:)      => null() !< momentum exchange coefficient over ice (m/s)

    ! others
    real(kind=kind_phys), pointer :: z01d(:)         => null() !< perturbation of momentum roughness length
    real(kind=kind_phys), pointer :: zt1d(:)         => null() !< perturbation of heat to momentum roughness length ratio
    logical,              pointer :: flag_guess(:)   => null() !< flag for guess run
    real(kind=kind_phys), pointer :: rb(:)           => null() !< bulk Richardson number at the surface
    real(kind=kind_phys), pointer :: fh2(:)          => null() !< Monin-Obukhov similarity parameter for heat at 2m
    real(kind=kind_phys), pointer :: fm10(:)         => null() !< Monin-Obukhov similarity parameter for momentum at 10m
    real(kind=kind_phys), pointer :: cdq(:)          => null() !< surface exchange coeff heat & moisture
    real(kind=kind_phys), pointer :: cd(:)           => null() !< surface exchange coeff for momentum
    real(kind=kind_phys), pointer :: hffac(:)        => null() !< surface upward sensible heat flux reduction factor from canopy heat storage
    real(kind=kind_phys), pointer :: stress(:)       => null() !< surface wind stress
    real(kind=kind_phys), pointer :: gflx(:)         => null() !< soil heat flux
    real(kind=kind_phys), pointer :: ep1d(:)         => null() !< surface upward potential latent heat flux
    contains
      procedure :: create => interstitial_create !< allocate array data
      procedure :: phys_reset  => interstitial_phys_reset !<   reset array data for physics
  end type MED_interstitial_type

!! \section arg_table_MED_control_type
!! \htmlinclude MED_control_type.html
!!
  type MED_control_type
    logical                       :: lseaspray                 !< flag for sea spray parameterization
    logical                       :: use_med_flux              !< flag for using atmosphere-ocean fluxes form mediator
    integer                       :: ivegsrc                   !< land use dataset choice 0 => USGS, 1 => IGBP, 2 => UMD
    integer                       :: lsm                       !< flag for land surface model
    integer                       :: lsm_noahmp                !< flag for NOAH MP land surface model
    logical                       :: redrag                    !< flag for reduced drag coeff. over sea
    integer                       :: sfc_z0_type               !< surface roughness options over water
    integer                       :: icplocn2atm               !< flag controlling whether to consider ocean current in air-sea flux calculation 
    logical                       :: thsfc_loc                 !< flag for reference pressure in theta calculation
    integer                       :: nstf_name(5)              !< NSSTM flag: off/uncoupled/coupled=0/1/2
    integer                       :: lkm                       !< 0 = no lake model, 1 = lake model, 2 = lake & nsst on lake points
    logical                       :: first_time_step           !< flag signaling first time step for time integration routine
    logical                       :: frac_grid                 !< flag for fractional grid
    logical                       :: cplwav2atm                !< default no wav->atm coupling
    logical                       :: restart                   !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
    logical                       :: cplice                    !< default no cplice collection (used together with cplflx)
    logical                       :: cplflx                    !< flag controlling cplflx collection (default off)
    logical                       :: cpl_fire                  !< flag controlling fire behavior collection (default off)
    integer                       :: kdt                       !< current forecast iteration
    real(kind=kind_phys)          :: min_lakeice               !< minimum lake ice value
    real(kind=kind_phys)          :: min_seaice                !< minimum sea ice value
    real(kind=kind_phys)          :: huge                      !< definition of NetCDF float FillValue
    logical                       :: lheatstrg                 !< flag for canopy heat storage parameterization
    real(kind=kind_phys)          :: h0facu                    !< canopy heat storage factor for sensible heat flux in unstable surface layer
    real(kind=kind_phys)          :: h0facs                    !< canopy heat storage factor for sensible heat flux in stable surface layer
    integer                       :: lsoil                     !< number of soil layers
    integer                       :: kice                      !< vertical loop extent for ice levels, start at 1
    integer                       :: lsm_ruc                   !< flag for RUC land surface model

    ! Lake variables
    logical                       :: frac_ice = .false.        !< flag for fractional ice when fractional grid is not in use
    logical                       :: use_lake2m = .false.      !< use 2m T & Q calculated by the lake model
    integer                       :: iopt_lake = 1             !< =1 flake, =2 clm lake
    integer                       :: iopt_lake_flake = 1
    integer                       :: iopt_lake_clm = 2

    logical                       :: diag_flux                 !< flag for flux method of 2-m diagnostics
    logical                       :: diag_log                  !< flag for log 2-m diagnostics
    contains
      procedure :: init  => control_initialize
  end type MED_control_type

!! \section arg_table_MED_coupling_type
!! \htmlinclude MED_coupling_type.html
!!
  type MED_coupling_type
    real(kind=kind_phys), pointer :: dtsfcin_med(:) => null() !< sfc latent heat flux over ocean
    real(kind=kind_phys), pointer :: dqsfcin_med(:) => null() !< sfc sensible heat flux over ocean
    contains
      procedure :: create  => coupling_create !< allocate array data
  end type MED_coupling_type

!! \section arg_table_MED_grid_type
!! \htmlinclude MED_grid_type.html
!!
  type MED_grid_type
    real(kind=kind_phys), pointer :: area(:)         => null() !< area of the grid cell
    real(kind=kind_phys), pointer :: xlat_d(:)       => null() !< latitude in degrees
    real(kind=kind_phys), pointer :: xlon_d(:)       => null() !< longtitude in degrees
    contains
      procedure :: create  => grid_create !< allocate array data
  end type MED_grid_type

!! \section arg_table_MED_sfcprop_type
!! \htmlinclude MED_sfcprop_type.html
!!
  type MED_sfcprop_type
    real(kind=kind_phys), pointer :: zorlw(:)        => null()  !< surface roughness length over water (cm)
    integer,              pointer :: vtype(:)        => null()  !< vegetation type
    real(kind=kind_phys), pointer :: shdmax(:)       => null()  !< max fractional coverage of green vegetation
    real(kind=kind_phys), pointer :: zorll(:)        => null()  !< surface roughness length over land (cm)
    real(kind=kind_phys), pointer :: zorli(:)        => null()  !< surface roughness length over ice (cm)
    real(kind=kind_phys), pointer :: zorlwav(:)      => null()  !< surface roughness length from wave model (cm)
    real(kind=kind_phys), pointer :: zorl(:)         => null()  !< surface roughness length (cm)
    real(kind=kind_phys), pointer :: slmsk(:)        => null()  !< sea/land mask array (sea:0,land:1,sea-ice:2)
    real(kind=kind_phys), pointer :: lakefrac(:)     => null()  !< lake fraction [0:1]
    real(kind=kind_phys), pointer :: lakedepth(:)    => null()  !< lake depth (m)
    real(kind=kind_phys), pointer :: landfrac(:)     => null()  !< fraction of horizontal grid area occupied by land
    real(kind=kind_phys), pointer :: snowd(:)        => null()  !< snow depth water equivalent in mm ; same as snwdph
    real(kind=kind_phys), pointer :: weasd(:)        => null()  !< water equiv of acc snow depth over land and sea ice
    real(kind=kind_phys), pointer :: tprcp(:)        => null()  !< total precipitation amount in each time step
    real(kind=kind_phys), pointer :: oceanfrac(:)    => null()  !< ocean fraction [0:1]
    real(kind=kind_phys), pointer :: fice(:)         => null()  !< ice fraction over open water
    real(kind=kind_phys), pointer :: hice(:)         => null()  !< sea ice thickness (m)
    real(kind=kind_phys), pointer :: tsfco(:)        => null()  !< sea surface temperature
    real(kind=kind_phys), pointer :: usfco(:)        => null()  !< sea surface ocean current (zonal) 
    real(kind=kind_phys), pointer :: vsfco(:)        => null()  !< sea surface ocean current (merdional)
    real(kind=kind_phys), pointer :: uustar(:)       => null()  !< boundary layer parameter
    real(kind=kind_phys), pointer :: tsfc(:)         => null()  !< surface skin temperature
    real(kind=kind_phys), pointer :: snodi(:)        => null()  !< water equivalent snow depth over ice (mm)
    real(kind=kind_phys), pointer :: snodl(:)        => null()  !< water equivalent snow depth over land (mm)
    real(kind=kind_phys), pointer :: qss(:)          => null()  !< surface air saturation specific humidity (kg/kg)
    real(kind=kind_phys), pointer :: weasdi(:)       => null()  !< water equiv of acc snow depth over ice (mm)
    real(kind=kind_phys), pointer :: weasdl(:)       => null()  !< water equiv of acc snow depth over land (mm)
    real(kind=kind_phys), pointer :: ffhh(:)         => null()  !< Monin-Obukhov similarity function for heat
    real(kind=kind_phys), pointer :: ffmm(:)         => null()  !< Monin-Obukhov similarity function for momentum
    real(kind=kind_phys), pointer :: evap(:)         => null()  !< kinematic surface upward latent heat flux (kg kg-1 m s-1)
    real(kind=kind_phys), pointer :: evap_fire(:)    => null()  !< kinematic surface upward latent heat flux of fire (kg kg-1 m s-1)
    real(kind=kind_phys), pointer :: hflx(:)         => null()  !< kinematic surface upward sensible heat flux (K m/s)
    real(kind=kind_phys), pointer :: hflx_fire(:)    => null()  !< kinematic surface upward sensible heat flux of fire (K m/s)
    real(kind=kind_phys), pointer :: tiice(:,:)      => null()  !< sea ice internal temperature
    real(kind=kind_phys), pointer :: t2m(:)          => null()  !< temperature at 2 m
    real(kind=kind_phys), pointer :: q2m(:)          => null()  !< specific humidity at 2 m
    real(kind=kind_phys), pointer :: f10m(:)         => null()  !< ratio of sigma level 1 wind and 10m wind
    contains
      procedure :: create  => sfcprop_create !< allocate array data
  end type MED_sfcprop_type

!! \section arg_table_MED_diag_type
!! \htmlinclude MED_diag_type.html
!!
  type MED_diag_type
    real(kind=kind_phys), pointer :: chh(:)          => null()  !< thermal exchange coefficient (kg m-2 s-1)
    real(kind=kind_phys), pointer :: cmm(:)          => null()  !< momentum exchange coefficient (m/s)
    real(kind=kind_phys), pointer :: dpt2m(:)        => null()  !< 2-m dewpoint (K)
    contains
      procedure :: create  => diag_create    !< allocate array data
  end type MED_diag_type

  public MED_init_type
  public MED_statein_type
  public MED_coupling_type
  public MED_control_type
  public MED_interstitial_type
  public MED_grid_type
  public MED_sfcprop_type
  public MED_diag_type

  contains

  subroutine statein_create(statein, im, model)
    implicit none
    class(MED_statein_type) :: statein
    integer, intent(in)     :: im
    type(MED_control_type), intent(in) :: model

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
    allocate(statein%stc(im,model%lsoil))
    statein%stc = clear_val

  end subroutine statein_create

  subroutine stateout_create(stateout, im)
    implicit none
    class(MED_stateout_type) :: stateout
    integer, intent(in)     :: im

    allocate(stateout%gu0(im))
    stateout%gu0 = clear_val
    allocate(stateout%gv0(im))
    stateout%gv0 = clear_val
    allocate(stateout%gt0(im))
    stateout%gt0 = clear_val
    allocate(stateout%gq0(im))
    stateout%gq0 = clear_val

  end subroutine stateout_create

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
    allocate(interstitial%use_lake_model(im))
    interstitial%use_lake_model = 0
    allocate(interstitial%lake_t2m(im))
    interstitial%lake_t2m=-9999
    allocate(interstitial%lake_q2m(im))
    interstitial%lake_q2m=-9999
    allocate(interstitial%wind(im))
    interstitial%wind = huge
    allocate(interstitial%flag_iter(im))
    interstitial%flag_iter = .true.
    allocate(interstitial%flag_lakefreeze(im))
    interstitial%flag_lakefreeze = .false.
    allocate(interstitial%qss_water(im))
    interstitial%qss_water = huge
    allocate(interstitial%cmm_ice(im))
    interstitial%cmm_ice = huge
    allocate(interstitial%cmm_land(im))
    interstitial%cmm_land = huge
    allocate(interstitial%cmm_water(im))
    interstitial%cmm_water = huge
    allocate(interstitial%chh_ice(im))
    interstitial%chh_ice = huge
    allocate(interstitial%chh_land(im))
    interstitial%chh_land = huge
    allocate(interstitial%chh_water(im))
    interstitial%chh_water = huge
    allocate(interstitial%gflx_water(im))
    interstitial%gflx_water = clear_val
    allocate(interstitial%evap_water(im))
    interstitial%evap_water = huge
    allocate(interstitial%hflx_water(im))
    interstitial%hflx_water = huge
    allocate(interstitial%hflx_land(im))
    interstitial%hflx_land = huge
    allocate(interstitial%hflx_ice(im))
    interstitial%hflx_ice = huge
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
    allocate(interstitial%lake(im))
    interstitial%lake = .false.
    allocate(interstitial%tprcp_water(im))
    interstitial%tprcp_water = huge

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
    allocate(interstitial%frland(im))
    interstitial%frland = clear_val
    allocate(interstitial%tprcp_land(im))
    interstitial%tprcp_land = huge
    allocate(interstitial%qss_land(im))
    interstitial%qss_land = huge
    allocate(interstitial%evap_land(im))
    interstitial%evap_land = huge
    allocate(interstitial%hflxq(im))
    interstitial%hflxq = clear_val
    allocate(interstitial%ep1d_land(im))
    interstitial%ep1d_land = huge
    allocate(interstitial%gflx_land(im))
    interstitial%gflx_land = clear_val

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
    allocate(interstitial%flag_cice(im))
    interstitial%flag_cice = .false.
    allocate(interstitial%tprcp_ice(im))
    interstitial%tprcp_ice = huge
    allocate(interstitial%islmsk(im))
    interstitial%islmsk = 0
    allocate(interstitial%islmsk_cice(im))
    interstitial%islmsk_cice = 0
    allocate(interstitial%qss_ice(im))
    interstitial%qss_ice = huge
    allocate(interstitial%ep1d_ice(im))
    interstitial%ep1d_ice = huge
    allocate(interstitial%gflx_ice(im))
    interstitial%gflx_ice = clear_val
    allocate(interstitial%evap_ice(im))
    interstitial%evap_ice = huge

    ! others
    allocate(interstitial%z01d(im))
    interstitial%z01d = clear_val
    allocate(interstitial%zt1d(im))
    interstitial%zt1d = clear_val
    allocate(interstitial%flag_guess(im))
    interstitial%flag_guess = .false.
    allocate(interstitial%rb(im))
    interstitial%rb = clear_val
    allocate(interstitial%fh2(im))
    interstitial%fh2 = clear_val
    allocate(interstitial%fm10(im))
    interstitial%fm10 = clear_val
    allocate(interstitial%cdq(im))
    interstitial%cdq_water = clear_val
    allocate(interstitial%cd(im))
    interstitial%cd = clear_val
    allocate(interstitial%ep1d(im))
    interstitial%ep1d = clear_val
    allocate(interstitial%hffac(im))
    interstitial%hffac = clear_val
    allocate(interstitial%stress(im))
    interstitial%stress = clear_val
    allocate(interstitial%gflx(im))
    interstitial%gflx = clear_val

  end subroutine interstitial_create

  subroutine interstitial_phys_reset(interstitial)
    implicit none
    class(MED_interstitial_type) :: interstitial

    interstitial%cd = clear_val
    interstitial%cd_ice = huge
    interstitial%cd_land = huge
    interstitial%cd_water = huge
    interstitial%cdq = clear_val
    interstitial%cdq_ice = huge
    interstitial%cdq_land = huge
    interstitial%cdq_water = huge
    interstitial%chh_ice = huge
    interstitial%chh_land = huge
    interstitial%chh_water = huge
    interstitial%cmm_ice = huge
    interstitial%cmm_land = huge
    interstitial%cmm_water = huge
    interstitial%dry = .false.
    interstitial%ep1d = clear_val
    interstitial%ep1d_ice = huge
    interstitial%ep1d_land = huge
    interstitial%ep1d_water = huge
    interstitial%evap_water = huge
    interstitial%evap_land = huge
    interstitial%evap_ice = huge
    interstitial%ffhh_ice = huge
    interstitial%ffhh_land = huge
    interstitial%ffhh_water = huge
    interstitial%ffmm_ice = huge
    interstitial%ffmm_land = huge
    interstitial%ffmm_water = huge
    Interstitial%fh2 = clear_val
    interstitial%fh2_ice = huge
    interstitial%fh2_land = huge
    interstitial%fh2_water = huge
    Interstitial%fm10 = clear_val
    interstitial%flag_cice = .false.
    interstitial%flag_guess = .false.
    interstitial%flag_iter = .true.
    interstitial%flag_lakefreeze = .false.
    interstitial%fm10_ice = huge
    interstitial%fm10_land = huge
    interstitial%fm10_water = huge
    interstitial%frland = clear_val
    interstitial%gflx = clear_val
    interstitial%gflx_ice = clear_val
    interstitial%gflx_land = clear_val
    interstitial%gflx_water = clear_val
    interstitial%hffac = clear_val
    interstitial%hflx_ice = huge
    interstitial%hflx_land = huge
    interstitial%hflx_water = huge
    interstitial%hflxq = clear_val
    interstitial%icy = .false.
    interstitial%islmsk = 0
    interstitial%islmsk_cice = 0
    interstitial%lake = .false.
    interstitial%prslki = clear_val
    interstitial%rb = clear_val
    interstitial%qss_ice = huge
    interstitial%qss_land = huge
    interstitial%qss_water = huge
    interstitial%rb_ice = huge
    interstitial%rb_land = huge
    interstitial%rb_water = huge
    interstitial%sigmaf = clear_val
    interstitial%stress = clear_val
    interstitial%stress_ice = huge
    interstitial%stress_land = huge
    interstitial%stress_water = huge
    interstitial%tisfc = clear_val
    interstitial%tprcp_water = huge
    interstitial%tprcp_land = huge
    interstitial%tprcp_ice = huge
    interstitial%tsfc_water = huge
    interstitial%tsfcl = clear_val
    interstitial%tsurf_ice = huge
    interstitial%tsurf_land = huge
    interstitial%tsurf_water = huge
    interstitial%use_lake_model = 0
    interstitial%lake_t2m = -9999
    interstitial%lake_q2m = -9999
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
    model%icplocn2atm = 0
    model%thsfc_loc = .true.
    model%lsm = 1
    model%lsm_noahmp = 2
    model%nstf_name = (/0,0,1,0,5/)
    model%lkm = 0
    model%first_time_step = .true.
    model%frac_grid = .false.
    model%cplwav2atm = .false.
    model%restart = .false.
    model%cplice = .false.
    model%cplflx = .false.
    model%cpl_fire = .false.
    model%kdt = 0 ! nint(Model%fhour*con_hr/Model%dtp)
    model%min_lakeice = 0.15d0
    model%min_seaice = 1.0d-11
    model%huge = 9.9692099683868690e36
    model%lheatstrg = .false.
    model%h0facu = 0.25
    model%h0facs = 1.0
    model%lsoil = 4
    model%kice = 2
    model%lsm_ruc = 3
    model%frac_ice = .false.
    model%use_lake2m = .false.
    model%iopt_lake = 1
    model%iopt_lake_flake = 1
    model%iopt_lake_clm = 2
    model%diag_flux = .false.
    model%diag_log = .false.

  end subroutine control_initialize

  subroutine coupling_create(coupling, im)
    implicit none
    class(MED_coupling_type) :: coupling
    integer, intent(in)      :: im

    allocate(coupling%dtsfcin_med(im))
    coupling%dtsfcin_med = clear_val
    allocate(coupling%dqsfcin_med(im))
    coupling%dqsfcin_med = clear_val

  end subroutine coupling_create

  subroutine grid_create(grid, im)
    implicit none
    class(MED_grid_type) :: grid
    integer, intent(in)  :: im

    allocate(grid%area(im))
    grid%area = clear_val
    allocate(grid%xlat_d(im))
    grid%xlat_d = clear_val
    allocate(grid%xlon_d(im))
    grid%xlon_d = clear_val

  end subroutine grid_create

  subroutine sfcprop_create(sfcprop, im, model)
    implicit none
    class(MED_sfcprop_type) :: sfcprop
    integer, intent(in)  :: im
    type(MED_control_type), intent(in) :: model

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
    allocate(sfcprop%slmsk(im))
    sfcprop%slmsk = clear_val
    allocate(sfcprop%lakefrac(im))
    sfcprop%lakefrac = clear_val
    allocate(sfcprop%lakedepth(im))
    sfcprop%lakedepth = clear_val
    allocate(sfcprop%landfrac(im))
    sfcprop%landfrac = clear_val
    allocate(sfcprop%snowd(im))
    sfcprop%snowd = clear_val
    allocate(sfcprop%weasd(im))
    sfcprop%weasd = clear_val
    allocate(sfcprop%tprcp(im))
    sfcprop%tprcp = clear_val
    allocate(sfcprop%oceanfrac(im))
    sfcprop%oceanfrac = clear_val
    allocate(sfcprop%fice(im))
    sfcprop%fice = clear_val
    allocate(sfcprop%hice(im))
    sfcprop%hice = clear_val
    allocate(sfcprop%tsfco(im))
    sfcprop%tsfco = clear_val
    allocate(sfcprop%usfco(im))
    sfcprop%usfco = clear_val
    allocate(sfcprop%vsfco(im))
    sfcprop%vsfco = clear_val
    allocate(sfcprop%uustar(im))
    sfcprop%uustar = clear_val
    allocate(sfcprop%tsfc(im))
    sfcprop%tsfc = clear_val
    allocate(sfcprop%snodi(im))
    sfcprop%snodi = clear_val
    allocate(sfcprop%snodl(im))
    sfcprop%snodl = clear_val
    allocate(sfcprop%qss(im))
    sfcprop%qss = clear_val
    allocate(sfcprop%weasdi(im))
    sfcprop%weasdi = clear_val
    allocate(sfcprop%weasdl(im))
    sfcprop%weasdl = clear_val
    allocate(sfcprop%ffhh(im))
    sfcprop%ffhh = clear_val
    allocate(sfcprop%ffmm(im))
    sfcprop%ffmm = clear_val
    allocate(sfcprop%evap(im))
    sfcprop%evap = clear_val
    allocate(sfcprop%evap_fire(im))
    sfcprop%evap_fire = clear_val
    allocate(sfcprop%hflx(im))
    sfcprop%hflx = clear_val
    allocate(sfcprop%hflx_fire(im))
    sfcprop%hflx_fire = clear_val
    allocate(sfcprop%tiice(im,model%kice))
    sfcprop%tiice = clear_val
    allocate(sfcprop%t2m(im))
    sfcprop%t2m = clear_val
    allocate(sfcprop%q2m(im))
    sfcprop%q2m = clear_val
    allocate(sfcprop%f10m(im))
    sfcprop%f10m = clear_val

  end subroutine sfcprop_create

  subroutine diag_create(diag, im)
    implicit none
    class(MED_diag_type) :: diag
    integer, intent(in)  :: im

    allocate(diag%chh(im))
    diag%chh = clear_val
    allocate(diag%cmm(im))
    diag%cmm = clear_val
    allocate(diag%dpt2m(im))
    diag%dpt2m = clear_val

  end subroutine diag_create

end module MED_typedefs
