module flux_atmocn_ccpp_mod

  use machine  , only: kp => kind_phys
  use funcphys , only: gpvs, fpvs, fpvsx
  use physcons , only: eps => con_eps
  use physcons , only: epsm1 => con_epsm1
  use physcons , only: grav => con_g
  use physcons , only: rvrdm1 => con_fvirt
  use physcons , only: cappa => con_rocp
  use physcons , only: hvap => con_hvap
  use physcons , only: cp => con_cp
  use physcons , only: rd => con_rd
  use physcons , only: rv => con_rv
  use physcons , only: hfus => con_hfus
  use physcons , only: p0 => con_p0
  use physcons , only: tice => con_tice
  use physcons , only: sbc => con_sbc
  use sfc_diff , only: sfc_diff_run
  use sfc_ocean, only: sfc_ocean_run
  use GFS_surface_composites_pre    , only: GFS_surface_composites_pre_run
  use GFS_surface_composites_post   , only: GFS_surface_composites_post_run
  use GFS_surface_loop_control_part1, only: GFS_surface_loop_control_part1_run
  use GFS_surface_loop_control_part2, only: GFS_surface_loop_control_part2_run
  use ufs_kind_mod
  use ufs_const_mod

  implicit none

  private ! default private

  public :: flux_atmOcn_ccpp ! computes atm/ocn fluxes

  !--- rename kinds for local readability only ---
  integer,parameter :: r8 = SHR_KIND_R8  ! 8 byte real

  !--- variables that need to carried through the iterations ---
  real(kp), allocatable, dimension(:) :: z0rl        , z0rl_wav  ,            &
                                         z0rl_wat    , z0rl_lnd  , z0rl_ice , &
                                         ustar       , fm        , fh       , &
                                         fm10        , hflx      , evap 

!===============================================================================
contains
!===============================================================================

  subroutine flux_atmOcn_ccpp(nMax, mask, psfc, pbot, tbot, qbot, zbot, &
             garea, ubot, usfc, vbot, vsfc, rbot, ts, lwdn, sen, lat, &
             lwup, evp, taux, tauy, qref, missval)

    implicit none

    !--- input arguments --------------------------------
    integer , intent(in)  :: nMax        ! data vector length
    integer , intent(in)  :: mask (nMax) ! ocn domain mask
    real(r8), intent(in)  :: psfc(nMax)  ! atm P (surface)                (Pa)
    real(r8), intent(in)  :: pbot(nMax)  ! atm P (bottom)                 (Pa)
    real(r8), intent(in)  :: tbot(nMax)  ! atm T (bottom)                 (K)
    real(r8), intent(in)  :: qbot(nMax)  ! atm specific humidity (bottom) (kg/kg)
    real(r8), intent(in)  :: zbot(nMax)  ! atm level height               (m)
    real(r8), intent(in)  :: garea(nMax) ! grid area                      (m^2)
    real(r8), intent(in)  :: ubot(nMax)  ! atm u wind (bottom)            (m/s)
    real(r8), intent(in)  :: usfc(nMax)  ! atm u wind (surface)           (m/s)
    real(r8), intent(in)  :: vbot(nMax)  ! atm v wind (bottom)            (m/s)    
    real(r8), intent(in)  :: vsfc(nMax)  ! atm v wind (surface)           (m/s)    
    real(r8), intent(in)  :: rbot(nMax)  ! atm density                    (kg/m^3)    
    real(r8), intent(in)  :: lwdn(nMax)  ! atm lw downward                (W/m^2)
    real(r8), intent(in)  :: ts(nMax)    ! ocn surface temperature        (K)
    real(r8), intent(in), optional :: missval ! masked value

    !--- output arguments -------------------------------
    real(r8), intent(out) :: sen(nMax)   ! heat flux: sensible            (W/m^2)
    real(r8), intent(out) :: lat(nMax)   ! heat flux: latent              (W/m^2)
    real(r8), intent(out) :: lwup(nMax)  ! heat flux: lw upward           (W/m^2)
    real(r8), intent(out) :: evp(nMax)   ! heat flux: evap                ((kg/s)/m^2)
    real(r8), intent(out) :: taux(nMax)  ! surface stress, zonal          (N)
    real(r8), intent(out) :: tauy(nMax)  ! surface stress, maridional     (N)
    real(r8), intent(out) :: qref(nMax)  ! diag: 2m ref humidity          (kg/kg)

    !--- local variables --------------------------------
    integer                     :: n           , iter      , ivegsrc   , &
                                   sfc_z0_type , errflg    , nstf_name1, &
                                   lkm         , nthreads  , kice      , &
                                   km          , lsm       , lsm_noahmp, &
                                   lsm_ruc
    real(kp)                    :: spval       , cpinv     , hvapi     , &
                                   elocp       , rch       , tem       , &
                                   min_lakeice , min_seaice, tgice     , &
                                   h0facu      , h0facs
    logical                     :: redrag      , thsfc_loc , lseaspray , &
                                   flag_restart, frac_grid , cplflx    , &
                                   cplice      , cplwav2atm, lheatstrg , &
                                   use_med_flux
    character(len=1024)         :: errmsg
    integer, dimension(nMax)    :: vegtype     , islmsk    , islmsk_cice 
    real(kp), dimension(nMax)   :: prsl1       , prslki    , prsik1    , &
                                   prslk1      , wind      , sigmaf    , &
                                   shdmax      , z0pert    , ztpert    , &
                                   tsurf_wat   , tsurf_lnd , tsurf_ice , &
                                   zvfun       , cm        , cm_wat    , &
                                   cm_lnd      , cm_ice    , ch        , &
                                   ch_wat      , ch_lnd    , ch_ice    , &
                                   rb          , rb_wat    , rb_lnd    , &
                                   rb_ice      , stress    ,             &
                                   stress_wat  , stress_lnd, stress_ice, &
                                   ztmax_wat   , ztmax_lnd , ztmax_ice , &
                                   landfrac    , lakefrac  , lakedepth , &
                                   oceanfrac   , frland    , hice      , &
                                   cice        , snowd     , snowd_lnd , &
                                   snowd_ice   , tprcp     , tprcp_wat , &
                                   tprcp_lnd   , tprcp_ice , weasd     , &
                                   weasd_lnd   , weasd_ice , hflxq     , &
                                   tsfco       , tsfcl     , tisfc     , &
                                   slmsk       , hffac     , vfrac     , &
                                   qss         ,                         &
                                   qss_wat     , qss_lnd   , qss_ice   , &
                                   tskin       ,                         &
                                   tskin_wat   , tskin_lnd , tskin_ice , &
                                   ustar_wat   , ustar_lnd , ustar_ice , &
                                   fm_wat      , fm_lnd    , fm_ice    , &
                                   fh_wat      , fh_lnd    , fh_ice    , &
                                   fm10_wat    , fm10_lnd  , fm10_ice  , &
                                   fh2         ,                         & 
                                   fh2_wat     , fh2_lnd   , fh2_ice   , &
                                   cmm         ,                         &
                                   cmm_wat     , cmm_lnd   , cmm_ice   , &
                                   chh         ,                         &
                                   chh_wat     , chh_lnd   , chh_ice   , &
                                   gflx        ,                         &
                                   gflx_wat    , gflx_lnd  , gflx_ice  , &
                                   ep1d        ,                         &
                                   ep1d_wat    , ep1d_lnd  , ep1d_ice  , &
                                   evap_wat    , evap_lnd  , evap_ice  , &
                                   hflx_wat    , hflx_lnd  , hflx_ice  , &
                                   tsfc        ,                         &
                                   tsfc_wat    , tsfc_lnd  , tsfc_ice  , &
                                   semis_rad   , emis_lnd  , emis_ice  , &
                                   semis_wat   , semis_lnd , semis_ice , &
                                   dqsfc       , dtsfc
    real(kp), dimension(nMax,1) :: tiice       , stc
    !integer                     :: naux2d
    !real(kp), dimension(nMax,2) :: aux2d
    logical, dimension(nMax)    :: flag_iter   , flag_guess, use_flake , &
                                   wet         , dry       , icy       , &
                                   flag_cice   , lake

    !--- local variables that are carried out -----------
    logical, save               :: flag_init = .true.
    integer, save               :: kdt = 0

    !--- parameters -------------------------------------
    real(kp), parameter         :: huge = 9.9692099683868690E36
    real(kp), parameter         :: zero = 0.0_kp
    real(kp), parameter         :: clear_val = zero

    !--- missing value --- 
    if (present(missval)) then
       spval = missval
    else
       spval = shr_const_spval
    endif

    !--- addtional constants ---
    cpinv = 1.0_kp/cp
    hvapi = 1.0_kp/hvap
    elocp = hvap/cp
 
    !--- compute some needed quantities ---
    wind(:) = sqrt(ubot(:)**2+vbot(:)**2)

    !--- compute dimensionless exner function ---
    prslk1(:) = (pbot(:)/p0)**cappa ! dimensionless_exner_function_at_surface_adjacent_layer
    prsik1(:) = (psfc(:)/p0)**cappa ! surface_dimensionless_exner_function
    prslki(:) = prsik1(:)/prslk1(:) ! ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer

    !--- initialization of variables ---
    kice          = 1              ! vertical_dimension_of_sea_ice
    km            = 1              ! vertical_dimension_of_soil
    tiice(:,:)    = 0.0_kp         ! temperature_in_ice_layer
    lheatstrg     = .true.         ! flag_for_canopy_heat_storage_in_land_surface_scheme
    h0facu        = 0.25_kp        ! multiplicative_tuning_parameter_for_reduced_surface_heat_fluxes_due_to_canopy_heat_storage
    h0facs        = 1.0            ! multiplicative_tuning_parameter_for_reduced_latent_heat_flux_due_to_canopy_heat_storage
    hflxq(:)      = 0.0_kp         ! kinematic_surface_upward_sensible_heat_flux_reduced_by_surface_roughness_and_vegetation
    hffac(:)      = 0.0_kp         ! surface_upward_sensible_heat_flux_reduction_factor
    stc(:,:)      = 0.0_kp         ! soil_temperature

    flag_restart  = .true.         ! flag_for_restart, restart run
    lkm           = 0              ! control_for_lake_surface_scheme
    frac_grid     = .true.         ! flag_for_fractional_landmask
    flag_cice(:)  = .true.         ! flag_for_cice
    cplflx        = .true.         ! flag_for_surface_flux_coupling
    cplice        = .true.         ! flag_for_sea_ice_coupling
    cplwav2atm    = .false.        ! flag_for_one_way_ocean_wave_coupling_to_atmosphere
    where (mask(:) /= 0)
    landfrac(:)   = 0.0_kp         ! land_area_fraction
    elsewhere
    landfrac(:)   = 1.0_kp         ! land_area_fraction
    end where 
    lakefrac(:)   = 0.0_kp         ! lake_area_fraction
    lakedepth(:)  = 0.0_kp         ! lake_depth
    where (mask(:) /= 0)
    oceanfrac(:)  = 1.0_kp         ! sea_area_fraction
    elsewhere
    oceanfrac(:)  = 0.0_kp         ! sea_area_fraction
    end where 
    frland(:)     = 0.0_kp         ! land_area_fraction_for_microphysics
    dry(:)        = .false.        ! flag_nonzero_land_surface_fraction, no land
    icy(:)        = .false.        ! flag_nonzero_sea_ice_surface_fraction, no sea-ice
    lake(:)       = .false.        ! flag_nonzero_lake_surface_fraction
    use_flake(:)  = .false.        ! flag_for_using_flake
    wet(:)        = .false.        ! flag_nonzero_wet_surface_fraction
    hice(:)       = 0.0_kp         ! sea_ice_thickness
    cice(:)       = 0.0_kp         ! sea_ice_area_fraction_of_sea_area_fraction

    if (flag_init) then
       allocate(z0rl(nMax))
       z0rl(:)    = 0.0_kp         ! surface_roughness_length
       allocate(z0rl_wat(nMax))
       z0rl_wat(:) = 0.0_kp        ! surface_roughness_length_over_water
       allocate(z0rl_lnd(nMax))
       z0rl_lnd(:) = 0.0_kp        ! surface_roughness_length_over_land
       allocate(z0rl_ice(nMax))
       z0rl_ice(:) = 0.0_kp        ! surface_roughness_length_over_ice
       allocate(z0rl_wav(nMax))
       z0rl_wav(:) = 0.0_kp        ! surface_roughness_length_from_wave_model
    end if

    snowd(:)      = 0.0_kp         ! lwe_surface_snow
    snowd_lnd(:)  = 0.0_kp         ! surface_snow_thickness_water_equivalent_over_land
    snowd_ice(:)  = 0.0_kp         ! surface_snow_thickness_water_equivalent_over_ice
    tprcp(:)      = 0.0_kp         ! nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep
    tprcp_wat(:)  = 0.0_kp         ! nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_water
    tprcp_lnd(:)  = 0.0_kp         ! nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_land
    tprcp_ice(:)  = 0.0_kp         ! nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_ice
    
    if (flag_init) then
       allocate(ustar(nMax))
       ustar(:)   = 0.0_kp         ! surface_friction_velocity
    end if

    ustar_wat(:)  = 0.0_kp         ! surface_friction_velocity_over_water
    ustar_lnd(:)  = 0.0_kp         ! surface_friction_velocity_over_land
    ustar_ice(:)  = 0.0_kp         ! surface_friction_velocity_over_ice
    weasd(:)      = 0.0_kp         ! lwe_thickness_of_surface_snow_amount
    weasd_lnd(:)  = 0.0_kp         ! water_equivalent_accumulated_snow_depth_over_land
    weasd_ice(:)  = 0.0_kp         ! water_equivalent_accumulated_snow_depth_over_ice
    tskin(:)      = 0.0_kp         ! surface_skin_temperature
    tskin_wat(:)  = 0.0_kp         ! surface_skin_temperature_over_water 
    tskin_lnd(:)  = 0.0_kp         ! surface_skin_temperature_over_land
    tskin_ice(:)  = 0.0_kp         ! surface_skin_temperature_over_ice
    tsfc(:)       = 0.0_kp         ! surface_skin_temperature
    tsfc_wat(:)   = 0.0_kp         ! surface_skin_temperature_over_water_interstitial
    tsfc_lnd(:)   = 0.0_kp         ! surface_skin_temperature_over_land_interstitial
    tsfc_ice(:)   = 0.0_kp         ! surface_skin_temperature_over_ice_interstitial
    tsfco(:)      = ts(:)          ! sea_surface_temperature
    tsurf_wat(:)  = 0.0_kp         ! surface_skin_temperature_after_iteration_over_water
    tsurf_lnd(:)  = 0.0_kp         ! surface_skin_temperature_after_iteration_over_land
    tsurf_ice(:)  = 0.0_kp         ! surface_skin_temperature_after_iteration_over_ice
    tisfc(:)      = 0.0_kp         ! sea_ice_temperature
    tgice         = tice           ! freezing_point_temperature_of_seawater
    islmsk(:)     = 0              ! sea_land_ice_mask, all sea 
    islmsk_cice(:) = 0             ! sea_land_ice_mask_cice, all sea
    slmsk(:)      = 0              ! area_type, all sea
    qss(:)        = qbot(:)        ! surface_specific_humidity ? not the lowest level
    qss_wat(:)    = qss(:)         ! surface_specific_humidity_over_water
    qss_lnd(:)    = 0.0_kp         ! surface_specific_humidity_over_land
    qss_ice(:)    = 0.0_kp         ! surface_specific_humidity_over_ice
    min_lakeice   = 0.15_kp        ! min_lake_ice_area_fraction
    min_seaice    = 1.0e-11_kp     ! min_sea_ice_area_fraction
    kdt           = kdt+1          ! index_of_timestep

    sigmaf(:)     = 0.0_kp         ! bounded_vegetation_area_fraction, no veg
    vegtype(:)    = 0              ! vegetation_type_classification
    shdmax(:)     = 0.0_kp         ! max_vegetation_area_fraction
    ivegsrc       = 1              ! control_for_vegetation_dataset, IGBP
    z0pert(:)     = 0.0_kp         ! perturbation_of_momentum_roughness_length
    ztpert(:)     = 0.0_kp         ! perturbation_of_heat_to_momentum_roughness_length_ratio
    flag_iter(:)  = .true.         ! flag_for_iteration
    redrag        = .true.         ! flag_for_limited_surface_roughness_length_over_ocean, redrag in input.nml
    sfc_z0_type   = 0              ! flag_for_surface_roughness_option_over_water, no change
    thsfc_loc     = .true.         ! flag_for_reference_pressure_theta
    cm(:)         = 0.0_kp         ! surface_drag_coefficient_for_momentum
    cm_wat(:)     = 0.0_kp         ! surface_drag_coefficient_for_momentum_in_air_over_water
    cm_lnd(:)     = 0.0_kp         ! surface_drag_coefficient_for_momentum_in_air_over_land
    cm_ice(:)     = 0.0_kp         ! surface_drag_coefficient_for_momentum_in_air_over_ice
    ch(:)         = 0.0_kp         ! surface_drag_coefficient_for_heat_and_moisture
    ch_wat(:)     = 0.0_kp         ! surface_drag_coefficient_for_heat_and_moisture_in_air_over_water
    ch_lnd(:)     = 0.0_kp         ! surface_drag_coefficient_for_heat_and_moisture_in_air_over_land
    ch_ice(:)     = 0.0_kp         ! surface_drag_coefficient_for_heat_and_moisture_in_air_over_ice
    rb(:)         = 0.0_kp         ! bulk_richardson_number_at_lowest_model_level
    rb_wat(:)     = 0.0_kp         ! bulk_richardson_number_at_lowest_model_level_over_water
    rb_lnd(:)     = 0.0_kp         ! bulk_richardson_number_at_lowest_model_level_over_land
    rb_ice(:)     = 0.0_kp         ! bulk_richardson_number_at_lowest_model_level_over_ice
    stress(:)     = 0.0_kp         ! surface_wind_stress
    stress_wat(:) = 0.0_kp         ! surface_wind_stress_over_water
    stress_lnd(:) = 0.0_kp         ! surface_wind_stress_over_land
    stress_ice(:) = 0.0_kp         ! surface_wind_stress_over_ice

    if (flag_init) then
       allocate(fm(nMax))
       fm(:)      = 0.0_kp        ! Monin_Obukhov_similarity_function_for_momentum
    end if

    fm_wat(:)     = 0.0_kp        ! Monin_Obukhov_similarity_function_for_momentum_over_water
    fm_lnd(:)     = 0.0_kp        ! Monin_Obukhov_similarity_function_for_momentum_over_land
    fm_ice(:)     = 0.0_kp        ! Monin_Obukhov_similarity_function_for_momentum_over_ice

    if (flag_init) then
       allocate(fh(nMax))
       fh(:)      = 0.0_kp        ! Monin_Obukhov_similarity_function_for_heat
    end if

    fh_wat(:)     = 0.0_kp        ! Monin_Obukhov_similarity_function_for_heat_over_water
    fh_lnd(:)     = 0.0_kp        ! Monin_Obukhov_similarity_function_for_heat_over_land
    fh_ice(:)     = 0.0_kp        ! Monin_Obukhov_similarity_function_for_heat_over_ice

    if (flag_init) then
       allocate(fm10(nMax))
       fm10(:)    = 0.0_kp        ! Monin_Obukhov_similarity_function_for_momentum
    end if

    fm10_wat(:)   = 0.0_kp        ! Monin_Obukhov_similarity_function_for_momentum_at_10m_over_water
    fm10_lnd(:)   = 0.0_kp        ! Monin_Obukhov_similarity_function_for_momentum_at_10m_over_land
    fm10_ice(:)   = 0.0_kp        ! Monin_Obukhov_similarity_function_for_momentum_at_10m_over_ice
    fh2(:)       = 0.0_kp         ! Monin_Obukhov_similarity_function_for_heat
    fh2_wat(:)   = 0.0_kp         ! Monin_Obukhov_similarity_function_for_heat_at_2m_over_water
    fh2_lnd(:)   = 0.0_kp         ! Monin_Obukhov_similarity_function_for_heat_at_2m_over_land
    fh2_ice(:)   = 0.0_kp         ! Monin_Obukhov_similarity_function_for_heat_at_2m_over_ice
    ztmax_wat(:) = 0.0_kp         ! bounded_surface_roughness_length_for_heat_over_water
    ztmax_lnd(:) = 0.0_kp         ! bounded_surface_roughness_length_for_heat_over_land
    ztmax_ice(:) = 0.0_kp         ! bounded_surface_roughness_length_for_heat_over_ice
    zvfun(:)     = 0.0_kp         ! function_of_surface_roughness_length_and_green_vegetation_fraction

    lseaspray    = .true.         ! flag_for_sea_spray
    cmm(:)       = 0.0_kp         ! surface_drag_wind_speed_for_momentum         
    cmm_wat(:)   = 0.0_kp         ! surface_drag_wind_speed_for_momentum_in_air_over_water
    cmm_lnd(:)   = 0.0_kp         ! surface_drag_wind_speed_for_momentum_in_air_over_land
    cmm_ice(:)   = 0.0_kp         ! surface_drag_wind_speed_for_momentum_in_air_over_ice
    chh(:)       = 0.0_kp         ! surface_drag_mass_flux_for_heat_and_moisture 
    chh_wat(:)   = 0.0_kp         ! surface_drag_mass_flux_for_heat_and_moisture_in_air_over_water
    chh_lnd(:)   = 0.0_kp         ! surface_drag_mass_flux_for_heat_and_moisture_in_air_over_land
    chh_ice(:)   = 0.0_kp         ! surface_drag_mass_flux_for_heat_and_moisture_in_air_over_ice
    gflx(:)      = 0.0_kp         ! upward_heat_flux_in_soil
    gflx_wat(:)  = 0.0_kp         ! upward_heat_flux_in_soil_over_water
    gflx_lnd(:)  = 0.0_kp         ! upward_heat_flux_in_soil_over_lnd
    gflx_ice(:)  = 0.0_kp         ! upward_heat_flux_in_soil_over_ice
    use_med_flux = .false.        ! flag_for_mediator_atmosphere_ocean_fluxes
    dqsfc(:)     = 0.0_kp         ! surface_upward_latent_heat_flux_over_ocean_from_coupled_process
    dtsfc(:)     = 0.0_kp         ! surface_upward_sensible_heat_flux_over_ocean_from_coupled_process

    if (flag_init) then
       allocate(evap(nMax))
       evap(:)   = 0.0_kp         ! kinematic_surface_upward_latent_heat_flux
    end if

    evap_wat(:)  = 0.0_kp         ! kinematic_surface_upward_latent_heat_flux_over_water
    evap_lnd(:)  = 0.0_kp         ! kinematic_surface_upward_latent_heat_flux_over_land
    evap_ice(:)  = 0.0_kp         ! kinematic_surface_upward_latent_heat_flux_over_ice

    if (flag_init) then
       allocate(hflx(nMax))
       hflx(:)   = 0.0_kp         ! kinematic_surface_upward_sensible_heat_flux
    end if

    hflx_wat(:)  = 0.0_kp         ! kinematic_surface_upward_sensible_heat_flux_over_water
    hflx_lnd(:)  = 0.0_kp         ! kinematic_surface_upward_sensible_heat_flux_over_land
    hflx_ice(:)  = 0.0_kp         ! kinematic_surface_upward_sensible_heat_flux_over_ice

    ep1d(:)      = 0.0_kp         ! surface_upward_potential_latent_heat_flux
    ep1d_wat(:)  = 0.0_kp         ! surface_upward_potential_latent_heat_flux_over_water
    ep1d_lnd(:)  = 0.0_kp         ! surface_upward_potential_latent_heat_flux_over_land
    ep1d_ice(:)  = 0.0_kp         ! surface_upward_potential_latent_heat_flux_over_ice

    lsm          = 2              ! control_for_land_surface_scheme 
    lsm_noahmp   = 2              ! identifier_for_noahmp_land_surface_scheme
    lsm_ruc      = 3              ! identifier_for_ruc_land_surface_scheme
    semis_rad(:) = 0.0_kp         ! surface_longwave_emissivity
    semis_lnd(:) = 0.0_kp         ! surface_longwave_emissivity_over_land_interstitial
    semis_ice(:) = 0.0_kp         ! surface_longwave_emissivity_over_ice_interstitial
    semis_wat(:) = 0.0_kp         ! surface_longwave_emissivity_over_water_interstitial
    emis_lnd(:)  = 0.0_kp         ! surface_longwave_emissivity_over_land
    emis_ice(:)  = 0.0_kp         ! surface_longwave_emissivity_over_ice

    !--- set up surface emissivity for lw radiation ---
    !--- semis_wat is constant and set to 0.97 in setemis() call --- 
    semis_wat(:) = 0.97

    !--- GFS surface scheme pre ---
    call GFS_surface_composites_pre_run( &
         nMax      , flag_init  , flag_restart, &
         lkm       , lsm        , lsm_noahmp  , &
         lsm_ruc   , frac_grid  , flag_cice   , &
         cplflx    , cplice     , cplwav2atm  , &
         landfrac  , lakefrac   , lakedepth   , &
         oceanfrac , frland     , dry         , &
         icy       , lake       , use_flake   , &
         wet       , hice       , cice        , &
         z0rl_wat  , z0rl_lnd   , z0rl_ice    , &
         snowd     , snowd_lnd  , snowd_ice   , &
         tprcp     ,                            &
         tprcp_wat , tprcp_lnd  , tprcp_ice   , &
         ustar     ,                            &
         ustar_wat , ustar_lnd  , ustar_ice   , &
         weasd     , weasd_lnd  , weasd_ice   , &
         ep1d_ice  , tskin      , tsfco       , &
         tskin_lnd , tskin_wat  , tskin_ice   , &
         tisfc     , tsurf_wat  , tsurf_lnd   , &
         tsurf_ice , gflx_ice   , tgice       , &
         islmsk    , islmsk_cice, slmsk       , &
         semis_rad , semis_wat  , semis_lnd   , &
         semis_ice , emis_lnd   , emis_ice    , &
         qss       , qss_wat    , qss_lnd     , &
         qss_ice   , min_lakeice, min_seaice  , &
         kdt       , errmsg     , errflg)

    !--- surface iteration loop ---
    do iter = 1, 2
       !--- calculate stability parameters ---
       call sfc_diff_run( &
            nMax      , rvrdm1     , eps         , &
            epsm1     , grav       , psfc        , &
            tbot      , qbot       , zbot        , &
            garea     , wind       , pbot        , &
            prslki    , prsik1     , prslk1      , &
            sigmaf    , vegtype    , shdmax      , &
            ivegsrc   , z0pert     , ztpert      , &
            flag_iter , redrag     , usfc        , &
            vsfc      , sfc_z0_type, wet         , &
            dry       , icy        , thsfc_loc   , &
            tskin_wat , tskin_lnd  , tskin_ice   , &
            tsurf_wat , tsurf_lnd  , tsurf_ice   , &
            z0rl_wat  , z0rl_lnd   , z0rl_ice    , &
            z0rl_wav  ,                            &
            ustar_wat , ustar_lnd  , ustar_ice   , &
            cm_wat    , cm_lnd     , cm_ice      , &
            ch_wat    , ch_lnd     , ch_ice      , &
            rb_wat    , rb_lnd     , rb_ice      , &
            stress_wat, stress_lnd , stress_ice  , &
            fm_wat    , fm_lnd     , fm_ice      , &
            fh_wat    , fh_lnd     , fh_ice      , &
            fm10_wat  , fm10_lnd   , fm10_ice    , &
            fh2_wat   , fh2_lnd    , fh2_ice     , &
            ztmax_wat , ztmax_lnd  , ztmax_ice   , &
            zvfun     , errmsg     , errflg)

       !--- update flag_guess ---
       call GFS_surface_loop_control_part1_run( &
            nMax       , iter      , wind        , &
            flag_guess , errmsg    , errflg)

       !--- calculate heat fluxes ---
       call sfc_ocean_run( &
            nMax        , hvap      , cp          , &
            rd          , eps       , epsm1       , &
            rvrdm1      , psfc      , ubot        , &
            vbot        , tbot      , qbot        , &
            tskin_wat   , cm_wat    , ch_wat      , &
            lseaspray   , fm_wat    , fm10_wat    , &
            pbot        , prslki    , wet         , &
            use_flake   , wind      , flag_iter   , &
            use_med_flux, dqsfc     , dtsfc       , &
            qss_wat     , cmm_wat   , chh_wat     , &
            gflx_wat    , evap_wat  , hflx_wat    , &
            ep1d_wat    , errmsg    , errflg)

       !--- update flag_guess and flag_iter ---
       call GFS_surface_loop_control_part2_run( &
            nMax       , iter      , wind        , &
            flag_guess , flag_iter , dry         , &
            wet        , icy       , nstf_name1  , &
            errmsg     , errflg)
    end do

    !--- GFS surface scheme post ---
    call GFS_surface_composites_post_run( &
         nMax      , kice       , km          , &
         rd        , rvrdm1     , cplflx      , &
         cplwav2atm, frac_grid  , flag_cice   , &
         thsfc_loc , islmsk     , dry         , &
         wet       , icy        , wind        , &
         tbot      , qbot       , pbot        , &
         landfrac  , lakefrac   , oceanfrac   , &
         z0rl      , z0rl_wat   , z0rl_lnd    , &
         z0rl_ice  , garea      , cm          , &
         cm_wat    , cm_lnd     , cm_ice      , &
         ch        , ch_wat     , ch_lnd      , &
         ch_ice    , rb         , rb_wat      , &
         rb_lnd    , rb_ice     , stress      , &
         stress_wat, stress_lnd , stress_ice  , &
         fm        , fm_wat     , fm_lnd      , &
         fm_ice    , fh         , fh_wat      , &
         fh_lnd    , fh_ice     , ustar       , &
         ustar_wat , ustar_lnd  , ustar_ice   , &
         fm10      , fm10_wat   , fm10_lnd    , &
         fm10_ice  , fh2        , fh2_wat     , &
         fh2_lnd   , fh2_ice    , tsurf_wat   , &
         tsurf_lnd , tsurf_ice  , cmm         , &
         cmm_wat   , cmm_lnd    , cmm_ice     , &
         chh       , chh_wat    , chh_lnd     , &
         chh_ice   , gflx       , gflx_wat    , &
         gflx_lnd  , gflx_ice   , ep1d        , &
         ep1d_wat  , ep1d_lnd   , ep1d_ice    , &
         weasd     , weasd_lnd  , weasd_ice   , &
         snowd     , snowd_lnd  , snowd_ice   , &
         tprcp     , tprcp_wat  , tprcp_lnd   , &
         tprcp_ice , evap       , evap_wat    , &
         evap_lnd  , evap_ice   , hflx        , &
         hflx_wat  , hflx_lnd   , hflx_ice    , &
         qss       , qss_wat    , qss_lnd     , &
         qss_ice   , tskin      , tsfco       , &
         tskin_lnd , tskin_wat  , tskin_ice   , &
         tisfc     , hice       , cice        , & 
         min_seaice,                            &
         tiice     , sigmaf     , zvfun       , &
         lheatstrg , h0facu     , h0facs      , &
         hflxq     , hffac      , stc         , &
         grav      , prsik1     , prslk1      , &
         prslki    , zbot       , ztmax_wat   , &
         ztmax_lnd , ztmax_ice  ,               &
         errmsg    , errflg)

    !--- unit conversion ---
    do n = 1, nMax
       if (mask(n) /= 0) then
          sen(n)  = -1.0_kp*hflx_wat(n)*rbot(n)*cp
          lat(n)  = -1.0_kp*evap_wat(n)*rbot(n)*hvap
          lwup(n) = semis_wat(n)*sbc*ts(n)**4+(1.0_r8-semis_wat(n))*lwdn(n)
          evp(n)  = lat(n)/hvap
          taux(n) = -1.0_kp*rbot(n)*stress(n)*ubot(n)/wind(n) 
          tauy(n) = -1.0_kp*rbot(n)*stress(n)*vbot(n)/wind(n)
          qref(n) = qss_wat(n)
       else
          sen(n)  = spval
          lat(n)  = spval
          lwup(n) = spval
          evap(n) = spval
          taux(n) = spval
          tauy(n) = spval
          qref(n) = spval
       end if
    end do

    flag_init = .false.

  end subroutine flux_atmOcn_ccpp

end module flux_atmocn_ccpp_mod
