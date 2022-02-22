module flux_atmocn_ccpp_mod

  use med_kind_mod,    only : R8=>SHR_KIND_R8
  use physcons,        only : p0 => con_p0
  use physcons,        only : cappa => con_rocp
  use physcons,        only : cp => con_cp
  use physcons,        only : hvap => con_hvap
  use physcons,        only : sbc => con_sbc
  use MED_data,        only : physics 
  use med_ccpp_driver, only : med_ccpp_driver_init
  use med_ccpp_driver, only : med_ccpp_driver_run
  use med_ccpp_driver, only : med_ccpp_driver_finalize
  use ufs_const_mod

  implicit none

  private ! default private

  public :: flux_atmOcn_ccpp ! computes atm/ocn fluxes

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
    integer :: n
    real(r8) :: spval, semis_water
    logical, save :: first_call = .true.
    character(len=*),parameter :: subname=' (flux_atmOcn_ccpp) '
    !---------------------------------------

    ! missing value
    if (present(missval)) then
       spval = missval
    else
       spval = shr_const_spval
    endif

    ! set up surface emissivity for lw radiation
    ! semis_wat is constant and set to 0.97 in setemis() call
    ! TODO: This could be a part of CCPP suite or provided by ESMF config
    semis_water = 0.97

    if (first_call) then
       ! allocate and initalize data structures
       call physics%statein%create(nMax,physics%model)
       call physics%interstitial%create(nMax)
       call physics%coupling%create(nMax)
       call physics%grid%create(nMax)
       call physics%sfcprop%create(nMax,physics%model)
       call physics%diag%create(nMax)

       ! initalize dimension 
       physics%init%im = nMax

       ! initalize model related parameters
       ! TODO: part of these need to be ingested from FV3 input.nml or configured through ESMF config file
       call physics%model%init()

       ! run CCPP init
       ! TODO: suite name need to be provided by ESMF config file
       call med_ccpp_driver_init('FV3_sfc_ocean')
       first_call = .false.
    end if

    ! fill in atmospheric forcing
    physics%statein%pgr(:)   = psfc(:)
    physics%statein%ugrs(:)  = ubot(:)
    physics%statein%vgrs(:)  = vbot(:)
    physics%statein%tgrs(:)  = tbot(:)
    physics%statein%qgrs(:)  = qbot(:)
    physics%statein%prsl(:)  = pbot(:)
    physics%statein%zlvl(:)  = zbot(:)
    physics%statein%prsik(:) = (psfc(:)/p0)**cappa
    physics%statein%prslk(:) = (pbot(:)/p0)**cappa
    physics%statein%u10m(:)  = usfc(:)
    physics%statein%v10m(:)  = vsfc(:)

    ! fill in grid related variables
    physics%grid%area(:) = garea(:)

    ! customization of host model options to calculate the fluxes
    ! TODO: this needs to be provided by config
    physics%model%lseaspray = .true.
    physics%model%ivegsrc = 1
    physics%model%redrag = .true.
    physics%model%lsm = 2
    physics%model%frac_grid = .true.
    physics%model%restart = .true.
    physics%model%cplice = .true.
    physics%model%cplflx = .true.
    physics%model%kdt = physics%model%kdt+1
    physics%model%lheatstrg = .true.

    ! reset physics variables
    call physics%interstitial%phys_reset()

    ! fill in required interstitial variables
    where (mask(:) /= 0)
       physics%interstitial%wet = .true.
    end where
    physics%interstitial%wind = sqrt(ubot(:)**2+vbot(:)**2)
    physics%interstitial%prslki = physics%statein%prsik(:)/physics%statein%prslk(:)
    physics%interstitial%tsurf_water = ts
    physics%interstitial%tsfc_water = ts
    physics%interstitial%qss_water = qbot

    ! fill in required sfcprop variables
    where (mask(:) /= 0)
       physics%sfcprop%oceanfrac = 1.0d0
    elsewhere
       physics%sfcprop%oceanfrac = 0.0d0
    end where
    physics%sfcprop%tsfco = ts
    physics%sfcprop%qss = qbot

    ! run CCPP physics
    ! TODO: suite name need to be provided by ESMF config file
    call med_ccpp_driver_run('FV3_sfc_ocean', 'physics')

    ! unit and sign conversion to be consistent with other flux scheme (CESM)
    do n = 1, nMax
       if (mask(n) /= 0) then
          sen(n)  = -1.0_r8*physics%interstitial%hflx_water(n)*rbot(n)*cp
          lat(n)  = -1.0_r8*physics%interstitial%evap_water(n)*rbot(n)*hvap
          lwup(n) = -1.0_r8*(semis_water*sbc*ts(n)**4+(1.0_r8-semis_water)*lwdn(n))
          evp(n)  = lat(n)/hvap
          taux(n) = rbot(n)*physics%interstitial%stress_water(n)*ubot(n)/physics%interstitial%wind(n)
          tauy(n) = rbot(n)*physics%interstitial%stress_water(n)*vbot(n)/physics%interstitial%wind(n)
          qref(n) = physics%interstitial%qss_water(n)
       else
          sen(n)  = spval
          lat(n)  = spval
          lwup(n) = spval
          evp(n)  = spval
          taux(n) = spval
          tauy(n) = spval
          qref(n) = spval
       end if
    end do

  end subroutine flux_atmOcn_ccpp

end module flux_atmocn_ccpp_mod
