module flux_atmocn_ccpp_mod

  use ESMF,            only : ESMF_GridComp, ESMF_SUCCESS
  use NUOPC,           only : NUOPC_CompAttributeGet

  use med_kind_mod,    only : R8=>SHR_KIND_R8, CS=>SHR_KIND_CS
  use physcons,        only : p0 => con_p0
  use physcons,        only : cappa => con_rocp
  use physcons,        only : cp => con_cp
  use physcons,        only : hvap => con_hvap
  use physcons,        only : sbc => con_sbc
  use MED_data,        only : physics 
  use med_utils_mod,   only : chkerr       => med_utils_chkerr
  use med_ccpp_driver, only : med_ccpp_driver_init
  use med_ccpp_driver, only : med_ccpp_driver_run
  use med_ccpp_driver, only : med_ccpp_driver_finalize
  use ufs_const_mod
  use med_internalstate_mod, only : aoflux_ccpp_suite

  implicit none

  private ! default private

  public :: flux_atmOcn_ccpp ! computes atm/ocn fluxes

  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine flux_atmOcn_ccpp(gcomp, mastertask, logunit, nMax, mask, psfc, pbot, &
             tbot, qbot, zbot, garea, ubot, usfc, vbot, vsfc, rbot, ts, lwdn, sen, lat, &
             lwup, evp, taux, tauy, qref, duu10n, missval)

    implicit none

    !--- input arguments --------------------------------
    type(ESMF_GridComp), intent(in) :: gcomp       ! gridded component
    logical , intent(in)  :: mastertask  ! master task
    integer , intent(in)  :: logunit     ! log file unit number
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
    real(r8), intent(out) :: sen(nMax)    ! heat flux: sensible            (W/m^2)
    real(r8), intent(out) :: lat(nMax)    ! heat flux: latent              (W/m^2)
    real(r8), intent(out) :: lwup(nMax)   ! heat flux: lw upward           (W/m^2)
    real(r8), intent(out) :: evp(nMax)    ! heat flux: evap                ((kg/s)/m^2)
    real(r8), intent(out) :: taux(nMax)   ! surface stress, zonal          (N)
    real(r8), intent(out) :: tauy(nMax)   ! surface stress, maridional     (N)
    real(r8), intent(out) :: qref(nMax)   ! diag: 2m ref humidity          (kg/kg)
    real(r8), intent(out) :: duu10n(nMax) ! diag: 10m wind speed squared (m/s)^2

    !--- local variables --------------------------------
    integer           :: n, rc
    real(r8)          :: spval
    logical           :: isPresent, isSet
    character(len=cs) :: cvalue
    real(r8), save    :: semis_water
    logical, save     :: first_call = .true.
    character(len=*), parameter :: subname=' (flux_atmOcn_ccpp) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! missing value
    if (present(missval)) then
       spval = missval
    else
       spval = shr_const_spval
    endif

    ! init CCPP and setup/allocate variables
    if (first_call) then
       ! determine CCPP/physics specific options
       ! semis_water, surface emissivity for lw radiation
       ! semis_wat is constant and set to 0.97 in setemis() call
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_semis_water", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       semis_water = 0.97
       if (isPresent .and. isSet) then
          read(cvalue,*) semis_water
       end if
       ! lseaspray
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_lseaspray", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%lseaspray = .true.
       if (isPresent .and. isSet) then
          if (trim(cvalue) .eq. '.false.' .or. trim(cvalue) .eq. 'false') physics%model%lseaspray = .false.
       end if
       ! ivegsrc
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_ivegsrc", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%ivegsrc = 1
       if (isPresent .and. isSet) then
          read(cvalue,*) physics%model%ivegsrc
       end if
       ! redrag 
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_redrag", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%redrag = .true.
       if (isPresent .and. isSet) then
          if (trim(cvalue) .eq. '.false.' .or. trim(cvalue) .eq. 'false') physics%model%redrag = .false.
       end if
       ! lsm
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_lsm", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%lsm = 1
       if (isPresent .and. isSet) then
          read(cvalue,*) physics%model%lsm
       end if
       ! frac_grid 
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_frac_grid", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%frac_grid = .true.
       if (isPresent .and. isSet) then
          if (trim(cvalue) .eq. '.false.' .or. trim(cvalue) .eq. 'false') physics%model%frac_grid = .false.
       end if
       ! restart
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_restart", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%restart = .true.
       if (isPresent .and. isSet) then
          if (trim(cvalue) .eq. '.false.' .or. trim(cvalue) .eq. 'false') physics%model%restart = .false.
       end if
       ! cplice
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_cplice", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%cplice = .true.
       if (isPresent .and. isSet) then
          if (trim(cvalue) .eq. '.false.' .or. trim(cvalue) .eq. 'false') physics%model%cplice = .false.
       end if
       ! cplflx
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_cplflx", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%cplflx = .true.
       if (isPresent .and. isSet) then
          if (trim(cvalue) .eq. '.false.' .or. trim(cvalue) .eq. 'false') physics%model%cplflx = .false.
       end if
       ! lheatstrg 
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_phy_lheatstrg", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       physics%model%lheatstrg = .true.
       if (isPresent .and. isSet) then
          if (trim(cvalue) .eq. '.false.' .or. trim(cvalue) .eq. 'false') physics%model%lheatstrg = .false.
       end if

       if (mastertask) then
          write(logunit,*) '========================================================'
          write(logunit,'(a,f5.2)') trim(subname)//' ccpp_phy_semis_water = ', semis_water
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_lseaspray   = ', physics%model%lseaspray
          write(logunit,'(a,i)')    trim(subname)//' ccpp_phy_ivegsrc     = ', physics%model%ivegsrc
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_redrag      = ', physics%model%redrag
          write(logunit,'(a,i)')    trim(subname)//' ccpp_phy_lsm         = ', physics%model%lsm
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_frac_grid   = ', physics%model%frac_grid
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_restart     = ', physics%model%restart
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_cplice      = ', physics%model%cplice
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_cplflx      = ', physics%model%cplflx
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_lheatstrg   = ', physics%model%lheatstrg
          write(logunit,*) '========================================================'
       end if

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
       call med_ccpp_driver_init(trim(aoflux_ccpp_suite))
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

    ! set counter
    physics%model%kdt = physics%model%kdt+1

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
    call med_ccpp_driver_run(trim(aoflux_ccpp_suite), 'physics')

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
          duu10n(n) = physics%interstitial%wind(n)*physics%interstitial%wind(n)
       else
          sen(n)  = spval
          lat(n)  = spval
          lwup(n) = spval
          evp(n)  = spval
          taux(n) = spval
          tauy(n) = spval
          qref(n) = spval
          duu10n(n) = spval
       end if
    end do

  end subroutine flux_atmOcn_ccpp

end module flux_atmocn_ccpp_mod
