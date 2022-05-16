module flux_atmocn_ccpp_mod

  use ESMF,            only : operator(-), operator(/)
  use ESMF,            only : ESMF_GridComp, ESMF_Time, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF,            only : ESMF_Clock, ESMF_TimeInterval, ESMF_ClockGet
  use ESMF,            only : ESMF_GridCompGetInternalState, ESMF_LOGMSG_INFO
  use ESMF,            only : ESMF_LogWrite
  use NUOPC,           only : NUOPC_CompAttributeGet
  use NUOPC_Mediator,  only : NUOPC_MediatorGet

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
  use ufs_io_mod,      only : read_initial, read_restart, write_restart
  use med_kind_mod,    only : R8=>SHR_KIND_R8, CS=>SHR_KIND_CS
  use med_kind_mod,    only : CL=>SHR_KIND_CL
  use med_utils_mod,   only : chkerr => med_utils_chkerr
  use med_internalstate_mod, only : aoflux_ccpp_suite, logunit
  use med_internalstate_mod, only : InternalState, mastertask
  use med_constants_mod,     only : dbug_flag => med_constants_dbug_flag

  implicit none

  private ! default private

  public :: flux_atmOcn_ccpp ! computes atm/ocn fluxes

  integer, save           :: restart_freq
  integer, save           :: layout(2)
  real(r8), save          :: semis_water
  character(len=cs), save :: starttype
  character(len=cl), save :: ini_file
  character(len=cl), save :: rst_file
  character(len=cl), save :: mosaic_file
  character(len=cl), save :: input_dir
  character(len=1) , save :: listDel  = ":"

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
    type(ESMF_Clock)        :: mclock
    type(ESMF_Time)         :: currtime, starttime
    type(ESMF_TimeInterval) :: timeStep
    type(InternalState)     :: is_local
    integer                 :: n, rc
    real(r8)                :: spval
    logical                 :: isPresent, isSet
    character(len=cs)       :: cvalue, cname
    logical, save           :: first_call = .true.
    character(len=*), parameter :: subname=' (flux_atmOcn_ccpp) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! missing value
    if (present(missval)) then
       spval = missval
    else
       spval = shr_const_spval
    endif

    !----------------------
    ! Determine clock, starttime and currtime
    !----------------------

    call NUOPC_MediatorGet(gcomp, mediatorClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet(mclock, currtime=currTime, starttime=startTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! init CCPP and setup/allocate variables
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
       call physics%model%init()

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

       ! determine CCPP/host model specific options
       ! restart interval, set it to < 0 for no restart
       call NUOPC_CompAttributeGet(gcomp, name="ccpp_restart_interval", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          read(cvalue,*) restart_freq
       else
          restart_freq = 3600 ! write restart file every hour
       end if

       ! file name for restart
       call NUOPC_CompAttributeGet(gcomp, name='ccpp_restart_file', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          rst_file = trim(cvalue)
       else
          rst_file = 'unset'
       end if

       ! file name for initial conditions
       call NUOPC_CompAttributeGet(gcomp, name='ccpp_ini_file_prefix', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          ini_file = trim(cvalue)
       else
          ini_file = 'INPUT/sfc_data.tile'
       end if

       ! name of mosaic file that will be used to read tiled files
       call NUOPC_CompAttributeGet(gcomp, name='ccpp_ini_mosaic_file', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (isPresent .and. isSet) then
          mosaic_file = trim(cvalue)
       else
          if (trim(rst_file) == 'unset') then
             call ESMF_LogWrite(trim(subname)//': ccpp_ini_mosaic_file is required to read tiled initial condition!', ESMF_LOGMSG_INFO)
             rc = ESMF_FAILURE
             return
          end if
       end if

       ! input directory for tiled CS grid files 
       call NUOPC_CompAttributeGet(gcomp, name='ccpp_input_dir', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (isPresent .and. isSet) then
          input_dir = trim(cvalue)
       else
          input_dir = "INPUT/"
       end if

       ! layout to to read tiled CS grid files
       call NUOPC_CompAttributeGet(gcomp, name='ccpp_ini_layout', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          do n = 1, 2
             call string_listGetName(cvalue, n, cname, rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             read(cname,*) layout(n)
          end do
       else
          if (trim(rst_file) == 'unset') then
             call ESMF_LogWrite(trim(subname)//': ccpp_ini_layout is required to read tiled initial condition!', ESMF_LOGMSG_INFO)
             rc = ESMF_FAILURE
             return
          end if
       end if

       if (mastertask) then
          write(logunit,*) '========================================================'
          write(logunit,'(a,f5.2)') trim(subname)//' ccpp_phy_semis_water  = ', semis_water
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_lseaspray    = ', physics%model%lseaspray
          write(logunit,'(a,i)')    trim(subname)//' ccpp_phy_ivegsrc      = ', physics%model%ivegsrc
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_redrag       = ', physics%model%redrag
          write(logunit,'(a,i)')    trim(subname)//' ccpp_phy_lsm          = ', physics%model%lsm
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_frac_grid    = ', physics%model%frac_grid
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_restart      = ', physics%model%restart
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_cplice       = ', physics%model%cplice
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_cplflx       = ', physics%model%cplflx
          write(logunit,'(a,l)')    trim(subname)//' ccpp_phy_lheatstrg    = ', physics%model%lheatstrg
          write(logunit,'(a,i)')    trim(subname)//' ccpp_restart_interval = ', restart_freq
          write(logunit,'(a)')      trim(subname)//' ccpp_ini_file_prefix  = ', trim(ini_file)
          write(logunit,'(a)')      trim(subname)//' ccpp_ini_mosaic_file  = ', trim(mosaic_file)
          write(logunit,'(a)')      trim(subname)//' ccpp_input_dir        = ', trim(input_dir)
          write(logunit,'(a)')      trim(subname)//' ccpp_restart_file     = ', trim(rst_file)
          do n = 1, 2
             write(logunit,'(a,i,a,i2)') trim(subname)//' ccpp_ini_layout(',n,') = ', layout(n)
          end do 
          write(logunit,*) '========================================================'
       end if

       ! read initial condition/restart
       call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) starttype
       if (trim(starttype) == trim('startup')) then
          call read_initial(gcomp, ini_file, mosaic_file, input_dir, layout, rc)
       else
          call read_restart(gcomp, rst_file, rc)
       end if

       ! run CCPP init
       ! TODO: suite name need to be provided by ESMF config file
       call med_ccpp_driver_init(trim(aoflux_ccpp_suite))
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
    physics%model%kdt = ((currTime-StartTime)/timeStep)+1
    if (mastertask .and. dbug_flag > 5) then
       write(logunit,'(a,i)') 'kdt = ', physics%model%kdt
    end if

    ! reset physics variables, mimic GFS_suite_interstitial_phys_reset
    call physics%interstitial%phys_reset()

    ! set required variables to mimic GFS_surface_generic_pre
    ! TODO: the wind calculation in GFS_surface_generic_pre has cnvwind adjustment
    physics%interstitial%wind = sqrt(ubot(:)*ubot(:)+vbot(:)*vbot(:))
    physics%interstitial%prslki = physics%statein%prsik(:)/physics%statein%prslk(:)

    ! set required variables to mimic GFS_surface_composites_pre (assumes no ice) 
    physics%interstitial%uustar_water(:) = physics%sfcprop%uustar(:) 
    physics%sfcprop%tsfco(:) = ts(:)
    physics%sfcprop%tsfc(:) = ts(:)
    physics%interstitial%tsfc_water(:) = physics%sfcprop%tsfc(:)
    physics%interstitial%tsurf_water(:) = physics%sfcprop%tsfc(:)
    physics%sfcprop%zorlw(:) = physics%sfcprop%zorl(:)
    do n = 1, nMax
       physics%sfcprop%zorlw(n) = max(1.0e-5, min(1.0d0, physics%sfcprop%zorlw(n)))
    end do

    ! other variables
    if (.not. first_call) physics%sfcprop%qss(:) = qbot(:)
    physics%interstitial%qss_water(:)  = physics%sfcprop%qss(:)

    ! calculate wet flag and ocean fraction based on masking, assumes full oceean
    where (mask(:) /= 0)
       physics%interstitial%wet = .true.
       physics%sfcprop%oceanfrac = 1.0d0
    elsewhere
       physics%sfcprop%oceanfrac = 0.0d0
    end where

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

    ! write restart file
    call write_restart(gcomp, restart_freq, rc)

    ! set first call flag
    first_call = .false.

  end subroutine flux_atmOcn_ccpp

  !===============================================================================
  subroutine string_listGetName(list, k, name, rc)

    ! ----------------------------------------------
    ! Get name of k-th field in list
    ! It is adapted from CDEPS, shr_string_listGetName
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*)     , intent(in)  :: list    ! list/string
    integer          , intent(in)  :: k       ! index of field
    character(*)     , intent(out) :: name    ! k-th name in list
    integer          , intent(out) :: rc

    ! local variables
    integer :: i,n     ! generic indecies
    integer :: kFlds   ! number of fields in list
    integer :: i0,i1   ! name = list(i0:i1)
    character(*), parameter :: subName = '(shr_string_listGetName)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    !--- check that this is a valid index ---
    kFlds = string_listGetNum(list)
    if (k < 1 .or. kFlds < k) then
      call ESMF_LogWrite(trim(subname)//": ERROR invalid index ", ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
    end if

    !--- start with whole list, then remove fields before and after desired
    !field ---
    i0 = 1
    i1 = len_trim(list)

    !--- remove field names before desired field ---
    do n=2,k
       i = index(list(i0:i1),listDel)
       i0 = i0 + i
    end do

    !--- remove field names after desired field ---
    if ( k < kFlds ) then
       i = index(list(i0:i1),listDel)
       i1 = i0 + i - 2
    end if

    !--- copy result into output variable ---
    name = list(i0:i1)//"   "

  end subroutine string_listGetName

  !===============================================================================
  integer function string_listGetNum(str)

    ! ----------------------------------------------
    ! Get number of fields in a string list
    ! It is adapted from CDEPS, string_listGetNum
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*), intent(in) :: str   ! string to search

    ! local variables
    integer :: count ! counts occurances of char
    character(*), parameter :: subName = '(string_listGetNum)'
    ! ----------------------------------------------

    string_listGetNum = 0

    if (len_trim(str) > 0) then
       count = string_countChar(str,listDel)
       string_listGetNum = count + 1
    endif

  end function string_listGetNum

  !===============================================================================
  integer function string_countChar(str,char,rc)

    ! ----------------------------------------------
    ! Count number of occurances of a character
    ! It is adapted from CDEPS, string_countChar
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*), intent(in)       :: str   ! string to search
    character(1), intent(in)       :: char  ! char to search for
    integer, intent(out), optional :: rc    ! return code

    ! local variables
    integer :: count    ! counts occurances of char
    integer :: n        ! generic index
    character(*), parameter :: subName = '(string_countChar)'
    ! ----------------------------------------------

    count = 0
    do n = 1, len_trim(str)
      if (str(n:n) == char) count = count + 1
    end do
    string_countChar = count

  end function string_countChar
end module flux_atmocn_ccpp_mod
