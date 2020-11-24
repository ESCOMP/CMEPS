module med_phases_prep_glc_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing glc export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use NUOPC_Model           , only : NUOPC_ModelGet
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_VM, ESMF_VMGet, ESMF_VMAllReduce, ESMF_REDUCE_SUM, ESMF_REDUCE_MAX
  use ESMF                  , only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockGetAlarm, ESMF_ClockAdvance, ESMF_ClockGet
  use ESMF                  , only : ESMF_Time, ESMF_TimeGet
  use ESMF                  , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmSet, ESMF_AlarmGet
  use ESMF                  , only : ESMF_AlarmIsRinging, ESMF_AlarmRingerOff
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleIsCreated
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF                  , only : ESMF_Mesh, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
  use esmFlds               , only : compglc, complnd, mapbilnr, mapconsd, compname
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use med_map_mod           , only : med_map_routehandles_init, med_map_rh_is_created
  use med_map_mod           , only : med_map_field_normalized, med_map_field
  use med_methods_mod       , only : FB_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_reset     => med_methods_FB_reset
  use med_utils_mod         , only : chkerr       => med_utils_ChkErr
  use med_time_mod          , only : med_time_alarmInit
  use glc_elevclass_mod     , only : glc_get_num_elevation_classes
  use glc_elevclass_mod     , only : glc_get_elevation_classes, glc_get_elevation_class
  use glc_elevclass_mod     , only : glc_get_fractional_icecov
  use perf_mod              , only : t_startf, t_stopf
  use shr_const_mod         , only : shr_const_pi
  use shr_mpi_mod           , only : shr_mpi_sum

  implicit none
  private

  public  :: med_phases_prep_glc_init
  public  :: med_phases_prep_glc_accum
  public  :: med_phases_prep_glc_avg

  private :: map_lnd2glc
  private :: med_phases_prep_glc_renormalize_smb

  ! glc fields with multiple elevation classes: lnd->glc
  ! - fields sent from lnd->med to glc    ARE     IN multiple elevation classes
  ! - fields sent from med->glc from land ARE NOT IN multiple elevation classes
  ! Need to keep track of the lnd->med fields destined for glc in the FBlndAccum field bundle.

  ! Needed for standard lnd->glc mapping
  type(ESMF_FieldBundle) :: FBlndAccum_lnd
  type(ESMF_FieldBundle) :: FBlndAccum_glc
  integer                :: FBlndAccumCnt
  character(len=14)      :: fldnames_fr_lnd(3) = (/'Flgl_qice_elev','Sl_tsrf_elev  ','Sl_topo_elev  '/)
  character(len=14)      :: fldnames_to_glc(2) = (/'Flgl_qice     ','Sl_tsrf       '/)

  ! Whether to renormalize the SMB for conservation.
  ! Should be set to true for 2-way coupled runs with evolving ice sheets.
  ! Does not need to be true for 1-way coupling.
  logical :: smb_renormalize

  ! Needed if renormalize SMB
  type(ESMF_Field)      :: field_icemask_g
  type(ESMF_Field)      :: field_icemask_l
  type(ESMF_Field)      :: field_frac_g
  type(ESMF_Field)      :: field_frac_l
  type(ESMF_Field)      :: field_frac_g_ec
  type(ESMF_Field)      :: field_frac_l_ec
  type(ESMF_Field)      :: field_lfrac_g

  character(len=*), parameter :: qice_fieldname       = 'Flgl_qice' ! Name of flux field giving surface mass balance
  character(len=*), parameter :: Sg_frac_fieldname    = 'Sg_ice_covered'
  character(len=*), parameter :: Sg_topo_fieldname    = 'Sg_topo'
  character(len=*), parameter :: Sg_icemask_fieldname = 'Sg_icemask'

  integer          :: ungriddedCount ! this equals the number of elevation classes + 1 (for bare land)
  logical          :: init_prep_glc = .false.
  type(ESMF_Clock) :: prepglc_clock
  character(*), parameter :: u_FILE_u  = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_prep_glc_init(gcomp, rc)

    !---------------------------------------
    ! Create land accumulation field bundles on and and glc grid and initialize accumulation count
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)       :: is_local
    integer                   :: i,n,ncnt
    type(ESMF_Mesh)           :: lmesh_glc
    type(ESMF_Mesh)           :: lmesh_lnd
    type(ESMF_Field)          :: lfield
    real(r8), pointer         :: data2d_in(:,:) => null()
    real(r8), pointer         :: data2d_out(:,:) => null()
    real(r8), pointer         :: dataptr1d(:) => null()
    character(len=CS)         :: glc_renormalize_smb
    logical                   :: glc_coupled_fluxes
    integer                   :: lsize
    logical                   :: isPresent
    integer                   :: fieldCount
    type(ESMF_Field), pointer :: fieldlist(:) => null()
    integer                   :: ungriddedUBound_output(1) ! currently the size must equal 1 for rank 2 fieldds
    type(ESMF_Clock)          :: med_clock
    type(ESMF_ALARM)          :: glc_avg_alarm
    logical                   :: glc_present
    character(len=CS)         :: glc_avg_period
    integer                   :: glc_cpl_dt
    character(len=CS)         :: cvalue
    character(len=*),parameter  :: subname=' (med_phases_prep_glc_init) '
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! First prepglc_clock from mclock - THIS CALL DOES NOT COPY ALARMS
    call NUOPC_ModelGet(gcomp, modelClock=med_clock,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    prepglc_clock = ESMF_ClockCreate(med_clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set alarm glc averaging interval
    call NUOPC_CompAttributeGet(gcomp, name="glc_avg_period", value=glc_avg_period, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(glc_avg_period) == 'yearly') then
       call med_time_alarmInit(prepglc_clock, glc_avg_alarm, 'nyears', opt_n=1, alarmname='alarm_glc_avg', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) then
          write(logunit,'(a,i10)') trim(subname)//&
               ' created alarm with averaging period for export to glc is yearly'
       end if
    else if (trim(glc_avg_period) == 'glc_coupling_period') then
       call NUOPC_CompAttributeGet(gcomp, name="glc_cpl_dt", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) glc_cpl_dt
       call med_time_alarmInit(prepglc_clock, glc_avg_alarm, 'nseconds', opt_n=glc_cpl_dt, alarmname='alarm_glc_avg', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) then
          write(logunit,'(a,i10)') trim(subname)//&
               ' created alarm with averaging period for export to glc (in seconds) ',glc_cpl_dt
       end if
    else
       call ESMF_LogWrite(trim(subname)// ": ERROR glc_avg_period = "//trim(glc_avg_period)//" not supported", &
            ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       RETURN
    end if
    call ESMF_AlarmSet(glc_avg_alarm, clock=prepglc_clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Count the number of fields outside of scalar data
    !---------------------------------------

    if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(complnd,complnd))) then
       ncnt = 0
       call ESMF_LogWrite(trim(subname)//": FBImp(complnd,complnd) is not created", ESMF_LOGMSG_INFO)
    else
       ! The scalar field has been removed from all mediator field bundles - so determine ncnt for below
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldCount=ncnt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! Create field bundles for the fldnames_fr_lnd that have an
    ! undistributed dimension corresponding to elevation classes
    !---------------------------------------

    if (ncnt > 0) then

       ! Create accumulation field bundle from land on the land grid (including bare land)
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fldnames_fr_lnd(1), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound_output, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ungriddedCount = ungriddedUBound_output(1)

       ! TODO: check that ungriddedCount = glc_nec+1
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(fieldlist(fieldcount))
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldlist=fieldlist, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(fieldlist(1), mesh=lmesh_lnd, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       deallocate(fieldlist)

       FBlndAccum_lnd = ESMF_FieldBundleCreate(name='FBlndAccum_lnd', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(fldnames_fr_lnd)
          lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name=fldnames_fr_lnd(n), &
               meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBlndAccum_lnd, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite(trim(subname)//' adding field '//trim(fldnames_fr_lnd(n))//' to FBLndAccum_lnd', &
               ESMF_LOGMSG_INFO)
       end do
       call FB_reset(FBlndAccum_lnd, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Create accumulation field bundle from land on the glc grid
       ! Determine glc mesh from the mesh from the first export field to glc
       ! However FBlndAccum_glc has the fields fldnames_fr_lnd BUT ON the glc grid
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(fieldlist(fieldcount))
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldlist=fieldlist, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(fieldlist(1), mesh=lmesh_glc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       deallocate(fieldlist)

       FBlndAccum_glc = ESMF_FieldBundleCreate(name='FBlndAccum_glc', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(fldnames_fr_lnd)
          lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name=fldnames_fr_lnd(n), &
               meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBlndAccum_glc, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do
       call FB_reset(FBlndAccum_glc, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Create land fraction field on glc mesh (this is just needed for normalization mapping)
       field_lfrac_g = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Create route handle if it has not been created
       if (.not. med_map_RH_is_created(is_local%wrap%RH(complnd,compglc,:),mapbilnr,rc=rc)) then
          call ESMF_LogWrite(trim(subname)//" mapbilnr is not created for lnd->glc mapping", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       end if

       ! Determine if renormalize smb
       call NUOPC_CompAttributeGet(gcomp, name='glc_renormalize_smb', value=glc_renormalize_smb, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! TODO: talk to Bill Sacks to determine if this is the correct logic
       glc_coupled_fluxes = is_local%wrap%med_coupling_active(compglc,complnd)
       ! Note glc_coupled_fluxes should be false in the no_evolve cases
       ! Goes back to the zero-gcm fluxes variable - if zero-gcm fluxes is true than do not renormalize
       ! The user can set this to true in an evolve cases

       select case (glc_renormalize_smb)
       case ('on')
          smb_renormalize = .true.
       case ('off')
          smb_renormalize = .false.
       case ('on_if_glc_coupled_fluxes')
          if (.not. glc_coupled_fluxes) then
             ! Do not renormalize if med_coupling_active is not true for compglc->complnd
             ! In this case, conservation is not important
             smb_renormalize = .false.
          else
             smb_renormalize = .true.
          end if
       case default
          write(logunit,*) subname,' ERROR: unknown value for glc_renormalize_smb: ', trim(glc_renormalize_smb)
          call ESMF_LogWrite(trim(subname)//' ERROR: unknown value for glc_renormalize_smb: '// trim(glc_renormalize_smb), &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
          rc = ESMF_FAILURE
          return
       end select
       if (mastertask) then
          write(logunit,'(a,l4)') trim(subname)//' smb_renormalize is ',smb_renormalize
       end if

       ! -------------------------------
       ! If smb will be renormalized then...
       ! -------------------------------
       if (smb_renormalize) then

          ! ice mask without elevation classes on glc and lnd
          field_icemask_g = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          field_icemask_l = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! ice fraction without multiple elevation classes on glc and lnd
          field_frac_g = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          field_frac_l = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! ice fraction in multiple elevation classes on glc and lnd - NOTE that this includes bare land
          field_frac_g_ec = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          field_frac_l_ec = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Create route handle if it has not been created - this will be needed to map the fractions
          if (.not. med_map_RH_is_created(is_local%wrap%RH(compglc,complnd,:),mapconsd,rc=rc)) then
             call med_map_routehandles_init( compglc, complnd, &
                  FBSrc=is_local%wrap%FBImp(compglc,compglc), &
                  FBDst=is_local%wrap%FBImp(compglc,complnd), &
                  mapindex=mapconsd, &
                  RouteHandle=is_local%wrap%RH, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       end if
    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_init

  !================================================================================================
  subroutine med_phases_prep_glc_accum(gcomp, rc)

    !---------------------------------------
    ! Carry out accumulation for the land-ice (glc) component
    ! Accumulation and averaging is done on the land input field to glc on the land grid
    ! Mapping from the land to the glc grid is then done after the accumulated fields have been
    ! time averaged
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: lfield
    type(ESMF_Alarm)    :: alarm
    integer             :: i,n,ncnt
    real(r8), pointer   :: data2d_in(:,:) => null()
    real(r8), pointer   :: data2d_out(:,:) => null()
    character(len=*),parameter  :: subname=' (med_phases_prep_glc_accum) '
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    if (.not. init_prep_glc) then
       call med_phases_prep_glc_init(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       init_prep_glc = .true.
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! accumulator land input to glc on land grid
    !---------------------------------------

    ! Advance prepglc_clock - this will make the prepglc_clock in sync with the mediator clock
    call ESMF_ClockAdvance(prepglc_clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Accumulate fields
    do n = 1, size(fldnames_fr_lnd)
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldname=trim(fldnames_fr_lnd(n)), &
            field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=data2d_in, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(FBlndAccum_lnd, fieldname=fldnames_fr_lnd(n), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=data2d_out, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do i = 1,size(data2d_out, dim=2)
          data2d_out(:,i) = data2d_out(:,i) + data2d_in(:,i)
       end do
    end do

    ! Increment accumulation counter
    FBlndAccumCnt = FBlndAccumCnt + 1

    if (dbug_flag > 1) then
       call FB_diagnose(FBlndAccum_lnd, string=trim(subname)// ' FBlndAccum_lnd ',  rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_accum

  !================================================================================================
  subroutine med_phases_prep_glc_avg(gcomp, rc)

    !---------------------------------------
    ! Prepare the GLC export Fields from the mediator
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: med_clock
    type(ESMF_Time)     :: med_currtime
    type(ESMF_Time)     :: prepglc_currtime
    integer             :: yr_med, mon_med, day_med, sec_med
    integer             :: yr_prepglc, mon_prepglc, day_prepglc, sec_prepglc
    type(ESMF_Alarm)    :: alarm
    type(ESMF_Field)    :: lfield
    integer             :: i, n, ncnt            ! counters
    real(r8), pointer   :: data2d(:,:) => null()
    real(r8), pointer   :: data2d_import(:,:) => null()
    character(len=*) , parameter   :: subname=' (med_phases_prep_glc_avg) '
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    if (.not. init_prep_glc) then
       call med_phases_prep_glc_init(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       init_prep_glc = .true.
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Determine if avg alarm is ringing - and if not ringing then return
    !---------------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=med_clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet(med_clock, currtime=med_currtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(med_currtime,yy=yr_med, mm=mon_med, dd=day_med, s=sec_med, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(prepglc_clock, currtime=prepglc_currtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(prepglc_currtime,yy=yr_prepglc, mm=mon_prepglc, dd=day_prepglc, s=sec_prepglc, rc=rc)

    write(logunit,'(a,4(i8,2x))') trim(subname)//'med clock yr, mon, day, sec = ',&
         yr_med,mon_med,day_med,sec_med
    write(logunit,'(a,4(i8,2x))') trim(subname)//'prep glc clock yr, mon, day, sec = ',&
         yr_prepglc,mon_prepglc,day_prepglc,sec_prepglc

    ! If the is ringing - turn it off and continue - otherwise reset field bundle to zero and return
    call ESMF_ClockGetAlarm(prepglc_clock, alarmname='alarm_glc_avg', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       ! Turn off the alarm
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Average import from accumulated land import FB
       call ESMF_LogWrite(trim(subname)//": glc_avg alarm is ringing - averaging input from lnd to glc", ESMF_LOGMSG_INFO)
       if (mastertask) then
          write(logunit,'(a)') trim(subname)//"glc_avg alarm is ringing - averaging input from lnd to glc"
       end if
       do n = 1, size(fldnames_fr_lnd)
          call ESMF_FieldBundleGet(FBlndAccum_lnd, fieldname=fldnames_fr_lnd(n), field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayptr=data2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (FBlndAccumCnt > 0) then
             ! If accumulation count is greater than 0, do the averaging
             data2d(:,:) = data2d(:,:) / real(FBlndAccumCnt)
          else
             ! If accumulation count is 0, then simply set the averaged field bundle values from the land
             ! to the import field bundle values
             call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldname=fldnames_fr_lnd(n), &
                  field=lfield, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield, farrayptr=data2d_import, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             data2d(:,:) = data2d_import(:,:)
          end if
       end do
       if (dbug_flag > 1) then
          call FB_diagnose(FBlndAccum_lnd, string=trim(subname)//' FBlndAccum for after avg for field bundle ', rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Initialize accumulated field bundle on the glc grid to zero before doing the mapping
       call FB_reset(FBlndAccum_glc, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Map accumulated field bundle from land grid (with elevation
       ! classes) to glc grid (without elevation classes) and set FBExp(compglc) data
       call map_lnd2glc(gcomp, rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       ! zero accumulator and accumulated field bundles on both land grid
       FBlndAccumCnt = 0
       call FB_reset(FBlndAccum_lnd, value=0.0_r8, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(compglc), string=trim(subname)//' FBexp(compglc) ', rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       endif

    else

       call ESMF_LogWrite(trim(subname)//": glc_avg alarm is not ringing - returning", ESMF_LOGMSG_INFO)
       ! Reset export field bundle to zero
       ! call FB_reset(is_local%wrap%FBExp(compglc), value=0.0_r8, rc=rc)
       ! if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! turn on stop timer and return

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_avg

  !================================================================================================

  subroutine map_lnd2glc(gcomp, rc)

    !---------------------------------------
    ! map accumulated land fields from the land to the glc mesh
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)          , intent(inout) :: gcomp
    integer                      , intent(out)   :: rc

    ! local variables
    type(InternalState) :: is_local
    real(r8), pointer   :: topolnd_g_ec(:,:) => null()     ! topo in elevation classes
    real(r8), pointer   :: dataptr_g(:) => null()          ! temporary data pointer for one elevation class
    real(r8), pointer   :: topoglc_g(:) => null()          ! ice topographic height on the glc grid extracted from glc import
    real(r8), pointer   :: data_ice_covered_g(:) => null() ! data for ice-covered regions on the GLC grid
    real(r8), pointer   :: ice_covered_g(:) => null()      ! if points on the glc grid is ice-covered (1) or ice-free (0)
    integer , pointer   :: elevclass_g(:) => null()        ! elevation classes glc grid
    real(r8), pointer   :: dataexp_g(:) => null()          ! pointer into
    real(r8), pointer   :: dataptr2d(:,:) => null()
    real(r8), pointer   :: dataptr1d(:) => null()
    real(r8)            :: elev_l, elev_u                  ! lower and upper elevations in interpolation range
    real(r8)            :: d_elev                          ! elev_u - elev_l
    integer             :: nfld, ec
    integer             :: i,j,n,g,lsize_g
    integer             :: ungriddedUBound_output(1)
    integer             :: fieldCount
    type(ESMF_Field)    :: lfield
    type(ESMF_Field)    :: field_lfrac_l
    type(ESMF_Field), pointer :: fieldlist_lnd(:) => null()
    type(ESMF_Field), pointer :: fieldlist_glc(:) => null()
    character(len=*) , parameter   :: subname=' (map_lnd2glc) '
    !---------------------------------------

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Map the accumulate land field from the land grid (in multiple elevation classes)
    ! to the glc grid (in multiple elevation classes) using bilinear interpolation
    ! ------------------------------------------------------------------------

    ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
    ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
    ! current glint implementation, which sets acab and artm to 0 over ocean (although
    ! notes that this could lead to a loss of conservation). Figure out how to handle
    ! this case.

    ! get fieldlist from FBlndAccum_lnd
    call ESMF_FieldBundleGet(FBlndAccum_lnd, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist_lnd(fieldcount))
    call ESMF_FieldBundleGet(FBlndAccum_lnd, fieldlist=fieldlist_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist_glc(fieldcount))
    call ESMF_FieldBundleGet(FBlndAccum_glc, fieldlist=fieldlist_glc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get land fraction field on land mesh
    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(complnd), fieldname='lfrac', field=field_lfrac_l, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: is this needed?
    call FB_reset(FBlndAccum_glc, value=0.0_r8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! map accumlated land fields and normalize by the land fraction
    do n = 1,fieldcount
       call med_map_field_normalized(  &
         field_src=fieldlist_lnd(n), &
         field_dst=fieldlist_glc(n), &
         routehandles=is_local%wrap%RH(complnd,compglc,:), &
         maptype=mapbilnr, &
         field_normsrc=field_lfrac_l, &
         field_normdst=field_lfrac_g, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

    if (dbug_flag > 1) then
       call FB_diagnose(FBlndAccum_lnd, string=trim(subname)//' FBlndAccum_lnd ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call FB_diagnose(is_local%wrap%FBfrac(complnd), string=trim(subname)//' FBFrac ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call FB_diagnose(FBlndAccum_glc, string=trim(subname)//' FBlndAccum_glc ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc grid gridcell (elevclass_g)
    ! ------------------------------------------------------------------------

    if (dbug_flag > 1) then
       call FB_diagnose(is_local%wrap%FBImp(compglc,compglc), &
            string=trim(subname)//' FBImp(compglc,compglc) ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc,compglc), fieldname=trim(Sg_frac_fieldname), &
         field=lfield, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=ice_covered_g, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc,compglc), fieldname=trim(Sg_topo_fieldname), &
         field=lfield, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=topoglc_g, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! get elevation classes including bare land
    ! for grid cells that are ice-free, the elevation class is set to 0.
    lsize_g = size(ice_covered_g)
    allocate(elevclass_g(lsize_g))
    call glc_get_elevation_classes(ice_covered_g, topoglc_g, elevclass_g, logunit)

    ! ------------------------------------------------------------------------
    ! Determine topo field in multiple elevation classes on the glc grid
    ! ------------------------------------------------------------------------

    call ESMF_FieldBundleGet(FBlndAccum_glc, fieldname='Sl_topo_elev', field=lfield, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=topolnd_g_ec, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Loop over fields in export field bundle to GLC
    ! ------------------------------------------------------------------------

    ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
    ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
    ! current glint implementation, which sets acab and artm to 0 over ocean (although
    ! notes that this could lead to a loss of conservation). Figure out how to handle this case.

    allocate(data_ice_covered_g(lsize_g))
    do nfld = 1, size(fldnames_to_glc)

       ! ------------------------------------------------------------------------
       ! Perform vertical interpolation of data onto ice sheet topography
       ! This maps all of the input elevation classes into an export to glc without elevation classes
       ! ------------------------------------------------------------------------

       ! Get a pointer to the land data in multiple elevation classes on the glc grid
       call ESMF_FieldBundleGet(FBlndAccum_glc, fieldname=trim(fldnames_fr_lnd(nfld)), &
            field=lfield, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=dataptr2d, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       ! Get a pointer to the data for the field that will be sent to glc (without elevation classes)
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldname=trim(fldnames_to_glc(nfld)), &
            field=lfield, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=dataexp_g, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       ! First set data_ice_covered_g to bare land everywehre
       data_ice_covered_g(:) = 0._r8

       ! Now overwrite with valid values
       do n = 1, lsize_g

          ! For each ice sheet point, find bounding EC values...
          if (topoglc_g(n) < topolnd_g_ec(2,n)) then

             ! lower than lowest mean EC elevation value
             data_ice_covered_g(n) = dataptr2d(2,n)

          else if (topoglc_g(n) >= topolnd_g_ec(ungriddedCount, n)) then

             ! higher than highest mean EC elevation value
             data_ice_covered_g(n) = dataptr2d(ungriddedCount,n)

          else

             ! do linear interpolation of data in the vertical
             do ec = 3, ungriddedCount
                if (topoglc_g(n) < topolnd_g_EC(ec,n)) then
                   elev_l = topolnd_g_EC(ec-1,n)
                   elev_u = topolnd_g_EC(ec  ,n)
                   d_elev = elev_u - elev_l
                   if (d_elev <= 0) then
                      ! This shouldn't happen, but handle it in case it does. In this case,
                      ! let's arbitrarily use the mean of the two elevation classes, rather
                      ! than the weighted mean.
                      write(logunit,*) subname//' WARNING: topo diff between elevation classes <= 0'
                      write(logunit,*) 'n, ec, elev_l, elev_u = ', n, ec, elev_l, elev_u
                      write(logunit,*) 'Simply using mean of the two elevation classes,'
                      write(logunit,*) 'rather than the weighted mean.'
                      data_ice_covered_g(n) = dataptr2d(ec-1,n) * 0.5_r8 &
                                            + dataptr2d(ec  ,n) * 0.5_r8
                   else
                      data_ice_covered_g(n) =  dataptr2d(ec-1,n) * (elev_u - topoglc_g(n)) / d_elev  &
                                             + dataptr2d(ec  ,n) * (topoglc_g(n) - elev_l) / d_elev
                   end if
                   exit
                end if
             end do
          end if  ! topoglc_g(n)

          if (elevclass_g(n) /= 0) then
             ! ice-covered cells have interpolated values
             dataexp_g(n) = data_ice_covered_g(n)
          else
             ! non ice-covered cells have bare land value
             dataexp_g(n) = real(dataptr2d(1,n))
          end if
       end do  ! lsize_g

    end do  ! end of loop over fields

    ! ------------------------------------------------------------------------
    ! Renormalize surface mass balance (smb, here named dataexp_g) so that the global
    ! integral on the glc grid is equal to the global integral on the land grid.
    ! ------------------------------------------------------------------------

    if (smb_renormalize) then
       call med_phases_prep_glc_renormalize_smb(gcomp, rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! clean up memory
    deallocate(elevclass_g)
    deallocate(data_ice_covered_g)

  end subroutine map_lnd2glc

  !================================================================================================

  subroutine med_phases_prep_glc_renormalize_smb(gcomp, rc)

    !------------------
    ! Renormalizes surface mass balance (smb, here named qice_g) so that the global
    ! integral on the glc grid is equal to the global integral on the land grid.
    !
    ! (1) Map Sg_icemask from the glc grid to the land grid.
    !     Because of coupler lags, the current Sg_icemask_l might not be up to date with Sg_icemask_g.
    ! (2) Map Sg_ice_covered from the glc grid (no elevation classes) to the
    !     land grid (multiple elevation classes)
    !
    ! This is required for conservation - although conservation is only necessary if we
    ! are running with a fully-interactive, two-way-coupled glc.
    !
    ! Note that, for a case with full two-way coupling, we will only conserve if the
    ! actual land cover used over the course of the year matches these currently-remapped
    ! values. This should generally be the case with the current coupling setup.
    !
    ! One could argue that it would be safer (for conservation purposes) if LND sent its
    ! grid cell average SMB values, or if it sent its own notion of the area in each
    ! elevation class for the purpose of creating grid cell average SMB values here. But
    ! these options cause problems if we're not doing full two-way coupling (e.g., in a TG
    ! case with dlnd, or in the common case where GLC is a diagnostic component that
    ! doesn't cause updates in the glacier areas in LND). In these cases without full
    ! two-way coupling, if we use the LND's notion of the area in each elevation class,
    ! then the conservation corrections would end up correcting for discrepancies in
    ! elevation class areas between LND and GLC, rather than just correcting for
    ! discrepancies arising from the remapping of SMB. (And before you get worried: It
    ! doesn't matter that we are not conserving in these cases without full two-way
    ! coupling, because GLC isn't connected with the rest of the system in terms of energy
    ! and mass in these cases. So in these cases, it's okay that the LND integral computed
    ! here differs from the integral that LND itself would compute.)
    !
    ! Note: Sg_icemask defines where the ice sheet model can receive a
    ! nonzero SMB from the land model.
    !
    ! For high-level design, see:
    ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit
    !------------------

    ! input/output variables
    type(ESMF_GridComp)    :: gcomp
    integer  , intent(out) :: rc          ! return error code

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_VM)       :: vm
    type(ESMF_Field)    :: lfield
    real(r8) , pointer  :: qice_g(:) => null()      ! SMB (Flgl_qice) on glc grid without elev classes
    real(r8) , pointer  :: qice_l_ec(:,:) => null() ! SMB (Flgl_qice) on land grid with elev classes
    real(r8) , pointer  :: topo_g(:) => null()      ! ice topographic height on the glc grid cell
    real(r8) , pointer  :: frac_g(:) => null()      ! total ice fraction in each glc cell
    real(r8) , pointer  :: frac_g_ec(:,:) => null() ! total ice fraction in each glc cell
    real(r8) , pointer  :: frac_l_ec(:,:) => null() ! EC fractions (Sg_ice_covered) on land grid
    real(r8) , pointer  :: icemask_g(:) => null()   ! icemask on glc grid
    real(r8) , pointer  :: icemask_l(:) => null()   ! icemask on land grid
    real(r8) , pointer  :: lfrac(:) => null()       ! land fraction on land grid
    real(r8) , pointer  :: dataptr1d(:) => null()   ! temporary 1d pointer
    real(r8) , pointer  :: dataptr2d(:,:) => null() ! temporary 2d pointer
    integer             :: ec                       ! loop index over elevation classes
    integer             :: n

    ! local and global sums of accumulation and ablation; used to compute renormalization factors
    real(r8) :: local_accum_lnd(1), global_accum_lnd(1)
    real(r8) :: local_accum_glc(1), global_accum_glc(1)
    real(r8) :: local_ablat_lnd(1), global_ablat_lnd(1)
    real(r8) :: local_ablat_glc(1), global_ablat_glc(1)

    ! renormalization factors (should be close to 1, e.g. in range 0.95 to 1.05)
    real(r8) :: accum_renorm_factor ! ratio between global accumulation on the two grids
    real(r8) :: ablat_renorm_factor ! ratio between global ablation on the two grids
    real(r8) :: effective_area      ! grid cell area multiplied by min(lfrac,icemask_l).
    real(r8), pointer :: area_g(:)  ! areas on glc grid
    character(len=*), parameter  :: subname=' (renormalize_smb) '
    !---------------------------------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Map icemask_g from the glc grid to the land grid.
    !---------------------------------------

    ! determine icemask_g and set as contents of field_icemask_g
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc,compglc), fieldname=trim(Sg_icemask_fieldname), &
         field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_icemask_g, farrayptr=icemask_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    icemask_g(:) = dataptr1d(:)

    ! map ice mask from glc to lnd with no normalization
    ! BUG(wjs, 2017-05-11, #1516) I think we actually want norm = .false. here, but this needs more thought
    ! Below the implementation is without normalization - this should be checked moving forwards
    call med_map_field(  &
         field_src=field_icemask_g, &
         field_dst=field_icemask_l, &
         routehandles=is_local%wrap%RH(compglc,complnd,:), &
         maptype=mapconsd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Map frac_field on glc grid without elevation classes to frac_field on land grid with elevation classes
    ! ------------------------------------------------------------------------

    ! get topo_g(:), the topographic height of each glc gridcell
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc,compglc), fieldname=trim(Sg_topo_fieldname), &
         field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=topo_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get frac_g(:), the total ice fraction in each glc gridcell
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc,compglc), fieldname=trim(Sg_frac_fieldname), &
         field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_frac_g, farrayptr=frac_g, rc=rc) ! module field
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    frac_g(:) = dataptr1d(:)

    ! get  frac_g_ec - the glc_elevclass gives the elevation class of each
    ! glc grid cell, assuming that the grid cell is ice-covered, spans [1 -> ungriddedcount]
    call ESMF_FieldGet(field_frac_g_ec, farrayptr=frac_g_ec, rc=rc) ! module field
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call glc_get_fractional_icecov(ungriddedCount-1, topo_g, frac_g, frac_g_ec, logunit)

    ! map the fraction in each elevation class from the glc grid to
    ! the land grid and normalize by the icemask on the glc grid
    call med_map_field_normalized(  &
         field_src=field_frac_g_ec, &
         field_dst=field_frac_l_ec, &
         routehandles=is_local%wrap%RH(compglc,complnd,:), &
         maptype=mapconsd, &
         field_normsrc=field_icemask_g, &
         field_normdst=field_icemask_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Sum qice_l_ec over all elevation classes for each local land grid cell then do a global sum
    !---------------------------------------

    ! get fractional ice coverage for each elevation class on the land grid, frac_l_ec(:,:)
    call ESMF_FieldGet(field_frac_l_ec, farrayptr=frac_l_ec, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! determine fraction on land grid, lfrac(:)
    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(complnd), fieldname='lfrac', field=lfield, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=lfrac, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! get icemask_l
    call ESMF_FieldGet(field_icemask_l, farrayptr=icemask_l, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! get qice_l_ec
    call ESMF_FieldBundleGet(FBlndAccum_lnd, trim(qice_fieldname)//'_elev', field=lfield, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=qice_l_ec, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    local_accum_lnd(1) = 0.0_r8
    local_ablat_lnd(1) = 0.0_r8
    do n = 1, size(lfrac)
       ! Calculate effective area for sum -  need the mapped icemask_l
       effective_area = min(lfrac(n), icemask_l(n)) * is_local%wrap%mesh_info(complnd)%areas(n)
       if (effective_area > 0.0_r8) then
          do ec = 1, ungriddedCount
             if (qice_l_ec(ec,n) >= 0.0_r8) then
                local_accum_lnd(1) = local_accum_lnd(1) + effective_area * frac_l_ec(ec,n) * qice_l_ec(ec,n)
             else
                local_ablat_lnd(1) = local_ablat_lnd(1) + effective_area * frac_l_ec(ec,n) * qice_l_ec(ec,n)
             endif
          end do ! ec
       end if ! if landmaks > 0
    enddo  ! n

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMAllreduce(vm, senddata=local_accum_lnd, recvdata=global_accum_lnd, count=1, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMAllreduce(vm, senddata=local_ablat_lnd, recvdata=global_ablat_lnd, count=1, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       write(logunit,'(a,d13.5)') trim(subname)//'global_accum_lnd = ', global_accum_lnd
       write(logunit,'(a,d13.5)') trim(subname)//'global_ablat_lnd = ', global_ablat_lnd
    endif

    !---------------------------------------
    ! Sum qice_g over local glc grid cells.
    !---------------------------------------

    ! determine qice_g
    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldname=trim(qice_fieldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=qice_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc,compglc), fieldname='Sg_area', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=area_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    local_accum_glc(1) = 0.0_r8
    local_ablat_glc(1) = 0.0_r8
    do n = 1, size(qice_g)
       if (qice_g(n) >= 0.0_r8) then
          local_accum_glc(1) = local_accum_glc(1) + icemask_g(n) * area_g(n) * qice_g(n)
       else
          local_ablat_glc(1) = local_ablat_glc(1) + icemask_g(n) * area_g(n) * qice_g(n)
       endif
    enddo  ! n
    call ESMF_VMAllreduce(vm, senddata=local_accum_glc, recvdata=global_accum_glc, count=1, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    call ESMF_VMAllreduce(vm, senddata=local_ablat_glc, recvdata=global_ablat_glc, count=1, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (mastertask) then
       write(logunit,'(a,d13.5)') trim(subname)//'global_accum_glc = ', global_accum_glc
       write(logunit,'(a,d13.5)') trim(subname)//'global_ablat_glc = ', global_ablat_glc
    endif

    ! Renormalize
    if (global_accum_glc(1) > 0.0_r8) then
       accum_renorm_factor = global_accum_lnd(1) / global_accum_glc(1)
    else
       accum_renorm_factor = 0.0_r8
    endif
    if (global_ablat_glc(1) < 0.0_r8) then  ! negative by definition
       ablat_renorm_factor = global_ablat_lnd(1) / global_ablat_glc(1)
    else
       ablat_renorm_factor = 0.0_r8
    endif
    if (mastertask) then
       write(logunit,'(a,d13.5)') trim(subname)//'accum_renorm_factor = ', accum_renorm_factor
       write(logunit,'(a,d13.5)') trim(subname)//'ablat_renorm_factor = ', ablat_renorm_factor
    endif

    do n = 1, size(qice_g)
       if (qice_g(n) >= 0.0_r8) then
          qice_g(n) = qice_g(n) * accum_renorm_factor
       else
          qice_g(n) = qice_g(n) * ablat_renorm_factor
       endif
    enddo

    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_renormalize_smb

end module med_phases_prep_glc_mod
