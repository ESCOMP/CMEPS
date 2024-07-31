module med_phases_prep_glc_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing glc export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use NUOPC_Model           , only : NUOPC_ModelGet
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_VM, ESMF_VMGet, ESMF_VMAllReduce, ESMF_REDUCE_SUM, ESMF_REDUCE_MAX
  use ESMF                  , only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockIsCreated
  use ESMF                  , only : ESMF_ClockGetAlarm, ESMF_ClockAdvance, ESMF_ClockGet
  use ESMF                  , only : ESMF_Time, ESMF_TimeGet
  use ESMF                  , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmSet, ESMF_AlarmGet
  use ESMF                  , only : ESMF_AlarmIsRinging, ESMF_AlarmRingerOff
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleIsCreated
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF                  , only : ESMF_Mesh, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8, ESMF_KIND_R8
  use ESMF                  , only : ESMF_DYNAMICMASK, ESMF_DynamicMaskSetR8R8R8, ESMF_DYNAMICMASKELEMENTR8R8R8
  use ESMF                  , only : ESMF_FieldRegrid, ESMF_REGION_EMPTY
  use med_internalstate_mod , only : complnd, compocn,  mapbilnr, mapconsd, compname, compglc
  use med_internalstate_mod , only : InternalState, maintask, logunit, map_fracname_lnd2glc
  use med_map_mod           , only : med_map_routehandles_init, med_map_rh_is_created
  use med_map_mod           , only : med_map_field_normalized, med_map_field
  use med_constants_mod     , only : dbug_flag        => med_constants_dbug_flag
  use med_constants_mod     , only : czero            => med_constants_czero
  use med_constants_mod     , only : shr_const_pi, shr_const_spval
  use med_methods_mod       , only : fldbun_getmesh   => med_methods_FB_getmesh
  use med_methods_mod       , only : fldbun_getdata2d => med_methods_FB_getdata2d
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_methods_mod       , only : fldbun_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : fldbun_reset     => med_methods_FB_reset
  use med_methods_mod       , only : fldbun_init      => med_methods_FB_init
  use med_methods_mod       , only : FB_check_for_nans => med_methods_FB_check_for_nans
  use med_methods_mod       , only : field_getdata2d  => med_methods_Field_getdata2d
  use med_methods_mod       , only : field_getdata1d  => med_methods_Field_getdata1d
  use med_methods_mod       , only : fldchk           => med_methods_FB_FldChk
  use med_utils_mod         , only : chkerr           => med_utils_ChkErr
  use med_time_mod          , only : med_time_alarmInit
  use glc_elevclass_mod     , only : glc_get_num_elevation_classes
  use glc_elevclass_mod     , only : glc_get_elevation_classes
  use glc_elevclass_mod     , only : glc_get_fractional_icecov
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_glc_init       ! called from med.F90
  public  :: med_phases_prep_glc_accum_lnd  ! called from med_phases_post_lnd_mod.F90
  public  :: med_phases_prep_glc_accum_ocn  ! called from med_phases_post_ocn_mod.F90
  public  :: med_phases_prep_glc_avg        ! called either from med_phases_post_lnd_mod.F90 or med_phases_prep_glc
  public  :: med_phases_prep_glc            ! called from nuopc run sequence

  private :: med_phases_prep_glc_map_lnd2glc
  private :: med_phases_prep_glc_renormalize_smb

  ! -----------------
  ! lnd -> glc
  ! -----------------

  ! glc fields with multiple elevation classes: lnd->glc
  ! - fields sent from lnd->med to glc    ARE     IN multiple elevation classes
  ! - fields sent from med->glc from land ARE NOT IN multiple elevation classes
  ! Need to keep track of the lnd->med fields destined for glc in the FBlndAccum field bundle.

  ! Whether to renormalize the SMB for conservation.
  ! Should be set to true for 2-way coupled runs with evolving ice sheets.
  ! Does not need to be true for 1-way coupling.
  logical :: smb_renormalize

  type(ESMF_FieldBundle), public :: FBlndAccum2glc_l
  integer               , public :: lndAccum2glc_cnt

  character(len=14)              :: fldnames_fr_lnd(3) = (/'Flgl_qice_elev','Sl_tsrf_elev  ','Sl_topo_elev  '/)
  character(len=14)              :: fldnames_to_glc(2) = (/'Flgl_qice     ','Sl_tsrf       '/)

  type, public :: toglc_frlnd_type
     character(CS)          :: name
     type(ESMF_FieldBundle) :: FBlndAccum2glc_g
     type(ESMF_Field)       :: field_icemask_g
     type(ESMF_Field)       :: field_frac_g
     type(ESMF_Field)       :: field_frac_g_ec
     type(ESMF_Field)       :: field_lfrac_g
     type(ESMF_Mesh)        :: mesh_g
  end type toglc_frlnd_type
  type(toglc_frlnd_type), allocatable :: toglc_frlnd(:)

  type(ESMF_Field)               :: field_normdst_l
  type(ESMF_Field)               :: field_icemask_l
  type(ESMF_Field)               :: field_frac_l
  type(ESMF_Field)               :: field_frac_l_ec

  character(len=*), parameter    :: qice_fieldname       = 'Flgl_qice' ! Name of flux field giving surface mass balance
  character(len=*), parameter    :: Sg_frac_fieldname    = 'Sg_ice_covered'
  character(len=*), parameter    :: Sg_topo_fieldname    = 'Sg_topo'
  character(len=*), parameter    :: Sg_icemask_fieldname = 'Sg_icemask'
  integer                        :: ungriddedCount ! this equals the number of elevation classes + 1 (for bare land)

  ! -----------------
  ! ocn -> glc
  ! -----------------

  type(ESMF_FieldBundle), public :: FBocnAccum2glc_o
  integer               , public :: ocnAccum2glc_cnt
  character(len=14)              :: fldnames_fr_ocn(2) = (/'So_t_depth','So_s_depth'/)  ! TODO: what else needs to be added here
  type(ESMF_DynamicMask)         :: dynamicOcnMask
  integer, parameter             :: num_ocndepths = 30

  type(ESMF_Clock)        :: prepglc_clock
  character(*), parameter :: u_FILE_u  = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_prep_glc_init(gcomp, rc)

    !---------------------------------------
    ! Create lnd accumulation field bundles on lnd and and glc mesh and initialize accumulation count
    ! Create ocn accumulation field bundle on onc mesh and initialize accumulation count
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,ns,nf
    type(ESMF_Mesh)     :: mesh_l
    type(ESMF_Mesh)     :: mesh_o
    type(ESMF_Field)    :: lfield
    character(len=CS)   :: glc_renormalize_smb
    logical             :: glc_coupled_fluxes
    integer             :: ungriddedUBound_output(1) ! currently the size must equal 1 for rank 2 fieldds
    character(len=*),parameter  :: subname=' (med_phases_prep_glc_init) '
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate module variables
    allocate(toglc_frlnd(is_local%wrap%num_icesheets))

    ! -------------------------------
    ! If will accumulate lnd2glc input on land grid
    ! -------------------------------

    if (is_local%wrap%accum_lnd2glc) then
       ! Create field bundles for the fldnames_fr_lnd that have an
       ! undistributed dimension corresponding to elevation classes (including bare land)
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fldnames_fr_lnd(1), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound_output, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ungriddedCount = ungriddedUBound_output(1)
       ! TODO: check that ungriddedCount = glc_nec+1

       ! Get land mesh
       call fldbun_getmesh(is_local%wrap%FBImp(complnd,complnd), mesh_l, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       FBlndAccum2glc_l = ESMF_FieldBundleCreate(name='FBlndAccum2glc_l', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(fldnames_fr_lnd)
          lfield = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, name=fldnames_fr_lnd(n), &
               meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBlndAccum2glc_l, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite(trim(subname)//' adding field '//trim(fldnames_fr_lnd(n))//' to FBLndAccum_l', &
               ESMF_LOGMSG_INFO)
       end do
       call fldbun_reset(FBlndAccum2glc_l, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! -------------------------------
    ! If lnd->glc couplng is active
    ! -------------------------------

    if (is_local%wrap%lnd2glc_coupling) then
       ! Create accumulation field bundles from land on each glc ice sheet mesh
       ! Determine glc mesh from the mesh from the first export field to glc
       ! However FBlndAccum2glc_g has the fields fldnames_fr_lnd BUT ON the glc grid
       do ns = 1,is_local%wrap%num_icesheets
          ! get mesh on glc grid
          call fldbun_getmesh(is_local%wrap%FBExp(compglc(ns)), toglc_frlnd(ns)%mesh_g, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! create accumulation field bundle on glc grid
          toglc_frlnd(ns)%FBlndAccum2glc_g = ESMF_FieldBundleCreate(rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do nf = 1,size(fldnames_fr_lnd)
             lfield = ESMF_FieldCreate(toglc_frlnd(ns)%mesh_g, ESMF_TYPEKIND_R8, name=fldnames_fr_lnd(nf), &
                  meshloc=ESMF_MESHLOC_ELEMENT, &
                  ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldBundleAdd(toglc_frlnd(ns)%FBlndAccum2glc_g, (/lfield/), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          end do
          call fldbun_reset(toglc_frlnd(ns)%FBlndAccum2glc_g, value=0.0_r8, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! create land fraction field on glc mesh (this is just needed for normalization mapping)
          toglc_frlnd(ns)%field_lfrac_g = ESMF_FieldCreate(toglc_frlnd(ns)%mesh_g, ESMF_TYPEKIND_R8, &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! create route handle if it has not been created
          if (.not. med_map_RH_is_created(is_local%wrap%RH(complnd,compglc(ns),:),mapbilnr,rc=rc)) then
             call ESMF_LogWrite(trim(subname)//" mapbilnr is not created for lnd->glc mapping", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          end if
       end do

       ! Determine if renormalize smb
       call NUOPC_CompAttributeGet(gcomp, name='glc_renormalize_smb', value=glc_renormalize_smb, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! TODO: talk to Bill Sacks to determine if this is the correct logic
       glc_coupled_fluxes = is_local%wrap%med_coupling_active(compglc(1),complnd)
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
       if (maintask) then
          write(logunit,'(a,l4)') trim(subname)//' smb_renormalize is ',smb_renormalize
       end if

       if (smb_renormalize) then
          ! field used in the normalization for the mapping
          field_normdst_l = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! ice mask without elevation classes on lnd
          field_icemask_l = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! ice fraction without multiple elevation classes on lnd
          field_frac_l = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! ice fraction in multiple elevation classes on lnd - NOTE that this includes bare land
          field_frac_l_ec = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Loop over ice sheets
          do ns = 1,is_local%wrap%num_icesheets
             ! ice mask without elevation classes on glc
             toglc_frlnd(ns)%field_icemask_g = ESMF_FieldCreate(toglc_frlnd(ns)%mesh_g, &
                  ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! ice fraction without multiple elevation classes on glc
             toglc_frlnd(ns)%field_frac_g = ESMF_FieldCreate(toglc_frlnd(ns)%mesh_g, &
                  ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! ice fraction in multiple elevation classes on glc - NOTE that this includes bare land
             toglc_frlnd(ns)%field_frac_g_ec = ESMF_FieldCreate(toglc_frlnd(ns)%mesh_g, &
                  ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
                  ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Create route handle if it has not been created - this will be needed to map the fractions
             if (.not. med_map_RH_is_created(is_local%wrap%RH(compglc(ns),complnd,:),mapconsd, rc=rc)) then
                if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compglc(ns),complnd))) then
                 call fldbun_init(is_local%wrap%FBImp(compglc(ns),complnd), is_local%wrap%flds_scalar_name, &
                      STgeom=is_local%wrap%NStateImp(complnd), &
                      STflds=is_local%wrap%NStateImp(compglc(ns)), &
                      name='FBImp'//trim(compname(compglc(ns)))//'_'//trim(compname(complnd)), rc=rc)
                end if
                call med_map_routehandles_init( compglc(ns), complnd, &
                     FBSrc=is_local%wrap%FBImp(compglc(ns),compglc(ns)), &
                     FBDst=is_local%wrap%FBImp(compglc(ns),complnd), &
                     mapindex=mapconsd, &
                     RouteHandle=is_local%wrap%RH, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end do
       end if
    end if

    ! -------------------------------
    ! If ocn->glc coupling is active
    ! -------------------------------

    if (is_local%wrap%ocn2glc_coupling) then
       ! Get ocean mesh
       call fldbun_getmesh(is_local%wrap%FBImp(compocn,compocn), mesh_o, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       FBocnAccum2glc_o = ESMF_FieldBundleCreate(name='FBocnAccum2glc_o', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(fldnames_fr_ocn)
          lfield = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, name=fldnames_fr_ocn(n), &
               meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/num_ocndepths/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBocnAccum2glc_o, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite(trim(subname)//' adding field '//trim(fldnames_fr_ocn(n))//' to FBOcnAccum2glc_o', &
               ESMF_LOGMSG_INFO)
       end do
       call fldbun_reset(FBocnAccum2glc_o, value=czero, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! create route handle if it has not been created
       do ns = 1,is_local%wrap%num_icesheets
          if (.not. med_map_RH_is_created(is_local%wrap%RH(compocn,compglc(ns),:),mapbilnr,rc=rc)) then
             call ESMF_LogWrite(trim(subname)//" mapbilnr is not created for ocn->glc mapping", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          end if
       end do

       ! Create a dynamic mask object
       ! The dynamic mask object further holds a pointer to the routine that will be called in order to
       ! handle dynamically masked elements - in this case its DynOcnMaskProc (see below)
       call ESMF_DynamicMaskSetR8R8R8(dynamicOcnMask, dynamicMaskRoutine=DynOcnMaskProc, &
            dynamicSrcMaskValue=1.e30_r8,  handleAllElements=.true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_init

  !================================================================================================
  subroutine med_phases_prep_glc_accum_lnd(gcomp, rc)

    !---------------------------------------
    ! Carry out accumulation for the lnd->glc and ocn->glc
    ! Accumulation and averaging is done on
    !  - on the land mesh for land input
    !  - on the ocean mesh for ocean input
    ! Mapping from the land to the glc grid and from the ocean to the glc grid
    ! is then done after the accumulated fields have been time averaged
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: i,n
    real(r8), pointer   :: data2d_in(:,:)
    real(r8), pointer   :: data2d_out(:,:)
    character(len=*),parameter :: subname=' (med_phases_prep_glc_accum) '
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! Accumulate fields from land on land mesh that will be sent to glc
    do n = 1, size(fldnames_fr_lnd)
       call fldbun_getdata2d(is_local%wrap%FBImp(complnd,complnd), fldnames_fr_lnd(n), data2d_in, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata2d(FBlndAccum2glc_l, fldnames_fr_lnd(n), data2d_out, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do i = 1,size(data2d_out, dim=2)
          data2d_out(:,i) = data2d_out(:,i) + data2d_in(:,i)
       end do
    end do
    lndAccum2glc_cnt = lndAccum2glc_cnt + 1
    if (dbug_flag > 1) then
       call fldbun_diagnose(FBlndAccum2glc_l, string=trim(subname)// ' FBlndAccum2glc_l ',  rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_accum_lnd

  !================================================================================================
  subroutine med_phases_prep_glc_accum_ocn(gcomp, rc)

    !---------------------------------------
    ! Carry out accumulation for ocn->glc
    ! Accumulation and averaging is done on
    !  - on the ocean mesh for ocean input
    ! Mapping from from the ocean to the glc grid is then done after
    ! the accumulated fields have been time averaged
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: i,n
    real(r8), pointer   :: data2d_in(:,:)
    real(r8), pointer   :: data2d_out(:,:)
    character(len=*),parameter  :: subname=' (med_phases_prep_glc_accum) '
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! Accumulate fields from ocean on ocean mesh that will be sent to glc
    do n = 1, size(fldnames_fr_ocn)
       call fldbun_getdata2d(is_local%wrap%FBImp(compocn,compocn), fldnames_fr_ocn(n), data2d_in, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata2d(FBocnAccum2glc_o, fldnames_fr_ocn(n), data2d_out, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do i = 1,size(data2d_out, dim=2)
          data2d_out(:,i) = data2d_out(:,i) + data2d_in(:,i)
       end do
    end do
    ocnAccum2glc_cnt = ocnAccum2glc_cnt + 1
    if (dbug_flag > 1) then
       call fldbun_diagnose(FBocnAccum2glc_o, string=trim(subname)// ' FBocnAccum2glc_o ',  rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_accum_ocn

  !================================================================================================
  subroutine med_phases_prep_glc_avg(gcomp, rc)

    !---------------------------------------
    ! Create module clock (prepglc_clock)
    ! Prepare the GLC export Fields from the mediator
    !---------------------------------------

    use med_phases_history_mod, only :  med_phases_history_write_lnd2glc

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: lfield_src
    type(ESMF_Field)    :: lfield_dst
    type(ESMF_Clock)    :: med_clock
    type(ESMF_Time)     :: med_currtime
    type(ESMF_Time)     :: prepglc_currtime
    type(ESMF_ALARM)    :: glc_avg_alarm
    character(len=CS)   :: glc_avg_period
    integer             :: glc_cpl_dt
    integer             :: yr_med, mon_med, day_med, sec_med
    integer             :: yr_prepglc, mon_prepglc, day_prepglc, sec_prepglc
    type(ESMF_Alarm)    :: alarm
    integer             :: n, ns
    real(r8), pointer   :: data2d(:,:)
    real(r8), pointer   :: data2d_import(:,:)
    character(len=CS)   :: cvalue
    logical             :: do_avg
    logical             :: isPresent, isSet
    logical             :: write_histaux_l2x1yrg
    character(len=*) , parameter   :: subname=' (med_phases_prep_glc) '

    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if
    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    if (.not. ESMF_ClockIsCreated(prepglc_clock)) then
       ! Initialize prepglc_clock from mclock - THIS CALL DOES NOT COPY ALARMS
       call NUOPC_ModelGet(gcomp, modelClock=med_clock,  rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       prepglc_clock = ESMF_ClockCreate(med_clock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Set alarm glc averaging interval
       call NUOPC_CompAttributeGet(gcomp, name="glc_avg_period", value=glc_avg_period, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (trim(glc_avg_period) == 'yearly') then
          call med_time_alarmInit(prepglc_clock, glc_avg_alarm, 'yearly', alarmname='alarm_glc_avg', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (maintask) then
             write(logunit,'(a,i10)') trim(subname)//&
                  ' created alarm with averaging period for export to glc is yearly'
          end if
       else if (trim(glc_avg_period) == 'glc_coupling_period') then
          call NUOPC_CompAttributeGet(gcomp, name="glc_cpl_dt", value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) glc_cpl_dt
          call med_time_alarmInit(prepglc_clock, glc_avg_alarm, 'nseconds', opt_n=glc_cpl_dt, alarmname='alarm_glc_avg', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (maintask) then
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
    end if

    ! Advance prepglc_clock - this will make the prepglc_clock in sync with the mediator clock
    call ESMF_ClockAdvance(prepglc_clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Check time
    if (dbug_flag > 5) then
       if (maintask) then
          call NUOPC_ModelGet(gcomp, modelClock=med_clock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockGet(med_clock, currtime=med_currtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(med_currtime,yy=yr_med, mm=mon_med, dd=day_med, s=sec_med, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockGet(prepglc_clock, currtime=prepglc_currtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(prepglc_currtime,yy=yr_prepglc, mm=mon_prepglc, dd=day_prepglc, s=sec_prepglc, rc=rc)
          if (maintask) then
             write(logunit,'(a,4(i8,2x))') trim(subname)//'med clock yr, mon, day, sec      = ',&
                  yr_med,mon_med,day_med,sec_med
             write(logunit,'(a,4(i8,2x))') trim(subname)//'prep glc clock yr, mon, day, sec = ',&
                  yr_prepglc,mon_prepglc,day_prepglc,sec_prepglc
          end if
       end if
    end if

    ! Determine if the alarm is ringing
    call ESMF_ClockGetAlarm(prepglc_clock, alarmname='alarm_glc_avg', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       do_avg = .true.
       call ESMF_LogWrite(trim(subname)//": glc_avg alarm is ringing - average input from lnd and ocn to glc", &
            ESMF_LOGMSG_INFO)
       if (maintask) then
          write(logunit,'(a)') trim(subname)//"glc_avg alarm is ringing - averaging input from lnd and ocn to glc"
       end if
       ! Turn off the alarm
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       do_avg = .false.
       call ESMF_LogWrite(trim(subname)//": glc_avg alarm is not ringing - returning", ESMF_LOGMSG_INFO)
    end if

    ! Average and map data from land (and possibly ocean)
    if (do_avg) then
       ! Always average import from accumulated land import data
       do n = 1, size(fldnames_fr_lnd)
          if (fldchk(FBlndAccum2glc_l, fldnames_fr_lnd(n), rc=rc)) then
             call fldbun_getdata2d(FBlndAccum2glc_l, fldnames_fr_lnd(n), data2d, rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (lndAccum2glc_cnt > 0) then
                ! If accumulation count is greater than 0, do the averaging
                data2d(:,:) = data2d(:,:) / real(lndAccum2glc_cnt)
             else
                ! If accumulation count is 0, then simply set the averaged field bundle values from the land
                ! to the import field bundle values
                call fldbun_getdata2d(is_local%wrap%FBImp(complnd,complnd), fldnames_fr_lnd(n), data2d_import, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                data2d(:,:) = data2d_import(:,:)
             end if
          end if
       end do

       if (is_local%wrap%ocn2glc_coupling) then
          ! Average import from accumulated ocn import data
          do n = 1, size(fldnames_fr_ocn)
             call fldbun_getdata2d(FBocnAccum2glc_o, fldnames_fr_ocn(n), data2d, rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (ocnAccum2glc_cnt > 0) then
                ! If accumulation count is greater than 0, do the averaging
                data2d(:,:) = data2d(:,:) / real(ocnAccum2glc_cnt)
             else
                ! If accumulation count is 0, then simply set the averaged field bundle values from the ocn
                ! to the import field bundle values
                call fldbun_getdata2d(is_local%wrap%FBImp(compocn,compocn), fldnames_fr_ocn(n), data2d_import, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                data2d(:,:) = data2d_import(:,:)
             end if
          end do
          if (dbug_flag > 1) then
             call fldbun_diagnose(FBocnAccum2glc_o, string=trim(subname)//' FBocnAccum for after avg for field bundle ', rc=rc)
             if (chkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Map accumulated ocean field from ocean mesh to land mesh and set FBExp(compglc(ns)) data
          ! Zero land accumulator and accumulated field bundles on ocean grid
          do n = 1,size(fldnames_fr_ocn)
             call ESMF_FieldBundleGet(FBocnAccum2glc_o, fldnames_fr_ocn(n), field=lfield_src, rc=rc)
             if (chkErr(rc,__LINE__,u_FILE_u)) return
             do ns = 1,is_local%wrap%num_icesheets
                call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc(ns)), fldnames_fr_ocn(n), field=lfield_dst, rc=rc)
                if (chkErr(rc,__LINE__,u_FILE_u)) return
                ! Do mapping of ocn to glc with dynamic masking
                call ESMF_FieldRegrid(lfield_src, lfield_dst, &
                     routehandle=is_local%wrap%RH(compocn,compglc(ns),mapbilnr), dynamicMask=dynamicOcnMask, &
                     zeroregion=ESMF_REGION_EMPTY, rc=rc)
                if (chkErr(rc,__LINE__,u_FILE_u)) return
                call fldbun_getdata2d(is_local%wrap%FBExp(compglc(ns)), fldnames_fr_ocn(n), data2d, rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                ! reset values of 0 to spval
                where (data2d == 0._r8) data2d = shr_const_spval
             end do
          end do
          ocnAccum2glc_cnt = 0
          call fldbun_reset(FBocnAccum2glc_o, value=czero, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Determine if auxiliary file will be written
       write_histaux_l2x1yrg = .false.
       if (lndAccum2glc_cnt > 0) then
          call NUOPC_CompAttributeGet(gcomp, name="histaux_l2x1yrg", value=cvalue, &
               isPresent=isPresent, isSet=isSet, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (isPresent .and. isSet) then
             read(cvalue,*) write_histaux_l2x1yrg
          end if
       end if

       ! Write auxiliary history file if flag is set and accumulation is being done
       if (is_local%wrap%lnd2glc_coupling) then
          ! Map accumulated field bundle from land grid (with elevation classes) to glc grid (without elevation classes)
          ! and set FBExp(compglc(ns)) data
          ! Zero land accumulator and accumulated field bundles on land grid
          call med_phases_prep_glc_map_lnd2glc(gcomp, rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return

          if (write_histaux_l2x1yrg) then
             call med_phases_history_write_lnd2glc(gcomp, FBlndAccum2glc_l, &
                  fldbun_glc=is_local%wrap%FBExp(compglc(:)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          lndAccum2glc_cnt = 0
          call fldbun_reset(FBlndAccum2glc_l, value=czero, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       else
          if (write_histaux_l2x1yrg) then
             call med_phases_history_write_lnd2glc(gcomp, FBlndAccum2glc_l, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if

       if (dbug_flag > 1) then
          do ns = 1,is_local%wrap%num_icesheets
             call fldbun_diagnose(is_local%wrap%FBExp(compglc(ns)), string=trim(subname)//' FBexp(compglc) ', rc=rc)
             if (chkErr(rc,__LINE__,u_FILE_u)) return
          end do
       endif
    end if

    ! Check for nans in fields export to glc
    do ns = 1,is_local%wrap%num_icesheets
       call FB_check_for_nans(is_local%wrap%FBExp(compglc(ns)), maintask, logunit, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_avg

  !================================================================================================
  subroutine med_phases_prep_glc(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    call med_phases_prep_glc_avg(gcomp, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_prep_glc

  !================================================================================================
  subroutine med_phases_prep_glc_map_lnd2glc(gcomp, rc)

    !---------------------------------------
    ! map accumulated land fields from the land to the glc mesh
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)          , intent(inout) :: gcomp
    integer                      , intent(out)   :: rc

    ! local variables
    type(InternalState) :: is_local
    real(r8), pointer   :: topolnd_g_ec(:,:)      ! topo in elevation classes
    real(r8), pointer   :: topoglc_g(:)           ! ice topographic height on the glc grid extracted from glc import
    real(r8), pointer   :: data_ice_covered_g(:)  ! data for ice-covered regions on the GLC grid
    real(r8), pointer   :: ice_covered_g(:)       ! if points on the glc grid is ice-covered (1) or ice-free (0)
    integer , pointer   :: elevclass_g(:)         ! elevation classes glc grid
    real(r8), pointer   :: dataexp_g(:)           ! pointer into
    real(r8), pointer   :: dataptr2d(:,:)
    real(r8)            :: elev_l, elev_u         ! lower and upper elevations in interpolation range
    real(r8)            :: d_elev                 ! elev_u - elev_l
    integer             :: nfld, ec
    integer             :: n,lsize_g,ns
    type(ESMF_Field)    :: field_lfrac_l
    integer             :: fieldCount
    character(len=3)    :: cnum
    type(ESMF_Field), pointer :: fieldlist_lnd(:)
    type(ESMF_Field), pointer :: fieldlist_glc(:)
    character(len=*) , parameter   :: subname=' (med_phases_prep_glc_map_lnd2glc) '
    !---------------------------------------

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Map the accumulate land field from the land grid (in multiple elevation classes)
    ! to the glc grid (in multiple elevation classes) using bilinear interpolation
    ! ------------------------------------------------------------------------

    ! Initialize accumulated field bundle on the glc grid to zero before doing the mapping
    do ns = 1,is_local%wrap%num_icesheets
       call fldbun_reset(toglc_frlnd(ns)%FBlndAccum2glc_g, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

    ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
    ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
    ! current glint implementation, which sets acab and artm to 0 over ocean (although
    ! notes that this could lead to a loss of conservation). Figure out how to handle
    ! this case.

    ! get fieldlist from FBlndAccum2glc_l
    call ESMF_FieldBundleGet(FBlndAccum2glc_l, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist_lnd(fieldcount))
    allocate(fieldlist_glc(fieldcount))
    call ESMF_FieldBundleGet(FBlndAccum2glc_l, fieldlist=fieldlist_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get land fraction field on land mesh
    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(complnd), fieldName=map_fracname_lnd2glc, field=field_lfrac_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! map accumlated land fields to each ice sheet (normalize by the land fraction in the mapping)
    do ns = 1,is_local%wrap%num_icesheets
       call fldbun_reset(toglc_frlnd(ns)%FBlndAccum2glc_g, value=0.0_r8, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    end do
    do ns = 1,is_local%wrap%num_icesheets
       call ESMF_FieldBundleGet(toglc_frlnd(ns)%FBlndAccum2glc_g, fieldlist=fieldlist_glc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do nfld = 1,fieldcount
          call med_map_field_normalized(  &
               field_src=fieldlist_lnd(nfld), &
               field_dst=fieldlist_glc(nfld), &
               routehandles=is_local%wrap%RH(complnd,compglc(ns),:), &
               maptype=mapbilnr, &
               field_normsrc=field_lfrac_l, &
               field_normdst=toglc_frlnd(ns)%field_lfrac_g, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do
    end do

    deallocate(fieldlist_lnd)
    deallocate(fieldlist_glc)

    if (dbug_flag > 1) then
       call fldbun_diagnose(FBlndAccum2glc_l, string=trim(subname)//' FBlndAccum2glc_l ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call fldbun_diagnose(is_local%wrap%FBfrac(complnd), string=trim(subname)//' FBFrac ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       do ns = 1,is_local%wrap%num_icesheets
          call fldbun_diagnose(toglc_frlnd(ns)%FBlndAccum2glc_g, string=trim(subname)//&
               ' FBlndAccum2glc_glc '//compname(compglc(ns)), rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end do
    endif

    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc grid gridcell (elevclass_g)
    ! ------------------------------------------------------------------------

    ! Loop over ice sheets
    do ns = 1,is_local%wrap%num_icesheets
       if (dbug_flag > 1) then
          write(cnum,'(a3)') ns
          call fldbun_diagnose(is_local%wrap%FBImp(compglc(ns),compglc(ns)), &
               string=trim(subname)//' FBImp(compglc,compglc) '//' for ice sheet '//trim(cnum), rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

       call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_frac_fieldname, ice_covered_g, rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_topo_fieldname, topoglc_g, rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       ! get elevation classes including bare land
       ! for grid cells that are ice-free, the elevation class is set to 0.
       lsize_g = size(ice_covered_g)
       allocate(elevclass_g(lsize_g))
       call glc_get_elevation_classes(ice_covered_g, topoglc_g, elevclass_g, logunit)

       ! Determine topo field in multiple elevation classes on the glc grid
       call fldbun_getdata2d(toglc_frlnd(ns)%FBlndAccum2glc_g, 'Sl_topo_elev', topolnd_g_ec, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       ! ------------------------------------------------------------------------
       ! Loop over fields in export field bundle to glc for ice sheet ns and
       ! perform vertical interpolation of data onto ice sheet topography
       ! This maps all of the input elevation classes into an export to glc without elevation classes
       ! ------------------------------------------------------------------------

       ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
       ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
       ! current glint implementation, which sets acab and artm to 0 over ocean (although
       ! notes that this could lead to a loss of conservation). Figure out how to handle this case.

       allocate(data_ice_covered_g(lsize_g))
       do nfld = 1, size(fldnames_to_glc)

          ! Get a pointer to the land data in multiple elevation classes on the glc grid
          call fldbun_getdata2d(toglc_frlnd(ns)%FBlndAccum2glc_g, fldnames_fr_lnd(nfld), dataptr2d, rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return

          ! Get a pointer to the data for the field that will be sent to glc (without elevation classes)
          call fldbun_getdata1d(is_local%wrap%FBExp(compglc(ns)), fldnames_to_glc(nfld), dataexp_g, rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return

          ! First set data_ice_covered_g to bare land everywehre
          data_ice_covered_g(:) = 0._r8

          ! Loop over land points and overwrite with valid values
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

          end do  ! end of loop over land points

       end do ! end loop over fields (nflds)

       ! ------------------------------------------------------------------------
       ! Renormalize surface mass balance (smb, here named dataexp_g) so that the global
       ! integral on the glc grid is equal to the global integral on the land grid.
       ! ------------------------------------------------------------------------

       ! No longer need to make a preemptive adjustment to qice_g to account for area differences
       ! between CISM and the coupler. In NUOPC, the area correction is done in! the cap not in the
       ! mediator, so to preserve the bilinear mapping values, do not need to do any area correction
       ! scaling in the CISM NUOPC cap

       if (smb_renormalize) then
          call med_phases_prep_glc_renormalize_smb(gcomp, ns, rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! clean up memory that is ice sheet dependent
       deallocate(elevclass_g)
       deallocate(data_ice_covered_g)

    end do ! end of loop over ice sheets

  end subroutine med_phases_prep_glc_map_lnd2glc

  !================================================================================================
  subroutine med_phases_prep_glc_renormalize_smb(gcomp, ns, rc)

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
    ! A note on the areas used in the global sums, with this CMEPS implementation: For
    ! glc, we use the internal model areas (sent from glc to the mediator), whereas for
    ! lnd, we use the mesh areas (which are the same as the areas used in mapping). The
    ! reason for this difference is that for lnd (as for most components), we plan to do
    ! area correction (correcting for the discrepancy between internal areas and mesh
    ! areas) in the cap, so that the fluxes received by the mediator are expressed
    ! per-unit-area according to the mesh areas. However, we are currently *not* planning
    ! to do this area correction for glc, because there are some fields that shouldn't be
    ! area corrected. Thus, for fluxes sent to glc, they are specified per-unit-area
    ! according to the internal areas, so we need to use glc's internal areas in the
    ! global sums. (If we did the area correction for qice sent to glc, we would probably
    ! want to make a preemptive adjustment to qice, multiplying by the inverse of the area
    ! corrections, as we did in MCT / CPL7. This ends up giving the same answer at the
    ! cost of additional complexity, and still requires sending glc's internal areas to
    ! the mediator so that it can do this preemptive adjustment.)
    !
    ! Note: Sg_icemask defines where the ice sheet model can receive a
    ! nonzero SMB from the land model.
    !
    ! For high-level design, see:
    ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit
    !------------------

    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    integer , intent(in)  :: ns          ! icesheet instance index
    integer , intent(out) :: rc          ! return error code

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_VM)       :: vm
    real(r8) , pointer  :: qice_g(:)       ! SMB (Flgl_qice) on glc grid without elev classes
    real(r8) , pointer  :: qice_l_ec(:,:)  ! SMB (Flgl_qice) on land grid with elev classes
    real(r8) , pointer  :: topo_g(:)       ! ice topographic height on the glc grid cell
    real(r8) , pointer  :: frac_g(:)       ! total ice fraction in each glc cell
    real(r8) , pointer  :: frac_g_ec(:,:)  ! total ice fraction in each glc cell
    real(r8) , pointer  :: frac_l_ec(:,:)  ! EC fractions (Sg_ice_covered) on land grid
    real(r8) , pointer  :: icemask_g(:)    ! icemask on glc grid
    real(r8) , pointer  :: icemask_l(:)    ! icemask on land grid
    real(r8) , pointer  :: lndfrac(:)      ! land fraction on land grid
    real(r8) , pointer  :: dataptr1d(:)    ! temporary 1d pointer
    integer             :: ec              ! loop index over elevation classes
    integer             :: n

    ! local and global sums of accumulation and ablation; used to compute renormalization factors
    real(r8) :: local_accum_lnd(1), global_accum_lnd(1)
    real(r8) :: local_accum_glc(1), global_accum_glc(1)
    real(r8) :: local_ablat_lnd(1), global_ablat_lnd(1)
    real(r8) :: local_ablat_glc(1), global_ablat_glc(1)

    ! renormalization factors (should be close to 1, e.g. in range 0.95 to 1.05)
    real(r8) :: accum_renorm_factor ! ratio between global accumulation on the two grids
    real(r8) :: ablat_renorm_factor ! ratio between global ablation on the two grids
    real(r8) :: effective_area      ! grid cell area multiplied by min(lndfrac,icemask_l).
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
    call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_icemask_fieldname, dataptr1d, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call field_getdata1d(toglc_frlnd(ns)%field_icemask_g, icemask_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    icemask_g(:) = dataptr1d(:)

    ! map ice mask from glc to lnd with no normalization
    call med_map_field(  &
         field_src=toglc_frlnd(ns)%field_icemask_g, &
         field_dst=field_icemask_l, &
         routehandles=is_local%wrap%RH(compglc(ns),complnd,:), &
         maptype=mapconsd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get icemask_l
    call field_getdata1d(field_icemask_l, icemask_l, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Map frac_field on glc grid without elevation classes to frac_field on land grid with elevation classes
    ! ------------------------------------------------------------------------

    ! get topo_g(:), the topographic height of each glc gridcell
    call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_topo_fieldname, topo_g, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get frac_g(:), the total ice fraction in each glc gridcell
    call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_frac_fieldname, frac_g, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get  frac_g_ec - the glc_elevclass gives the elevation class of each
    ! glc grid cell, assuming that the grid cell is ice-covered, spans [1 -> ungriddedcount]
    call field_getdata2d(toglc_frlnd(ns)%field_frac_g_ec, frac_g_ec, rc=rc) ! module field
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call glc_get_fractional_icecov(ungriddedCount-1, topo_g, frac_g, frac_g_ec, logunit)

    ! map fraction in each elevation class from the glc grid to the land grid and normalize by the icemask on the
    ! glc grid
    call med_map_field_normalized(  &
         field_src=toglc_frlnd(ns)%field_frac_g_ec, &
         field_dst=field_frac_l_ec, &
         routehandles=is_local%wrap%RH(compglc(ns),complnd,:), &
         maptype=mapconsd, &
         field_normsrc=toglc_frlnd(ns)%field_icemask_g, &
         field_normdst=field_normdst_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Sum qice_l_ec over all elevation classes for each local land grid cell then do a global sum
    !---------------------------------------

    ! get fractional ice coverage for each elevation class on the land grid, frac_l_ec(:,:)
    call field_getdata2d(field_frac_l_ec, frac_l_ec, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! determine fraction on land grid, lndfrac(:)
    call fldbun_getdata1d(is_local%wrap%FBFrac(complnd), map_fracname_lnd2glc, lndfrac, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! get qice_l_ec
    call fldbun_getdata2d(FBlndAccum2glc_l, trim(qice_fieldname)//'_elev', qice_l_ec, rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    local_accum_lnd(1) = 0.0_r8
    local_ablat_lnd(1) = 0.0_r8
    do n = 1, size(lndfrac)
       ! Calculate effective area for sum -  need the mapped icemask_l
       effective_area = min(lndfrac(n), icemask_l(n)) * is_local%wrap%mesh_info(complnd)%areas(n)
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
    if (maintask) then
       write(logunit,'(a,d21.10)') trim(subname)//'global_accum_lnd = ', global_accum_lnd
       write(logunit,'(a,d21.10)') trim(subname)//'global_ablat_lnd = ', global_ablat_lnd
    endif

    !---------------------------------------
    ! Sum qice_g over local glc grid cells.
    !---------------------------------------

    ! determine qice_g
    call fldbun_getdata1d(is_local%wrap%FBExp(compglc(ns)), qice_fieldname, qice_g, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get areas internal to glc grid
    call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), 'Sg_area', area_g, rc)
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
    if (maintask) then
       write(logunit,'(a,d21.10)') trim(subname)//'global_accum_glc = ', global_accum_glc
       write(logunit,'(a,d21.10)') trim(subname)//'global_ablat_glc = ', global_ablat_glc
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
    if (maintask) then
       write(logunit,'(a,d21.10)') trim(subname)//'accum_renorm_factor = ', accum_renorm_factor
       write(logunit,'(a,d21.10)') trim(subname)//'ablat_renorm_factor = ', ablat_renorm_factor
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

  !================================================================================================
  subroutine dynOcnMaskProc(dynamicMaskList, dynamicSrcMaskValue, dynamicDstMaskValue, rc)

    use ESMF, only : ESMF_RC_ARG_BAD

    ! input/output arguments
    type(ESMF_DynamicMaskElementR8R8R8) , pointer :: dynamicMaskList(:)
    real(ESMF_KIND_R8), intent(in), optional  :: dynamicSrcMaskValue
    real(ESMF_KIND_R8), intent(in), optional  :: dynamicDstMaskValue
    integer           , intent(out)           :: rc

    ! local variables
    integer  :: no, ni
    real(ESMF_KIND_R8)  :: renorm
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Below - ONLY if you do NOT have the source masked out then do
    ! the regridding (which is done explicitly here)

    if (associated(dynamicMaskList)) then
       do no = 1, size(dynamicMaskList)
          dynamicMaskList(no)%dstElement = czero ! set to zero
          renorm = 0.d0 ! reset
          do ni = 1, size(dynamicMaskList(no)%factor)
             ! Need to multiply by .90 to handle averaging of input fields before remapping is called
             if ( dynamicMaskList(no)%srcElement(ni) > 0.d0 .and. &
                  dynamicMaskList(no)%srcElement(ni) < dynamicSrcMaskValue*.90) then
                dynamicMaskList(no)%dstElement = dynamicMaskList(no)%dstElement + &
                     (dynamicMaskList(no)%factor(ni) * dynamicMaskList(no)%srcElement(ni))
                renorm = renorm + dynamicMaskList(no)%factor(ni)
             endif
          enddo
          if (renorm > 0.d0) then
             dynamicMaskList(no)%dstElement = dynamicMaskList(no)%dstElement / renorm
          else if (present(dynamicSrcMaskValue)) then
             dynamicMaskList(no)%dstElement = dynamicSrcMaskValue
          else
             rc = ESMF_RC_ARG_BAD  ! error detected
             return
          endif
       enddo
    endif

  end subroutine DynOcnMaskProc

end module med_phases_prep_glc_mod
