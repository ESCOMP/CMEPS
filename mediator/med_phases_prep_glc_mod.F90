module med_phases_prep_glc_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing glc export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_VM, ESMF_VMGet, ESMF_VMAllReduce, ESMF_REDUCE_SUM
  use ESMF                  , only : ESMF_FieldBundle
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_RouteHandleIsCreated
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleIsCreated
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF                  , only : ESMF_Array, ESMF_ArrayGet, ESMF_ArrayCreate, ESMF_ArrayDestroy
  use ESMF                  , only : ESMF_DistGrid
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet, ESMF_MESHLOC_ELEMENT
  use ESMF                  , only : ESMF_TYPEKIND_R8
  use esmFlds               , only : compglc, complnd, mapbilnr, mapconsf, compname
  use esmFlds               , only : med_fldlist_type
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_map_mod           , only : med_map_FB_Regrid_Norm
  use med_map_mod           , only : med_map_Fractions_Init
  use med_methods_mod       , only : FB_Init      => med_methods_FB_init
  use med_methods_mod       , only : FB_getFldPtr => med_methods_FB_getFldPtr
  use med_methods_mod       , only : FB_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_reset     => med_methods_FB_reset
  use med_utils_mod         , only : chkerr       => med_utils_ChkErr
  use glc_elevclass_mod     , only : glc_get_num_elevation_classes
  use glc_elevclass_mod     , only : glc_get_elevation_classes
  use glc_elevclass_mod     , only : glc_get_fractional_icecov
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_glc_init
  public  :: med_phases_prep_glc_accum
  public  :: med_phases_prep_glc_avg

  private :: med_phases_prep_glc_map_lnd2glc
  private :: med_phases_prep_glc_renormalize_smb

  ! glc fields with multiple elevation classes: lnd->glc
  ! - fields sent from lnd->med to glc    ARE     IN multiple elevation classes
  ! - fields sent from med->glc from land ARE NOT IN multiple elevation classes
  ! Need to keep track of the lnd->med fields destined for glc in the FBlndAccum field bundle.

  ! Needed for standard lnd->glc mapping
  type(ESMF_FieldBundle) :: FBlndAccum_lnd
  type(ESMF_FieldBundle) :: FBlndAccum_glc
  integer                :: FBlndAccumCnt
  type(med_fldlist_type) :: fldlist_lnd2glc
  character(len=14)      :: fldnames_fr_lnd(3) = (/'Flgl_qice_elev','Sl_tsrf_elev  ','Sl_topo_elev  '/)
  character(len=14)      :: fldnames_to_glc(2) = (/'Flgl_qice     ','Sl_tsrf       '/)

  ! Whether to renormalize the SMB for conservation.
  ! Should be set to true for 2-way coupled runs with evolving ice sheets.
  ! Does not need to be true for 1-way coupling.
  logical :: smb_renormalize

  ! Needed if renormalize SMB
  type(ESMF_FieldBundle)      :: FBglc_icemask
  type(ESMF_FieldBundle)      :: FBlnd_icemask
  type(med_fldlist_type)      :: fldlist_glc2lnd_icemask
  type(ESMF_FieldBundle)      :: FBglc_frac
  type(ESMF_FieldBundle)      :: FBlnd_frac
  type(med_fldlist_type)      :: fldlist_glc2lnd_frac
  real(r8) , pointer          :: aream_l(:)         ! cell areas on land grid, for mapping
  real(r8) , pointer          :: aream_g(:)         ! cell areas on glc grid, for mapping
  character(len=*), parameter :: qice_fieldname   = 'Flgl_qice' ! Name of flux field giving surface mass balance
  character(len=*), parameter :: Sg_frac_field    = 'Sg_ice_covered'
  character(len=*), parameter :: Sg_topo_field    = 'Sg_topo'
  character(len=*), parameter :: Sg_icemask_field = 'Sg_icemask'

  ! Size of undistributed dimension from land
  integer :: ungriddedCount ! this equals the number of elevation classes + 1 (for bare land)

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
    type(InternalState) :: is_local
    integer             :: i,n,ncnt
    type(ESMF_Mesh)     :: lmesh_glc
    type(ESMF_Mesh)     :: lmesh_lnd
    type(ESMF_Field)    :: lfield
    real(r8), pointer   :: data2d_in(:,:)
    real(r8), pointer   :: data2d_out(:,:)
    real(r8), pointer   :: dataptr1d(:)
    character(len=CS)   :: glc_renormalize_smb
    logical             :: glc_coupled_fluxes
    type(ESMF_Array)    :: larray
    type(ESMF_DistGrid) :: ldistgrid
    integer             :: lsize
    logical             :: isPresent
    integer             :: ungriddedUBound_output(1) ! currently the size must equal 1 for rank 2 fieldds
    character(len=*),parameter  :: subname='(med_phases_prep_glc_init)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

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

       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), fldnames_fr_lnd(1), fldptr2=data2d_in, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, mesh=lmesh_lnd, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       FBlndAccum_lnd = ESMF_FieldBundleCreate(name='FBlndAccum_lnd', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(fldnames_fr_lnd)
          lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name=fldnames_fr_lnd(n), &
               meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBlndAccum_lnd, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do
       call FB_reset(FBlndAccum_lnd, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Create accumulation field bundle from land on the glc grid
       ! Determine glc mesh from the mesh from the first export field to glc
       ! However FBlndAccum_glc has the fields fldnames_fr_lnd BUT ON the glc grid

       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fldnames_to_glc(1), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, mesh=lmesh_glc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

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

       FBlndAccumCnt = 0

       allocate(fldlist_lnd2glc%flds(3))
       fldlist_lnd2glc%flds(1)%shortname = trim(fldnames_fr_lnd(1))
       fldlist_lnd2glc%flds(2)%shortname = trim(fldnames_fr_lnd(2))
       fldlist_lnd2glc%flds(3)%shortname = trim(fldnames_fr_lnd(3))
       fldlist_lnd2glc%flds(1)%mapindex(compglc) = mapbilnr
       fldlist_lnd2glc%flds(2)%mapindex(compglc) = mapbilnr
       fldlist_lnd2glc%flds(3)%mapindex(compglc) = mapbilnr
       fldlist_lnd2glc%flds(1)%mapnorm(compglc) = 'lfrac'
       fldlist_lnd2glc%flds(2)%mapnorm(compglc) = 'lfrac'
       fldlist_lnd2glc%flds(3)%mapnorm(compglc) = 'lfrac'

       ! Create route handle if it has not been created
       if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(complnd,compglc,mapbilnr))) then
          call med_map_Fractions_init( gcomp, complnd, compglc, &
               FBSrc=FBlndAccum_lnd, &
               FBDst=FBlndAccum_glc, &
               RouteHandle=is_local%wrap%RH(complnd,compglc,mapbilnr), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
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

       ! -------------------------------
       ! If smb will be renormalized then...
       ! -------------------------------
       if (smb_renormalize) then

          ! -------------------------------
          ! get land and glc meshes and determine areas on corresponding meshes
          ! -------------------------------
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), mesh=lmesh_lnd, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), mesh=lmesh_glc, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! determine areas on land mesh
          call ESMF_MeshGet(lmesh_lnd, numOwnedElements=lsize, elementDistGrid=ldistgrid, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
          allocate(aream_l(lsize), dataptr1d(lsize))
          lArray = ESMF_ArrayCreate(ldistgrid, dataptr1d, rc=rc)
          call ESMF_MeshGet(lmesh_lnd, elemMaskArray=lArray, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
          aream_l(:) = dataptr1d(:)
          call ESMF_ArrayDestroy(larray, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
          deallocate(dataptr1d)

          ! determine areas on glc mesh
          call ESMF_MeshGet(lmesh_glc, numOwnedElements=lsize, elementDistGrid=ldistgrid, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
          allocate(aream_g(lsize), dataptr1d(lsize))
          lArray = ESMF_ArrayCreate(ldistgrid, dataptr1d, rc=rc)
          call ESMF_MeshGet(lmesh_glc, elemMaskArray=lArray, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
          aream_g(:) = dataptr1d(:)
          call ESMF_ArrayDestroy(larray, rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
          deallocate(dataptr1d)

          ! -------------------------------
          ! ice mask without elevation classes on glc and lnd
          ! -------------------------------
          call FB_init(FBglc_icemask, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBExp(compglc), fieldnameList=(/Sg_icemask_field/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return

          call FB_init(FBlnd_icemask, is_local%wrap%flds_scalar_name, &
               FBgeom=is_local%wrap%FBImp(complnd,complnd), fieldNameList=(/Sg_icemask_field/), rc=rc)
          if (chkerr(rc,__line__,u_file_u)) return

          allocate(fldlist_glc2lnd_icemask%flds(1))
          fldlist_glc2lnd_icemask%flds(1)%shortname = Sg_icemask_field
          fldlist_glc2lnd_icemask%flds(1)%mapindex(complnd) = mapconsf
          fldlist_glc2lnd_icemask%flds(1)%mapnorm(complnd) = 'none'

          ! -------------------------------
          ! ice fraction in multiple elevation classes on glc and lnd - NOTE that this includes bare land
          ! -------------------------------
          FBglc_frac = ESMF_FieldBundleCreate(rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name=trim(Sg_frac_field), meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBglc_frac, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          FBlnd_frac = ESMF_FieldBundleCreate(rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name=trim(Sg_frac_field), meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(FBlnd_frac, (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          allocate(fldlist_glc2lnd_frac%flds(1))
          fldlist_glc2lnd_frac%flds(1)%shortname = Sg_frac_field
          fldlist_glc2lnd_frac%flds(1)%mapindex(complnd) = mapconsf
          fldlist_glc2lnd_frac%flds(1)%mapnorm(complnd) = trim(Sg_icemask_field)  ! will use FBglc_icemask

          ! Create route handle if it has not been created
          if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compglc,complnd,mapconsf))) then
             call med_map_Fractions_init( gcomp, compglc, complnd, &
                  FBSrc=FBlndAccum_glc, &
                  FBDst=FBlndAccum_lnd, &
                  RouteHandle=is_local%wrap%RH(compglc,complnd,mapconsf), rc=rc)
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
    ! Accumulation and averaging is done on the land input to the river component on the land grid
    ! Mapping from the land to the glc grid is then done with the time averaged fields
    !---------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: i,n,ncnt
    real(r8), pointer   :: data2d_in(:,:)
    real(r8), pointer   :: data2d_out(:,:)
    logical             :: first_time = .true.
    character(len=*),parameter  :: subname='(med_phases_prep_glc_accum)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

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
    ! accumulator land input to glc on land grid
    !---------------------------------------

    if (ncnt > 0) then

       ! Initialize module variables needed to accumulate input to glc
       if (first_time) then
          call  med_phases_prep_glc_init(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          first_time = .false.
       end if

       do n = 1, size(fldnames_fr_lnd)
          call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), &
               fldnames_fr_lnd(n), fldptr2=data2d_in, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(FBlndAccum_lnd, fldnames_fr_lnd(n), fldptr2=data2d_out, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do i = 1,size(data2d_out, dim=2)
             data2d_out(:,i) = data2d_out(:,i) + data2d_in(:,i)
          end do
       end do

       FBlndAccumCnt = FBlndAccumCnt + 1

       if (dbug_flag > 1) then
          call FB_diagnose(FBlndAccum_lnd, string=trim(subname)// ' FBlndAccum_lnd ',  rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

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
    type(InternalState)    :: is_local
    integer                :: i, n, ncnt            ! counters
    real(r8), pointer      :: data2d(:,:)
    character(len=*) , parameter   :: subname='(med_phases_prep_glc_avg)'
    !---------------------------------------

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
    ! Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldCount=ncnt, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt ==  0) then
       call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBExp(compglc), returning", &
            ESMF_LOGMSG_INFO)

    else

       !---------------------------------------
       ! Average import from accumulated land import FB
       !---------------------------------------

       do n = 1, size(fldnames_fr_lnd)
          call FB_GetFldPtr(FBlndAccum_lnd, fldnames_fr_lnd(n), fldptr2=data2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          data2d(:,:) = data2d(:,:) / real(FBlndAccumCnt)
       end do

       if (dbug_flag > 1) then
          call FB_diagnose(FBlndAccum_lnd, string=trim(subname)//&
               ' FBlndAccum for after avg for field bundle ', rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       ! Map accumulated field bundle from land grid (with elevation classes) to glc grid (without elevation classes)
       ! and set FBExp(compglc) data
       !---------------------------------------

       call FB_reset(FBlndAccum_glc, value=0.0_r8, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call med_phases_prep_glc_map_lnd2glc(gcomp, rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(compglc), string=trim(subname)//' FBexp(compglc) ', rc=rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       endif

       !---------------------------------------
       ! zero accumulator and accumulated field bundles
       !---------------------------------------

       FBlndAccumCnt = 0

       call FB_reset(FBlndAccum_lnd, value=0.0_r8, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call FB_reset(FBlndAccum_glc, value=0.0_r8, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       ! update local scalar data - set valid input flag to .true.  TODO:
       !---------------------------------------

       !call med_infodata_set_valid_glc_input(.true., med_infodata, rc=rc) ! TODO: fix this
       !if (chkErr(rc,__LINE__,u_FILE_u)) return

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_avg

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
    real(r8), pointer   :: topolnd_g_ec(:,:)     ! topo in elevation classes
    real(r8), pointer   :: dataptr_g(:)          ! temporary data pointer for one elevation class
    real(r8), pointer   :: topoglc_g(:)          ! ice topographic height on the glc grid extracted from glc import
    real(r8), pointer   :: data_ice_covered_g(:) ! data for ice-covered regions on the GLC grid
    real(r8), pointer   :: glc_ice_covered(:)    ! if points on the glc grid is ice-covered (1) or ice-free (0)
    integer , pointer   :: glc_elevclass(:)      ! elevation classes glc grid
    real(r8), pointer   :: dataexp_g(:)          ! pointer into
    real(r8), pointer   :: dataptr2d(:,:)
    real(r8), pointer   :: dataptr1d(:)
    real(r8)            :: elev_l, elev_u        ! lower and upper elevations in interpolation range
    real(r8)            :: d_elev                ! elev_u - elev_l
    integer             :: nfld, ec
    integer             :: i,j,n,g,ncnt, lsize_g
    character(len=*) , parameter   :: subname='(med_phases_prep_glc_mod:med_phases_prep_glc_map_lnd2glc)'
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

    call FB_reset(FBlndAccum_glc, value=0.0_r8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call med_map_FB_Regrid_Norm( &
         fldsSrc=fldList_lnd2glc%flds, &
         srccomp=complnd, &
         destcomp=compglc, &
         FBSrc=FBlndAccum_lnd, &
         FBDst=FBlndAccum_glc, &
         FBFracSrc=is_local%wrap%FBFrac(complnd), &
         FBNormOne=is_local%wrap%FBNormOne(complnd,compglc,:), &
         RouteHandles=is_local%wrap%RH(complnd,compglc,:), &
         string=trim(compname(complnd))//'2'//trim(compname(compglc)), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call FB_diagnose(FBlndAccum_lnd, string=trim(subname)//' FBlndAccum_lnd ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call FB_diagnose(is_local%wrap%FBfrac(complnd), string=trim(subname)//' FBFrac ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
       call FB_diagnose(FBlndAccum_glc, string=trim(subname)//' FBlndAccum_glc ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc point on glc grid (output is topoglc_g)
    ! ------------------------------------------------------------------------

    if (dbug_flag > 1) then
       call FB_diagnose(is_local%wrap%FBImp(compglc,compglc), &
            string=trim(subname)//' FBImp(compglc,compglc) ', rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), Sg_frac_field, fldptr1=glc_ice_covered, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), Sg_topo_field, fldptr1=topoglc_g, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! get elevation classes with bare land
    ! for grid cells that are ice-free, the elevation class is set to 0.
    lsize_g = size(glc_ice_covered)
    allocate(glc_elevclass(lsize_g))
    call glc_get_elevation_classes(glc_ice_covered, topoglc_g, glc_elevclass, logunit)

    ! ------------------------------------------------------------------------
    ! Determine topo field in multiple elevation classes on the glc grid
    ! ------------------------------------------------------------------------

    call FB_getFldPtr(FBlndAccum_glc, 'Sl_topo_elev', fldptr2=topolnd_g_ec, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Loop over fields in export field bundle to GLC
    ! ------------------------------------------------------------------------

    ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
    ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
    ! current glint implementation, which sets acab and artm to 0 over ocean (although
    ! notes that this could lead to a loss of conservation). Figure out how to handle this case.

    allocate(data_ice_covered_g (lsize_g))
    do nfld = 1, size(fldnames_to_glc)

       ! ------------------------------------------------------------------------
       ! Perform vertical interpolation of data onto ice sheet topography
       ! This maps all of the input elevation classes into an export to glc without elevation classes
       ! ------------------------------------------------------------------------

       ! Get a pointer to the land data in multiple elevation classes on the glc grid
       call FB_getFldPtr(FBlndAccum_glc, fldnames_fr_lnd(nfld), fldptr2=dataptr2d, rc=rc)
       if (chkErr(rc,__LINE__,u_FILE_u)) return

       ! Get a pointer to the data for the field that will be sent to glc (without elevation classes)
       call FB_getFldPtr(is_local%wrap%FBExp(compglc), fldnames_to_glc(nfld), fldptr1=dataexp_g, rc=rc)
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

          if (glc_elevclass(n) /= 0) then
             ! ice-covered cells have interpolated values
             dataexp_g(n) = data_ice_covered_g(n)
          else
             ! non ice-covered cells have bare land value
             dataexp_g(n) = real(dataptr2d(1,n))
          end if
       end do  ! lsize_g

       ! ------------------------------------------------------------------------
       ! Renormalize surface mass balance (smb, here named dataexp_g) so that the global
       ! integral on the glc grid is equal to the global integral on the land grid.
       ! ------------------------------------------------------------------------

       ! No longer need to make a preemptive adjustment to qice_g to account for area differences
       ! between CISM and the coupler. In NUOPC, the area correction is done in! the cap not in the
       ! mediator, so to preserve the bilinear mapping values, do not need to do any area correction
       ! scaling in the CISM NUOPC cap

       if (smb_renormalize) then
          call med_phases_prep_glc_renormalize_smb(gcomp, rc)
          if (chkErr(rc,__LINE__,u_FILE_u)) return
       end if

    end do  ! end of loop over fields

    ! clean up memory
    deallocate(glc_elevclass)
    deallocate(data_ice_covered_g)

  end subroutine med_phases_prep_glc_map_lnd2glc

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
    ! For high-level design, see:
    ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit
    !------------------

    ! input/output variables
    type(ESMF_GridComp)    :: gcomp
    integer  , intent(out) :: rc          ! return error code

    ! local variables
    ! Note: Sg_icemask defines where the ice sheet model can receive a nonzero SMB from the land model.
    type(InternalState) :: is_local
    type(ESMF_VM)       :: vm
    type(ESMF_Mesh)     :: lmesh
    type(ESMF_Field)    :: lfield
    real(r8) , pointer  :: qice_l(:,:)        ! SMB (Flgl_qice) on land grid with elev classes
    real(r8) , pointer  :: qice_g(:)          ! SMB (Flgl_qice) on glc grid without elev classes
    real(r8) , pointer  :: glc_topo_g(:)      ! ice topographic height on the glc grid cell
    real(r8) , pointer  :: glc_frac_g(:)      ! total ice fraction in each glc cell
    real(r8) , pointer  :: glc_frac_g_ec(:,:) ! total ice fraction in each glc cell
    real(r8) , pointer  :: glc_frac_l_ec(:,:) ! EC fractions (Sg_ice_covered) on land grid
    real(r8) , pointer  :: Sg_icemask_g(:)    ! icemask on glc grid
    real(r8) , pointer  :: Sg_icemask_l(:)    ! icemask on land grid
    real(r8) , pointer  :: lfrac(:)           ! land fraction on land grid
    real(r8) , pointer  :: dataptr1d(:)       ! temporary 1d pointer
    real(r8) , pointer  :: dataptr2d(:,:)     ! temporary 2d pointer
    integer             :: ec                 ! loop index over elevation classes
    integer             :: n

    ! local and global sums of accumulation and ablation; used to compute renormalization factors
    real(r8), target :: local_accum_lnd(1), global_accum_lnd(1)
    real(r8), target :: local_accum_glc(1), global_accum_glc(1)
    real(r8), target :: local_ablat_lnd(1), global_ablat_lnd(1)
    real(r8), target :: local_ablat_glc(1), global_ablat_glc(1)

    ! renormalization factors (should be close to 1, e.g. in range 0.95 to 1.05)
    real(r8) :: accum_renorm_factor ! ratio between global accumulation on the two grids
    real(r8) :: ablat_renorm_factor ! ratio between global ablation on the two grids

    real(r8) :: effective_area      ! grid cell area multiplied by min(lfrac,Sg_icemask_l).
    ! This is the area that can contribute SMB to the ice sheet model.
    character(len=*), parameter  :: subname='(med_phases_prep_glc_renormalize_smb)'
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
    ! Get vm
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Map Sg_icemask_g from the glc grid to the land grid.
    !---------------------------------------

    ! determine Sg_icemask_g and set as contents of FBglc_icemask
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_icemask_field), fldptr1=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBglc_icemask, trim(Sg_icemask_field), fldptr1=Sg_icemask_g,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    Sg_icemask_g(:) = dataptr1d(:)

    ! BUG(wjs, 2017-05-11, #1516) I think we actually want norm = .false. here, but this
    ! requires some more thought
    call med_map_FB_Regrid_Norm( &
         fldsSrc=fldList_glc2lnd_icemask%flds, &
         srccomp=compglc, &
         destcomp=complnd, &
         FBSrc=FBglc_icemask, &
         FBDst=FBlnd_icemask, &
         FBFracSrc=FBglc_icemask, &  ! this will not be used since are asking for 'none' normalization
         FBNormOne=is_local%wrap%FBNormOne(compglc,complnd,:), &
         RouteHandles=is_local%wrap%RH(compglc,complnd,:), &
         string='mapping Sg_imask_g to Sg_imask_l (from glc to land)', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Map frac_field on glc grid without elevation classes to frac_field on land grid with elevation classes
    ! ------------------------------------------------------------------------

    ! set FBglc_frac on the glc grid (fractional ice coverage per elevation class)
    ! glc_topo_g(:) is the topographic height of each glc gridcell
    ! glc_frac_g(:) is the total ice fraction in each glc gridcell
    ! glc_frac_g_ec(:,:) are the glc fractions on the glc grid for each elevation class (inner dimension)
    ! setting glc_frac_g_ec (in the call to glc_get_fractional_icecov) sets the contents of FBglc_frac
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_topo_field), fldptr1=glc_topo_g,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_frac_field), fldptr1=glc_frac_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBglc_frac, Sg_frac_field, fldptr2=glc_frac_g_ec, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! note that nec = ungriddedCount - 1
    call glc_get_fractional_icecov(ungriddedCount-1, glc_topo_g, glc_frac_g, glc_frac_g_ec, logunit)

    ! map fraction in each elevation class from the glc grid to the land grid
    call med_map_FB_Regrid_Norm( &
         fldsSrc=fldList_glc2lnd_frac%flds, &
         srccomp=compglc, &
         destcomp=complnd, &
         FBSrc=FBglc_frac, & ! this has multiple elevation classes
         FBDst=FBlnd_frac, & ! this has multiple elvation classes
         FBFracSrc=FBglc_icemask, &  ! this is used with  a mapnorm of Sg_icemask_field
         FBNormOne=is_local%wrap%FBNormOne(compglc,complnd,:), & ! this will not be used
         RouteHandles=is_local%wrap%RH(compglc,complnd,:), &
         string='mapping elevation class fractions from glc to land ', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get fractional ice coverage for each elevation class on the land grid, glc_frac_l_ec(:,:)
    call FB_getFldPtr(FBlnd_frac, trim(Sg_frac_field), fldptr2=glc_frac_l_ec, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! determine fraction on land grid, lfrac(:)
    call FB_getFldPtr(is_local%wrap%FBFrac(complnd), 'lfrac', fldptr1=lfrac, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! get Sg_icemask_l(:)
    call FB_getFldPtr(FBlnd_icemask, Sg_icemask_field, fldptr1=Sg_icemask_l, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    ! determine qice_l
    call FB_getFldPtr(FBlndAccum_lnd, trim(qice_fieldname)//'_elev', fldptr2=qice_l, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Sum qice_l over all elevation classes for each local land grid cell then do a global sum
    !---------------------------------------

    local_accum_lnd(1) = 0.0_r8
    local_ablat_lnd(1) = 0.0_r8
    do n = 1, size(lfrac)
       ! Calculate effective area for sum -  need the mapped Sg_icemask_l
       effective_area = min(lfrac(n), Sg_icemask_l(n)) * aream_l(n)

       do ec = 1, ungriddedCount
          if (qice_l(ec,n) >= 0.0_r8) then
             local_accum_lnd(1) = local_accum_lnd(1) + effective_area * glc_frac_l_ec(ec,n) * qice_l(ec,n)
          else
             local_ablat_lnd(1) = local_ablat_lnd(1) + effective_area * glc_frac_l_ec(ec,n) * qice_l(ec,n)
          endif
       enddo  ! ec
    enddo  ! n
    call ESMF_VMAllreduce(vm, senddata=local_accum_lnd, recvdata=global_accum_lnd, count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    call ESMF_VMAllreduce(vm, senddata=local_ablat_lnd, recvdata=global_ablat_lnd, count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)

    !---------------------------------------
    ! Sum qice_g over local glc grid cells.
    !---------------------------------------

    ! TODO: is the following a problem
    ! Note: This sum uses the coupler areas (aream_g), which differ from the native CISM areas.
    ! But since the original qice_g (from bilinear remapping) has been multiplied by
    ! area_g/aream_g above, this calculation is equivalent to multiplying the original qice_g
    ! by the native CISM areas (area_g).
    ! If Flgl_qice were changed to a state (and not included in seq_flds_x2g_fluxes),
    ! then it would be appropriate to use the native CISM areas in this sum.

    ! determine qice_g
    call FB_getFldPtr(is_local%wrap%FBExp(compglc), qice_fieldname, fldptr1=qice_g, rc=rc)
    if (chkErr(rc,__LINE__,u_FILE_u)) return

    local_accum_glc(1) = 0.0_r8
    local_ablat_glc(1) = 0.0_r8
    do n = 1, size(qice_g)
       if (qice_g(n) >= 0.0_r8) then
          local_accum_glc(1) = local_accum_glc(1) + Sg_icemask_g(n) * aream_g(n) * qice_g(n)
       else
          local_ablat_glc(1) = local_ablat_glc(1) + Sg_icemask_g(n) * aream_g(n) * qice_g(n)
       endif
    enddo  ! n
    call ESMF_VMAllreduce(vm, senddata=local_accum_glc, recvdata=global_accum_glc, count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    call ESMF_VMAllreduce(vm, senddata=local_ablat_glc, recvdata=global_ablat_glc, count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)

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
       write(logunit,*) 'accum_renorm_factor = ', accum_renorm_factor
       write(logunit,*) 'ablat_renorm_factor = ', ablat_renorm_factor
    endif

    do n = 1, size(qice_g)
       if (qice_g(n) >= 0.0_r8) then
          qice_g(n) = qice_g(n) * accum_renorm_factor
       else
          qice_g(n) = qice_g(n) * ablat_renorm_factor
       endif
    enddo

  end subroutine med_phases_prep_glc_renormalize_smb

end module med_phases_prep_glc_mod
