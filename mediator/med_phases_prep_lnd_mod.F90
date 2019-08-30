module med_phases_prep_lnd_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing land export from mediator
  !-----------------------------------------------------------------------------

  use ESMF                  , only : operator(/=)
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_RouteHandle, ESMF_RouteHandleIsCreated
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_StateGet, ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshLoc, ESMF_MESHLOC_ELEMENT
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF                  , only : ESMF_TYPEKIND_R8
  use esmFlds               , only : complnd, compatm, compglc, ncomps, compname, mapconsf
  use esmFlds               , only : fldListFr, fldListTo
  use esmFlds               , only : shr_nuopc_fldlist_type
  use shr_nuopc_methods_mod , only : FB_getFieldN    => shr_nuopc_methods_FB_getFieldN
  use shr_nuopc_methods_mod , only : FB_getFldPtr    => shr_nuopc_methods_FB_getFldPtr
  use shr_nuopc_methods_mod , only : FB_getNumFlds   => shr_nuopc_methods_FB_getNumFlds
  use shr_nuopc_methods_mod , only : FB_init         => shr_nuopc_methods_FB_init
  use shr_nuopc_methods_mod , only : FB_diagnose     => shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : FB_FldChk       => shr_nuopc_methods_FB_FldChk
  use shr_nuopc_methods_mod , only : State_GetScalar => shr_nuopc_methods_State_GetScalar
  use shr_nuopc_methods_mod , only : State_SetScalar => shr_nuopc_methods_State_SetScalar
  use shr_nuopc_utils_mod   , only : chkerr          => shr_nuopc_utils_ChkErr
  use med_constants_mod     , only : R8, CS
  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, logunit
  use med_map_mod           , only : med_map_FB_Regrid_Norm
  use med_map_mod           , only : med_map_Fractions_Init
  use med_merge_mod         , only : med_merge_auto
  use glc_elevclass_mod     , only : glc_mean_elevation_virtual
  use glc_elevclass_mod     , only : glc_get_num_elevation_classes
  use glc_elevclass_mod     , only : glc_get_fractional_icecov
  use perf_mod              , only : t_startf, t_stopf
  use shr_sys_mod           , only : shr_sys_abort

  implicit none
  private

  public  :: med_phases_prep_lnd

  private :: med_map_glc2lnd_init
  private :: med_map_glc2lnd

  ! private module variables
  character(len =*), parameter :: Sg_icemask        = 'Sg_icemask'
  character(len =*), parameter :: Sg_frac           = 'Sg_ice_covered'
  character(len =*), parameter :: Sg_frac_x_icemask = 'Sg_frac_times_icemask'
  character(len =*), parameter :: Sg_topo           = 'Sg_topo'
  character(len =*), parameter :: Flgg_hflx         = 'Flgg_hflx'

  type(ESMF_FieldBundle) :: FBglc_icemask        ! no elevation classes
  type(ESMF_FieldBundle) :: FBglc_frac_x_icemask ! elevation classes
  type(ESMF_FieldBundle) :: FBlnd_frac_x_icemask ! elevation classes
  type(ESMF_FieldBundle) :: FBglc_ec
  type(ESMF_FieldBundle) :: FBlnd_ec

  type(shr_nuopc_fldlist_type) :: fldlist_glc2lnd_frac_x_icemask
  type(shr_nuopc_fldlist_type) :: fldlist_glc2lnd_norm_icemask
  type(shr_nuopc_fldlist_type) :: fldlist_glc2lnd_norm_none

  integer :: nec ! number of elevation classes (excluding bare land)

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_prep_lnd(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_StateItem_Flag) :: itemType
    type(InternalState)       :: is_local
    integer                   :: n1,ncnt
    real(r8)                  :: nextsw_cday
    logical                   :: first_call = .true.
    character(len=*), parameter :: subname='(med_phases_prep_lnd)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call FB_getNumFlds(is_local%wrap%FBExp(complnd), trim(subname)//"FBexp(complnd)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- map to create FBimp(:,complnd)
       !---------------------------------------

       do n1 = 1,ncomps
          if (is_local%wrap%med_coupling_active(n1,complnd)) then
             ! The following will map all atm->lnd, rof->lnd, Sg_icemask_field and Sg_topo_field
             call med_map_FB_Regrid_Norm( &
                  fldsSrc=fldListFr(n1)%flds, &
                  srccomp=n1, destcomp=complnd, &
                  FBSrc=is_local%wrap%FBImp(n1,n1), &
                  FBDst=is_local%wrap%FBImp(n1,complnd), &
                  FBFracSrc=is_local%wrap%FBFrac(n1), &
                  FBNormOne=is_local%wrap%FBNormOne(n1,complnd,:), &
                  RouteHandles=is_local%wrap%RH(n1,complnd,:), &
                  string=trim(compname(n1))//'2'//trim(compname(complnd)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       enddo

       !---------------------------------------
       !--- auto merges to create FBExp(complnd)
       !---------------------------------------

       ! The following will merge all fields in fldsSrc (for glc these
       ! are Sg_icemask and Sg_topo)
       call med_merge_auto(trim(compname(complnd)), &
            is_local%wrap%FBExp(complnd), &
            is_local%wrap%FBFrac(complnd), &
            is_local%wrap%FBImp(:,complnd), &
            fldListTo(complnd), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

       if (first_call) then
          call med_map_glc2lnd_init(gcomp, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          first_call = .false.
       end if

       ! The will following will map and merge Sg_frac and Sg_topo
       ! (and in the future Flgg_hflx)
       call med_map_glc2lnd(gcomp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       !--- update scalar data
       !---------------------------------------

       call ESMF_StateGet(is_local%wrap%NStateImp(compatm), trim(is_local%wrap%flds_scalar_name), itemType, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (itemType /= ESMF_STATEITEM_NOTFOUND) then
          ! send nextsw_cday to land - first obtain it from atm import
          call State_GetScalar(&
               scalar_value=nextsw_cday, &
               scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
               state=is_local%wrap%NstateImp(compatm), &
               flds_scalar_name=is_local%wrap%flds_scalar_name, &
               flds_scalar_num=is_local%wrap%flds_scalar_num, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call State_SetScalar(&
               scalar_value=nextsw_cday, &
               scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
               state=is_local%wrap%NstateExp(complnd), &
               flds_scalar_name=is_local%wrap%flds_scalar_name, &
               flds_scalar_num=is_local%wrap%flds_scalar_num, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- diagnose
       !---------------------------------------

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(complnd), &
               string=trim(subname)//' FBexp(complnd) ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- clean up
       !---------------------------------------

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_lnd

  !================================================================================================

  subroutine med_map_glc2lnd_init(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: lfield
    type(ESMF_Mesh)     :: lmesh_lnd
    type(ESMF_Mesh)     :: lmesh_glc
    character(len=*) , parameter   :: subname='(med_map_glc2lnd_mod:med_map_glc2lnd_init)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Set the module variable for the number of elevation classes
    !---------------------------------------

    nec = glc_get_num_elevation_classes()

    !---------------------------------------
    ! Get the glc and land meshes
    !---------------------------------------

    call FB_getFieldN(is_local%wrap%FBExp(complnd), 1, field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh_lnd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_getFieldN(is_local%wrap%FBImp(compglc,compglc), 1, field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh_glc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------
    ! Create module field bundlel ice masks on glc grid (no elevation classes)
    ! -------------------------------

    FBglc_icemask = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name=trim(Sg_icemask), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBglc_icemask, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------
    ! Create module field bundle for icefrac times icemask in multiple elevation classes
    ! on glc and land (fraction in this elevation class) x (icemask) for a given elevation class
    ! -------------------------------

    FBglc_frac_x_icemask = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name=trim(Sg_frac_x_icemask), meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/nec+1/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBglc_frac_x_icemask, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    FBlnd_frac_x_icemask = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name=trim(Sg_frac_x_icemask), meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/nec+1/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBlnd_frac_x_icemask, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(fldlist_glc2lnd_frac_x_icemask%flds(1))
    fldlist_glc2lnd_frac_x_icemask%flds(1)%shortname = trim(Sg_frac_x_icemask)
    fldlist_glc2lnd_frac_x_icemask%flds(1)%mapindex(complnd) = mapconsf
    fldlist_glc2lnd_frac_x_icemask%flds(1)%mapnorm(complnd) = 'none'

    ! -------------------------------
    ! Create generic multiple elevation module field bundles that can be re-used
    ! -------------------------------

    FBglc_ec = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name='field_ec', meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/nec+1/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBglc_ec, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    FBlnd_ec = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name='field_ec', meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/nec+1/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBlnd_ec, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! determine different mapping types and normalization for generic field bundles above
    allocate(fldlist_glc2lnd_norm_icemask%flds(1))
    fldlist_glc2lnd_norm_icemask%flds(1)%shortname = 'field_ec'
    fldlist_glc2lnd_norm_icemask%flds(1)%mapindex(complnd) = mapconsf
    fldlist_glc2lnd_norm_icemask%flds(1)%mapnorm(complnd) = trim(Sg_icemask)

    allocate(fldlist_glc2lnd_norm_none%flds(1))
    fldlist_glc2lnd_norm_none%flds(1)%shortname = 'field_ec'
    fldlist_glc2lnd_norm_none%flds(1)%mapindex(complnd) = mapconsf
    fldlist_glc2lnd_norm_none%flds(1)%mapnorm(complnd) = 'none'

    ! -------------------------------
    ! Create route handle if it has not been created
    ! -------------------------------

    if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compglc,complnd,mapconsf))) then
       call med_map_Fractions_init( gcomp, compglc, complnd, &
            FBSrc=FBglc_ec, FBDst=FBlnd_ec, &
            RouteHandle=is_local%wrap%RH(compglc,complnd,mapconsf), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! -------------------------------
    ! Currently cannot map hflx in multiple elevation classes from glc to land
    ! -------------------------------

    if (FB_fldchk(is_local%wrap%FBExp(complnd), trim(Flgg_hflx), rc=rc)) then
       call shr_sys_abort('ERROR: Flgg_hflx to land has not been implemented yet')
    end if

  end subroutine med_map_glc2lnd_init

  !================================================================================================

  subroutine med_map_glc2lnd( gcomp, rc)

    !------------------
    ! Maps fields from the GLC grid to the LND grid.
    ! On the GLC grid the fields will not have elevation classes.
    ! On the LND grid they will have elevation classes.
    !------------------

    ! input/output variables
    type(ESMF_GridComp)    , intent(inout) :: gcomp
    integer                , intent(out)   :: rc

    ! local variables
    type(InternalState)   :: is_local
    integer               :: ec, nfld, n, l
    real(r8)              :: topo_virtual
    real(r8), pointer     :: icemask_g(:)             ! glc ice mask field on glc grid
    real(r8), pointer     :: frac_g(:)                ! total ice fraction in each glc cell
    real(r8), pointer     :: frac_g_ec(:,:)           ! glc fractions on the glc grid
    real(r8), pointer     :: frac_l_ec(:,:)           ! glc fractions on the land grid
    real(r8), pointer     :: frac_x_icemask_g_ec(:,:) ! (glc fraction) x (icemask), on the glc grid
    real(r8), pointer     :: topo_g(:)                ! topographic height of each glc cell (no elevation classes)
    real(r8), pointer     :: topo_l_ec(:,:)           ! topographic height in each land gridcell for each elevation class
    real(r8), pointer     :: topo_x_icemask_g(:,:)
    real(r8), pointer     :: topo_x_icemask_l(:,:)
    real(r8), pointer     :: frac_x_icemask_g(:,:)
    real(r8), pointer     :: frac_x_icemask_l(:,:)
    real(r8), pointer     :: dataptr1d(:)
    real(r8), pointer     :: dataptr2d_exp(:,:)
    character(len=*), parameter :: subname = 'med_map_glc2lnd'
    !-----------------------------------------------------------------------

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
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------
    ! Map fractional ice coverage to the land grid (multiple elevation classes)
    !---------------------------------

    ! Set contents of FBglc_ec to contain frac_g_ec
    ! (fractional ice coverage for each elevation class on the glc grid)

    ! set FBglc_ec on the glc grid (fractional ice coverage per elevation class)
    ! topo_g(:) is the topographic height of each glc gridcell
    ! frac_g(:) is the total ice fraction in each glc gridcell
    ! frac_g_ec(:,:) are the glc fractions on the glc grid for each elevation class (inner dimension)
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_topo), fldptr1=topo_g,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_frac), fldptr1=frac_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBglc_ec, 'field_ec', fldptr2=frac_g_ec, rc=rc) 
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! compute frac_g_ec
    call glc_get_fractional_icecov(nec, topo_g, frac_g, frac_g_ec, logunit)

    ! Set the contents of FBglc_icemask
    call FB_getFldPtr(FBglc_icemask, trim(Sg_icemask), fldptr1=icemask_g,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_icemask), fldptr1=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    icemask_g(:) = dataptr1d(:)

    ! Only include grid cells that are both (a) within the icemask and (b) in this elevation class
    call FB_getFldPtr(FBglc_frac_x_icemask, trim(Sg_frac_x_icemask), fldptr2=frac_x_icemask_g_ec, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do ec = 1, nec+1
       frac_x_icemask_g_ec(ec,:) = frac_g_ec(ec,:) * icemask_g(:)
    end do

    ! map to lnd and normalize by Sg_icemask_field
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": calling mapping elevation class fractions from glc to land", ESMF_LOGMSG_INFO)
    end if
    call med_map_FB_Regrid_Norm( &
         fldsSrc=fldList_glc2lnd_norm_icemask%flds, &
         srccomp=compglc, &
         destcomp=complnd, &
         FBSrc=FBglc_ec, & ! this has multiple elevation classes
         FBDst=FBlnd_ec, & ! this has multiple elvation classes
         FBFracSrc=FBglc_icemask, &  ! this is used with  a mapnorm of Sg_icemask_field
         FBNormOne=is_local%wrap%FBNormOne(compglc,complnd,:), & ! this will not be used
         RouteHandles=is_local%wrap%RH(compglc,complnd,:), &
         string='mapping elevation class fractions from glc to land ', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(FBlnd_ec, 'field_ec', fldptr2=frac_l_ec, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBExp(complnd), trim(Sg_frac)//'_elev', fldptr2=dataptr2d_exp, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr2d_exp(:,:) = frac_l_ec(:,:)

    !---------------------------------
    ! Map topo to the land grid (multiple elevation classes)
    !---------------------------------

    ! Note that all topo values in FBimp(compglc,compglc) do not have elevation class dependence
    ! Normalize by frac_x_icemask_g_ec - this is what introduces
    ! elevation class information from the glc grid (without elevation classes) to the
    ! land grid (with elevation classes)
    ! Note that bare land values are mapped in the same way as ice-covered values

    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_topo), fldptr1=topo_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBglc_ec, 'field_ec', fldptr2=topo_x_icemask_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ec = 1,nec+1
       do l = 1,size(topo_g)
          topo_x_icemask_g(:,:) = topo_g(l) * frac_x_icemask_g_ec(ec,l)
       end do
    end do

    ! map FBglc_topo_x_icemask from glc to land (with multiple elevation classes) - no normalization
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": calling mapping of topo from glc to land", ESMF_LOGMSG_INFO)
    end if
    call med_map_FB_Regrid_Norm( &
         fldsSrc=fldList_glc2lnd_norm_none%flds, &
         srccomp=compglc, &
         destcomp=complnd, &
         FBSrc=FBglc_ec, &
         FBDst=FBlnd_ec, &
         FBFracSrc=is_local%wrap%FBFrac(compglc), & ! this will not be used
         FBNormOne=is_local%wrap%FBNormOne(compglc,complnd,:), & ! this will not be used
         RouteHandles=is_local%wrap%RH(compglc,complnd,:), &
         string='mapping topo from glc to land (with elevation classes)', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBlnd_ec, 'field_ec', fldptr2=topo_l_ec , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! map FBglc_frac_x_icemask from glc to land (with multiple elevation classes) - no normalization
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": calling mapping of frac_x_icemask from glc to land", ESMF_LOGMSG_INFO)
    end if
    call med_map_FB_Regrid_Norm( &
         fldsSrc=fldList_glc2lnd_frac_x_icemask%flds, &
         srccomp=compglc, &
         destcomp=complnd, &
         FBSrc=FBglc_frac_x_icemask, &
         FBDst=FBlnd_frac_x_icemask, &
         FBFracSrc=is_local%wrap%FBFrac(compglc), & ! this will not be used
         FBNormOne=is_local%wrap%FBNormOne(compglc,complnd,:), & ! this will not be used
         RouteHandles=is_local%wrap%RH(compglc,complnd,:), &
         string='mapping frac_x_icemask from glc to land (with elevation classes)', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBlnd_ec, 'field_ec', fldptr2=frac_x_icemask_l , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! set Sg_topo values in export state to land (in multiple elevation classes)
    ! also set the topo field for virtual columns, in a given elevation class.
    ! This is needed because virtual columns (i.e., elevation classes that have no
    ! contributing glc grid cells) won't have any topographic information mapped onto
    ! them, so would otherwise end up with an elevation of 0.
    call FB_getFldPtr(is_local%wrap%FBExp(complnd), trim(Sg_topo)//'_elev', fldptr2=dataptr2d_exp, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ec = 1,nec+1
       !write(6,*)'DEBUG: ec = ',ec
       topo_virtual = glc_mean_elevation_virtual(ec)
       !write(6,*)'DEBUG: ec, topo_virtual = ',ec,topo_virtual
       do l = 1,size(frac_x_icemask_l, dim=2)
          if (frac_l_ec(ec,l) <= 0._r8) then
             dataptr2d_exp(ec,l) = topo_virtual
          else
             if (frac_x_icemask_l(ec,l) == 0.0_r8) then
                dataptr2d_exp(ec,l) = 0.0_r8
             else
                dataptr2d_exp(ec,l) = topo_l_ec(ec,l) / frac_x_icemask_l(ec,l)
             end if
          end if
       end do
    end do

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
       call t_stopf('MED:'//subname)
    end if

  end subroutine med_map_glc2lnd

end module med_phases_prep_lnd_mod
