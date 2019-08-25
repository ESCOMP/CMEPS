module med_phases_prep_lnd_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing land export from mediator
  !-----------------------------------------------------------------------------

  use ESMF                  , only : operator(/=)
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_RouteHandle
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
  use shr_nuopc_methods_mod , only : State_GetScalar => shr_nuopc_methods_State_GetScalar
  use shr_nuopc_methods_mod , only : State_SetScalar => shr_nuopc_methods_State_SetScalar
  use shr_nuopc_utils_mod   , only : chkerr          => shr_nuopc_utils_ChkErr
  use med_constants_mod     , only : R8, CS
  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, logunit
  use med_map_mod           , only : med_map_FB_Regrid_Norm
  use med_merge_mod         , only : med_merge_auto
  use glc_elevclass_mod     , only : glc_mean_elevation_virtual
  use glc_elevclass_mod     , only : glc_get_num_elevation_classes
  use glc_elevclass_mod     , only : glc_get_fractional_icecov
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  ! private module variables
  character(len=*), parameter :: Sg_icemask_field = 'Sg_icemask'
  character(len=*), parameter :: Sg_frac_field    = 'Sg_ice_covered'
  character(len=*), parameter :: Sg_topo_field    = 'Sg_topo'
  character(len=*), parameter :: Flgg_hflx_field  = 'Flgg_hflx'

  type(ESMF_FieldBundle) :: FBglc_icemask ! no elevation classes
  type(ESMF_FieldBundle) :: FBlnd_icemask ! no elevation classes
  type(ESMF_FieldBundle) :: FBglc_frac
  type(ESMF_FieldBundle) :: FBlnd_frac
  type(ESMF_FieldBundle) :: FBglc_frac_times_icemask
  type(ESMF_FieldBundle) :: FBlnd_frac_times_icemask
  type(ESMF_FieldBundle) :: FBglc_topo
  type(ESMF_FieldBundle) :: FBlnd_topo
  type(ESMF_FieldBundle) :: FBglc_hflx
  type(ESMF_FieldBundle) :: FBlnd_hflx

  type(shr_nuopc_fldlist_type) :: fldlist_glc2lnd_frac
  type(shr_nuopc_fldlist_type) :: fldlist_glc2lnd_topo
  type(shr_nuopc_fldlist_type) :: fldlist_glc2lnd_hflx

  integer :: num_elevation_classes
  integer :: nec

  character(*) , parameter :: u_FILE_u = &
       __FILE__

  public  :: med_phases_prep_lnd
  private :: med_map_glc2lnd_init
  private :: med_map_glc2lnd

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
             if (n1 == compglc) then
                if (first_call) then
                   call med_map_glc2lnd_init(gcomp, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   first_call = .false.
                end if
                call med_map_glc2lnd(gcomp, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
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
          endif
       enddo

       !---------------------------------------
       !--- auto merges to create FBExp(complnd)
       !---------------------------------------

       call med_merge_auto(trim(compname(complnd)), &
            is_local%wrap%FBExp(complnd), &
            is_local%wrap%FBFrac(complnd), &
            is_local%wrap%FBImp(:,complnd), &
            fldListTo(complnd), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(complnd), &
               string=trim(subname)//' FBexp(complnd) ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

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
    integer             :: ubound, ec
    type(ESMF_Field)    :: lfield 
    type(ESMF_Mesh)     :: lmesh_lnd
    type(ESMF_Mesh)     :: lmesh_glc
    real(r8), pointer   :: dataptr1d(:)
    real(r8), pointer   :: dataptr2d(:,:)
    character(len=*) , parameter   :: subname='(med_map_glc2lnd_mod:med_map_glc2lnd_init)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    num_elevation_classes = glc_get_num_elevation_classes()
    ubound = num_elevation_classes 
    
    !---------------------------------------
    ! Get the glc and lnd meshes
    !---------------------------------------

    call FB_getFieldN(is_local%wrap%FBImp(complnd,complnd), 1, field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh_lnd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_getFieldN(is_local%wrap%FBImp(compglc,compglc), 1, field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh_glc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------
    ! mapped ice mask on land (no elevation classes)
    ! -------------------------------

    FBglc_icemask = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name=trim(Sg_icemask_field), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBglc_icemask, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------
    ! ice fraction in multiple elevation classes on glc and lnd
    ! -------------------------------

    FBglc_frac = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name=trim(Sg_frac_field), meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/0/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBglc_frac, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    FBlnd_frac = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name=trim(Sg_frac_field), meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/0/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBlnd_frac, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(fldlist_glc2lnd_frac%flds(1))
    fldlist_glc2lnd_frac%flds(1)%shortname = Sg_frac_field
    fldlist_glc2lnd_frac%flds(1)%mapindex(complnd) = mapconsf
    fldlist_glc2lnd_frac%flds(1)%mapnorm(complnd) = trim(Sg_icemask_field)

    ! -------------------------------
    ! icefrac time icemask  in multiple elevation classes on glc and land
    ! (fraction in this elevation class) x (icemask) for a given elevation class
    ! -------------------------------

    FBglc_frac_times_icemask = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name='Sg_frac_times_icemask', meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBglc_frac_times_icemask, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    FBlnd_frac_times_icemask = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name='Sg_frac_times_icemask', meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBlnd_frac_times_icemask, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------
    ! topo in multiple elevation classes on glc grid and land grid
    ! -------------------------------

    ! Note that all topo values in FBglc_topo will be the same for each elevation class - the elevation class distinctio
    ! will arise due to the scalaing by glc_frac_times_icemask in the normalization of the mapping

    FBglc_topo = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name=trim(Sg_topo_field), meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBglc_topo, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    FBlnd_topo = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name=trim(Sg_topo_field), meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBlnd_topo, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_topo_field), fldptr1=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBglc_topo, trim(Sg_topo_field), fldptr2=dataptr2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ec = 0, num_elevation_classes
       dataptr2d(ec+1,:) = dataptr1d(:)
    end do

    allocate(fldlist_glc2lnd_topo%flds(1))
    fldlist_glc2lnd_topo%flds(1)%shortname = trim(Sg_topo_field)
    fldlist_glc2lnd_topo%flds(1)%mapindex(complnd) = mapconsf
    fldlist_glc2lnd_topo%flds(1)%mapnorm(complnd) = 'Sg_frac_times_icemask'

    ! -------------------------------
    ! hflx glc and lnd in multiple elevation classes
    ! -------------------------------

    FBglc_hflx  = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_glc, ESMF_TYPEKIND_R8, name=trim(Flgg_hflx_field), meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBglc_hflx, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Flgg_hflx_field), fldptr1=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBglc_topo, trim(Flgg_hflx_field), fldptr2=dataptr2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ec = 0, num_elevation_classes
       dataptr2d(ec+1,:) = dataptr1d(:)
    end do

    FBlnd_hflx  = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh_lnd, ESMF_TYPEKIND_R8, name=trim(Flgg_hflx_field), meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBlnd_hflx, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(fldlist_glc2lnd_hflx%flds(1))
    fldlist_glc2lnd_hflx%flds(1)%shortname = trim(Flgg_hflx_field)
    fldlist_glc2lnd_hflx%flds(1)%mapindex(complnd) = mapconsf
    fldlist_glc2lnd_hflx%flds(1)%mapnorm(complnd) = 'Sg_frac_times_icemask'

  end subroutine med_map_glc2lnd_init

  !================================================================================================

  subroutine med_map_glc2lnd( gcomp, rc)

    !------------------
    ! Maps fields from the GLC grid to the LND grid.
    ! On the GLC grid the fields will not have elevation classes.
    ! On the LND grid they will have elevation classes.
    !
    ! Maps frac_field, topo_field, plus all fields defined in extra_fields (for now just hflx_field)
    ! without the '_elev' suffix.
    !
    ! Assumes that 
    ! - FBglc_icemask contains icemask_field (NOT mapped here, but needed as an input to the mapping)
    ! - FBglc_frac  and FBlnd_frac  contain frac_field
    !
    ! Assumes that FBlnd (on land grid) contains:
    ! - frac_field (with multiple elevation classes in undistributed dimension)
    ! - topo_field (with multiple elevation classes in undistributed dimension)
    ! - And similarly for each field in extra_fields
    !------------------

    ! input/output variables
    type(ESMF_GridComp)    , intent(inout) :: gcomp
    integer                , intent(out)   :: rc

    ! local variables
    type(InternalState)    :: is_local
    integer                :: ec, nfld, n, l
    real(r8)               :: topo_virtual
    integer  , allocatable :: glc_elevclass(:)                   ! elevation class of each glc cell (assuming cell is ice-covered)
    real(r8) , pointer     :: glc_icemask_g(:)                   ! glc ice mask field on glc grid
    real(r8) , pointer     :: glc_frac_g(:)                      ! total ice fraction in each glc cell
    real(r8) , pointer     :: glc_frac_g_ec(:,:)                 ! glc fractions on the glc grid
    real(r8) , pointer     :: glc_frac_l_ec(:,:)                 ! glc fractions on the land grid
    real(r8) , pointer     :: glc_frac_times_icemask_g_ec(:,:)   ! (glc fraction) x (icemask), on the glc grid
    real(r8) , pointer     :: glc_topo_g(:)                      ! topographic height of each glc cell
    real(r8) , pointer     :: dataptr1d(:)
    real(r8) , pointer     :: dataptr2d(:,:)
    real(r8) , pointer     :: dataptr2d_exp(:,:)
    integer                :: lsize_g 
    character(len=3)       :: cvalue
    character(len=*), parameter :: subname = 'map_glc2lnd_elevclass'
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
    ! Set FBglc_frac (fractional ice coverage for each elevation class on the glc grid)
    !---------------------------------

    nec = glc_get_num_elevation_classes()

    ! set FBglc_frac on the glc grid (fractional ice coverage per elevation class)
    ! glc_topo_g(:) is the topographic height of each glc gridcell
    ! glc_frac_g(:) is the total ice fraction in each glc gridcell
    ! glc_frac_g_ec(:,:) are the glc fractions on the glc grid for each elevation class (inner dimension)
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_topo_field), fldptr1=glc_topo_g,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_frac_field), fldptr1=glc_frac_g, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBglc_frac, Sg_frac_field, fldptr2=glc_frac_g_ec, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call glc_get_fractional_icecov(nec, glc_topo_g, glc_frac_g, glc_frac_g_ec, logunit)

    !---------------------------------------
    ! Set FBglc_icemask and  FBglc_frac_times_icemask on the glc grid
    !---------------------------------------

    call FB_getFldPtr(FBglc_icemask, trim(Sg_icemask_field), &
         fldptr1=glc_icemask_g,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), trim(Sg_icemask_field), &
         fldptr1=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    glc_icemask_g(:) = dataptr1d(:)

    ! Only include grid cells that are both (a) within the icemask and (b) in this elevation class
    call FB_getFldPtr(FBglc_frac_times_icemask, 'Sg_frac_times_icemask', &
         fldptr2=glc_frac_times_icemask_g_ec, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do ec = 0, num_elevation_classes
       glc_frac_times_icemask_g_ec(ec+1,:) = glc_frac_g_ec(ec+1,:) * glc_icemask_g(:)
    end do

    !---------------------------------
    ! Map fractional ice coverage to the land grid (multiple elevation classes)
    !---------------------------------

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

    call FB_getFldPtr(FBlnd_frac, trim(Sg_frac_field), &
         fldptr2=dataptr2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBExp(complnd), trim(Sg_frac_field)//'_elev', &
         fldptr2=dataptr2d_exp, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr2d_exp(:,:) = dataptr2d(:,:)

    !---------------------------------
    ! Map topo to the land grid (multiple elevation classes)
    !---------------------------------

    call med_map_FB_Regrid_Norm( &
         fldsSrc=fldList_glc2lnd_topo%flds, &
         srccomp=compglc, &
         destcomp=complnd, &
         FBSrc=FBglc_topo, &
         FBDst=FBlnd_topo, &
         FBFracSrc=FBglc_frac_times_icemask, & ! this will be used with a mapnorm of 'Sg_frac_times_icemask'
         FBNormOne=is_local%wrap%FBNormOne(compglc,complnd,:), & ! this will not be used
         RouteHandles=is_local%wrap%RH(compglc,complnd,:), &
         string='mapping topo from glc to land (with elevation classes)', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(FBlnd_topo, trim(Sg_topo_field), &
         fldptr2=dataptr2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBExp(complnd), trim(Sg_topo_field)//'_elev', &
         fldptr2=dataptr2d_exp, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr2d_exp(:,:) = dataptr2d(:,:)

    ! set the topo field for virtual columns, in a given elevation class.
    ! This is needed because virtual columns (i.e., elevation classes that have no
    ! contributing glc grid cells) won't have any topographic information mapped onto
    ! them, so would otherwise end up with an elevation of 0.
       
    call FB_getFldPtr(is_local%wrap%FBExp(complnd), trim(Sg_topo_field)//'_elev', fldptr2=dataptr2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(FBlnd_frac, trim(Sg_frac_field), fldptr2=glc_frac_l_ec, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do ec = 0, num_elevation_classes
       topo_virtual = glc_mean_elevation_virtual(ec)
       do n = 1,size(dataptr2d, dim=2)
          if (glc_frac_l_ec(ec+1,n) <= 0._r8) then
             dataptr2d(ec+1,n) = topo_virtual
          end if
       end do
    end do  

    !---------------------------------
    ! Map hflx field to the land grid (multiple elvation lcasses)
    !---------------------------------

    ! Normalizing by glc_frac_times_icemask_g_ec - this is what introduces
    ! elevation class information from the glc grid (without elevation classes) to the 
    ! land grid (with elevation classes)
    ! Note that bare land values are mapped in the same way as ice-covered values

    call med_map_FB_Regrid_Norm( &
         fldsSrc=fldList_glc2lnd_hflx%flds, &
         srccomp=compglc, &
         destcomp=complnd, &
         FBSrc=FBglc_hflx, &
         FBDst=FBlnd_hflx, &
         FBFracSrc=FBglc_frac_times_icemask, & ! this will be used with a mapnorm of 'Sg_frac_times_icemask'
         FBNormOne=is_local%wrap%FBNormOne(compglc,complnd,:), & ! this will not be used
         RouteHandles=is_local%wrap%RH(compglc,complnd,:), &
         string='mapping hflx from glc to land (with elevation classes)', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(FBlnd_topo, trim(Flgg_hflx_field), &
         fldptr2=dataptr2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBExp(complnd), trim(Flgg_hflx_field)//'_elev', &
         fldptr2=dataptr2d_exp, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr2d_exp(:,:) = dataptr2d(:,:)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
       call t_stopf('MED:'//subname)
    end if

  end subroutine med_map_glc2lnd

end module med_phases_prep_lnd_mod
