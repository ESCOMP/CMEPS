module med_phases_prep_lnd_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing land export from mediator
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use ESMF                  , only : operator(/=)
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_StateGet, ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshLoc, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF                  , only : ESMF_RouteHandle, ESMF_RouteHandleIsCreated
  use esmFlds               , only : complnd, compatm, comprof, ncomps, compname
  use esmFlds               , only : max_icesheets, num_icesheets, compglc
  use esmFlds               , only : mapbilnr, mapconsd, mapconsf, compname
  use esmFlds               , only : fldListTo
  use med_methods_mod       , only : FB_diagnose     => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_FldChk       => med_methods_FB_fldchk
  use med_methods_mod       , only : State_GetScalar => med_methods_State_GetScalar
  use med_methods_mod       , only : State_SetScalar => med_methods_State_SetScalar
  use med_utils_mod         , only : chkerr          => med_utils_ChkErr
  use med_constants_mod     , only : dbug_flag       => med_constants_dbug_flag
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_map_mod           , only : med_map_rh_is_created, med_map_routehandles_init
  use med_map_mod           , only : med_map_field_packed, med_map_field_normalized, med_map_field
  use med_merge_mod         , only : med_merge_auto
  use glc_elevclass_mod     , only : glc_get_num_elevation_classes
  use glc_elevclass_mod     , only : glc_mean_elevation_virtual
  use glc_elevclass_mod     , only : glc_get_fractional_icecov
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_lnd

  private :: map_glc2lnd_init
  private :: map_glc2lnd

  ! private module variables
  character(len =*), parameter :: Sg_icemask = 'Sg_icemask'
  character(len =*), parameter :: Sg_icemask_coupled_fluxes = 'Sg_icemask_coupled_fluxes'
  character(len =*), parameter :: Sg_frac = 'Sg_ice_covered'
  character(len =*), parameter :: Sg_frac_x_icemask = 'Sg_frac_times_icemask'
  character(len =*), parameter :: Sg_topo = 'Sg_topo'
  character(len =*), parameter :: Flgg_hflx = 'Flgg_hflx'

  type, public :: ice_sheet_tolnd_type
     character(CS)    :: name
     logical          :: is_active
     type(ESMF_Field) :: field_icemask_g           ! no elevation classes
     type(ESMF_Field) :: field_frac_g_ec           ! elevation classes
     type(ESMF_Field) :: field_frac_x_icemask_g_ec ! elevation classes
     type(ESMF_Field) :: field_topo_x_icemask_g_ec ! elevation classes
     type(ESMF_Mesh)  :: mesh_g
  end type ice_sheet_tolnd_type
  type(ice_sheet_tolnd_type) :: ice_sheet_tolnd(max_icesheets)

  type(ESMF_field) :: field_icemask_l                ! no elevation classes
  type(ESMF_Field) :: field_frac_l_ec                ! elevation classes
  type(ESMF_Field) :: field_frac_x_icemask_l_ec      ! elevation classes
  type(ESMF_Field) :: field_topo_x_icemask_l_ec      ! elevation classes

  ! the number of elevation classes (excluding bare land) = ungriddedCount - 1
  integer :: ungriddedCount ! this equals the number of elevation classes + 1 (for bare land)

  logical :: glc2lnd_coupling = .false.
  logical :: cism_evolve = .false.

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
    integer                   :: n1,ncnt,ns
    real(r8)                  :: nextsw_cday
    logical                   :: first_call = .true.
    logical                   :: isPresent
    character(CL)             :: cvalue
    character(len=*), parameter :: subname='(med_phases_prep_lnd)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine if coupling to CISM is 2-way
    if (first_call) then
       call NUOPC_CompAttributeGet(gcomp, name="cism_evolve", isPresent=isPresent, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (isPresent) then
          call NUOPC_CompAttributeGet(gcomp, name="cism_evolve", value=cvalue, isPresent=isPresent, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          read (cvalue,*) cism_evolve
          if (mastertask) then
             write(logunit,'(a,l7)') trim(subname)//' cism_evolve = ',cism_evolve
          end if
       end if
    end if

    ! Count the number of fields outside of scalar data, if zero, then return
    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), fieldCount=ncnt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       ! map to create FBimp(:,complnd) - other than glc->lnd
       !---------------------------------------

       call t_startf('MED:'//trim(subname)//' map')
       do n1 = 1,ncomps
          ! Skip glc here and handle it below - only components that are allowed to couple to
          ! land are atm, rof and glc
          if (is_local%wrap%med_coupling_active(n1,complnd)) then
             if (n1 == comprof .or. n1 == compatm) then
                call med_map_field_packed( &
                     FBSrc=is_local%wrap%FBImp(n1,n1), &
                     FBDst=is_local%wrap%FBImp(n1,complnd), &
                     FBFracSrc=is_local%wrap%FBFrac(n1), &
                     field_normOne=is_local%wrap%field_normOne(n1,complnd,:), &
                     packed_data=is_local%wrap%packed_data(n1,complnd,:), &
                     routehandles=is_local%wrap%RH(n1,complnd,:), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if
       end do
       call t_stopf('MED:'//trim(subname)//' map')

       !---------------------------------------
       ! auto merges to create FBExp(complnd) - other than glc->lnd
       !---------------------------------------

       ! The following will merge all fields in fldsSrc
       ! Note that there is no merge for glc Sg_icemask and Sg_icemask_coupled_fluxes)
       ! these are mapped and set directly in map_glc2lnd

       call t_startf('MED:'//trim(subname)//' merge')
       call med_merge_auto(complnd, &
            is_local%wrap%med_coupling_active(:,complnd), &
            is_local%wrap%FBExp(complnd), &
            is_local%wrap%FBFrac(complnd), &
            is_local%wrap%FBImp(:,complnd), &
            fldListTo(complnd), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' merge')

       !---------------------------------------
       ! mapping and merging for glc->lnd
       !---------------------------------------

       ! determine if there will be any glc to lnd coupling
       if (first_call) then
          do ns = 1,num_icesheets
             if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
                glc2lnd_coupling = .true.
                exit
             end if
          end do
       end if

       ! the following is only done if glc->lnd coupling is active
       if (glc2lnd_coupling) then
          if (first_call) then
             call t_startf('MED:'//trim(subname)//' glc2lnd init')
             call map_glc2lnd_init(gcomp, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call t_stopf('MED:'//trim(subname)//' glc2lnd init')
          end if
          ! The will following will map and merge Sg_frac and Sg_topo (and in the future Flgg_hflx)
          call t_startf('MED:'//trim(subname)//' glc2lnd')
          call map_glc2lnd(gcomp, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call t_stopf('MED:'//trim(subname)//' glc2lnd')
       end if

       !---------------------------------------
       ! update scalar data
       !---------------------------------------

       call ESMF_StateGet(is_local%wrap%NStateImp(compatm), trim(is_local%wrap%flds_scalar_name), itemType, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (itemType /= ESMF_STATEITEM_NOTFOUND) then
          call t_startf('MED:'//trim(subname)//' nextsw_cday')
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
          call t_stopf('MED:'//trim(subname)//' nextsw_cday')
       end if

       ! diagnose
       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(complnd), &
               string=trim(subname)//' FBexp(complnd) ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Reset first call logical
    first_call = .false.

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_lnd

  !================================================================================================

  subroutine map_glc2lnd_init(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState)       :: is_local
    type(ESMF_Field)          :: lfield_l
    type(ESMF_Mesh)           :: mesh_l
    integer                   :: ungriddedUBound_output(1)
    integer                   :: fieldCount
    integer                   :: ns,n
    type(ESMF_Field), pointer :: fieldlist(:) => null()
    character(len=*) , parameter   :: subname='(map_glc2lnd_mod:map_glc2lnd_init)'
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

    ! Determine number of elevation classes by querying a field that has elevation classes in it
    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), fieldname='Sg_topo_elev', field=lfield_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield_l, ungriddedUBound=ungriddedUBound_output, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ungriddedCount = ungriddedUBound_output(1)
    ! TODO: check that ungriddedCount = glc_nec+1

    ! -------------------------------
    ! Create module fields on land mesh
    ! -------------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), fieldlist=fieldlist, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(fieldlist(1), mesh=mesh_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    deallocate(fieldlist)

    field_icemask_l = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    field_frac_l_ec = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    field_frac_x_icemask_l_ec = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    field_topo_x_icemask_l_ec = ESMF_FieldCreate(mesh_l, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! create module fields on glc mesh
    !---------------------------------------

    do ns = 1,max_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then

          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc(ns),compglc(ns)), &
               fieldCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldlist(fieldcount))
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc(ns),compglc(ns)), &
               fieldlist=fieldlist, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(fieldlist(1), mesh=ice_sheet_tolnd(ns)%mesh_g, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          deallocate(fieldlist)

          ice_sheet_tolnd(ns)%field_icemask_g = ESMF_FieldCreate(ice_sheet_tolnd(ns)%mesh_g, &
               ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ice_sheet_tolnd(ns)%field_frac_g_ec = ESMF_FieldCreate(ice_sheet_tolnd(ns)%mesh_g, &
               ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ice_sheet_tolnd(ns)%field_frac_x_icemask_g_ec = ESMF_FieldCreate(ice_sheet_tolnd(ns)%mesh_g, &
               ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ice_sheet_tolnd(ns)%field_topo_x_icemask_g_ec = ESMF_FieldCreate(ice_sheet_tolnd(ns)%mesh_g, &
               ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
               ungriddedLbound=(/1/), ungriddedUbound=(/ungriddedCount/), gridToFieldMap=(/2/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Create route handle if it has not been created
          if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compglc(ns),complnd,mapconsd), rc=rc)) then
             call med_map_routehandles_init( compglc(ns), complnd, &
                  ice_sheet_tolnd(ns)%field_icemask_g, lfield_l, &
                  mapindex=mapconsd, &
                  routehandles=is_local%wrap%rh(compglc(ns),complnd,:), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

    ! Currently cannot map hflx in multiple elevation classes from glc to land
    if (FB_fldchk(is_local%wrap%FBExp(complnd), trim(Flgg_hflx), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//'ERROR: Flgg_hflx to land has not been implemented yet', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       rc = ESMF_FAILURE
       return
    end if

  end subroutine map_glc2lnd_init

  !================================================================================================
  subroutine map_glc2lnd( gcomp, rc)

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
    type(ESMF_Field)      :: lfield
    type(ESMF_Field)      :: lfield_src
    type(ESMF_Field)      :: lfield_dst
    integer               :: ec, l, g, ns, n
    real(r8)              :: topo_virtual
    real(r8), pointer     :: icemask_g(:) => null()             ! glc ice mask field on glc grid
    real(r8), pointer     :: frac_g(:) => null()                ! total ice fraction in each glc cell
    real(r8), pointer     :: frac_g_ec(:,:) => null()           ! glc fractions on the glc grid
    real(r8), pointer     :: frac_l_ec(:,:) => null()           ! glc fractions on the land grid
    real(r8), pointer     :: topo_g(:) => null()                ! topo height of each glc cell (no elev classes)
    real(r8), pointer     :: topo_l_ec(:,:) => null()           ! topo height in each land gridcell for each elev class
    real(r8), pointer     :: frac_x_icemask_g_ec(:,:) => null() ! (glc fraction) x (icemask), on the glc grid
    real(r8), pointer     :: frac_x_icemask_l_ec(:,:) => null()
    real(r8), pointer     :: topo_x_icemask_g_ec(:,:) => null()
    real(r8), pointer     :: topo_x_icemask_l_ec(:,:) => null()
    real(r8), pointer     :: dataptr1d(:) => null()
    real(r8), pointer     :: dataptr2d(:,:) => null()
    real(r8), pointer     :: frac_l_ec_sum(:,:) => null()
    real(r8), pointer     :: topo_l_ec_sum(:,:) => null()
    real(r8), pointer     :: dataptr1d_src(:) => null()
    real(r8), pointer     :: dataptr1d_dst(:) => null()
    character(len=*), parameter :: subname = 'map_glc2lnd'
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
    ! Get pointers into land export field bundle (this is summed over all ice sheets)
    !---------------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), fieldname=trim(Sg_frac)//'_elev', &
         field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_fieldGet(lfield, farrayptr=frac_l_ec_sum, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    frac_l_ec_sum(:,:) = 0._r8

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), fieldname=trim(Sg_topo)//'_elev', &
         field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayptr=topo_l_ec_sum, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    topo_l_ec_sum(:,:) = 0._r8

    !---------------------------------
    ! Map fractional ice coverage to the land grid (multiple elevation classes)
    !---------------------------------

    ! Map Sg_icemask and Sg_icemask_coupled_fluxes (no elevation classes)
    do ns = 1,num_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
          call t_startf('MED:'//trim(subname)//' glc2lnd ')
          call med_map_field_packed( &
               FBSrc=is_local%wrap%FBImp(compglc(ns),compglc(ns)), &
               FBDst=is_local%wrap%FBImp(compglc(ns),complnd), &
               FBFracSrc=is_local%wrap%FBFrac(compglc(ns)), &
               field_normOne=is_local%wrap%field_normOne(compglc(ns),complnd,:), &
               packed_data=is_local%wrap%packed_data(compglc(ns),complnd,:), &
               routehandles=is_local%wrap%RH(compglc(ns),complnd,:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call t_stopf('MED:'//trim(subname)//' glc2lnd')
       end if
    end do

    ! Get Sg_icemask on land as sum of all ice sheets (no elevation classes)
    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), Sg_icemask, &
         field=lfield_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield_dst, farrayPtr=dataptr1d_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr1d_dst(:) = 0._r8
    do ns = 1,num_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc(ns),complnd), Sg_icemask, &
               field=lfield_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield_src, farrayPtr=dataptr1d_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d_dst(:) = dataptr1d_dst(:) + dataptr1d_src(:) 
       end if
    end do

    ! Get Sg_icemask_coupled_fluxes on land as sum of all ice sheets (no elevation classes)
    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), Sg_icemask_coupled_fluxes, &
         field=lfield_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield_dst, farrayPtr=dataptr1d_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr1d_dst(:) = 0._r8
    do ns = 1,num_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc(ns),complnd), Sg_icemask_coupled_fluxes, &
               field=lfield_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield_src, farrayPtr=dataptr1d_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d_dst(:) = dataptr1d_dst(:) + dataptr1d_src(:) 
       end if
    end do

    do ns = 1,num_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
          ! Set contents of FBglc_ec to contain frac_g_ec
          ! (fractional ice coverage for each elevation class on the glc grid)

          ! set FBglc_ec on the glc grid (fractional ice coverage per elevation class)
          ! topo_g(:) is the topographic height of each glc gridcell
          ! frac_g(:) is the total ice fraction in each glc gridcell
          ! frac_g_ec(:,:) are the glc fractions on the glc grid for each elevation class (inner dimension)
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc(ns),compglc(ns)), fieldname=trim(Sg_topo), &
               field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayptr=topo_g, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! compute frac_g_ec
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc(ns),compglc(ns)), fieldname=trim(Sg_frac), &
               field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayptr=frac_g, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(ice_sheet_tolnd(ns)%field_frac_g_ec, farrayptr=frac_g_ec, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call glc_get_fractional_icecov(ungriddedCount-1, topo_g, frac_g, frac_g_ec, logunit)

          ! compute icemask_g
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc(ns),compglc(ns)), fieldname=trim(Sg_icemask), &
               field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(ice_sheet_tolnd(ns)%field_icemask_g, farrayptr=icemask_g,  rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          icemask_g(:) = dataptr1d(:)

          ! compute frac_x_icemask_g_ec
          ! only include grid cells that are both (a) within the icemask and (b) in this elevation class
          call ESMF_FieldGet(ice_sheet_tolnd(ns)%field_frac_x_icemask_g_ec, farrayptr=frac_x_icemask_g_ec, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do ec = 1, ungriddedCount
             frac_x_icemask_g_ec(ec,:) = frac_g_ec(ec,:) * icemask_g(:)
          end do

          ! map frac_g_ec to frac_l_ec and normalize by icemask_g
          if (dbug_flag > 1) then
             call ESMF_LogWrite(trim(subname)//": calling mapping elevation class fractions from glc to land", &
                  ESMF_LOGMSG_INFO)
          end if
          call med_map_field_normalized(  &
               field_src=ice_sheet_tolnd(ns)%field_frac_g_ec, &
               field_dst=field_frac_l_ec, &
               routehandles=is_local%wrap%RH(compglc(ns),complnd,:), &
               maptype=mapconsd, &
               field_normsrc=ice_sheet_tolnd(ns)%field_icemask_g, &
               field_normdst=field_icemask_l, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! now set values in land export state for Sg_frac_elev (this is summed over all ice sheets)
          call ESMF_fieldGet(field_frac_l_ec, farrayptr=frac_l_ec, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          frac_l_ec_sum(:,:) = frac_l_ec_sum(:,:) + frac_l_ec(:,:)

          !---------------------------------
          ! Map topo to the land grid (multiple elevation classes)
          !---------------------------------

          ! Note that all topo values in FBimp(compglc(ns),compglc(ns)) do not have elevation class dependence
          ! Normalize by frac_x_icemask_g_ec - this is what introduces
          ! elevation class information from the glc grid (without elevation classes) to the
          ! land grid (with elevation classes)
          ! Note that bare land values are mapped in the same way as ice-covered values

          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compglc(ns),compglc(ns)), fieldname=trim(Sg_topo), &
               field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_fieldGet(lfield, farrayptr=topo_g, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(ice_sheet_tolnd(ns)%field_topo_x_icemask_g_ec, farrayptr=topo_x_icemask_g_ec, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do ec = 1,ungriddedCount
             do l = 1,size(topo_g)
                topo_x_icemask_g_ec(ec,l) = topo_g(l) * frac_x_icemask_g_ec(ec,l)
             end do
          end do

          ! map field_topo_x_icemask_g_ec from glc to land (with multiple elevation classes) - no normalization
          if (dbug_flag > 1) then
             call ESMF_LogWrite(trim(subname)//": calling mapping of topo from glc to land", ESMF_LOGMSG_INFO)
          end if
          call med_map_field(  &
               field_src=ice_sheet_tolnd(ns)%field_topo_x_icemask_g_ec, &
               field_dst=field_topo_x_icemask_l_ec, &
               routehandles=is_local%wrap%RH(compglc(ns),complnd,:), &
               maptype=mapconsd, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field_topo_x_icemask_l_ec, farrayptr=topo_l_ec, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! map FBglc_frac_x_icemask from glc to land (with multiple elevation classes) - no normalization
          if (dbug_flag > 1) then
             call ESMF_LogWrite(trim(subname)//": calling mapping of frac_x_icemask from glc to land", ESMF_LOGMSG_INFO)
          end if
          call med_map_field(  &
               field_src=ice_sheet_tolnd(ns)%field_frac_x_icemask_g_ec, &
               field_dst=field_frac_x_icemask_l_ec, &
               routehandles=is_local%wrap%RH(compglc(ns),complnd,:), &
               maptype=mapconsd, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field_frac_x_icemask_l_ec, farrayptr=frac_x_icemask_l_ec, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! set Sg_topo values in export state to land (in multiple elevation classes)
          ! also set the topo field for virtual columns, in a given elevation class.
          ! This is needed because virtual columns (i.e., elevation classes that have no
          ! contributing glc grid cells) won't have any topographic information mapped onto
          ! them, so would otherwise end up with an elevation of 0.
          do ec = 1,ungriddedCount
             topo_virtual = glc_mean_elevation_virtual(ec-1) ! glc_mean_elevation_virtual uses 0:glc_nec
             do l = 1,size(frac_x_icemask_l_ec, dim=2)
                if (frac_l_ec_sum(ec,l) <= 0._r8) then
                   topo_l_ec_sum(ec,l) = topo_l_ec_sum(ec,l) + topo_virtual
                else
                   if (frac_x_icemask_l_ec(ec,l) /= 0.0_r8) then
                      topo_l_ec_sum(ec,l) = topo_l_ec_sum(ec,l) + topo_l_ec(ec,l) / frac_x_icemask_l_ec(ec,l)
                   end if
                end if
             end do
          end do
       end if
    end do

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine map_glc2lnd

end module med_phases_prep_lnd_mod
