module med_phases_post_glc_mod

  !-----------------------------------------------------------------------------
  ! Mediator phase for mapping glc->lnd and glc->ocn after the receive of glc
  ! ASSUMES that multiple ice sheets do not overlap
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_StateGet, ESMF_StateItem_Flag
  use ESMF                  , only : ESMF_Mesh, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF                  , only : ESMF_RouteHandle, ESMF_RouteHandleIsCreated
  use med_internalstate_mod , only : compatm, compice, complnd, comprof, compocn, compname, compglc
  use med_internalstate_mod , only : mapbilnr, mapconsd, compname
  use med_internalstate_mod , only : InternalState, maintask, logunit
  use med_methods_mod       , only : fldbun_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : fldbun_fldchk    => med_methods_FB_fldchk
  use med_methods_mod       , only : fldbun_getmesh   => med_methods_FB_getmesh
  use med_methods_mod       , only : fldbun_getdata1d => med_methods_FB_getdata1d
  use med_methods_mod       , only : fldbun_getdata2d => med_methods_FB_getdata2d
  use med_methods_mod       , only : field_getdata1d  => med_methods_Field_getdata1d
  use med_methods_mod       , only : field_getdata2d  => med_methods_Field_getdata2d
  use med_utils_mod         , only : chkerr           => med_utils_ChkErr
  use med_constants_mod     , only : dbug_flag        => med_constants_dbug_flag
  use med_map_mod           , only : med_map_rh_is_created, med_map_routehandles_init
  use med_map_mod           , only : med_map_field_packed, med_map_field_normalized, med_map_field
  use glc_elevclass_mod     , only : glc_mean_elevation_virtual, glc_get_fractional_icecov
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_post_glc

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
  type(ice_sheet_tolnd_type), allocatable :: ice_sheet_tolnd(:)

  type(ESMF_field) :: field_icemask_l                ! no elevation classes
  type(ESMF_Field) :: field_frac_l_ec                ! elevation classes
  type(ESMF_Field) :: field_frac_x_icemask_l_ec      ! elevation classes
  type(ESMF_Field) :: field_topo_x_icemask_l_ec      ! elevation classes

  ! the number of elevation classes (excluding bare land) = ungriddedCount - 1
  integer :: ungriddedCount ! this equals the number of elevation classes + 1 (for bare land)

  logical :: cism_evolve = .false.
  logical :: glc2lnd_coupling = .false.
  logical :: glc2rof_coupling = .false.
  logical :: glc2ice_coupling = .false.

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!================================================================================================
contains
!================================================================================================

  subroutine med_phases_post_glc(gcomp, rc)

    use NUOPC_Mediator        , only : NUOPC_MediatorGet
    use ESMF                  , only : ESMF_Clock, ESMF_ClockIsCreated
    use med_phases_history_mod, only : med_phases_history_write_comp

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)          :: dClock
    type(InternalState)       :: is_local
    integer                   :: ns
    logical                   :: first_call = .true.
    logical                   :: isPresent
    character(CL)             :: cvalue
    character(len=*), parameter :: subname='(med_phases_post_glc)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (first_call) then
       ! determine if there will be any glc to lnd coupling
       do ns = 1,is_local%wrap%num_icesheets
          if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
             glc2lnd_coupling = .true.
             exit
          end if
       end do
       ! determine if there will be any glc to ocn coupling
       do ns = 1,is_local%wrap%num_icesheets
          if (is_local%wrap%med_coupling_active(compglc(ns),comprof)) then
             glc2rof_coupling = .true.
             exit
          end if
       end do
       ! determine if there will be any glc to ice coupling
       do ns = 1,is_local%wrap%num_icesheets
          if (is_local%wrap%med_coupling_active(compglc(ns),compice)) then
             glc2ice_coupling = .true.
             exit
          end if
       end do
       if (maintask) then
          write(logunit,'(a,L1)') trim(subname) // 'glc2lnd_coupling is ',glc2lnd_coupling
          write(logunit,'(a,L1)') trim(subname) // 'glc2rof_coupling is ',glc2rof_coupling
          write(logunit,'(a,L1)') trim(subname) // 'glc2ice_coupling is ',glc2ice_coupling
       end if

       ! determine if coupling to CISM is 2-way
       call NUOPC_CompAttributeGet(gcomp, name="cism_evolve", isPresent=isPresent, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (isPresent) then
          call NUOPC_CompAttributeGet(gcomp, name="cism_evolve", value=cvalue, isPresent=isPresent, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          read (cvalue,*) cism_evolve
          if (maintask) then
             write(logunit,'(a,l7)') trim(subname)//' cism_evolve = ',cism_evolve
          end if
       end if
    end if

    !---------------------------------------
    ! glc->rof mapping
    !---------------------------------------

    if (glc2rof_coupling) then
       do ns = 1,is_local%wrap%num_icesheets
          if (is_local%wrap%med_coupling_active(compglc(ns),comprof)) then
             call med_map_field_packed( &
                  FBSrc=is_local%wrap%FBImp(compglc(ns),compglc(ns)), &
                  FBDst=is_local%wrap%FBImp(compglc(ns),comprof), &
                  FBFracSrc=is_local%wrap%FBFrac(compglc(ns)), &
                  field_normOne=is_local%wrap%field_normOne(compglc(ns),comprof,:), &
                  packed_data=is_local%wrap%packed_data(compglc(ns),comprof,:), &
                  routehandles=is_local%wrap%RH(compglc(ns),comprof,:), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
    end if

    !---------------------------------------
    ! glc->ice mapping
    !---------------------------------------
    if (glc2ice_coupling) then
       ! Fill this in
    end if

    !---------------------------------------
    ! glc->lnd mapping and custom merging of all ice sheets onto land mesh
    !---------------------------------------
    if (glc2lnd_coupling) then
       ! The will following will map and merge Sg_frac and Sg_topo (and in the future Flgg_hflx)
       call t_startf('MED:'//trim(subname)//' glc2lnd ')
       do ns = 1,is_local%wrap%num_icesheets
          if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
             call med_map_field_packed( &
                  FBSrc=is_local%wrap%FBImp(compglc(ns),compglc(ns)), &
                  FBDst=is_local%wrap%FBImp(compglc(ns),complnd), &
                  FBFracSrc=is_local%wrap%FBFrac(compglc(ns)), &
                  field_normOne=is_local%wrap%field_normOne(compglc(ns),complnd,:), &
                  packed_data=is_local%wrap%packed_data(compglc(ns),complnd,:), &
                  routehandles=is_local%wrap%RH(compglc(ns),complnd,:), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
       call t_stopf('MED:'//trim(subname)//' glc2lnd')

       ! The following is only done if glc->lnd coupling is active
       if (first_call) then
          call map_glc2lnd_init(gcomp, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! The will following will map and merge Sg_frac and Sg_topo (and in the future Flgg_hflx)
       call map_glc2lnd(gcomp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Reset first call logical
    first_call = .false.

    ! Write glc inst, avg or aux if requested in mediator attributes
    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_ClockIsCreated(dclock)) then
       do ns = 1,is_local%wrap%num_icesheets
          call med_phases_history_write_comp(gcomp, compglc(ns), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end do
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_glc

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
    integer                   :: ns
    character(len=*) , parameter   :: subname='(map_glc2lnd_init)'
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
    call ESMF_FieldBundleGet(is_local%wrap%FBExp(complnd), 'Sg_topo_elev', field=lfield_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield_l, ungriddedUBound=ungriddedUBound_output, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ungriddedCount = ungriddedUBound_output(1)
    ! TODO: check that ungriddedCount = glc_nec+1

    ! -------------------------------
    ! Create module fields on land mesh
    ! -------------------------------

    call fldbun_getmesh(is_local%wrap%FBExp(complnd), mesh_l, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

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

    ! allocate module variable
    allocate(ice_sheet_tolnd(is_local%wrap%num_icesheets))

    do ns = 1,is_local%wrap%num_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then

          call fldbun_getmesh(is_local%wrap%FBImp(compglc(ns),compglc(ns)), ice_sheet_tolnd(ns)%mesh_g, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
                  ice_sheet_tolnd(ns)%field_icemask_g, field_icemask_l, &
                  mapindex=mapconsd, &
                  routehandles=is_local%wrap%rh(compglc(ns),complnd,:), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

    ! Currently cannot map hflx in multiple elevation classes from glc to land
    if (fldbun_fldchk(is_local%wrap%FBExp(complnd), trim(Flgg_hflx), rc=rc)) then
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
    integer               :: ec, l, ns
    real(r8)              :: topo_virtual
    real(r8), pointer     :: icemask_g(:)              ! glc ice mask field on glc grid
    real(r8), pointer     :: frac_g(:)                 ! total ice fraction in each glc cell
    real(r8), pointer     :: frac_g_ec(:,:)            ! glc fractions on the glc grid
    real(r8), pointer     :: frac_l_ec(:,:)            ! glc fractions on the land grid
    real(r8), pointer     :: topo_g(:)                 ! topo height of each glc cell (no elev classes)
    real(r8), pointer     :: topo_l_ec(:,:)            ! topo height in each land gridcell for each elev class
    real(r8), pointer     :: frac_x_icemask_g_ec(:,:)  ! (glc fraction) x (icemask), on the glc grid
    real(r8), pointer     :: frac_x_icemask_l_ec(:,:)
    real(r8), pointer     :: topo_x_icemask_g_ec(:,:)
    real(r8), pointer     :: dataptr1d(:)
    real(r8), pointer     :: frac_l_ec_sum(:,:)
    real(r8), pointer     :: topo_l_ec_sum(:,:)
    real(r8), pointer     :: dataptr1d_src(:)
    real(r8), pointer     :: dataptr1d_dst(:)
    real(r8), pointer     :: icemask_l(:)
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

    call fldbun_getdata2d(is_local%wrap%FBExp(complnd), trim(Sg_frac)//'_elev', frac_l_ec_sum, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    frac_l_ec_sum(:,:) = 0._r8

    call fldbun_getdata2d(is_local%wrap%FBExp(complnd), trim(Sg_topo)//'_elev', topo_l_ec_sum, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    topo_l_ec_sum(:,:) = 0._r8

    !---------------------------------
    ! Map fractional ice coverage to the land grid (multiple elevation classes)
    !---------------------------------

    ! Map Sg_icemask and Sg_icemask_coupled_fluxes (no elevation classes)
    do ns = 1,is_local%wrap%num_icesheets
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
    call fldbun_getdata1d(is_local%wrap%FBExp(complnd), Sg_icemask, dataptr1d_dst, rc)
    dataptr1d_dst(:) = 0._r8
    do ns = 1,is_local%wrap%num_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
          call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),complnd), Sg_icemask, dataptr1d_src, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d_dst(:) = dataptr1d_dst(:) + dataptr1d_src(:)
       end if
    end do

    ! Get Sg_icemask_coupled_fluxes on land as sum of all ice sheets (no elevation classes)
    call fldbun_getdata1d(is_local%wrap%FBExp(complnd), Sg_icemask_coupled_fluxes, dataptr1d_dst, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr1d_dst(:) = 0._r8
    do ns = 1,is_local%wrap%num_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then
          call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),complnd), Sg_icemask_coupled_fluxes, dataptr1d_src, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d_dst(:) = dataptr1d_dst(:) + dataptr1d_src(:)
       end if
    end do

    do ns = 1,is_local%wrap%num_icesheets
       if (is_local%wrap%med_coupling_active(compglc(ns),complnd)) then

          ! Set (fractional ice coverage for each elevation class on the glc grid)

          ! get topo_g(:) - the topographic height of each glc gridcell
          call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_topo, topo_g, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! get frac_g(:) - the total ice fraction in each glc gridcell
          call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_frac, frac_g, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! compute frac_g_ec(:,:) - the glc fractions on the glc grid for each elevation class (inner dimension)
          call field_getdata2d(ice_sheet_tolnd(ns)%field_frac_g_ec, frac_g_ec, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call glc_get_fractional_icecov(ungriddedCount-1, topo_g, frac_g, frac_g_ec, logunit)

          ! compute icemask_g
          call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_icemask, dataptr1d, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call field_getdata1d(ice_sheet_tolnd(ns)%field_icemask_g, icemask_g,  rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          icemask_g(:) = dataptr1d(:)

          ! compute frac_x_icemask_g_ec
          ! only include grid cells that are both (a) within the icemask and (b) in this elevation class
          call field_getdata2d(ice_sheet_tolnd(ns)%field_frac_x_icemask_g_ec, frac_x_icemask_g_ec, rc=rc)
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
          call field_getdata2d(field_frac_l_ec, frac_l_ec, rc=rc)
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

          call fldbun_getdata1d(is_local%wrap%FBImp(compglc(ns),compglc(ns)), Sg_topo, topo_g, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call field_getdata2d(ice_sheet_tolnd(ns)%field_topo_x_icemask_g_ec, topo_x_icemask_g_ec, rc=rc)
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
          call field_getdata2d(field_topo_x_icemask_l_ec, topo_l_ec, rc=rc)
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
          call field_getdata2d(field_frac_x_icemask_l_ec, frac_x_icemask_l_ec, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call field_getdata1d(field_icemask_l, icemask_l, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! set Sg_topo values in export state to land (in multiple elevation classes)
          ! also set the topo field for virtual columns, in a given elevation class.
          ! This is needed because virtual columns (i.e., elevation classes that have no
          ! contributing glc grid cells) won't have any topographic information mapped onto
          ! them, so would otherwise end up with an elevation of 0.
          ! ASSUME that multiple ice sheets do not overlap
          do ec = 1,ungriddedCount
             topo_virtual = glc_mean_elevation_virtual(ec-1) ! glc_mean_elevation_virtual uses 0:glc_nec
             do l = 1,size(frac_x_icemask_l_ec, dim=2)
                if (icemask_l(l) > 0._r8) then
                   ! We only do this where icemask_l > 0 to avoid adding topo_virtual
                   ! multiple times. If icemask_l == 0 for all ice sheets, then lnd should
                   ! ignore the topo values from glc, so it's safe to leave them unset; if
                   ! icemask_l is 0 for this ice sheet but > 0 for some other ice sheet,
                   ! then we'll get the appropriate topo setting from that other ice
                   ! sheet.
                   !
                   ! Note that frac_l_ec_sum is the sum over ice sheets we have handled so
                   ! far in the outer loop over ice sheets. At first glance, that could
                   ! seem wrong (because what if a later ice sheet causes this sum to
                   ! become greater than 0?), and it may be that we should rework this for
                   ! clarity. However, since icemask_l > 0 (which is the ice mask for this
                   ! ice sheet) and we assume that multiple ice sheets do not overlap, we
                   ! can be confident that no other ice sheet will contribute to
                   ! frac_l_ec_sum for this land point, so if it is <= 0 at this point,
                   ! it should remain <= 0.
                   if (frac_l_ec_sum(ec,l) <= 0._r8) then
                      ! This is formulated as an addition for consistency with other
                      ! additions to the *_sum variables, but in practice only one ice
                      ! sheet will contribute to any land point, given the assumption of
                      ! non-overlapping ice sheet domains. (If more than one ice sheet
                      ! contributed to a given land point, the following line would do the
                      ! wrong thing, since it would add topo_virtual multiple times.)
                      topo_l_ec_sum(ec,l) = topo_l_ec_sum(ec,l) + topo_virtual
                   else
                      if (frac_x_icemask_l_ec(ec,l) /= 0.0_r8) then
                         topo_l_ec_sum(ec,l) = topo_l_ec_sum(ec,l) + topo_l_ec(ec,l) / frac_x_icemask_l_ec(ec,l)
                      end if
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

end module med_phases_post_glc_mod
