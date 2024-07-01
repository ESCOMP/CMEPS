module med_phases_aofluxes_mod

  ! --------------------------------------------------------------------------
  ! Determine atm/ocn flux calculation in mediator - for one of 3 cases:
  ! if aoflux grid is ocn
  !  - map atm attributes of aoflux_in to ocn and map aoflux_out back to atm
  ! if aoflux grid is atm
  !  - map ocn attributes of oaflux_in to atm and map aoflux_out back to ocn
  ! if aoflux grid is exchange
  !  - map both atm and ocn attributes of aoflux_in to xgrid and then
  !    map aoflux_out from xgrid to both atm and ocn grid
  ! --------------------------------------------------------------------------

  use ESMF                  , only : operator(/=)
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_CoordSys_Flag
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate, ESMF_FieldIsCreated, ESMF_FieldDestroy
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldRegridGetArea
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_RouteHandle, ESMF_FieldRegrid, ESMF_FieldRegridStore
  use ESMF                  , only : ESMF_REGRIDMETHOD_CONSERVE_2ND, ESMF_REGRIDMETHOD_CONSERVE
  use ESMF                  , only : ESMF_REGRIDMETHOD_PATCH, ESMF_REGRIDMETHOD_BILINEAR, ESMF_COORDSYS_CART
  use ESMF                  , only : ESMF_TERMORDER_SRCSEQ, ESMF_REGION_TOTAL, ESMF_MESHLOC_ELEMENT, ESMF_MAXSTR
  use ESMF                  , only : ESMF_XGRIDSIDE_B, ESMF_XGRIDSIDE_A, ESMF_END_ABORT, ESMF_LOGERR_PASSTHRU
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet, ESMF_XGrid, ESMF_XGridCreate, ESMF_TYPEKIND_R8
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LOGMSG_ERROR, ESMF_FAILURE
  use ESMF                  , only : ESMF_Finalize, ESMF_LogFoundError
  use ESMF                  , only : ESMF_XGridGet, ESMF_MeshCreate, ESMF_MeshWrite, ESMF_KIND_R8
  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : InternalState, maintask, logunit
  use med_internalstate_mod , only : compatm, compocn, compwav, coupling_mode, aoflux_code, mapconsd, mapconsf, mapfcopy
  use med_constants_mod     , only : dbug_flag    => med_constants_dbug_flag
  use med_utils_mod         , only : memcheck     => med_memcheck
  use med_utils_mod         , only : chkerr       => med_utils_chkerr
  use perf_mod              , only : t_startf, t_stopf
#ifndef CESMCOUPLED
  use ufs_const_mod         , only : rearth => SHR_CONST_REARTH
  use ufs_const_mod         , only : pi => SHR_CONST_PI
#else
  use shr_const_mod         , only : rearth => SHR_CONST_REARTH
  use shr_const_mod         , only : pi => SHR_CONST_PI
#endif

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public routines
  !--------------------------------------------------------------------------

  public :: med_phases_aofluxes_init_fldbuns
  public :: med_phases_aofluxes_run
  public :: med_aofluxes_map_ogrid2agrid_output
  public :: med_aofluxes_map_xgrid2agrid_output
  public :: med_aofluxes_map_xgrid2ogrid_output
  public :: med_aofluxes_map_agrid2ogrid_output

  !--------------------------------------------------------------------------
  ! Private routines
  !--------------------------------------------------------------------------

  private :: med_aofluxes_init
  private :: med_aofluxes_init_ogrid
  private :: med_aofluxes_init_agrid
  private :: med_aofluxes_init_xgrid
  private :: med_aofluxes_map_ogrid2xgrid_input
  private :: med_aofluxes_map_agrid2xgrid_input
  private :: med_aofluxes_map_ogrid2agrid_input
  private :: med_aofluxes_update
  private :: set_aoflux_in_pointers
  private :: set_aoflux_out_pointers
  private :: fldbun_getfldptr

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  logical :: flds_wiso  ! use case
  logical :: compute_atm_dens
  logical :: compute_atm_thbot
  integer :: ocn_surface_flux_scheme ! use case
  logical :: add_gusts

  character(len=CS), pointer :: fldnames_ocn_in(:)
  character(len=CS), pointer :: fldnames_atm_in(:)
  character(len=CS), pointer :: fldnames_aof_out(:)

  ! following is needed for atm/ocn fluxes on atm grid
  type(ESMF_FieldBundle) :: FBocn_a ! ocean fields need for aoflux calc on atm grid

  ! following is needed for atm/ocn fluxes on the exchange grid
  type(ESMF_FieldBundle) :: FBocn_x               ! input ocn fields
  type(ESMF_FieldBundle) :: FBatm_x               ! input atm fields
  type(ESMF_FieldBundle) :: FBaof_x               ! output aoflux fields
  type(ESMF_RouteHandle) :: rh_ogrid2xgrid        ! ocn->xgrid mapping
  type(ESMF_RouteHandle) :: rh_agrid2xgrid        ! atm->xgrid mapping
  type(ESMF_RouteHandle) :: rh_xgrid2ogrid        ! xgrid->ocn mapping
  type(ESMF_RouteHandle) :: rh_xgrid2agrid        ! xgrid->atm mapping
  type(ESMF_RouteHandle) :: rh_agrid2xgrid_2ndord ! atm->xgrid mapping 2nd order conservative
  type(ESMF_RouteHandle) :: rh_agrid2xgrid_bilinr ! atm->xgrid mapping bilinear
  type(ESMF_RouteHandle) :: rh_agrid2xgrid_patch  ! atm->xgrid mapping patch
  type(ESMF_XGrid)       :: xgrid
  type(ESMF_Field)       :: field_o
  type(ESMF_Field)       :: field_x

  type aoflux_in_type
     ! input: ocn
     real(R8) , pointer :: uocn        (:) => null() ! ocn velocity, zonal
     real(R8) , pointer :: vocn        (:) => null() ! ocn velocity, meridional
     real(R8) , pointer :: tocn        (:) => null() ! ocean temperature
     real(R8) , pointer :: roce_16O    (:) => null() ! ocn H2O ratio
     real(R8) , pointer :: roce_HDO    (:) => null() ! ocn HDO ratio
     real(R8) , pointer :: roce_18O    (:) => null() ! ocn H218O ratio
     ! input: atm
     real(R8) , pointer :: zbot        (:) => null() ! atm level height
     real(R8) , pointer :: ubot        (:) => null() ! atm velocity, zonal
     real(R8) , pointer :: vbot        (:) => null() ! atm velocity, meridional
     real(R8) , pointer :: usfc        (:) => null() ! atm surface velocity, zonal
     real(R8) , pointer :: vsfc        (:) => null() ! atm surface velocity, meridional
     real(R8) , pointer :: thbot       (:) => null() ! atm potential T
     real(R8) , pointer :: shum        (:) => null() ! atm specific humidity
     real(R8) , pointer :: pbot        (:) => null() ! atm bottom pressure
     real(R8) , pointer :: psfc        (:) => null() ! atm surface pressure
     real(R8) , pointer :: dens        (:) => null() ! atm bottom density
     real(R8) , pointer :: tbot        (:) => null() ! atm bottom surface T
     real(R8) , pointer :: shum_16O    (:) => null() ! atm H2O tracer
     real(R8) , pointer :: shum_HDO    (:) => null() ! atm HDO tracer
     real(R8) , pointer :: shum_18O    (:) => null() ! atm H218O tracer
     real(R8) , pointer :: lwdn        (:) => null() ! atm downward longwave heat flux
     real(R8) , pointer :: rainc       (:) => null() ! convective rain flux
     ! local size and computational mask and area: on aoflux grid
     integer            :: lsize                     ! local size
     integer  , pointer :: mask        (:) => null() ! integer ocn domain mask: 0 <=> inactive cell
     real(R8) , pointer :: rmask       (:) => null() ! real    ocn domain mask: 0 <=> inactive cell
     real(R8) , pointer :: garea       (:) => null() ! atm grid area
  end type aoflux_in_type

  type aoflux_out_type
     real(R8) , pointer :: sen         (:) => null() ! heat flux: sensible
     real(R8) , pointer :: lat         (:) => null() ! heat flux: latent
     real(R8) , pointer :: lwup        (:) => null() ! lwup over ocean
     real(R8) , pointer :: evap        (:) => null() ! water flux: evaporation
     real(R8) , pointer :: evap_16O    (:) => null() ! H2O flux: evaporation
     real(R8) , pointer :: evap_HDO    (:) => null() ! HDO flux: evaporation
     real(R8) , pointer :: evap_18O    (:) => null() ! H218O flux: evaporation
     real(R8) , pointer :: taux        (:) => null() ! wind stress, zonal
     real(R8) , pointer :: tauy        (:) => null() ! wind stress, meridional
     real(R8) , pointer :: tref        (:) => null() ! diagnostic: 2m ref T
     real(R8) , pointer :: qref        (:) => null() ! diagnostic: 2m ref Q
     real(R8) , pointer :: u10         (:) => null() ! diagnostic: 10m wind speed
     real(R8) , pointer :: duu10n      (:) => null() ! diagnostic: 10m wind speed squared
     real(R8) , pointer :: ugust_out   (:) => null() ! diagnostic: gust wind added
     real(R8) , pointer :: u10_withGust(:) => null() ! diagnostic: gust wind added
     real(R8) , pointer :: u10res      (:) => null() ! diagnostic: no gust wind added
     real(R8) , pointer :: ustar       (:) => null() ! saved ustar
     real(R8) , pointer :: re          (:) => null() ! saved re
     real(R8) , pointer :: ssq         (:) => null() ! saved sq
  end type aoflux_out_type

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine med_phases_aofluxes_init_fldbuns(gcomp, rc)

    use ESMF            , only : ESMF_FieldBundleIsCreated
    use esmFlds         , only : med_fldList_GetNumFlds
    use esmFlds         , only : med_fldList_GetFldNames
    use esmFlds         , only : med_fldList_GetaofluxfldList
    use esmFlds         , only : med_fldList_type
    use med_methods_mod , only : FB_init => med_methods_FB_init
    use med_internalstate_mod, only : compname

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: n
    integer             :: fieldcount
    type(med_fldList_type), pointer :: fldListMed_aoflux
    type(InternalState) :: is_local
    character(len=*),parameter :: subname=' (med_phases_aofluxes_init_fldbuns) '
    !---------------------------------------

    ! Create field bundles for mediator ocean/atmosphere flux computation
    ! This is needed regardless of the grid on which the atm/ocn flux computation is done on
    fldListMed_aoflux => med_fldList_GetaofluxFldList()
    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Set module variable fldnames_aof_out
    fieldCount = med_fldList_GetNumFlds(fldListMed_aoflux)
    allocate(fldnames_aof_out(fieldCount))
    call med_fldList_getfldnames(fldListMed_aoflux%fields, fldnames_aof_out, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize FBMed_aoflux_a
    call FB_init(is_local%wrap%FBMed_aoflux_a, is_local%wrap%flds_scalar_name, &
         STgeom=is_local%wrap%NStateImp(compatm), fieldnamelist=fldnames_aof_out, name='FBMed_aoflux_a', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (maintask) then
       write(logunit,*)
       write(logunit,'(a)') trim(subname)//' initialized FB FBMed_aoflux_a'
    end if

    ! Initialize FBMed_aoflux_o
    call FB_init(is_local%wrap%FBMed_aoflux_o, is_local%wrap%flds_scalar_name, &
         STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames_aof_out, name='FBMed_aoflux_o', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (maintask) then
       write(logunit,'(a)') trim(subname)//' initialized FB FBMed_aoflux_o'
       write(logunit,'(a)') trim(subname)//' following are the fields in FBMed_aoflux_o and FBMed_aoflux_a'
       do n = 1,fieldcount
          write(logunit,'(a)')'   FBmed_aoflux fieldname = '//trim(fldnames_aof_out(n))
       end do
    end if

    ! Create required field bundles
    if (is_local%wrap%aoflux_grid == 'ogrid' .or. is_local%wrap%aoflux_grid == 'agrid') then

       ! Create the field bundle is_local%wrap%FBImp(compatm,compocn) if needed
       if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compatm,compocn), rc=rc)) then
          if (maintask) then
             write(logunit,'(a)') trim(subname)//' creating field bundle FBImp(compatm,compocn)'
          end if
          call FB_init(is_local%wrap%FBImp(compatm,compocn), is_local%wrap%flds_scalar_name, &
               STgeom=is_local%wrap%NStateImp(compocn), STflds=is_local%wrap%NStateImp(compatm), &
               name='FBImp'//trim(compname(compatm))//'_'//trim(compname(compocn)), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (maintask) then
          write(logunit,'(a)') trim(subname)//' initializing FB for '// &
               trim(compname(compatm))//'_'//trim(compname(compocn))
       end if

       ! Create the field bundle is_local%wrap%FBImp(compocn,compatm) if needed
       if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compocn,compatm), rc=rc)) then
          if (maintask) then
             write(logunit,'(a)') trim(subname)//' creating field bundle FBImp(compocn,compatm)'
          end if
          call FB_init(is_local%wrap%FBImp(compocn,compatm), is_local%wrap%flds_scalar_name, &
               STgeom=is_local%wrap%NStateImp(compatm), STflds=is_local%wrap%NStateImp(compocn), &
               name='FBImp'//trim(compname(compocn))//'_'//trim(compname(compatm)), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (maintask) then
          write(logunit,'(a)') trim(subname)//' initializing FB for '// &
               trim(compname(compocn))//'_'//trim(compname(compatm))
       end if

    end if

  end subroutine med_phases_aofluxes_init_fldbuns

  !================================================================================
  subroutine med_phases_aofluxes_run(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Compute atm/ocn fluxes
    !-----------------------------------------------------------------------

    use NUOPC           , only : NUOPC_CompAttributeGet
    use ESMF            , only : ESMF_FieldBundleIsCreated
    use med_methods_mod , only : FB_diagnose  => med_methods_FB_diagnose
    use med_phases_history_mod, only : med_phases_history_write_med

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)          :: is_local
    type(aoflux_in_type)  , save :: aoflux_in
    type(aoflux_out_type) , save :: aoflux_out
    logical               , save :: aoflux_created
    logical               , save :: first_call = .true.
    character(len=*),parameter :: subname=' (med_phases_aofluxes_run) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (first_call) then
       ! If field bundles have been created for the ocean/atmosphere flux computation
       if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a, rc=rc) .and. &
            ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o, rc=rc)) then

          ! Allocate memroy for the aoflux module data type (mediator atm/ocn field bundle on the ocean grid)
          call med_aofluxes_init(gcomp, aoflux_in, aoflux_out, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          aoflux_created = .true.
       else
          aoflux_created = .false.
       end if

       ! Now set first_call to .false.
       first_call = .false.
    end if

    ! Return if there is no aoflux has not been created
    if ( aoflux_created) then
       ! Start time timer
       call t_startf('MED:'//subname)
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
       endif
       call memcheck(subname, 5, maintask)

       ! Calculate atm/ocn fluxes on the destination grid
       call med_aofluxes_update(gcomp, aoflux_in, aoflux_out, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Write mediator aofluxes
       call med_phases_history_write_med(gcomp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBMed_aoflux_o, &
               string=trim(subname) //' FBAMed_aoflux_o' , rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       call t_stopf('MED:'//subname)
    end if

  end subroutine med_phases_aofluxes_run

  !================================================================================
  subroutine med_aofluxes_init(gcomp, aoflux_in, aoflux_out, rc)

    use NUOPC           , only : NUOPC_CompAttributeGet
    use ESMF            , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LogFoundError
    use ESMF            , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU
    use ESMF            , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF            , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldBundle
    use med_methods_mod , only : FB_fldchk    => med_methods_FB_FldChk
#ifdef CESMCOUPLED
    use shr_flux_mod    , only : shr_flux_adjust_constants
#else
    use flux_atmocn_mod , only : flux_adjust_constants
#endif

    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)   , intent(inout) :: gcomp
    type(aoflux_in_type)  , intent(inout) :: aoflux_in
    type(aoflux_out_type) , intent(inout) :: aoflux_out
    integer               , intent(out)   :: rc

    ! local variables
    type(InternalState) :: is_local
    character(CL)       :: cvalue
    real(R8)            :: flux_convergence        ! convergence criteria for implicit flux computation
    integer             :: flux_max_iteration      ! maximum number of iterations for convergence
    logical             :: coldair_outbreak_mod    ! cold air outbreak adjustment  (Mahrt & Sun 1995,MWR)
    logical             :: isPresent, isSet
    character(*),parameter  :: subName = '(med_aofluxes_init) '
    !-----------------------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS
    call memcheck(subname, 5, maintask)

    call t_startf('MED:'//subname)

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Initialize module variables
    !----------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_wiso
    else
       flds_wiso = .false.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ocn_surface_flux_scheme', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) ocn_surface_flux_scheme
    else
       ocn_surface_flux_scheme = 0
    end if
#ifdef CESMCOUPLED
    if (maintask) then
       write(logunit,*)
       write(logunit,'(a)') trim(subname)//' ocn_surface_flux_scheme is '//trim(cvalue)
    end if
#endif

    call NUOPC_CompAttributeGet(gcomp, name='add_gusts', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) add_gusts
    else
       add_gusts = .false.
    end if

    ! bottom level potential temperature and/or botom level density
    ! will need to be computed if not received from the atm
    if (FB_fldchk(is_local%Wrap%FBImp(Compatm,Compatm), 'Sa_ptem', rc=rc)) then
       compute_atm_thbot = .false.
    else
       compute_atm_thbot = .true.
    end if
    if (FB_fldchk(is_local%Wrap%FBImp(Compatm,Compatm), 'Sa_dens', rc=rc)) then
       compute_atm_dens = .false.
    else
       compute_atm_dens = .true.
    end if

    !----------------------------------
    ! Initialize aoflux
    !----------------------------------

    if (is_local%wrap%aoflux_grid == 'ogrid') then  ! aoflux_grid is ocn
       call med_aofluxes_init_ogrid(gcomp, aoflux_in, aoflux_out, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (is_local%wrap%aoflux_grid == 'agrid') then ! aoflux_grid is atm
       call med_aofluxes_init_agrid(gcomp, aoflux_in, aoflux_out, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (is_local%wrap%aoflux_grid == 'xgrid') then ! aoflux_grid is exchange grid
       call med_aofluxes_init_xgrid(gcomp, aoflux_in, aoflux_out, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------------------
    ! Initialize flux_adjust_constants
    !----------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='coldair_outbreak_mod', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) coldair_outbreak_mod
    else
       coldair_outbreak_mod = .false.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='flux_max_iteration', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flux_max_iteration
    else
       flux_max_iteration = 1
    end if
    call NUOPC_CompAttributeGet(gcomp, name='flux_convergence', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flux_convergence
    else
       flux_convergence = 0.0_r8
    end if

#ifdef CESMCOUPLED
    call shr_flux_adjust_constants(&
         flux_convergence_tolerance=flux_convergence, &
         flux_convergence_max_iteration=flux_max_iteration, &
         coldair_outbreak_mod=coldair_outbreak_mod)
#else
    call flux_adjust_constants(&
         flux_convergence_tolerance=flux_convergence, &
         flux_convergence_max_iteration=flux_max_iteration, &
         coldair_outbreak_mod=coldair_outbreak_mod)
#endif

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_aofluxes_init

  !===============================================================================
  subroutine med_aofluxes_init_ogrid(gcomp, aoflux_in, aoflux_out, rc)

    ! --------------------------------------------
    ! Initialize aoflux data type and compute mask
    ! for computations on ocn grid
    ! --------------------------------------------

    use ESMF        , only : ESMF_FieldBundleIsCreated
    use esmFlds     , only : med_fldlist_GetaofluxfldList
    use esmFlds     , only : med_fldList_type
    use med_map_mod , only : med_map_packed_field_create
    use med_methods_mod , only : FB_fldchk    => med_methods_FB_FldChk

    ! Arguments
    type(ESMF_GridComp)   , intent(inout) :: gcomp
    type(aoflux_in_type)  , intent(inout) :: aoflux_in
    type(aoflux_out_type) , intent(inout) :: aoflux_out
    integer               , intent(out)   :: rc
    !
    ! Local variables
    type(med_fldList_type), pointer :: FldListMed_aoflux
    type(InternalState) :: is_local
    character(len=CX)   :: tmpstr
    integer             :: lsize
    type(ESMF_Field)    :: lfield
    type(ESMF_Mesh)     :: lmesh
    real(R8), pointer   :: garea(:) => null()
    type(ESMF_CoordSys_Flag)   :: coordSys
    character(len=*),parameter :: subname=' (med_aofluxes_init_ocngrid) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    FldListMed_aoflux => med_fldlist_GetaofluxFldList()
    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! input fields from atm and ocn on aofluxgrid
    ! ------------------------
    call set_aoflux_in_pointers(is_local%wrap%FBImp(compatm,compocn), is_local%wrap%FBImp(compocn,compocn), &
         aoflux_in, lsize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! output fields from aoflux calculation
    ! ------------------------
    call set_aoflux_out_pointers(is_local%wrap%FBMed_aoflux_o, lsize, aoflux_out, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! set aoflux computational mask on ocn grid
    ! ------------------------
    ! default compute everywhere, then "turn off" gridcells
    allocate(aoflux_in%mask(lsize))
    aoflux_in%mask(:) = 1
    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux_in%rmask),sum(aoflux_in%mask)
    call ESMF_LogWrite(trim(subname)//" : maskA= "//trim(tmpstr), ESMF_LOGMSG_INFO)
    where (aoflux_in%rmask(:) == 0._R8) aoflux_in%mask(:) = 0   ! like nint
    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux_in%rmask),sum(aoflux_in%mask)
    call ESMF_LogWrite(trim(subname)//" : maskB= "//trim(tmpstr), ESMF_LOGMSG_INFO)

    ! ------------------------
    ! setup grid area
    ! ------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBArea(compocn), 'area', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(aoflux_in%garea(lsize))
    call ESMF_FieldGet(lfield, farrayPtr=garea, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    is_local%wrap%aoflux_mesh = ESMF_MeshCreate(lmesh, rc=rc)
    call ESMF_MeshGet(lmesh, coordSys=coordSys, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (coordSys /= ESMF_COORDSYS_CART) then
       ! Convert square radians to square meters
       aoflux_in%garea(:) = garea(:)*(rearth**2)
    else
       aoflux_in%garea(:) = garea(:)
    end if

    ! ------------------------
    ! create packed mapping from ocn->atm if aoflux_grid is ocn
    ! ------------------------
    if (is_local%wrap%aoflux_grid == 'ogrid') then
       if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o) .and. &
            ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a)) then
          call med_map_packed_field_create(destcomp=compatm, &
               flds_scalar_name=is_local%wrap%flds_scalar_name, &
               fieldsSrc=fldListMed_aoflux, &
               FBSrc=is_local%wrap%FBMed_aoflux_o, &
               FBDst=is_local%wrap%FBMed_aoflux_a, &
               packed_data=is_local%wrap%packed_data_aoflux_o2a(:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

  end subroutine med_aofluxes_init_ogrid

  !===============================================================================
  subroutine med_aofluxes_init_agrid(gcomp, aoflux_in, aoflux_out, rc)

    ! --------------------------------------------
    ! Initialize aoflux data type and compute mask for computations on atm grid
    ! - all aoflux fields are on the atm mesh
    ! - input atm aoflux attributes are just pointers into is_local%wrap%FBImp(compatm,compatm)
    ! - input ocn aoflux attributes are just pointers into is_local%wrap%FBImp(compocn,compatm)
    ! - output aoflux attributes are on the atm mesh
    ! --------------------------------------------

    use med_methods_mod, only : FB_init => med_methods_FB_init
    use med_map_mod    , only : med_map_rh_is_created, med_map_field

    ! Arguments
    type(ESMF_GridComp)   , intent(inout) :: gcomp
    type(aoflux_in_type)  , intent(inout) :: aoflux_in
    type(aoflux_out_type) , intent(inout) :: aoflux_out
    integer               , intent(out)   :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    integer             :: lsize,n
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    real(r8), pointer   :: dataptr1d(:)
    type(ESMF_Mesh)     :: mesh_src
    type(ESMF_Mesh)     :: mesh_dst
    integer             :: maptype
    type(ESMF_Field)    :: lfield
    type(ESMF_Mesh)     :: lmesh
    real(R8), pointer   :: garea(:) => null()
    type(ESMF_CoordSys_Flag)   :: coordSys
    character(len=*),parameter :: subname=' (med_aofluxes_init_atmgrid) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! input fields from atm and ocn on atm grid
    ! ------------------------

    if (flds_wiso) then
       allocate(fldnames_ocn_in(5))
       fldnames_ocn_in = (/'So_omask    ','So_t        ','So_u        ','So_v        ','So_roce_wiso' /)
    else
       allocate(fldnames_ocn_in(4))
       fldnames_ocn_in = (/'So_omask','So_t    ','So_u    ','So_v    '/)
    end if
    call FB_init(FBocn_a, is_local%wrap%flds_scalar_name, &
         FBgeom=is_local%wrap%FBImp(compatm,compatm), fieldnamelist=fldnames_ocn_in, name='FBocn_a', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call set_aoflux_in_pointers(is_local%wrap%FBImp(compatm,compatm), FBocn_a, aoflux_in, lsize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! output fields from aoflux calculation on atm grid
    ! ------------------------

    call set_aoflux_out_pointers(is_local%wrap%FBMed_aoflux_a, lsize, aoflux_out, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! Determine maptype for ocn->atm mapping
    ! ------------------------

    if (med_map_RH_is_created(is_local%wrap%RH(compocn,compatm,:), mapfcopy, rc=rc)) then
       maptype = mapfcopy
    else if (med_map_RH_is_created(is_local%wrap%RH(compocn,compatm,:), mapconsd, rc=rc)) then
       maptype = mapconsd
    else
       call ESMF_LogWrite(trim(subname)//&
            ": maptype for atm->ocn mapping of So_mask must be either mapfcopy or mapconsd", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    end if

    ! ------------------------
    ! set aoflux computational mask on atm grid
    ! ------------------------

    ! Compute mask is the ocean mask mapped to atm grid (conservatively without fractions)
    ! This computes So_omask in FBocn_a - but the assumption is that it already is there
    ! Compute mask is the ocean mask mapped to atm grid (conservatively without fractions)
    ! This computes So_omask in FBocn_a - but the assumption is that it already is there
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), 'So_omask', field=field_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBocn_a, 'So_omask', field=field_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call med_map_field( field_src=field_src, field_dst=field_dst, &
         routehandles=is_local%wrap%RH(compocn,compatm,:), maptype=maptype, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_dst, farrayptr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(aoflux_in%mask(lsize))
    do n = 1,lsize
       if (dataptr1d(n) == 0._r8) then
          aoflux_in%mask(n) = 0
       else
          aoflux_in%mask(n) = 1
       end if
    enddo

    ! ------------------------
    ! setup grid area
    ! ------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBArea(compatm), 'area', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(aoflux_in%garea(lsize))
    call ESMF_FieldGet(lfield, farrayPtr=garea, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    is_local%wrap%aoflux_mesh = ESMF_MeshCreate(lmesh, rc=rc)
    call ESMF_MeshGet(lmesh, coordSys=coordSys, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (coordSys /= ESMF_COORDSYS_CART) then
       ! Convert square radians to square meters
       aoflux_in%garea(:) = garea(:)*(rearth**2)
    else
       aoflux_in%garea(:) = garea(:)
    end if

    ! ------------------------
    ! set one normalization for ocn-atm mapping if needed
    ! ------------------------

    if (.not. ESMF_FieldIsCreated(is_local%wrap%field_NormOne(compocn,compatm,maptype))) then
       ! Get source mesh
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), 'So_omask', field=field_src, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_src, mesh=mesh_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       field_src = ESMF_FieldCreate(mesh_src, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_src, farrayptr=dataPtr1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       dataptr1d(:) = 1.0_R8

       ! Create field is_local%wrap%field_NormOne(compocn,compatm,maptype) and fill in its values
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compatm), 'So_omask', field=field_dst, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_dst, mesh=mesh_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       is_local%wrap%field_NormOne(compocn,compatm,maptype) = ESMF_FieldCreate(mesh_dst, &
            ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       call med_map_field( field_src=field_src, field_dst=is_local%wrap%field_NormOne(compocn,compatm,maptype), &
            routehandles=is_local%wrap%RH(compocn,compatm,:), maptype=maptype, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldDestroy(field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine med_aofluxes_init_agrid

  !===============================================================================
  subroutine med_aofluxes_init_xgrid(gcomp, aoflux_in, aoflux_out, rc)

    ! --------------------------------------------
    ! Initialize aoflux data type and compute mask
    ! for computations on exchange grid
    ! --------------------------------------------

    ! Arguments
    type(ESMF_GridComp)   , intent(inout) :: gcomp
    type(aoflux_in_type)  , intent(inout) :: aoflux_in
    type(aoflux_out_type) , intent(inout) :: aoflux_out
    integer               , intent(out)   :: rc

    ! Local variables
    integer              :: lsize
    type(InternalState)  :: is_local
    type(ESMF_Field)     :: field_a
    type(ESMF_Field)     :: field_o
    type(ESMF_Field)     :: lfield
    type(ESMF_Mesh)      :: lmesh
    type(ESMF_Mesh)      :: ocn_mesh
    type(ESMF_Mesh)      :: atm_mesh
    type(ESMF_Mesh)      :: xch_mesh
    real(r8), pointer    :: dataptr(:)
    integer              :: fieldcount
    integer              :: stp ! srcTermProcessing is declared inout and must have variable not constant
    type(ESMF_CoordSys_Flag)           :: coordSys
    real(ESMF_KIND_R8)    ,allocatable :: garea(:)
    character(len=*),parameter :: subname=' (med_aofluxes_init_xgrid) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! create the aoflux exchange grid
    ! ------------------------

    ! determine atm mesh
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), fieldname='Sa_z', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=atm_mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! determine ocn mesh
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), fieldname='So_t', field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=ocn_mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! create exchange grid - assume that atm mask is always 1
    xgrid = ESMF_XGridCreate(sideBMesh=(/ocn_mesh/), sideAMesh=(/atm_mesh/), sideBMaskValues=(/0/), &
         storeOverlay=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! write meshes for debug purpose
    if (dbug_flag > 20) then
       call ESMF_MeshWrite(atm_mesh, filename="atm_mesh", rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshWrite(ocn_mesh, filename="ocn_mesh", rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_XGridGet(xgrid, mesh=xch_mesh, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshWrite(xch_mesh, filename="xch_mesh", rc=rc)
    end if

    ! create module field on exchange grid and set its initial value to 1
    field_x = ESMF_FieldCreate(xgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_x, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 1.0_r8

    ! ------------------------
    ! input fields from atm and ocn on xgrid
    ! ------------------------

    ! Create FBatm_x and FBocn_x (module variables)
    FBatm_x = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    FBocn_x = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call set_aoflux_in_pointers(FBatm_x, FBocn_x, aoflux_in, lsize, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(FBatm_x, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fldnames_atm_in(fieldcount))
    call ESMF_FieldBundleGet(FBatm_x, fieldnamelist=fldnames_atm_in, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(FBocn_x, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fldnames_ocn_in(fieldcount))
    call ESMF_FieldBundleGet(FBocn_x, fieldnamelist=fldnames_ocn_in, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! output fields from aoflux calculation on exchange grid
    ! ------------------------

    FBaof_x = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call set_aoflux_out_pointers(FBaof_x, lsize, aoflux_out, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! create the routehandles atm->xgrid and xgrid->atm
    ! ------------------------

    ! create temporary field
    field_a = ESMF_FieldCreate(atm_mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_a, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 1.0_r8

    ! create agrid->xgrid route handles
    call ESMF_FieldRegridStore(xgrid, field_a, field_x, routehandle=rh_agrid2xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(xgrid, field_a, field_x, routehandle=rh_agrid2xgrid_2ndord, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE_2ND, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (trim(coupling_mode) == 'cesm') then
       stp = 1
       call ESMF_FieldRegridStore(field_a, field_x, routehandle=rh_agrid2xgrid_bilinr, &
            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, dstMaskValues=(/0/), srcTermProcessing=stp, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldRegridStore(field_a, field_x, routehandle=rh_agrid2xgrid_patch, &
            regridmethod=ESMF_REGRIDMETHOD_PATCH, dstMaskValues=(/0/), srcTermProcessing=stp, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! create xgrid->zgrid route handle
    call ESMF_FieldRegridStore(xgrid, field_x, field_a, routehandle=rh_xgrid2agrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! destroy temporary field
    call ESMF_FieldDestroy(field_a, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! create the routehandles ocn->xgrid and xgrid->ocn
    ! ------------------------

    ! TODO: the second order conservative route handle below error out in its creation

    field_o = ESMF_FieldCreate(ocn_mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 1.0_r8
    call ESMF_FieldRegridStore(xgrid, field_o, field_x, routehandle=rh_ogrid2xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(xgrid, field_x, field_o, routehandle=rh_xgrid2ogrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! call ESMF_FieldRegridStore(xgrid, field_o, field_x, routehandle=rh_ogrid2xgrid_2ndord, &
    !      regridmethod=ESMF_REGRIDMETHOD_CONSERVE_2ND, rc=rc)
    ! if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldDestroy(field_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! setup the compute mask - default compute everywhere for exchange grid
    ! ------------------------

    allocate(aoflux_in%mask(lsize))
    aoflux_in%mask(:) = 1

    ! ------------------------
    ! setup grid area
    ! ------------------------

    allocate(garea(lsize))
    allocate(aoflux_in%garea(lsize))
    call ESMF_XGridGet(xgrid, mesh=lmesh, coordSys=coordSys, area=garea, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    is_local%wrap%aoflux_mesh = ESMF_MeshCreate(lmesh, rc=rc)
    if (coordSys /= ESMF_COORDSYS_CART) then
       ! Convert square radians to square meters
       aoflux_in%garea(:) = garea(:)*(rearth**2)
    else
       aoflux_in%garea(:) = garea(:)
    end if
    deallocate(garea)

  end subroutine med_aofluxes_init_xgrid

  !===============================================================================
  subroutine med_aofluxes_update(gcomp, aoflux_in, aoflux_out, rc)

    !-----------------------------------------------------------------------
    ! Determine atm/ocn fluxes eother on atm, ocn or exchange grid
    ! The module arrays are set via pointers to the mediator internal states
    ! in med_ocnatm_init and are used below.
    !  1) Create input on aoflux grid
    !  2) Update atmosphere/ocean surface fluxes
    !  3) Map aoflux output to relevant atm/ocn grid(s)
    !-----------------------------------------------------------------------

    use ESMF           , only : ESMF_GridComp
    use ESMF           , only : ESMF_LogWrite, ESMF_LogMsg_Info, ESMF_SUCCESS
    use med_map_mod    , only : med_map_field_packed, med_map_rh_is_created
    use med_map_mod    , only : med_map_routehandles_init
    use med_methods_mod, only : FB_fldchk => med_methods_FB_fldchk
    use med_methods_mod, only : FB_diagnose  => med_methods_FB_diagnose
#ifdef CESMCOUPLED
    use shr_flux_mod   , only : flux_atmocn
#else
    use flux_atmocn_mod, only : flux_atmocn
#endif
#ifdef UFS_AOFLUX
    use flux_atmocn_ccpp_mod, only : flux_atmocn_ccpp
#endif

    ! Arguments
    type(ESMF_GridComp)                   :: gcomp
    type(aoflux_in_type)  , intent(inout) :: aoflux_in
    type(aoflux_out_type) , intent(inout) :: aoflux_out
    integer               , intent(out)   :: rc
    !
    ! Local variables
    type(InternalState)      :: is_local
    integer                  :: n                          ! indices
    real(r8), parameter      :: qmin = 1.0e-8_r8
    real(r8), parameter      :: p0 = 100000.0_r8           ! reference pressure in Pa
    real(r8), parameter      :: rcp = 0.286_r8             ! gas constant of air / specific heat capacity at a constant pressure
    real(r8), parameter      :: rdair = 287.058_r8         ! dry air gas constant in J/K/kg
    integer                  :: maptype
    type(ESMF_Field)         :: field_src
    type(ESMF_Field)         :: field_dst
    character(*),parameter   :: subName = '(med_aofluxes_update) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Create input on aoflux grid
    !----------------------------------

    if (is_local%wrap%aoflux_grid == 'ogrid') then

       ! Do nothing - mapping of input atm to ogrid is in med_phases_post_atm
       ! via the call to med_map_field_packed

    else if (is_local%wrap%aoflux_grid == 'agrid') then

       call med_aofluxes_map_ogrid2agrid_input(gcomp, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    else if (is_local%wrap%aoflux_grid == 'xgrid') then

       call med_aofluxes_map_agrid2xgrid_input(gcomp, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call med_aofluxes_map_ogrid2xgrid_input(gcomp, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    end if

    !----------------------------------
    ! Calculate quantities if they are not defined
    !----------------------------------

    ! Note pbot, tbot and shum have already been mapped or are available on the aoflux grid
    if (compute_atm_thbot) then
       do n = 1,aoflux_in%lsize
          if (aoflux_in%mask(n) /= 0.0_r8) then
             aoflux_in%thbot(n) = aoflux_in%tbot(n)*((p0/aoflux_in%pbot(n))**rcp)
          end if
       end do
    end if
    if (compute_atm_dens) then
       if (trim(aoflux_code) == 'ccpp' .and. &
          (trim(coupling_mode) == 'ufs.frac.aoflux')) then
          ! Add limiting factor to humidity to be consistent with UFS aoflux calculation
          do n = 1,aoflux_in%lsize
             if (aoflux_in%mask(n) /= 0.0_r8) then
                aoflux_in%shum(n) = max(aoflux_in%shum(n), qmin)
             end if
          end do
          ! Use pbot as psfc for the initial pass since psfc provided by UFS atm is zero
          if (maxval(aoflux_in%psfc, mask=(aoflux_in%mask/= 0.0_r8)) < 100.0_r8) then
             aoflux_in%psfc(:) = aoflux_in%pbot(:)
             call ESMF_LogWrite(trim(subname)//" : using pbot as psfc for initial pass!", ESMF_LOGMSG_INFO)
          end if
       end if
       do n = 1,aoflux_in%lsize
          if (aoflux_in%mask(n) /= 0.0_r8) then
             aoflux_in%dens(n) = aoflux_in%pbot(n)/(rdair*(1.0_r8 + 0.608_r8*aoflux_in%shum(n))*aoflux_in%tbot(n))
          end if
       end do
    end if

    !----------------------------------
    ! Update atmosphere/ocean surface fluxes
    !----------------------------------

#ifdef CESMCOUPLED
    call flux_atmocn (logunit=logunit, &
         nMax=aoflux_in%lsize, &
         zbot=aoflux_in%zbot, ubot=aoflux_in%ubot, vbot=aoflux_in%vbot, thbot=aoflux_in%thbot, qbot=aoflux_in%shum, &
         rainc=aoflux_in%rainc, &
         s16O=aoflux_in%shum_16O, sHDO=aoflux_in%shum_HDO, s18O=aoflux_in%shum_18O, rbot=aoflux_in%dens, &
         tbot=aoflux_in%tbot, us=aoflux_in%uocn, vs=aoflux_in%vocn, pslv=aoflux_in%psfc, ts=aoflux_in%tocn, &
         mask=aoflux_in%mask, seq_flux_atmocn_minwind=0.5_r8, &
         sen=aoflux_out%sen, lat=aoflux_out%lat, lwup=aoflux_out%lwup, &
         r16O=aoflux_in%roce_16O, rhdo=aoflux_in%roce_HDO, r18O=aoflux_in%roce_18O, &
         evap=aoflux_out%evap, evap_16O=aoflux_out%evap_16O, evap_HDO=aoflux_out%evap_HDO, evap_18O=aoflux_out%evap_18O, &
         taux=aoflux_out%taux, tauy=aoflux_out%tauy, tref=aoflux_out%tref, qref=aoflux_out%qref, &
         ocn_surface_flux_scheme=ocn_surface_flux_scheme, &
         add_gusts=add_gusts, &
         duu10n=aoflux_out%duu10n, &
         ugust_out = aoflux_out%ugust_out, &
         u10res = aoflux_out%u10res, &
         ustar_sv=aoflux_out%ustar, re_sv=aoflux_out%re, ssq_sv=aoflux_out%ssq, &
         missval=0.0_r8)

#else
#ifdef UFS_AOFLUX
     if (trim(aoflux_code) == 'ccpp') then
       call flux_atmocn_ccpp(gcomp=gcomp, maintask=maintask, logunit=logunit, &
            nMax=aoflux_in%lsize, psfc=aoflux_in%psfc, &
            pbot=aoflux_in%pbot, tbot=aoflux_in%tbot, qbot=aoflux_in%shum, lwdn=aoflux_in%lwdn, &
            zbot=aoflux_in%zbot, garea=aoflux_in%garea, ubot=aoflux_in%ubot, usfc=aoflux_in%usfc, vbot=aoflux_in%vbot, &
            vsfc=aoflux_in%vsfc, rbot=aoflux_in%dens, ts=aoflux_in%tocn, mask=aoflux_in%mask, &
            sen=aoflux_out%sen, lat=aoflux_out%lat, lwup=aoflux_out%lwup, evp=aoflux_out%evap, &
            taux=aoflux_out%taux, tauy=aoflux_out%tauy, tref=aoflux_out%tref, qref=aoflux_out%qref, &
            duu10n=aoflux_out%duu10n, ustar_sv=aoflux_out%ustar, re_sv=aoflux_out%re, ssq_sv=aoflux_out%ssq, &
            missval=0.0_r8)
     else
#endif
       call flux_atmocn (logunit=logunit, &
            nMax=aoflux_in%lsize, mask=aoflux_in%mask, &
            zbot=aoflux_in%zbot, ubot=aoflux_in%ubot, vbot=aoflux_in%vbot, thbot=aoflux_in%thbot, qbot=aoflux_in%shum, &
            rbot=aoflux_in%dens, tbot=aoflux_in%tbot, us=aoflux_in%uocn, vs=aoflux_in%vocn, ts=aoflux_in%tocn, &
            ocn_surface_flux_scheme=ocn_surface_flux_scheme, &
            sen=aoflux_out%sen, lat=aoflux_out%lat, lwup=aoflux_out%lwup, evap=aoflux_out%evap, &
            taux=aoflux_out%taux, tauy=aoflux_out%tauy, tref=aoflux_out%tref, qref=aoflux_out%qref, &
            duu10n=aoflux_out%duu10n, &
            missval=0.0_r8)
#ifdef UFS_AOFLUX
     end if
#endif

#endif

    do n = 1,aoflux_in%lsize
       if (aoflux_in%mask(n) /= 0) then
          aoflux_out%u10(n)          = aoflux_out%u10res(n)
          aoflux_out%u10_withGust(n) = sqrt(aoflux_out%duu10n(n))
       end if
    enddo

    !----------------------------------
    ! map aoflux output to relevant atm/ocn grid(s)
    !----------------------------------

    if (is_local%wrap%aoflux_grid == 'ogrid') then

       ! mapping aoflux from ogrid to agrid is done in med_phases_prep_atm
       ! which is called from med_phases_prep_atm (since need to use updated ocean fractions)

    else if (is_local%wrap%aoflux_grid == 'agrid') then

       if (is_local%wrap%med_coupling_active(compatm,compocn)) then
          call med_aofluxes_map_agrid2ogrid_output(gcomp, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

    else if (is_local%wrap%aoflux_grid == 'xgrid') then

       ! mapping aoflux from xgrid to agrid is done in med_aofluxes_map_xgrid2agrid_output
       ! which is called from med_phases_prep_atm (since need to use updated ocean fractions)
       call med_aofluxes_map_xgrid2ogrid_output(gcomp, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    end if

    ! map taux and tauy from ocean to wave grid if stresses are needed on the wave grid
    if ( FB_fldchk(is_local%wrap%FBExp(compwav), 'Fwxx_taux', rc=rc) .and. &
         FB_fldchk(is_local%wrap%FBExp(compwav), 'Fwxx_tauy', rc=rc)) then
       maptype = mapconsf
       if (.not. med_map_RH_is_created(is_local%wrap%RH(compocn,compwav,:), maptype, rc=rc)) then
          call med_map_routehandles_init( compocn, compwav, &
               FBSrc=is_local%wrap%FBImp(compocn,compocn), &
               FBDst=is_local%wrap%FBImp(compwav,compwav), &
               mapindex=maptype, RouteHandle=is_local%wrap%RH, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, 'Faox_taux', field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compwav), 'Fwxx_taux', field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldRegrid(field_src, field_dst, &
            routehandle=is_local%wrap%RH(compocn, compwav, maptype), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, 'Faox_tauy', field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(is_local%wrap%FBExp(compwav), 'Fwxx_tauy', field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldRegrid(field_src, field_dst, &
            routehandle=is_local%wrap%RH(compocn, compwav, maptype), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_aofluxes_update

  !================================================================================
  subroutine med_aofluxes_map_ogrid2agrid_input(gcomp, rc)

    ! aoflux is on agrid and this maps the ogrid input to the agrid

    use med_map_mod, only : med_map_RH_is_created

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    real(r8), pointer   :: data_normdst(:)
    real(r8), pointer   :: data_dst(:)
    integer             :: nf,n
    integer             :: maptype
    character(*),parameter  :: subName = '(med_aofluxes_map_ogrid2agrid_input) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Map input ocn to agrid
    do nf = 1,size(fldnames_ocn_in)
       ! Create source field
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), fldnames_ocn_in(nf), field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Create destination field
       call ESMF_FieldBundleGet(FBocn_a, fldnames_ocn_in(nf), field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Determine maptype from ocn->atm
       if (med_map_RH_is_created(is_local%wrap%RH(compocn,compatm,:), mapfcopy, rc=rc)) then
          maptype = mapfcopy
       else if (med_map_RH_is_created(is_local%wrap%RH(compocn,compatm,:), mapconsd, rc=rc)) then
          maptype = mapconsd
       else
          call ESMF_LogWrite(trim(subname)//&
               ": maptype for atm->ocn mapping of aofluxes from atm->ocn either mapfcopy or mapconsd", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       end if

       ! Map ocn->atm conservatively without fractions
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=is_local%wrap%RH(compocn,compatm, maptype), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Normalization of map by 'one'
       if (maptype /= mapfcopy) then
          call ESMF_FieldGet(is_local%wrap%field_normOne(compocn,compatm,maptype), farrayPtr=data_normdst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field_dst, farrayptr=data_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do n = 1,size(data_dst)
             if (data_normdst(n) == 0.0_r8) then
                data_dst(n) = 0.0_r8
             else
                data_dst(n) = data_dst(n)/data_normdst(n)
             end if
          end do
       end if
    end do

  end subroutine med_aofluxes_map_ogrid2agrid_input

  !================================================================================
  subroutine med_aofluxes_map_agrid2xgrid_input(gcomp, rc)

    ! Map input atm to xgrid

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    integer             :: nf
    character(*),parameter  :: subName = '(med_aofluxes_map_agrid2xgrid_input) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do nf = 1,size(fldnames_atm_in)
       ! Get the source field
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), fldnames_atm_in(nf), field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Get the destination field
       call ESMF_FieldBundleGet(FBatm_x, fldnames_atm_in(nf), field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Map atm->xgrid
       if (trim(fldnames_atm_in(nf)) == 'Sa_u' .or. (trim(fldnames_atm_in(nf)) == 'Sa_v')) then
          if (trim(coupling_mode) == 'cesm') then
             call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_agrid2xgrid_patch, &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
          else
             call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_agrid2xgrid_2ndord, &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
          end if
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          if (trim(coupling_mode) == 'cesm') then
             call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_agrid2xgrid_bilinr, &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
          else
             call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_agrid2xgrid, &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
          end if
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine med_aofluxes_map_agrid2xgrid_input

  !================================================================================
  subroutine med_aofluxes_map_ogrid2xgrid_input(gcomp, rc)

    ! Map input ocn to xgrid

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    integer             :: nf
    character(*),parameter  :: subName = '(med_aofluxes_map_ogrid2xgrid_input) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do nf = 1,size(fldnames_ocn_in)
       ! Create source field
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), fldnames_ocn_in(nf), field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! Create destination field
       call ESMF_FieldBundleGet(FBocn_x, fldnames_ocn_in(nf), field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! Map ocn->xgrid conservatively without fractions
       if (trim(fldnames_atm_in(nf)) == 'So_u' .or. (trim(fldnames_atm_in(nf)) == 'So_v')) then
          call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_ogrid2xgrid, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       else
          call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_ogrid2xgrid, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       end if
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

end subroutine med_aofluxes_map_ogrid2xgrid_input

  !================================================================================
  subroutine med_aofluxes_map_ogrid2agrid_output(gcomp, rc)

    use med_map_mod, only : med_map_field_packed

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    character(*),parameter  :: subName = '(med_aofluxes_map_ogrid2agrid_output) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call med_map_field_packed( &
         FBSrc=is_local%wrap%FBMed_aoflux_o, &
         FBDst=is_local%wrap%FBMed_aoflux_a, &
         FBFracSrc=is_local%wrap%FBFrac(compocn), &
         field_normOne=is_local%wrap%field_normOne(compocn,compatm,:), &
         packed_data=is_local%wrap%packed_data_aoflux_o2a(:), &
         routehandles=is_local%wrap%RH(compocn,compatm,:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_aofluxes_map_ogrid2agrid_output

  !================================================================================
  subroutine med_aofluxes_map_agrid2ogrid_output(gcomp, rc)

    ! map aoflux from agrid to ogrid
    use med_map_mod    , only : med_map_field_packed, med_map_rh_is_created

    ! Arguments
    type(ESMF_GridComp)                   :: gcomp
    integer               , intent(out)   :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    integer             :: nf                     ! indices
    integer             :: maptype
    character(*),parameter  :: subName = '(med_aofluxes_map_agrid2ogrid_output) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do nf = 1,size(fldnames_aof_out)
       ! Create source field
       call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_a, fldnames_aof_out(nf), field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! Create destination field
       call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fldnames_aof_out(nf), field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! Map atm->ocn conservatively WITHOUT fractions
       if (med_map_RH_is_created(is_local%wrap%RH(compatm,compocn,:), mapfcopy, rc=rc)) then
          maptype = mapfcopy
       else if (med_map_RH_is_created(is_local%wrap%RH(compatm,compocn,:), mapconsf, rc=rc)) then
          maptype = mapconsf
       else
          call ESMF_LogWrite(trim(subname)//&
               ": maptype for atm->ocn mapping of aofluxes from atm->ocn either mapfcopy or mapconsf", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_FieldRegrid(field_src, field_dst, &
            routehandle=is_local%wrap%RH(compatm, compocn, maptype), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine med_aofluxes_map_agrid2ogrid_output

!================================================================================
  subroutine med_aofluxes_map_xgrid2agrid_output(gcomp, rc)

    use ESMF, only : ESMF_FieldBundleIsCreated

    ! Arguments
    type(ESMF_GridComp)                   :: gcomp
    integer               , intent(out)   :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    type(ESMF_Field)    :: lfield
    integer             :: n,nf                     ! indices
    real(r8), pointer   :: data_src(:)
    real(r8), pointer   :: data_src_save(:)
    real(r8), pointer   :: data_dst(:)
    real(r8), pointer   :: ofrac_x(:)
    real(r8), pointer   :: ofrac_a(:)
    character(*),parameter  :: subName = '(med_aofluxes_map_xgrid2agrid_output) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. ESMF_FieldBundleIsCreated(FBaof_x)) then
       RETURN
    end if

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Map ocn fraction on ocn mesh to xgrid
    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), 'ofrac', field=field_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegrid(field_o, field_x, routehandle=rh_ogrid2xgrid, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_x, farrayptr=ofrac_x, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do nf = 1,size(fldnames_aof_out)

       ! Get the source field
       call ESMF_FieldBundleGet(FBaof_x, fldnames_aof_out(nf), field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! map aoflux from xgrid to agrid followed by normalization by 'one'
       call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_a, fldnames_aof_out(nf), field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_src, farrayptr=data_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(data_src_save(size(data_src)))
       data_src_save(:) = data_src(:)
       do n = 1,size(data_src)
          data_src(n) = data_src(n) * ofrac_x(n)
       end do
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_xgrid2agrid, &
            termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       data_src(:) = data_src_save(:)
       deallocate(data_src_save)
       call ESMF_FieldGet(field_dst, farrayptr=data_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! normalization by '1./ofrac_a'
       call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compatm), 'ofrac', field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayptr=ofrac_a, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(ofrac_a)
          if (ofrac_a(n) == 0.0_r8) then
             data_dst(n) = 0.0_r8
          else
             data_dst(n) = data_dst(n)/ofrac_a(n)
          end if
       end do

    end do

  end subroutine med_aofluxes_map_xgrid2agrid_output

!================================================================================
  subroutine med_aofluxes_map_xgrid2ogrid_output(gcomp, rc)

    ! map aoflx output from xgrid->ogrid

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    integer , intent(out) :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    integer             :: nf                     ! indices
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    character(*),parameter  :: subName = '(med_aofluxes_map_xgrid2ogrid_output) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do nf = 1,size(fldnames_aof_out)
       ! Get the source field
       call ESMF_FieldBundleGet(FBaof_x, fldnames_aof_out(nf), field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! map aoflx from xgrid->ogrid conservatively
       call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fldnames_aof_out(nf), field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_xgrid2ogrid, &
            termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine med_aofluxes_map_xgrid2ogrid_output

!================================================================================
  subroutine set_aoflux_in_pointers(fldbun_a, fldbun_o, aoflux_in, lsize, xgrid, rc)

    ! Set pointers for aoflux_in attributes
    ! Note that if computation is on the xgrid, fldbun_a and fldbun_o are both fldbun_x

    use med_methods_mod , only : FB_fldchk    => med_methods_FB_FldChk

    ! input/output variables
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun_a
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun_o
    type(aoflux_in_type)       , intent(inout) :: aoflux_in
    integer                    , intent(out)   :: lsize
    type(ESMF_Xgrid), optional , intent(inout) :: xgrid
    integer                    , intent(out)   :: rc
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! ------------------------
    ! input fields from atm on aoflux grid
    ! ------------------------

    ! Determine lsize from first field
    call fldbun_getfldptr(fldbun_a, 'Sa_z', aoflux_in%zbot, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(aoflux_in%zbot)
    aoflux_in%lsize = lsize

    ! bulk formula quantities for ufs non-frac with med-aoflux
    if (trim(coupling_mode) == 'ufs.nfrac.aoflux' .and. ocn_surface_flux_scheme == -1) then
       call fldbun_getfldptr(fldbun_a, 'Sa_u10m', aoflux_in%ubot, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_v10m', aoflux_in%vbot, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_t2m', aoflux_in%tbot, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_q2m', aoflux_in%shum, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call fldbun_getfldptr(fldbun_a, 'Sa_u', aoflux_in%ubot, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_v', aoflux_in%vbot, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_tbot', aoflux_in%tbot, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_shum', aoflux_in%shum, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (add_gusts) then
          call fldbun_getfldptr(fldbun_a, 'Faxa_rainc', aoflux_in%rainc, xgrid=xgrid, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! extra fields for ufs.frac.aoflux
    if (trim(coupling_mode) == 'ufs.frac.aoflux') then
       call fldbun_getfldptr(fldbun_a, 'Sa_u10m', aoflux_in%usfc, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_v10m', aoflux_in%vsfc, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Faxa_lwdn', aoflux_in%lwdn, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! bottom level potential temperature will need to be computed if not received from the atm
    if (compute_atm_thbot) then
       allocate(aoflux_in%thbot(lsize))
    else
       call fldbun_getfldptr(fldbun_a, 'Sa_ptem', aoflux_in%thbot, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! bottom level density will need to be computed if not received from the atm
    if (compute_atm_dens) then
       allocate(aoflux_in%dens(lsize))
    else
       call fldbun_getfldptr(fldbun_a, 'Sa_dens', aoflux_in%dens, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (FB_fldchk(fldbun_a, 'Sa_pslv', rc=rc)) then
       call fldbun_getfldptr(fldbun_a, 'Sa_pslv', aoflux_in%psfc, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! if either density or potential temperature are computed, will need bottom level pressure
    if (compute_atm_dens .or. compute_atm_thbot) then
       call fldbun_getfldptr(fldbun_a, 'Sa_pbot', aoflux_in%pbot, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (trim(coupling_mode) == 'ufs.frac.aoflux') then
          call fldbun_getfldptr(fldbun_a, 'Sa_pslv', aoflux_in%psfc, xgrid=xgrid, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    if (flds_wiso) then
       call fldbun_getfldptr(fldbun_a, 'Sa_shum_16O', aoflux_in%shum_16O, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_shum_18O', aoflux_in%shum_18O, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_a, 'Sa_shum_HDO', aoflux_in%shum_HDO, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux_in%shum_16O(lsize)); aoflux_in%shum_16O(:) = 0._R8
       allocate(aoflux_in%shum_18O(lsize)); aoflux_in%shum_18O(:) = 0._R8
       allocate(aoflux_in%shum_HDO(lsize)); aoflux_in%shum_HDO(:) = 0._R8
    end if

    ! ------------------------
    ! input fields from ocn on aoflux_grid
    ! ------------------------

    ! point directly into input field bundle from ocean on the ocean grid
    call fldbun_getfldptr(fldbun_o, 'So_omask', aoflux_in%rmask, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun_o, 'So_t', aoflux_in%tocn, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun_o, 'So_u', aoflux_in%uocn, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun_o, 'So_v', aoflux_in%vocn, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call fldbun_getfldptr(fldbun_o, 'So_roce_16O', aoflux_in%roce_16O, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_o, 'So_roce_18O', aoflux_in%roce_18O, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun_o, 'So_roce_HDO', aoflux_in%roce_HDO, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux_in%roce_16O(aoflux_in%lsize)); aoflux_in%roce_16O(:) = 0._R8
       allocate(aoflux_in%roce_18O(aoflux_in%lsize)); aoflux_in%roce_18O(:) = 0._R8
       allocate(aoflux_in%roce_HDO(aoflux_in%lsize)); aoflux_in%roce_HDO(:) = 0._R8
    end if

  end subroutine set_aoflux_in_pointers

  !================================================================================
  subroutine set_aoflux_out_pointers(fldbun, lsize, aoflux_out, xgrid, rc)

    ! input/output variables
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun
    integer                    , intent(in)    :: lsize
    type(aoflux_out_type)      , intent(inout) :: aoflux_out
    type(ESMF_Xgrid), optional , intent(inout) :: xgrid
    integer                    , intent(out)   :: rc

    rc = ESMF_SUCCESS
    !-----------------------------------------------------------------------

    call fldbun_getfldptr(fldbun, 'So_tref', aoflux_out%tref, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'So_qref', aoflux_out%qref, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'So_ustar', aoflux_out%ustar, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'So_re', aoflux_out%re, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'So_ssq', aoflux_out%ssq, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'So_u10', aoflux_out%u10, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'So_duu10n', aoflux_out%duu10n, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call fldbun_getfldptr(fldbun, 'So_ugustOut', aoflux_out%ugust_out, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'So_u10withGust', aoflux_out%u10_withGust, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'So_u10res', aoflux_out%u10res, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'Faox_taux', aoflux_out%taux, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'Faox_tauy', aoflux_out%tauy, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'Faox_lat', aoflux_out%lat, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'Faox_sen', aoflux_out%sen, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'Faox_evap', aoflux_out%evap, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call fldbun_getfldptr(fldbun, 'Faox_lwup', aoflux_out%lwup, xgrid=xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call fldbun_getfldptr(fldbun, 'Faox_evap_16O', aoflux_out%evap_16O, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun, 'Faox_evap_18O', aoflux_out%evap_18O, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call fldbun_getfldptr(fldbun, 'Faox_evap_HDO', aoflux_out%evap_HDO, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux_out%evap_16O(lsize)); aoflux_out%evap_16O(:) = 0._R8
       allocate(aoflux_out%evap_18O(lsize)); aoflux_out%evap_18O(:) = 0._R8
       allocate(aoflux_out%evap_HDO(lsize)); aoflux_out%evap_HDO(:) = 0._R8
    end if

    if (add_gusts) then
       call fldbun_getfldptr(fldbun, 'So_ugustOut', aoflux_out%ugust_out, xgrid=xgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux_out%ugust_out(lsize)); aoflux_out%ugust_out(:) = 0._R8
    end if

  end subroutine set_aoflux_out_pointers

  !================================================================================
  subroutine fldbun_getfldptr(fldbun, fldname, fldptr, xgrid, rc)

    ! input/output variables
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun
    character(len=*)           , intent(in)    :: fldname
    real(r8)                   , pointer       :: fldptr(:)
    type(ESMF_Xgrid), optional , intent(in)    :: xgrid
    integer                    , intent(out)   :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    !-----------------------------------------------------------------------
    rc = ESMF_SUCCESS

    if (present(xgrid)) then
       lfield = ESMF_FieldCreate(xgrid, typekind=ESMF_TYPEKIND_R8, name=trim(fldname), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(fldbun, (/lfield/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_FieldBundleGet(fldbun, trim(fldname), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine fldbun_getfldptr

end module med_phases_aofluxes_mod
