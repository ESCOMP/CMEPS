module med_phases_aofluxes_mod

  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate, ESMF_FieldIsCreated, ESMF_FieldDestroy
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_RouteHandle, ESMF_FieldRegrid, ESMF_FieldRegridStore
  use ESMF                  , only : ESMF_TERMORDER_SRCSEQ, ESMF_REGION_TOTAL, ESMF_MESHLOC_ELEMENT, ESMF_MAXSTR
  use ESMF                  , only : ESMF_XGRIDSIDE_B, ESMF_XGRIDSIDE_A, ESMF_END_ABORT, ESMF_LOGERR_PASSTHRU
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet, ESMF_XGrid, ESMF_XGridCreate, ESMF_TYPEKIND_R8
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF                  , only : ESMF_Finalize, ESMF_LogFoundError
  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : InternalState
  use med_internalstate_mod , only : mastertask, logunit
  use med_constants_mod     , only : dbug_flag    => med_constants_dbug_flag
  use med_utils_mod         , only : memcheck     => med_memcheck
  use med_utils_mod         , only : chkerr       => med_utils_chkerr
  use med_methods_mod       , only : FB_fldchk    => med_methods_FB_FldChk
  use med_methods_mod       , only : FB_GetFldPtr => med_methods_FB_GetFldPtr
  use med_methods_mod       , only : FB_diagnose  => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_init      => med_methods_FB_init
  use med_map_mod           , only : med_map_field
  use esmFlds               , only : compatm, compocn, coupling_mode, mapconsd, mapconsf
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public routines
  !--------------------------------------------------------------------------

  public :: med_phases_aofluxes_run

  !--------------------------------------------------------------------------
  ! Private routines
  !--------------------------------------------------------------------------

  private :: med_aofluxes_init
  private :: med_aofluxes_init_ogrid
  private :: med_aofluxes_init_agrid
  private :: med_aofluxes_init_xgrid
  private :: med_aofluxes_update

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  logical :: flds_wiso  ! use case
  logical :: compute_atm_dens
  logical :: compute_atm_thbot
  integer :: ocn_surface_flux_scheme ! use case

  character(len=CS), allocatable :: fldnames_ocn_in(:)
  character(len=CS), allocatable :: fldnames_atm_in(:)
  character(len=CS), allocatable :: fldnames_aof_out(:)

  ! following is needed for atm/ocn fluxes on atm grid
  type(ESMF_FieldBundle) :: FBocn_a ! ocean fields need for aoflux calc on atm grid

  ! following is needed for atm/ocn fluxes on the exchange grid
  type(ESMF_FieldBundle) :: FBocn_x        ! input ocn fields
  type(ESMF_FieldBundle) :: FBatm_x        ! input atm fields
  type(ESMF_FieldBundle) :: FBaof_x        ! output aoflux fields
  type(ESMF_RouteHandle) :: rh_ogrid2xgrid ! ocn->xgrid mapping
  type(ESMF_RouteHandle) :: rh_xgrid2ogrid ! xgrid->ocn mapping
  type(ESMF_RouteHandle) :: rh_agrid2xgrid ! atm->xgrid mapping
  type(ESMF_RouteHandle) :: rh_xgrid2agrid ! xgrid->atm mapping
  type(ESMF_Field)       :: field_ogrid2xgrid_normone
  type(ESMF_Field)       :: field_xgrid2agrid_normone

  type aoflux_type
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
     real(R8) , pointer :: thbot       (:) => null() ! atm potential T
     real(R8) , pointer :: shum        (:) => null() ! atm specific humidity
     real(R8) , pointer :: pbot        (:) => null() ! atm bottom pressure
     real(R8) , pointer :: dens        (:) => null() ! atm bottom density
     real(R8) , pointer :: tbot        (:) => null() ! atm bottom surface T
     real(R8) , pointer :: shum_16O    (:) => null() ! atm H2O tracer
     real(R8) , pointer :: shum_HDO    (:) => null() ! atm HDO tracer
     real(R8) , pointer :: shum_18O    (:) => null() ! atm H218O tracer

     ! output: on aoflux grid
     ! if aoflux grid is ocn - then need to map these to the atm
     ! if aoflux grid is atm - then need to map these to the ocn
     ! if aoflux grid is exchange - will map back to both the atm and ocn
     integer  , pointer :: mask        (:) => null() ! integer ocn domain mask: 0 <=> inactive cell
     real(R8) , pointer :: rmask       (:) => null() ! real    ocn domain mask: 0 <=> inactive cell
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
     real(R8) , pointer :: ustar       (:) => null() ! saved ustar
     real(R8) , pointer :: re          (:) => null() ! saved re
     real(R8) , pointer :: ssq         (:) => null() ! saved sq

     ! local size
     integer            :: lsize                     ! local size

     ! logical for creation of aoflux instance
     logical            :: created                   ! has this data type been created

     real(r8), pointer :: mapping_norm_one(:)

  end type aoflux_type

  character(len=CS) :: aoflux_grid

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine med_phases_aofluxes_run(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Compute atm/ocn fluxes
    !-----------------------------------------------------------------------

    use ESMF  , only : ESMF_FieldBundleIsCreated
    use NUOPC , only : NUOPC_CompAttributeGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    type(aoflux_type) , save   :: aoflux
    logical           , save   :: first_call = .true.
    character(len=*),parameter :: subname='(med_phases_aofluxes_run)'
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
          call med_aofluxes_init(gcomp, aoflux, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          aoflux%created = .true.
       else
          aoflux%created = .false.
       end if

       ! Now set first_call to .false.
       first_call = .false.
    end if

    ! Return if there is no aoflux has not been created
    if ( aoflux%created) then
       ! Start time timer
       call t_startf('MED:'//subname)
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
       endif
       call memcheck(subname, 5, mastertask)

       ! Calculate atm/ocn fluxes on the destination grid
       call med_aofluxes_update(gcomp, aoflux, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBMed_aoflux_o, &
               string=trim(subname) //' FBAMed_aoflux_o' , rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       call t_stopf('MED:'//subname)
    end if

  end subroutine med_phases_aofluxes_run

  !================================================================================
  subroutine med_aofluxes_init(gcomp, aoflux, rc)

    use ESMF         , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LogFoundError
    use ESMF         , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU
    use ESMF         , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF         , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldBundle
    use NUOPC        , only : NUOPC_CompAttributeGet
    use esmFlds      , only : coupling_mode
    use shr_flux_mod , only : shr_flux_adjust_constants

    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(aoflux_type)   , intent(inout) :: aoflux
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n
    character(CL)       :: cvalue
    character(len=CX)   :: tmpstr
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
    call memcheck(subname, 5, mastertask)

    call t_startf('MED:'//subname)

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Initialize aoflux
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

    if (is_local%wrap%aoflux_grid == 'ogrid') then  ! aoflux_grid is ocn
       call med_aofluxes_init_ogrid(gcomp, aoflux, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (is_local%wrap%aoflux_grid == 'agrid') then ! aoflux_grid is atm
       call med_aofluxes_init_agrid(gcomp, aoflux, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (is_local%wrap%aoflux_grid == 'xgrid') then ! aoflux_grid is exchange grid
       call med_aofluxes_init_xgrid(gcomp, aoflux, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------------------
    ! Initialize shr_flux_adjust_constants
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
    call shr_flux_adjust_constants(&
         flux_convergence_tolerance=flux_convergence, &
         flux_convergence_max_iteration=flux_max_iteration, &
         coldair_outbreak_mod=coldair_outbreak_mod)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_aofluxes_init

  !===============================================================================
  subroutine med_aofluxes_init_ogrid(gcomp, aoflux, rc)

    ! --------------------------------------------
    ! Initialize aoflux data type and compute mask
    ! for computations on ocn grid
    ! --------------------------------------------

    ! Arguments
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(aoflux_type)   , intent(inout) :: aoflux
    integer             , intent(out)   :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    character(len=CX)   :: tmpstr
    integer             :: lsize
    integer             :: fieldcount
    character(len=*),parameter :: subname='(med_aofluxes_init_ocngrid)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fldnames_aof_out(fieldcount))
    call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fieldNameList=fldnames_aof_out, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! input fields from ocn
    ! ------------------------

    ! point directly into input field bundle from ocean on the ocean grid
    call FB_GetFldPtr(Is_local%Wrap%FBImp(compocn,compocn), fldname='So_omask', fldptr1=aoflux%rmask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(aoflux%rmask)
    aoflux%lsize = lsize
    call FB_GetFldPtr(Is_local%Wrap%FBImp(compocn,compocn), fldname='So_t', fldptr1=aoflux%tocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(Is_local%Wrap%FBImp(compocn,compocn), fldname='So_u', fldptr1=aoflux%uocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(Is_local%Wrap%FBImp(compocn,compocn), fldname='So_v', fldptr1=aoflux%vocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call FB_GetFldPtr(Is_local%Wrap%FBImp(compocn,compocn), fldname='So_roce_16O', fldptr1=aoflux%roce_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(Is_local%Wrap%FBImp(compocn,compocn), fldname='So_roce_18O', fldptr1=aoflux%roce_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(Is_local%Wrap%FBImp(compocn,compocn), fldname='So_roce_HDO', fldptr1=aoflux%roce_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%roce_16O(lsize)); aoflux%roce_16O(:) = 0._R8
       allocate(aoflux%roce_18O(lsize)); aoflux%roce_18O(:) = 0._R8
       allocate(aoflux%roce_HDO(lsize)); aoflux%roce_HDO(:) = 0._R8
    end if

    ! ------------------------
    ! input fields from atm
    ! ------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_z', fldptr1=aoflux%zbot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! bulk formula quantities for nems_orig_data
    if (trim(coupling_mode) == 'nems_orig_data' .and. ocn_surface_flux_scheme == -1) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_u10m', fldptr1=aoflux%ubot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_v10m', fldptr1=aoflux%vbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_t2m', fldptr1=aoflux%tbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_q2m', fldptr1=aoflux%shum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_u', fldptr1=aoflux%ubot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_v', fldptr1=aoflux%vbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_tbot', fldptr1=aoflux%tbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_shum', fldptr1=aoflux%shum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! bottom level potential temperature will need to be computed if not received from the atm
    if (FB_fldchk(is_local%wrap%FBImp(compatm,compocn), 'Sa_ptem', rc=rc)) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_ptem', fldptr1=aoflux%thbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_thbot = .false.
    else
       allocate(aoflux%thbot(lsize))
       compute_atm_thbot = .true.
    end if

    ! bottom level density will need to be computed if not received from the atm
    if (FB_fldchk(is_local%wrap%FBImp(compatm,compocn), 'Sa_dens', rc=rc)) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_dens', fldptr1=aoflux%dens, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_dens = .false.
    else
       compute_atm_dens = .true.
       allocate(aoflux%dens(lsize))
    end if

    ! if either density or potential temperature are computed, will need bottom level pressure
    if (compute_atm_dens .or. compute_atm_thbot) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_pbot', fldptr1=aoflux%pbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (flds_wiso) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_shum_16O', fldptr1=aoflux%shum_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_shum_18O', fldptr1=aoflux%shum_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Sa_shum_HDO', fldptr1=aoflux%shum_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%shum_16O(lsize)); aoflux%shum_16O(:) = 0._R8
       allocate(aoflux%shum_18O(lsize)); aoflux%shum_18O(:) = 0._R8
       allocate(aoflux%shum_HDO(lsize)); aoflux%shum_HDO(:) = 0._R8
    end if

    ! ------------------------
    ! output fields from aoflux calculation
    ! ------------------------

    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='So_tref', fldptr1=aoflux%tref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='So_qref', fldptr1=aoflux%qref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='So_ustar', fldptr1=aoflux%ustar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='So_re', fldptr1=aoflux%re, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='So_ssq', fldptr1=aoflux%ssq, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='So_u10', fldptr1=aoflux%u10, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='So_duu10n', fldptr1=aoflux%duu10n, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_taux', fldptr1=aoflux%taux, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_tauy', fldptr1=aoflux%tauy, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_lat', fldptr1=aoflux%lat, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_sen', fldptr1=aoflux%sen, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_evap', fldptr1=aoflux%evap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_evap_16O', fldptr1=aoflux%evap_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_evap_18O', fldptr1=aoflux%evap_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_evap_HDO', fldptr1=aoflux%evap_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%evap_16O(lsize)); aoflux%evap_16O(:) = 0._R8
       allocate(aoflux%evap_18O(lsize)); aoflux%evap_18O(:) = 0._R8
       allocate(aoflux%evap_HDO(lsize)); aoflux%evap_HDO(:) = 0._R8
    end if
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_o, fldname='Faox_lwup', fldptr1=aoflux%lwup, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! compute mask
    ! ------------------------

    ! default compute everywhere, then "turn off" gridcells
    allocate(aoflux%mask(lsize))
    aoflux%mask(:) = 1
    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux%rmask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : maskA= "//trim(tmpstr), ESMF_LOGMSG_INFO)
    where (aoflux%rmask(:) == 0._R8) aoflux%mask(:) = 0   ! like nint
    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux%rmask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : maskB= "//trim(tmpstr), ESMF_LOGMSG_INFO)

  end subroutine med_aofluxes_init_ogrid

  !===============================================================================
  subroutine med_aofluxes_init_agrid(gcomp, aoflux, rc)

    ! --------------------------------------------
    ! Initialize aoflux data type and compute mask
    ! for computations on atm grid
    ! all aoflux fields are on the atm mesh
    ! - input atm aoflux attributes are just pointers into
    !   is_local%wrap%FBImp(compatm,compatm)
    ! - input ocn aoflux attributes are just pointers into
    !   is_local%wrap%FBImp(compocn,compatm)
    ! - output aoflux attributes are on the atm mesh
    ! --------------------------------------------

    ! Arguments
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(aoflux_type)   , intent(inout) :: aoflux
    integer             , intent(out)   :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    integer             :: lsize,n
    integer             :: fieldcount
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    real(r8), pointer   :: dataptr1d(:)
    type(ESMF_Mesh)     :: mesh_src
    type(ESMF_Mesh)     :: mesh_dst
    character(len=*),parameter :: subname='(med_aofluxes_init_atmgrid)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_a, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fldnames_aof_out(fieldcount))
    call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_a, fieldNameList=fldnames_aof_out, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! input ocn fields on atm grid
    ! ------------------------

    ! create field bundles FBocn_a for the above fieldlist
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

    ! point directly into FBocn_a for ocean fields on the atm grid
    call FB_GetFldPtr(FBocn_a, fldname='So_t', fldptr1=aoflux%tocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(aoflux%tocn)
    aoflux%lsize = lsize
    call FB_GetFldPtr(FBocn_a, fldname='So_u', fldptr1=aoflux%uocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBocn_a, fldname='So_v', fldptr1=aoflux%vocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call FB_GetFldPtr(FBocn_a, fldname='So_roce_16O', fldptr1=aoflux%roce_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBocn_a, fldname='So_roce_18O', fldptr1=aoflux%roce_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBocn_a, fldname='So_roce_HDO', fldptr1=aoflux%roce_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%roce_16O(lsize)); aoflux%roce_16O(:) = 0._R8
       allocate(aoflux%roce_18O(lsize)); aoflux%roce_18O(:) = 0._R8
       allocate(aoflux%roce_HDO(lsize)); aoflux%roce_HDO(:) = 0._R8
    end if

    ! ------------------------
    ! input atm fields on atm grid
    ! ------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_z', fldptr1=aoflux%zbot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! bulk formula quantities for nems_orig_data
    if (trim(coupling_mode) == 'nems_orig_data' .and. ocn_surface_flux_scheme == -1) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_u10m', fldptr1=aoflux%ubot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_v10m', fldptr1=aoflux%vbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_t2m', fldptr1=aoflux%tbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_q2m', fldptr1=aoflux%shum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_u', fldptr1=aoflux%ubot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_v', fldptr1=aoflux%vbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_tbot', fldptr1=aoflux%tbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_shum', fldptr1=aoflux%shum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! bottom level potential temperature will need to be computed if not received from the atm
    if (FB_fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc)) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_ptem', fldptr1=aoflux%thbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_thbot = .false.
    else
       allocate(aoflux%thbot(lsize))
       compute_atm_thbot = .true.
    end if

    ! bottom level density will need to be computed if not received from the atm
    if (FB_fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_dens', rc=rc)) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_dens', fldptr1=aoflux%dens, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_dens = .false.
    else
       compute_atm_dens = .true.
       allocate(aoflux%dens(lsize))
    end if

    ! if either density or potential temperature are computed, will need bottom level pressure
    if (compute_atm_dens .or. compute_atm_thbot) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_pbot', fldptr1=aoflux%pbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (flds_wiso) then
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_shum_16O', fldptr1=aoflux%shum_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_shum_18O', fldptr1=aoflux%shum_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), fldname='Sa_shum_HDO', fldptr1=aoflux%shum_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%shum_16O(lsize)); aoflux%shum_16O(:) = 0._R8
       allocate(aoflux%shum_18O(lsize)); aoflux%shum_18O(:) = 0._R8
       allocate(aoflux%shum_HDO(lsize)); aoflux%shum_HDO(:) = 0._R8
    end if

    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_tref', fldptr1=aoflux%tref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_qref', fldptr1=aoflux%qref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_ustar', fldptr1=aoflux%ustar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_re', fldptr1=aoflux%re, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_ssq', fldptr1=aoflux%ssq, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_u10', fldptr1=aoflux%u10, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_duu10n', fldptr1=aoflux%duu10n, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_taux', fldptr1=aoflux%taux, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_tauy', fldptr1=aoflux%tauy, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_lat', fldptr1=aoflux%lat, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_sen', fldptr1=aoflux%sen, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_evap', fldptr1=aoflux%evap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_evap_16O', fldptr1=aoflux%evap_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_evap_18O', fldptr1=aoflux%evap_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_evap_HDO', fldptr1=aoflux%evap_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%evap_16O(lsize)); aoflux%evap_16O(:) = 0._R8
       allocate(aoflux%evap_18O(lsize)); aoflux%evap_18O(:) = 0._R8
       allocate(aoflux%evap_HDO(lsize)); aoflux%evap_HDO(:) = 0._R8
    end if
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_lwup', fldptr1=aoflux%lwup, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! output fields from aoflux calculation on atm grid
    ! ------------------------

    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_tref', fldptr1=aoflux%tref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_qref', fldptr1=aoflux%qref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_ustar', fldptr1=aoflux%ustar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_re', fldptr1=aoflux%re, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_ssq', fldptr1=aoflux%ssq, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_u10', fldptr1=aoflux%u10, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='So_duu10n', fldptr1=aoflux%duu10n, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_taux', fldptr1=aoflux%taux, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_tauy', fldptr1=aoflux%tauy, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_lat', fldptr1=aoflux%lat, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_sen', fldptr1=aoflux%sen, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_evap', fldptr1=aoflux%evap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_evap_16O', fldptr1=aoflux%evap_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_evap_18O', fldptr1=aoflux%evap_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_evap_HDO', fldptr1=aoflux%evap_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%evap_16O(lsize)); aoflux%evap_16O(:) = 0._R8
       allocate(aoflux%evap_18O(lsize)); aoflux%evap_18O(:) = 0._R8
       allocate(aoflux%evap_HDO(lsize)); aoflux%evap_HDO(:) = 0._R8
    end if
    call FB_GetFldPtr(is_local%wrap%FBMed_aoflux_a, fldname='Faox_lwup', fldptr1=aoflux%lwup, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! setup the compute mask.
    ! ------------------------

    ! Compute mask is the ocean mask mapped to atm grid (conservatively without fractions)
    ! This computes So_omask in FBocn_a - but the assumption is that it already is there
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), 'So_omask', field=field_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBocn_a, 'So_omask', field=field_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call med_map_field( field_src=field_src, field_dst=field_dst, &
         routehandles=is_local%wrap%RH(compocn,compatm,:), maptype=mapconsd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_dst, farrayptr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(dataptr1d)
    allocate(aoflux%mask(lsize))
    do n = 1,lsize
       if (dataptr1d(n) == 0._r8) then
          aoflux%mask(n) = 0
       else
          aoflux%mask(n) = 1
       end if
    enddo

    if (.not. ESMF_FieldIsCreated(is_local%wrap%field_NormOne(compocn,compatm,mapconsd))) then
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

       ! Create field is_local%wrap%field_NormOne(compocn,compatm,mapconsd) and fill in its values
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compatm), 'So_omask', field=field_dst, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_dst, mesh=mesh_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       is_local%wrap%field_NormOne(compocn,compatm,mapconsd) = ESMF_FieldCreate(mesh_dst, &
            ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       call med_map_field( field_src=field_src, field_dst=is_local%wrap%field_NormOne(compocn,compatm,mapconsd), &
            routehandles=is_local%wrap%RH(compocn,compatm,:), maptype=mapconsd, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldDestroy(field_src, rc=rc, noGarbage=.true.)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine med_aofluxes_init_agrid

  !===============================================================================
  subroutine med_aofluxes_init_xgrid(gcomp, aoflux, rc)

    ! --------------------------------------------
    ! Initialize aoflux data type and compute mask
    ! for computations on exchange grid
    ! --------------------------------------------

    ! Arguments
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(aoflux_type)   , intent(inout) :: aoflux
    integer             , intent(out)   :: rc
    !
    ! Local variables
    integer              :: n
    integer              :: lsize
    type(InternalState)  :: is_local
    type(ESMF_Field)     :: lfield_a
    type(ESMF_Field)     :: lfield_o
    type(ESMF_Field)     :: lfield_x
    type(ESMF_Field)     :: lfield
    integer              :: elementCount
    type(ESMF_Mesh)      :: ocn_mesh
    type(ESMF_Mesh)      :: atm_mesh
    integer, allocatable :: ocn_mask(:)
    type(ESMF_XGrid)     :: xgrid
    type(ESMF_Field)     :: field_src  ! needed for normalization 
    type(ESMF_Field)     :: field_dst  ! needed for normalization 
    type(ESMF_Mesh)      :: mesh_src   ! needed for normalization 
    type(ESMF_Mesh)      :: mesh_dst   ! needed for normalization 
    real(r8), pointer    :: dataptr1d(:) 
    integer              :: fieldcount
    character(ESMF_MAXSTR),allocatable :: fieldNameList(:)
    character(len=*),parameter :: subname='(med_aofluxes_init_xgrid)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! Create exchange grid
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
    xgrid = ESMF_XGridCreate(sideBMesh=(/ocn_mesh/), sideAMesh=(/atm_mesh/), sideBMaskValues=(/0/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! input fields from ocn on xgrid
    ! ------------------------

    ! Create FBocn_x (module variable)
    FBocn_x = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call add_fld_to_xfldbun(xgrid, FBocn_x, 'So_t', aoflux%tocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBocn_x, 'So_u', aoflux%uocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBocn_x, 'So_v', aoflux%vocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(aoflux%tocn)
    aoflux%lsize = lsize
    if (flds_wiso) then
       call add_fld_to_xfldbun(xgrid, FBocn_x, 'So_roce_16O', aoflux%roce_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBocn_x, 'So_roce_18O', aoflux%roce_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBocn_x, 'So_roce_HDO', aoflux%roce_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%roce_16O(lsize)); aoflux%roce_16O(:) = 0._R8
       allocate(aoflux%roce_18O(lsize)); aoflux%roce_18O(:) = 0._R8
       allocate(aoflux%roce_HDO(lsize)); aoflux%roce_HDO(:) = 0._R8
    end if

    call ESMF_FieldBundleGet(FBocn_x, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fldnames_ocn_in(fieldcount))
    call ESMF_FieldBundleGet(FBocn_x, fieldnamelist=fldnames_ocn_in, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! input fields from atm on xgrid
    ! ------------------------

    FBatm_x = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_z', aoflux%zbot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (trim(coupling_mode) == 'nems_orig_data' .and. ocn_surface_flux_scheme == -1) then
       ! bulk formula quantities for nems_orig_data
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_u10m', aoflux%ubot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_v10m', aoflux%vbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_t2m' , aoflux%tbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_q2m' , aoflux%shum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_u'   , aoflux%ubot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_v'   , aoflux%vbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_tbot', aoflux%tbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_shum', aoflux%shum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! bottom level potential temperature will need to be computed if not received from the atm
    if (FB_fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc)) then
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_ptem', aoflux%thbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_thbot = .false.
    else
       allocate(aoflux%thbot(lsize))
       compute_atm_thbot = .true.
    end if

    ! bottom level density will need to be computed if not received from the atm
    if (FB_fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_dens', rc=rc)) then
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_dens', aoflux%dens, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_dens = .false.
    else
       compute_atm_dens = .true.
       allocate(aoflux%dens(lsize))
    end if

    ! if either density or potential temperature are computed, will need bottom level pressure
    if (compute_atm_dens .or. compute_atm_thbot) then
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_pbot', aoflux%pbot, rc=rc)
    end if

    if (flds_wiso) then
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_shum_16O', aoflux%shum_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_shum_18O', aoflux%shum_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call add_fld_to_xfldbun(xgrid, FBatm_x, 'Sa_shum_HDO', aoflux%shum_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%shum_16O(lsize)); aoflux%shum_16O(:) = 0._R8
       allocate(aoflux%shum_18O(lsize)); aoflux%shum_18O(:) = 0._R8
       allocate(aoflux%shum_HDO(lsize)); aoflux%shum_HDO(:) = 0._R8
    end if

    call ESMF_FieldBundleGet(FBatm_x, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fldnames_atm_in(fieldcount))
    call ESMF_FieldBundleGet(FBatm_x, fieldnamelist=fldnames_atm_in, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! output fields from aoflux calculation on exchange grid
    ! ------------------------

    FBaof_x = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call add_fld_to_xfldbun(xgrid, FBaof_x, 'So_tref'  , aoflux%tref   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'So_qref'  , aoflux%qref   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'So_ustar' , aoflux%ustar  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'So_re'    , aoflux%re     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'So_ssq'   , aoflux%ssq    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'So_u10'   , aoflux%u10    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'So_duu10n', aoflux%duu10n , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_taux', aoflux%taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_tauy', aoflux%tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_lat' , aoflux%lat    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_sen' , aoflux%sen    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_evap', aoflux%evap   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_lwup', aoflux%lwup   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_evap_16O', aoflux%evap_16O, rc=rc)
       call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_evap_18O', aoflux%evap_18O, rc=rc)
       call add_fld_to_xfldbun(xgrid, FBaof_x, 'Faox_evap_HDO', aoflux%evap_HDO, rc=rc)
    else
       allocate(aoflux%evap_16O(lsize)); aoflux%evap_16O(:) = 0._R8
       allocate(aoflux%evap_18O(lsize)); aoflux%evap_18O(:) = 0._R8
       allocate(aoflux%evap_HDO(lsize)); aoflux%evap_HDO(:) = 0._R8
    end if

    call ESMF_FieldBundleGet(FBaof_x, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fldnames_aof_out(fieldcount))
    call ESMF_FieldBundleGet(FBaof_x, fieldnamelist=fldnames_aof_out, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! create the routehandles atm->xgrid and xgrid->atm
    ! ------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), trim(fldnames_atm_in(1)), field=lfield_a, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBatm_x, trim(fldnames_atm_in(1)), field=lfield_x, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(xgrid, lfield_a, lfield_x, routehandle=rh_agrid2xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(xgrid, lfield_x, lfield_a, routehandle=rh_xgrid2agrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! create the routehandles ocn->xgrid and xgrid->ocn
    ! ------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), trim(fldnames_ocn_in(1)), field=lfield_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBocn_x, trim(fldnames_ocn_in(1)), field=lfield_x, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(xgrid, lfield_o, lfield_x, routehandle=rh_ogrid2xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(xgrid, lfield_x, lfield_o, routehandle=rh_xgrid2ogrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! setup the compute mask - default compute everywhere for exchange grid
    ! ------------------------

    allocate(aoflux%mask(lsize))
    aoflux%mask(:) = 1

    ! ------------------------
    ! Determine one normalization field for ocn->xgrid
    ! ------------------------

    ! Create temporary source field on ocn mesh and set its value to 1.
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), 'So_t', field=lfield_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield_o, mesh=ocn_mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield_o = ESMF_FieldCreate(ocn_mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield_o, farrayptr=dataPtr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr1d(:) = 1.0_R8

    ! Create field_ogrid2xgrid_normone (module variable)
    field_ogrid2xgrid_normone = ESMF_FieldCreate(xgrid, ESMF_TYPEKIND_R8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegrid(lfield_o, field_ogrid2xgrid_normone, routehandle=rh_ogrid2xgrid, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Destroy temporary field
    call ESMF_FieldDestroy(lfield_o, rc=rc, noGarbage=.true.)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! Determine one normalization field for xgrid->atm
    ! ------------------------
    ! Create temporary field on xgrid and set its value to 1.
    lfield_x = ESMF_FieldCreate(xgrid, typekind=ESMF_TYPEKIND_R8, name='Sa_z', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield_x, farrayptr=dataPtr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr1d(:) = 1.0_R8

    ! Create field_xgrid2agrid_normone (module variable) - on the atm mesh
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), 'Sa_z', field=lfield_a, rc=rc)    
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield_a, mesh=atm_mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_xgrid2agrid_normone = ESMF_FieldCreate(atm_mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegrid(lfield_x, field_xgrid2agrid_normone, routehandle=rh_xgrid2agrid, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Destroy temporary field on xgrid
    call ESMF_FieldDestroy(lfield_x, rc=rc, noGarbage=.true.)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_aofluxes_init_xgrid

  !===============================================================================
  subroutine med_aofluxes_update(gcomp, aoflux, rc)

    !-----------------------------------------------------------------------
    ! Determine atm/ocn fluxes eother on atm, ocn or exchange grid
    ! The module arrays are set via pointers to the mediator internal states
    ! in med_ocnatm_init and are used below.
    !  1) Create input on aoflux grid
    !  2) Update atmosphere/ocean surface fluxes
    !  3) Map aoflux output to relevant atm/ocn grid(s)
    !-----------------------------------------------------------------------

    use ESMF          , only : ESMF_GridComp
    use ESMF          , only : ESMF_LogWrite, ESMF_LogMsg_Info, ESMF_SUCCESS
    use med_map_mod   , only : med_map_field_packed
    use shr_flux_mod  , only : shr_flux_atmocn

    ! Arguments
    type(ESMF_GridComp)               :: gcomp
    type(aoflux_type) , intent(inout) :: aoflux
    integer           , intent(out)   :: rc
    !
    ! Local variables
    type(InternalState) :: is_local
    type(ESMF_Field)    :: field_src
    type(ESMF_Field)    :: field_dst
    integer             :: n,i,nf                     ! indices
    real(r8), pointer   :: data_normdst(:)
    real(r8), pointer   :: data_dst(:)
    character(*),parameter  :: subName = '(med_aofluxes_update) '
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

       ! Map input ocn to agrid
       do nf = 1,size(fldnames_ocn_in)
          ! Create source field
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), fldnames_ocn_in(nf), field=field_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Create destination field
          call ESMF_FieldBundleGet(FBocn_a, fldnames_ocn_in(nf), field=field_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Map ocn->atm conservatively without fractions
          call ESMF_FieldRegrid(field_src, field_dst, &
               routehandle=is_local%wrap%RH(compocn,compatm, mapconsd), &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)

          ! One normalization
          call ESMF_FieldGet(is_local%wrap%field_normOne(compocn,compatm,mapconsd), farrayPtr=data_normdst, rc=rc)
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
       end do


    else if (is_local%wrap%aoflux_grid == 'xgrid') then

       ! Map input atm to xgrid
       do nf = 1,size(fldnames_atm_in)
          ! Get the source field
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), fldnames_atm_in(nf), field=field_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Get the destination field
          call ESMF_FieldBundleGet(FBatm_x, fldnames_atm_in(nf), field=field_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Map atm->xgrid conservatively
          call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_agrid2xgrid, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       end do

       ! map input ocn to xgrid
       do nf = 1,size(fldnames_ocn_in)
          ! Create source field
          call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), fldnames_ocn_in(nf), field=field_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Create destination field
          call ESMF_FieldBundleGet(FBocn_x, fldnames_ocn_in(nf), field=field_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Map ocn->xgrid conservatively without fractions
          call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_ogrid2xgrid, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       end do

    end if

    !----------------------------------
    ! Calculate quantities if they are not defined
    !----------------------------------

    ! Note pbot, tbot and shum have already been mapped or are available on the aoflux grid
    if (compute_atm_thbot) then
       do n = 1,aoflux%lsize
          if (aoflux%mask(n) /= 0._r8) then
             aoflux%thbot(n) = aoflux%tbot(n)*((100000._R8/aoflux%pbot(n))**0.286_R8)
          end if
       end do
    end if
    if (compute_atm_dens) then
       do n = 1,aoflux%lsize
          if (aoflux%mask(n) /= 0._r8) then
             aoflux%dens(n) = aoflux%pbot(n)/(287.058_R8*(1._R8 + 0.608_R8*aoflux%shum(n))*aoflux%tbot(n))
          end if
       end do
    end if

    !----------------------------------
    ! Update atmosphere/ocean surface fluxes
    !----------------------------------

    call shr_flux_atmocn (&
         nMax=aoflux%lsize, zbot=aoflux%zbot, ubot=aoflux%ubot, vbot=aoflux%vbot, thbot=aoflux%thbot, &
         qbot=aoflux%shum, s16O=aoflux%shum_16O, sHDO=aoflux%shum_HDO, s18O=aoflux%shum_18O, rbot=aoflux%dens, &
         tbot=aoflux%tbot, us=aoflux%uocn, vs=aoflux%vocn, &
         ts=aoflux%tocn, mask=aoflux%mask, seq_flux_atmocn_minwind=0.5_r8, &
         sen=aoflux%sen, lat=aoflux%lat, lwup=aoflux%lwup, &
         r16O=aoflux%roce_16O, rhdo=aoflux%roce_HDO, r18O=aoflux%roce_18O, &
         evap=aoflux%evap, evap_16O=aoflux%evap_16O, evap_HDO=aoflux%evap_HDO, evap_18O=aoflux%evap_18O, &
         taux=aoflux%taux, tauy=aoflux%tauy, tref=aoflux%tref, qref=aoflux%qref, &
         ocn_surface_flux_scheme=ocn_surface_flux_scheme, &
         duu10n=aoflux%duu10n, ustar_sv=aoflux%ustar, re_sv=aoflux%re, ssq_sv=aoflux%ssq, &
         missval = 0.0_r8)

    do n = 1,aoflux%lsize
       if (aoflux%mask(n) /= 0) then
          aoflux%u10(n) = sqrt(aoflux%duu10n(n))
       end if
    enddo

    !----------------------------------
    ! map aoflux to output to relevant atm/ocn grid(s)
    !----------------------------------

    if (is_local%wrap%aoflux_grid == 'ogrid') then

       ! mapping aoflux from ogrid to agrid is done in med_phases_prep_atm using updated ocean fractions
       ! on the atm grid

    else if (is_local%wrap%aoflux_grid == 'agrid') then

       if (is_local%wrap%med_coupling_active(compatm,compocn)) then
          ! map aoflux from agrid to ogrid
          do nf = 1,size(fldnames_aof_out)
             ! Create source field
             call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_a, fldnames_aof_out(nf), field=field_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             ! Create destination field
             call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fldnames_aof_out(nf), field=field_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Map atm->ocn conservatively WITHOUT fractions
             call ESMF_FieldRegrid(field_src, field_dst, &
                                !routehandle=is_local%wrap%RH(compatm, compocn, mapconsd), &
                  routehandle=is_local%wrap%RH(compatm, compocn, mapconsf), &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)

             ! One normalization
             ! call ESMF_FieldGet(is_local%wrap%field_normOne(compatm,compocn,mapconsd), farrayPtr=data_normdst, rc=rc)
             ! if (chkerr(rc,__LINE__,u_FILE_u)) return
             ! call ESMF_FieldGet(field_dst, farrayptr=data_dst, rc=rc)
             ! if (chkerr(rc,__LINE__,u_FILE_u)) return
             ! do n = 1,size(data_dst)
             !    if (data_normdst(n) == 0.0_r8) then
             !       data_dst(n) = 0.0_r8
             !    else
             !       data_dst(n) = data_dst(n)/data_normdst(n)
             !    end if
             ! end do
          end do
       end if

    else if (is_local%wrap%aoflux_grid == 'xgrid') then

       ! map aoflux from xgrid to agrid and ogrid
       do nf = 1,size(fldnames_aof_out)
          ! Get the source field
          call ESMF_FieldBundleGet(FBaof_x, fldnames_aof_out(nf), field=field_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Get the destination field
          call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_a, fldnames_aof_out(nf), field=field_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Map xgrid->agrid conservatively
          call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_xgrid2agrid, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
          ! One normalization
          call ESMF_FieldGet(field_xgrid2agrid_normone, farrayPtr=data_normdst, rc=rc)
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

          ! Get the destination field
          call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fldnames_aof_out(nf), field=field_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Map xgrid->ogrid conservatively
          call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_xgrid2agrid, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       end do

    end if

    call t_stopf('MED:'//subname)

  end subroutine med_aofluxes_update

  !================================================================================
  subroutine add_fld_to_xfldbun(xgrid, fldbun, fldname, aoflux_dataptr, rc)

    ! input/output variables
    type(ESMF_Xgrid)       , intent(in)    :: xgrid
    type(ESMF_FieldBundle) , intent(inout) :: fldbun
    character(len=*)       , intent(in)    :: fldname
    real(r8)               , pointer       :: aoflux_dataptr(:)
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    !-----------------------------------------------------------------------
    rc = ESMF_SUCCESS

    lfield = ESMF_FieldCreate(xgrid, typekind=ESMF_TYPEKIND_R8, name=trim(fldname), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(fldbun, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=aoflux_dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine add_fld_to_xfldbun

  !================================================================================
  subroutine fldbun_getfldptr(fldbun, fldname, fldptr, rc)

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)  :: fldbun
    character(len=*)       , intent(in)  :: fldname
    real(r8)               , pointer     :: fldptr(:)
    integer                , intent(out) :: rc 

    ! local variables
    type(ESMF_Field) :: lfield
    !-----------------------------------------------------------------------
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fldbun, trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine fldbun_getfldptr

end module med_phases_aofluxes_mod
