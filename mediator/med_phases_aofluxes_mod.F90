module med_phases_aofluxes_mod

  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF                  , only : ESMF_RouteHandle, ESMF_FieldRegrid, ESMF_FieldRegridStore
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet, ESMF_XGrid, ESMF_XGridCreate, ESMF_TYPEKIND_R8
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
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
  use esmFlds               , only : compatm, compocn, coupling_mode
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
  character(len=CS), allocatable :: fldnames_aoflux_out(:)

  ! following is needed for atm/ocn fluxes on atm grid
  type(ESMF_FieldBundle) :: FBocn_a ! ocean fields need for aoflux calc on atm grid

  ! following is needed for atm/ocn fluxes on the exchange grid
  type(ESMF_FieldBundle) :: FBocn_x        ! input ocn fields on exchange grid
  type(ESMF_FieldBundle) :: FBatm_x        ! input atm fields on exchange grid
  type(ESMF_FieldBundle) :: FBaoflux_x     ! output atm/ocn fluxes on exchange grid
  type(ESMF_RouteHandle) :: rh_ogrid2xgrid ! ocn->xgrid mapping
  type(ESMF_RouteHandle) :: rh_xgrid2ogrid ! xgrid->ocn mapping
  type(ESMF_RouteHandle) :: rh_agrid2xgrid ! atm->xgrid mapping
  type(ESMF_RouteHandle) :: rh_xgrid2agrid ! xgrid->atm mapping

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

     ! input: mask on aoflux grid
     integer  , pointer :: mask        (:) => null() ! ocn domain mask: 0 <=> inactive cell
     real(R8) , pointer :: rmask       (:) => null() ! ocn domain mask: 0 <=> inactive cell

     ! output: on aoflux grid
     ! if aoflux grid is ocn - then need to map these to the atm
     ! if aoflux grid is atm - then need to map these to the ocn
     ! if aoflux grid is exchange - will map back to both the atm and ocn
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
    integer             :: fieldcount
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

    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_wiso
    else
       flds_wiso = .false.
    end if

    !----------------------------------
    ! Determine input ocn field names in aoflux
    !----------------------------------

    if (flds_wiso) then
       allocate(fldnames_ocn_in(5))
       fldnames_ocn_in = (/'So_omask   ','So_t       ','So_u       ','So_v       ','So_roce_wiso' /)
    else
       allocate(fldnames_ocn_in(4))
       fldnames_ocn_in = (/'So_omask   ' ,'So_t       ','So_u       ','So_v       '/)
    end if

    !----------------------------------
    ! Determine input atm field names in aoflux
    !----------------------------------

    if (trim(coupling_mode) == 'nems_orig_data' .and. ocn_surface_flux_scheme == -1) then
       ! TODO: fix this
       allocate(fldnames_atm_in(8))
       ! fldnames_atm_in = (/'Sa_z    ','Sa_u10m ','Sa_v10m ','Sa_t2m  ',&
       !                    'Sa_q2m  ','Sa_pbot ','Sa_dens ','Sa_ptem'/)
       ! fldnames_atm_in = (/'Sa_u10m','Sa_v10m','Sa_t2m ','Sa_q2m '/)
    else if (flds_wiso) then
       allocate(fldnames_atm_in(9))
       fldnames_atm_in = (/'Sa_z        ','Sa_u        ','Sa_v        ','Sa_tbot     ',&
                           'Sa_shum     ','Sa_pbot     ','Sa_dens     ','Sa_ptem     ','Sa_shum_wiso'/)
    else
       allocate(fldnames_atm_in(8))
       fldnames_atm_in = (/'Sa_z   ','Sa_u   ','Sa_v   ','Sa_tbot',&
                           'Sa_shum','Sa_pbot','Sa_dens','Sa_ptem'/)
    end if

    !----------------------------------
    ! Determine aoflux output field names (same for ocn and atm grid)
    !----------------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fldnames_aoflux_out(fieldcount))
    call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fieldNameList=fldnames_aoflux_out, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Initialize aoflux
    !----------------------------------

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
    character(len=*),parameter :: subname='(med_aofluxes_init_ocngrid)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
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
    character(len=CX)   :: tmpstr
    integer             :: lsize
    real(R8), pointer   :: ofrac(:) => null()
    real(R8), pointer   :: ifrac(:) => null()
    character(len=*),parameter :: subname='(med_aofluxes_init_atmgrid)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! create field bundles FBocn_a for the above fieldlist
    call FB_init(FBocn_a, is_local%wrap%flds_scalar_name, &
         FBgeom=is_local%wrap%FBImp(compatm,compatm), fieldnamelist=fldnames_ocn_in, name='FBocn_a', rc=rc)

    ! ------------------------
    ! input fields from ocn on atm grid
    ! ------------------------

    ! point directly into FBocn_a for ocean fields on the atm grid
    call FB_GetFldPtr(FBocn_a, fldname='So_omask', fldptr1=aoflux%rmask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(aoflux%rmask)
    aoflux%lsize = lsize
    call FB_GetFldPtr(FBocn_a, fldname='So_t', fldptr1=aoflux%tocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
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
    ! input fields from atm on atm grid
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

    ! setup the compute mask. - default compute everywhere, then
    ! "turn off" gridcells want the ocean mask ocean mask mapped to
    ! atm grid, but do not have access to the ocean mask mapped to
    ! the atm grid.  on the atm grid, the mask is all oness it's
    ! just all 1's so not very useful.  next look at ofrac+ifrac in
    ! fractions. want to compute on all non-land points.  using
    ! ofrac alone will exclude points that are currently all sea
    ! ice but that later could be less that 100% covered in ice.

    allocate(aoflux%mask(lsize))
    aoflux%mask(:) = 1

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux%rmask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : maskA= "//trim(tmpstr), ESMF_LOGMSG_INFO)

    where (aoflux%rmask(:) == 0._R8) aoflux%mask(:) = 0   ! like nint

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux%rmask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : maskB= "//trim(tmpstr), ESMF_LOGMSG_INFO)

    call FB_getFldPtr(is_local%wrap%FBFrac(compatm) , fldname='ofrac' , fldptr1=ofrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBFrac(compatm) , fldname='ifrac' , fldptr1=ifrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    where (ofrac(:) + ifrac(:) <= 0.0_R8) aoflux%mask(:) = 0

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
    character(len=CX)    :: tmpstr
    type(ESMF_Field)     :: lfield_a
    type(ESMF_Field)     :: lfield_o
    type(ESMF_Field)     :: lfield_x
    type(ESMF_Field)     :: lfield
    integer              :: elementCount
    type(ESMF_Mesh)      :: ocn_mesh
    type(ESMF_Mesh)      :: atm_mesh
    integer, allocatable :: ocn_mask(:)
    type(ESMF_XGrid)     :: aoflux_xgrid 
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

    ! determine atm mesh - assume that atm mask is always 1!
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), fieldname=fldnames_atm_in(1), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=atm_mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! determine ocn mesh and ocn mask
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), fieldname=fldnames_ocn_in(1), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=ocn_mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(ocn_mesh, elementCount=elementCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ocn_mask(elementCount))
    call ESMF_MeshGet(ocn_mesh, elementMask=ocn_mask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! create the exchange grid -
    ! first reset ocean mask to be 1 minus ocn mask - to determine what points to mask out for xgrid creation
    do n = 1,size(ocn_mask)
       ocn_mask(n) = 1 - ocn_mask(n)
    end do
    aoflux_xgrid = ESMF_XGridCreate(sideAMesh=(/ocn_mesh/), sideBMesh=(/atm_mesh/), sideAMaskValues=ocn_mask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! input fields from ocn on xchange grid
    ! ------------------------

    ! Create the FBocn_x
    FBocn_x = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,size(fldnames_ocn_in)
       lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name=trim(fldnames_ocn_in(n)), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBocn_x, (/lfield_x/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

    call FB_GetFldPtr(FBocn_x, fldname='So_omask', fldptr1=aoflux%rmask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lsize = size(aoflux%rmask)
    aoflux%lsize = lsize

    call FB_GetFldPtr(FBocn_x, fldname='So_t', fldptr1=aoflux%tocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBocn_x, fldname='So_u', fldptr1=aoflux%uocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBocn_x, fldname='So_v', fldptr1=aoflux%vocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call FB_GetFldPtr(FBocn_x, fldname='So_roce_16O', fldptr1=aoflux%roce_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBocn_x, fldname='So_roce_18O', fldptr1=aoflux%roce_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBocn_x, fldname='So_roce_HDO', fldptr1=aoflux%roce_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%roce_16O(lsize)); aoflux%roce_16O(:) = 0._R8
       allocate(aoflux%roce_18O(lsize)); aoflux%roce_18O(:) = 0._R8
       allocate(aoflux%roce_HDO(lsize)); aoflux%roce_HDO(:) = 0._R8
    end if

    ! ------------------------
    ! input fields from atm on exchange grid
    ! ------------------------

    ! Note - the mapping from FBImp(compatm,compatm) to FBatm_x is done in med_phase
    ! med_phases_postatm

    ! create FBatm_x
    FBatm_x = ESMF_FieldBundleCreate(rc=rc)
    do n = 1,size(fldnames_atm_in)
       lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name=trim(fldnames_atm_in(n)), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBatm_x, (/lfield_x/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

    call FB_GetFldPtr(FBatm_x, fldname='Sa_z', fldptr1=aoflux%zbot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! bulk formula quantities for nems_orig_data
    if (trim(coupling_mode) == 'nems_orig_data' .and. ocn_surface_flux_scheme == -1) then
       call FB_GetFldPtr(FBatm_x, fldname='Sa_u10m', fldptr1=aoflux%ubot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBatm_x, fldname='Sa_v10m', fldptr1=aoflux%vbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBatm_x, fldname='Sa_t2m', fldptr1=aoflux%tbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBatm_x, fldname='Sa_q2m', fldptr1=aoflux%shum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call FB_GetFldPtr(FBatm_x, fldname='Sa_u', fldptr1=aoflux%ubot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBatm_x, fldname='Sa_v', fldptr1=aoflux%vbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBatm_x, fldname='Sa_tbot', fldptr1=aoflux%tbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBatm_x, fldname='Sa_shum', fldptr1=aoflux%shum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! bottom level potential temperature will need to be computed if not received from the atm
    if (FB_fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc)) then
       call FB_GetFldPtr(FBatm_x, fldname='Sa_ptem', fldptr1=aoflux%thbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_thbot = .false.
    else
       allocate(aoflux%thbot(lsize))
       compute_atm_thbot = .true.
    end if

    ! bottom level density will need to be computed if not received from the atm
    if (FB_fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_dens', rc=rc)) then
       call FB_GetFldPtr(FBatm_x, fldname='Sa_dens', fldptr1=aoflux%dens, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_dens = .false.
    else
       compute_atm_dens = .true.
       allocate(aoflux%dens(lsize))
    end if

    ! if either density or potential temperature are computed, will need bottom level pressure
    if (compute_atm_dens .or. compute_atm_thbot) then
       call FB_GetFldPtr(FBatm_x, fldname='Sa_pbot', fldptr1=aoflux%pbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (flds_wiso) then
       call FB_GetFldPtr(FBatm_x, fldname='Sa_shum_16O', fldptr1=aoflux%shum_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBatm_x, fldname='Sa_shum_18O', fldptr1=aoflux%shum_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBatm_x, fldname='Sa_shum_HDO', fldptr1=aoflux%shum_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%shum_16O(lsize)); aoflux%shum_16O(:) = 0._R8
       allocate(aoflux%shum_18O(lsize)); aoflux%shum_18O(:) = 0._R8
       allocate(aoflux%shum_HDO(lsize)); aoflux%shum_HDO(:) = 0._R8
    end if

    ! ------------------------
    ! output fields from aoflux calculation on exchange grid
    ! ------------------------

    ! create FBaoflux_x and then set pointers into each field added

    FBaoflux_x = ESMF_FieldBundleCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='So_tref', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='So_tref', fldptr1=aoflux%tref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='So_qref', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='So_qref', fldptr1=aoflux%qref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='So_ustar', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='So_ustar', fldptr1=aoflux%ustar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='So_re', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='So_re', fldptr1=aoflux%re, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='So_ssq', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='So_ssq', fldptr1=aoflux%ssq, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='So_u10', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='So_u10', fldptr1=aoflux%u10, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='So_duu10n', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='So_duu10n', fldptr1=aoflux%duu10n, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_taux', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='Faox_taux', fldptr1=aoflux%taux, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_tauy', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='Faox_tauy', fldptr1=aoflux%tauy, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_lat', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='Faox_lat', fldptr1=aoflux%lat, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_sen', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='Faox_sen', fldptr1=aoflux%sen, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_evap', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='Faox_evap', fldptr1=aoflux%evap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_evap_16O', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBaoflux_x, fldname='Faox_evap_16O', fldptr1=aoflux%evap_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_evap_18O', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBaoflux_x, fldname='Faox_evap_18O', fldptr1=aoflux%evap_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_evap_HDO', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBaoflux_x, fldname='Faox_evap_HDO', fldptr1=aoflux%evap_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%evap_16O(lsize)); aoflux%evap_16O(:) = 0._R8
       allocate(aoflux%evap_18O(lsize)); aoflux%evap_18O(:) = 0._R8
       allocate(aoflux%evap_HDO(lsize)); aoflux%evap_HDO(:) = 0._R8
    end if

    lfield_x = ESMF_FieldCreate(aoflux_xgrid, typekind=ESMF_TYPEKIND_R8, name='Faox_lwup', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBaoflux_x, (/lfield_x/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBaoflux_x, fldname='Faox_lwup', fldptr1=aoflux%lwup, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! setup the compute mask - default compute everywhere for exchange grid
    ! ------------------------

    allocate(aoflux%mask(lsize))
    aoflux%mask(:) = 1

    ! ------------------------
    ! create the routehandles atm->xgrid and xgrid->atm
    ! ------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compatm,compatm), trim(fldnames_atm_in(1)), field=lfield_a, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBatm_x, trim(fldnames_atm_in(1)), field=lfield_x, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(aoflux_xgrid, lfield_a, lfield_x, routehandle=rh_agrid2xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(aoflux_xgrid, lfield_x, lfield_a, routehandle=rh_xgrid2agrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------
    ! create the routehandles ocn->xgrid and xgrid->ocn
    ! ------------------------

    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), trim(fldnames_ocn_in(1)), field=lfield_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBocn_x, trim(fldnames_ocn_in(1)), field=lfield_x, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(aoflux_xgrid, lfield_o, lfield_x, routehandle=rh_ogrid2xgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridStore(aoflux_xgrid, lfield_x, lfield_o, routehandle=rh_xgrid2ogrid, rc=rc)
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
    use ESMF          , only : ESMF_TERMORDER_SRCSEQ, ESMF_REGION_TOTAL
    use esmFlds       , only : mapconsd
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
               routehandle=is_local%wrap%RH(compocn,compatm,mapconsd), &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)

          ! One normalization
          call ESMF_FieldGet(is_local%wrap%field_normOne(compatm,compatm,mapconsd), farrayPtr=data_normdst, rc=rc)
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
    call t_stopf('MED:'//subname)

    !----------------------------------
    ! map aoflux to output to relevant atm/ocn grid(s)
    !----------------------------------

    if (is_local%wrap%aoflux_grid == 'ogrid') then

       ! map aoflux from ogrid to agrid
       ! fluxes are computed on the ocean grid - this is set up from the data in esmFldsExchange_xxx
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBMed_aoflux_o, &
            FBDst=is_local%wrap%FBMed_aoflux_a, &
            FBFracSrc=is_local%wrap%FBFrac(compocn), &
            field_normOne=is_local%wrap%field_normOne(compocn,compatm,:), &
            packed_data=is_local%wrap%packed_data_aoflux_o2a(:), &
            routehandles=is_local%wrap%RH(compocn,compatm,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    else if (is_local%wrap%aoflux_grid == 'agrid') then

       ! map aoflux from agrid to ogrid
       do nf = 1,size(fldnames_aoflux_out)
          ! Create source field
          call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_a, fldnames_aoflux_out(nf), field=field_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Create destination field
          call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fldnames_aoflux_out(nf), field=field_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Map atm->ocn conservatively WITHOUT fractions
          call ESMF_FieldRegrid(field_src, field_dst, &
               routehandle=is_local%wrap%RH(compatm, compocn, mapconsd), &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)

          ! One normalization
          call ESMF_FieldGet(is_local%wrap%field_normOne(compatm,compocn,mapconsd), farrayPtr=data_normdst, rc=rc)
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

       ! map aoflux from xgrid to agrid and ogrid
       do nf = 1,size(fldnames_aoflux_out)
          ! Get the source field
          call ESMF_FieldBundleGet(FBaoflux_x, fldnames_aoflux_out(nf), field=field_src, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Get the destination field
          call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_a, fldnames_aoflux_out(nf), field=field_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Map xgrid->agrid conservatively
          call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_xgrid2agrid, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)

          ! Get the destination field
          call ESMF_FieldBundleGet(is_local%wrap%FBMed_aoflux_o, fldnames_aoflux_out(nf), field=field_dst, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ! Map xgrid->ogrid conservatively
          call ESMF_FieldRegrid(field_src, field_dst, routehandle=rh_xgrid2agrid, &
               termorderflag=ESMF_TERMORDER_SRCSEQ, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       end do

    end if

  end subroutine med_aofluxes_update

end module med_phases_aofluxes_mod
