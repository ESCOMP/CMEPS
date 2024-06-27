module MED

  !-----------------------------------------------------------------------------
  ! Mediator Initialization
  !
  ! Note on time management:
  ! Each time loop has its own associated clock object. NUOPC manages
  ! these clock objects, i.e. their creation and destruction, as well as
  ! startTime, endTime, timeStep adjustments during the execution. The
  ! outer most time loop of the run sequence is a special case. It uses
  ! the driver clock itself. If a single outer most loop is defined in
  ! the run sequence provided by freeFormat, this loop becomes the driver
  ! loop level directly. Therefore, setting the timeStep or runDuration
  ! for the outer most time loop results in modifying the driver clock
  ! itself. However, for cases with cocnatenated loops on the upper level
  ! of the run sequence in freeFormat, a single outer loop is added
  ! automatically during ingestion, and the driver clock is used for this
  ! loop instead.
  !-----------------------------------------------------------------------------

  use ESMF                     , only : ESMF_VMLogMemInfo
  use NUOPC_Model              , only : SetVM
  use med_kind_mod             , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod        , only : dbug_flag          => med_constants_dbug_flag
  use med_constants_mod        , only : spval_init         => med_constants_spval_init
  use med_constants_mod        , only : spval              => med_constants_spval
  use med_constants_mod        , only : czero              => med_constants_czero
  use med_utils_mod            , only : chkerr             => med_utils_ChkErr
  use med_methods_mod          , only : Field_GeomPrint    => med_methods_Field_GeomPrint
  use med_methods_mod          , only : State_GeomPrint    => med_methods_State_GeomPrint
  use med_methods_mod          , only : State_reset        => med_methods_State_reset
  use med_methods_mod          , only : State_getNumFields => med_methods_State_getNumFields
  use med_methods_mod          , only : State_GetScalar    => med_methods_State_GetScalar
  use med_methods_mod          , only : FB_Init            => med_methods_FB_init
  use med_methods_mod          , only : FB_Init_pointer    => med_methods_FB_Init_pointer
  use med_methods_mod          , only : FB_Reset           => med_methods_FB_Reset
  use med_methods_mod          , only : FB_diagnose        => med_methods_FB_diagnose
  use med_methods_mod          , only : FB_getFieldN       => med_methods_FB_getFieldN
  use med_methods_mod          , only : clock_timeprint    => med_methods_clock_timeprint
  use med_utils_mod            , only : memcheck           => med_memcheck
  use med_time_mod             , only : med_time_alarmInit
  use med_internalstate_mod    , only : InternalState, med_internalstate_init, med_internalstate_coupling
  use med_internalstate_mod    , only : med_internalstate_defaultmasks, logunit, maintask
  use med_internalstate_mod    , only : ncomps, compname
  use med_internalstate_mod    , only : compmed, compatm, compocn, compice, complnd, comprof, compwav, compglc
  use med_internalstate_mod    , only : coupling_mode, aoflux_code, aoflux_ccpp_suite
  use esmFlds                  , only : med_fldList_GetocnalbfldList, med_fldList_type
  use esmFlds                  , only : med_fldList_GetNumFlds, med_fldList_GetFldNames, med_fldList_GetFldInfo
  use esmFlds                  , only : med_fldList_Document_Mapping, med_fldList_Document_Merging
  use esmFlds                  , only : med_fldList_GetfldListFr, med_fldList_GetfldListTo, med_fldList_Realize
  use esmFldsExchange_ufs_mod  , only : esmFldsExchange_ufs
  use esmFldsExchange_cesm_mod , only : esmFldsExchange_cesm
  use esmFldsExchange_hafs_mod , only : esmFldsExchange_hafs
  use med_phases_profile_mod   , only : med_phases_profile_finalize

  implicit none
  private

  public  SetServices
  public  SetVM
  private InitializeP0
  private AdvertiseFields ! advertise fields
  private RealizeFieldsWithTransferProvided ! realize connected Fields with transfer action "provide"
  private ModifyDecompofMesh ! optionally modify the decomp/distr of transferred Grid/Mesh
  private RealizeFieldsWithTransferAccept ! realize all Fields with transfer action "accept"
  private DataInitialize     ! finish initialization and resolve data dependencies
  private SetRunClock
  private med_meshinfo_create
  private med_grid_write
  private med_finalize

  character(len=*), parameter :: u_FILE_u  = &
       __FILE__

  logical :: profile_memory = .false.

  logical, allocatable :: compDone(:) ! component done flag

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)

    use ESMF                    , only: ESMF_SUCCESS, ESMF_GridCompSetEntryPoint
    use ESMF                    , only: ESMF_METHOD_INITIALIZE, ESMF_METHOD_RUN
    use ESMF                    , only: ESMF_GridComp, ESMF_MethodRemove
    use NUOPC                   , only: NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize, NUOPC_NoOP
    use NUOPC_Mediator          , only: mediator_routine_SS             => SetServices
    use NUOPC_Mediator          , only: mediator_routine_Run            => routine_Run
    use NUOPC_Mediator          , only: mediator_label_DataInitialize   => label_DataInitialize
    use NUOPC_Mediator          , only: mediator_label_Advance          => label_Advance
    use NUOPC_Mediator          , only: mediator_label_CheckImport      => label_CheckImport
    use NUOPC_Mediator          , only: mediator_label_TimestampExport  => label_TimestampExport
    use NUOPC_Mediator          , only: mediator_label_SetRunClock      => label_SetRunClock
    use NUOPC_Mediator          , only: mediator_label_Finalize         => label_Finalize
    use med_phases_history_mod  , only: med_phases_history_write
    use med_phases_restart_mod  , only: med_phases_restart_write
    use med_phases_prep_atm_mod , only: med_phases_prep_atm
    use med_phases_prep_ice_mod , only: med_phases_prep_ice
    use med_phases_prep_lnd_mod , only: med_phases_prep_lnd
    use med_phases_prep_wav_mod , only: med_phases_prep_wav_accum
    use med_phases_prep_wav_mod , only: med_phases_prep_wav_avg
    use med_phases_prep_glc_mod , only: med_phases_prep_glc
    use med_phases_prep_rof_mod , only: med_phases_prep_rof
    use med_phases_prep_ocn_mod , only: med_phases_prep_ocn_accum
    use med_phases_prep_ocn_mod , only: med_phases_prep_ocn_avg
    use med_phases_post_atm_mod , only: med_phases_post_atm
    use med_phases_post_ice_mod , only: med_phases_post_ice
    use med_phases_post_lnd_mod , only: med_phases_post_lnd
    use med_phases_post_glc_mod , only: med_phases_post_glc
    use med_phases_post_ocn_mod , only: med_phases_post_ocn
    use med_phases_post_rof_mod , only: med_phases_post_rof
    use med_phases_post_wav_mod , only: med_phases_post_wav
    use med_phases_ocnalb_mod   , only: med_phases_ocnalb_run
    use med_phases_aofluxes_mod , only: med_phases_aofluxes_run
    use med_diag_mod            , only: med_phases_diag_accum, med_phases_diag_print
    use med_diag_mod            , only: med_phases_diag_atm
    use med_diag_mod            , only: med_phases_diag_lnd
    use med_diag_mod            , only: med_phases_diag_rof
    use med_diag_mod            , only: med_phases_diag_glc
    use med_diag_mod            , only: med_phases_diag_ocn
    use med_diag_mod            , only: med_phases_diag_ice_ice2med, med_phases_diag_ice_med2ice
    use med_fraction_mod        , only: med_fraction_init, med_fraction_set
    use med_phases_profile_mod  , only: med_phases_profile
#ifdef CDEPS_INLINE
    use med_phases_cdeps_mod    , only: med_phases_cdeps_run
#endif

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname = '('//__FILE__//':SetServices)'
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    !------------------
    ! the NUOPC model component mediator_routine_SS will register the generic methods
    !------------------

    call NUOPC_CompDerive(gcomp, mediator_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! set entry point for methods that require specific implementation
    ! Provide InitializeP0 to switch from default IPDv00 to IPDv03
    !------------------

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p1: advertise Fields
    !------------------

    ! Mediator advertises its import and export Fields and sets the TransferOfferGeomObject Attribute.
    ! The TransferOfferGeomObject is a String value indicating a component's
    ! intention to transfer the underlying Grid or Mesh on which an advertised Field object is defined.
    ! The valid values are: [will provide, can provide, cannot provide]

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p1"/), userRoutine=AdvertiseFields, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p3: realize connected Fields with transfer action "provide"
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p3"/), userRoutine=RealizeFieldsWithTransferProvided, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p4: optionally modify the decomp/distr of transferred Grid/Mesh
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p4"/), userRoutine=ModifyDecompofMesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p5: realize all Fields with transfer action "accept"
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p5"/), userRoutine=RealizeFieldsWithTransferAccept, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! attach specializing method for DataInitialize
    !------------------

    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! setup mediator history phases for all output variables
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_history_write"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_history_write", specRoutine=med_phases_history_write, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! setup mediator restart phase
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_restart_write"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_restart_write", specRoutine=med_phases_restart_write, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_restart_write", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! setup mediator profile phase
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_profile"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_profile", specRoutine=med_phases_profile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_profile", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep and post routines for atm
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_atm"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_atm", specRoutine=med_phases_prep_atm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_post_atm"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_post_atm", specRoutine=med_phases_post_atm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep and post routines for ocn
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_ocn_accum"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_ocn_accum", specRoutine=med_phases_prep_ocn_accum, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_prep_ocn_accum", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_ocn_avg"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_ocn_avg", specRoutine=med_phases_prep_ocn_avg, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_post_ocn"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_post_ocn", specRoutine=med_phases_post_ocn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep and post routines for ice
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_ice"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_ice", specRoutine=med_phases_prep_ice, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! note that med_fraction_set is now called from post_ice
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_post_ice"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_post_ice", specRoutine=med_phases_post_ice, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep/post routines for lnd
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_lnd"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_lnd", specRoutine=med_phases_prep_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_post_lnd"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_post_lnd", specRoutine=med_phases_post_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep/post routines for rof
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_rof"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_rof", specRoutine=med_phases_prep_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! post routine for rof (mapping to lnd, ocn, ice)
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_post_rof"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_post_rof", specRoutine=med_phases_post_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep/post routines for wav
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_wav_accum"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_wav_accum", specRoutine=med_phases_prep_wav_accum, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_prep_wav_accum", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_wav_avg"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_wav_avg", specRoutine=med_phases_prep_wav_avg, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_post_wav"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_post_wav", specRoutine=med_phases_post_wav, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep/post routines for glc
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_glc"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_glc", specRoutine=med_phases_prep_glc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! post routine for glc (mapping to lnd, ocn, ice)
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_post_glc"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_post_glc", specRoutine=med_phases_post_glc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! phase routine for ocean albedo computation
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_ocnalb_run"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_ocnalb_run", specRoutine=med_phases_ocnalb_run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! phase routine for ocn/atm flux computation
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_aofluxes_run"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_aofluxes_run", specRoutine=med_phases_aofluxes_run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! phase routine for updating fractions
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_fraction_set"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_fraction_set", specRoutine=med_fraction_set, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_fraction_set", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! phase routines for budget diagnostics
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_atm"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaselabel="med_phases_diag_atm", specRoutine=med_phases_diag_atm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_atm", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_lnd"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaselabel="med_phases_diag_lnd", specRoutine=med_phases_diag_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_lnd", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_rof"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaselabel="med_phases_diag_rof", specRoutine=med_phases_diag_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_rof", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_ocn"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaselabel="med_phases_diag_ocn", specRoutine=med_phases_diag_ocn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_ocn", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_glc"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaselabel="med_phases_diag_glc", specRoutine=med_phases_diag_glc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_glc", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_ice_ice2med"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaselabel="med_phases_diag_ice_ice2med", specRoutine=med_phases_diag_ice_ice2med, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_ice_ice2med", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_ice_med2ice"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaselabel="med_phases_diag_ice_med2ice", specRoutine=med_phases_diag_ice_med2ice, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_ice_med2ice", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_accum"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaselabel="med_phases_diag_accum", specRoutine=med_phases_diag_accum, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_accum", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_diag_print"/), userRoutine=mediator_routine_Run, rc=rc)
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_diag_print", specRoutine=med_phases_diag_print, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_diag_print", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

#ifdef CDEPS_INLINE
    !------------------
    ! phase routine for cdeps inline capability
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_cdeps_run"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_cdeps_run", specRoutine=med_phases_cdeps_run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
#endif

    !------------------
    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
    !------------------

    call ESMF_MethodRemove(gcomp, mediator_label_CheckImport, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_CheckImport, specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
    !------------------
    ! This is called every time you enter a mediator phase

    call ESMF_MethodRemove(gcomp, mediator_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_SetRunClock, specRoutine=SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
    !------------------

    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Finalize, &
         specRoutine=med_finalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    use ESMF  , only : ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_VM, ESMF_SUCCESS
    use ESMF  , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_AttributeGet, ESMF_AttributeSet
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_METHOD_INITIALIZE
    use NUOPC , only : NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet
    use med_internalstate_mod, only : maintask, logunit, diagunit
#ifdef CESMCOUPLED
    use nuopc_shr_methods, only : set_component_logging
    use shr_log_mod, only : shr_log_unit
#endif
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    character(len=CL) :: cvalue
    integer           :: localPet
    integer           :: i
    logical           :: isPresent, isSet
    character(len=CX) :: msgString
    character(len=CX) :: diro
    character(len=CX) :: logfile
    character(len=CX) :: diagfile
    character(len=*), parameter :: subname = '('//__FILE__//':InitializeP0)'
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    maintask = .false.
    if (localPet == 0) maintask=.true.

    ! Determine mediator logunit
    if (maintask) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, isPresent=isPresent, isSet=isSet, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (.not. isPresent .and. .not. isSet) then
          diro = './'
       end if
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, isPresent=isPresent, isSet=isSet, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (.not. isPresent .and. .not. isSet) then
          logfile = 'mediator.log'
       end if
#ifdef CESMCOUPLED
       call set_component_logging(gcomp, maintask, logunit, shr_log_unit, rc)
#else
       open(newunit=logunit,file=trim(diro)//"/"//trim(logfile))
#endif
       call NUOPC_CompAttributeGet(gcomp, name="do_budgets", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          if (trim(cvalue) .eq. '.true.') then
             i = index(logfile, '.log')
             diagfile = "diags"//logfile(i:)
             open(newunit=diagunit, file=trim(diro)//"/"//trim(diagfile))
          endif
       end if
    else
       logUnit = 6
    endif

    ! Obtain verbosity level
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=cvalue, defaultValue="max", &
         convention="NUOPC", purpose="Instance", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (maintask) then
       write(logunit,'(a)')trim(subname)//": Mediator verbosity is set to "//trim(cvalue)
    end if

    ! Obtain profiling level
    call NUOPC_CompAttributeGet(gcomp, name="Profiling", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (maintask) then
          write(logunit,'(a)') trim(subname)//": Mediator profiling is set to "//trim(cvalue)
       end if
    end if

    ! Obtain dbug_flag setting if present; otherwise use default value in med_constants
    call NUOPC_CompAttributeGet(gcomp, name='dbug_flag', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
     read(cvalue,*) dbug_flag
    end if
    write(msgString,'(A,i6)') trim(subname)//': Mediator dbug_flag is ',dbug_flag
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine InitializeP0

  !-----------------------------------------------------------------------

  subroutine AdvertiseFields(gcomp, importState, exportState, clock, rc)

    ! Mediator advertises its import and export Fields and sets the
    ! TransferOfferGeomObject Attribute.

    use ESMF  , only : ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_SUCCESS, ESMF_LogFoundAllocError
    use ESMF  , only : ESMF_StateIsCreated
    use ESMF  , only : ESMF_LogMsg_Info, ESMF_LogWrite
    use ESMF  , only : ESMF_END_ABORT, ESMF_Finalize, ESMF_MAXSTR
    use NUOPC , only : NUOPC_AddNamespace, NUOPC_Advertise, NUOPC_AddNestedState
    use NUOPC , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use esmFlds, only : med_fldlist_init1, med_fld_GetFldInfo, med_fldList_entry_type
    use med_phases_history_mod, only : med_phases_history_init
    use med_methods_mod       , only : mediator_checkfornans

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=CS)   :: stdname, shortname
    integer             :: ncomp, ns
    logical             :: isPresent, isSet
    character(len=CS)   :: transferOffer
    character(len=CS)   :: cvalue
    character(len=8)    :: cnum
    type(InternalState) :: is_local
    type(med_fldlist_type), pointer :: fldListFr, fldListTo
    type(med_fldList_entry_type), pointer :: fld
    integer             :: stat
    character(len=*), parameter :: subname = '('//__FILE__//':AdvertiseFields)'
    !-----------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    !------------------
    ! Allocate memory for the internal state
    !------------------

    allocate(is_local%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Allocation of the internal state memory failed.", line=__LINE__, file=u_FILE_u)) then
       return  ! bail out
    end if

    call ESMF_GridCompSetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call med_internalstate_init(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! Allocate memory for history module variables
    !------------------
    call med_phases_history_init()

    !------------------
    ! add a namespace (i.e. nested state)  for each import and export component state in the mediator's InternalState
    !------------------

    ! Namespaces are implemented via nested states. This creates a nested state inside of
    ! state. The nested state is returned as nestedState. nestedStateName will be used to name the
    ! newly created nested state.

    call NUOPC_AddNamespace(importState, namespace="ATM", nestedStateName="AtmImp", &
         nestedState=is_local%wrap%NStateImp(compatm), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="ATM", nestedStateName="AtmExp", &
         nestedState=is_local%wrap%NStateExp(compatm), rc=rc)

    call NUOPC_AddNamespace(importState, namespace="OCN", nestedStateName="OcnImp", &
         nestedState=is_local%wrap%NStateImp(compocn), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="OCN", nestedStateName="OcnExp", &
         nestedState=is_local%wrap%NStateExp(compocn), rc=rc)

    call NUOPC_AddNamespace(importState, namespace="ICE", nestedStateName="IceImp", &
         nestedState=is_local%wrap%NStateImp(compice), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="ICE", nestedStateName="IceExp", &
         nestedState=is_local%wrap%NStateExp(compice), rc=rc)

    call NUOPC_AddNamespace(importState, namespace="LND", nestedStateName="LndImp", &
         nestedState=is_local%wrap%NStateImp(complnd), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="LND", nestedStateName="LndExp", &
         nestedState=is_local%wrap%NStateExp(complnd), rc=rc)

    call NUOPC_AddNamespace(importState, namespace="ROF", nestedStateName="RofImp", &
         nestedState=is_local%wrap%NStateImp(comprof), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="ROF", nestedStateName="RofExp", &
         nestedState=is_local%wrap%NStateExp(comprof), rc=rc)

    call NUOPC_AddNamespace(importState, namespace="WAV", nestedStateName="WavImp", &
         nestedState=is_local%wrap%NStateImp(compwav), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="WAV", nestedStateName="WavExp", &
         nestedState=is_local%wrap%NStateExp(compwav), rc=rc)

    ! Only create nested states for active land-ice sheets
    do ns = 1,is_local%wrap%num_icesheets
       write(cnum,'(i0)') ns
       call NUOPC_AddNestedState(importState, CplSet="GLC"//trim(cnum), &
            nestedState=is_local%wrap%NStateImp(compglc(ns)), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_AddNestedState(exportState, CplSet="GLC"//trim(cnum), &
            nestedState=is_local%wrap%NStateExp(compglc(ns)), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! Determine aoflux grid
    call NUOPC_CompAttributeGet(gcomp, name='aoflux_grid', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (.not. isPresent .and. .not. isSet) then
       cvalue = 'ogrid'
    end if
    is_local%wrap%aoflux_grid = trim(cvalue)

    ! Determine aoflux scheme that will be used to compute atmosphere-ocean fluxes [cesm|ccpp]
    ! TODO: If ccpp is not available it will be always run in cesm mode independent from aoflux_code option
    call NUOPC_CompAttributeGet(gcomp, name='aoflux_code', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (.not. isPresent .and. .not. isSet) then
       cvalue = 'cesm'
    end if
    aoflux_code = trim(cvalue)
    if (maintask) then
       write(logunit,*) '========================================================'
       write(logunit,'(a)')trim(subname)//' Mediator aoflux scheme is '//trim(aoflux_code)
       write(logunit,*) '========================================================'
    end if

    ! Determine CCPP suite if aoflux scheme set to 'ccpp'
    if (trim(aoflux_code) == 'ccpp') then
       call NUOPC_CompAttributeGet(gcomp, name='aoflux_ccpp_suite', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (.not. isPresent .and. .not. isSet) then
          call ESMF_LogWrite("aoflux_ccpp_suite need to be provided when aoflux_code is set to 'ccpp'", ESMF_LOGMSG_INFO)
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       aoflux_ccpp_suite = trim(cvalue)
       if (maintask) then
          write(logunit,*) '========================================================'
          write(logunit,'(a)')trim(subname)//' Mediator aoflux CCPP suite is '//trim(aoflux_ccpp_suite)
          write(logunit,*) '========================================================'
       end if
    end if

    !------------------
    ! Initialize mediator flds
    !------------------

    call NUOPC_CompAttributeGet(gcomp, name='coupling_mode', value=coupling_mode, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('coupling_mode = '// trim(coupling_mode), ESMF_LOGMSG_INFO)
    if (maintask) then
       write(logunit,*) '========================================================'
       write(logunit,'(a)')trim(subname)//' Mediator Coupling Mode is '//trim(coupling_mode)
       write(logunit,*) '========================================================'
       write(logunit,*)
    end if

    ! Initialize memory for fldlistTo and fldlistFr - this is need for the calls below for the
    ! advertise phase
    call med_fldlist_init1(ncomps)

    if (trim(coupling_mode) == 'cesm') then
       call esmFldsExchange_cesm(gcomp, phase='advertise', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (coupling_mode(1:3) == 'ufs') then
       call esmFldsExchange_ufs(gcomp, phase='advertise', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (coupling_mode(1:4) == 'hafs') then
       call esmFldsExchange_hafs(gcomp, phase='advertise', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
        call ESMF_LogWrite(trim(coupling_mode)//' is not a valid coupling_mode', ESMF_LOGMSG_INFO)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Set default masking for mapping
    call med_internalstate_defaultmasks(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! Determine component present indices
    !------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    is_local%wrap%flds_scalar_name = trim(cvalue)

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue, *) is_local%wrap%flds_scalar_num

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) is_local%wrap%flds_scalar_index_nx

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) is_local%wrap%flds_scalar_index_ny

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNTile", value=cvalue, &
         isPresent=isPresent, isSet=isSet,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) is_local%wrap%flds_scalar_index_ntile
    else
       is_local%wrap%flds_scalar_index_ntile = 0
    end if

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) is_local%wrap%flds_scalar_index_nextsw_cday
    else
       is_local%wrap%flds_scalar_index_nextsw_cday = 0
    end if

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxPrecipFactor", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) is_local%wrap%flds_scalar_index_precip_factor
    else
       is_local%wrap%flds_scalar_index_precip_factor = 0
    end if

    !------------------
    ! Advertise import/export mediator field names
    !------------------

    do ncomp = 1,ncomps
       if (ncomp /= compmed) then
          if (maintask) write(logunit,*)
          fldListFr => med_fldList_GetFldListFr(ncomp)
          fld => fldListFr%fields
          do while(associated(fld))
             call med_fld_GetFldInfo(fld, stdname=stdname, shortname=shortname)
             if (maintask) then
                write(logunit,'(a)') trim(subname)//':Fr_'//trim(compname(ncomp))//': '//trim(shortname)
             end if
             if (trim(shortname) == is_local%wrap%flds_scalar_name) then
                transferOffer = 'will provide'
             else
                transferOffer = 'cannot provide'
             end if
             call NUOPC_Advertise(is_local%wrap%NStateImp(ncomp), &
                  standardName=stdname, shortname=shortname, name=shortname, &
                  TransferOfferGeomObject=transferOffer, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_LogWrite(subname//':Fr_'//trim(compname(ncomp))//': '//trim(shortname), ESMF_LOGMSG_INFO)
             fld => fld%next
          end do

          fldListTo => med_fldList_GetFldListTo(ncomp)
          fld => fldListTo%fields
          do while(associated(fld))
             call med_fld_GetFldInfo(fld, stdname=stdname, shortname=shortname, rc=rc)
             if (maintask) then
                write(logunit,'(a)') trim(subname)//':To_'//trim(compname(ncomp))//': '//trim(shortname)
             end if
             if (trim(shortname) == is_local%wrap%flds_scalar_name) then
                transferOffer = 'will provide'
             else
                transferOffer = 'cannot provide'
             end if
             call NUOPC_Advertise(is_local%wrap%NStateExp(ncomp), &
                  standardName=stdname, shortname=shortname, name=shortname, &
                  TransferOfferGeomObject=transferOffer, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_LogWrite(subname//':To_'//trim(compname(ncomp))//': '//trim(shortname), ESMF_LOGMSG_INFO)
             fld => fld%next
          end do
       end if
    end do ! end of ncomps loop

    ! Should mediator check for NaNs?
    call NUOPC_CompAttributeGet(gcomp, name="check_for_nans", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(isPresent .and. isSet) then
       read(cvalue, *) mediator_checkfornans
    else
       mediator_checkfornans = .false.
    endif
    if(maintask) then
       write(logunit,*) ' check_for_nans is ',mediator_checkfornans
       if(mediator_checkfornans) then
          write(logunit,*) ' Fields will be checked for NaN values when passed from mediator to component'
       else
          write(logunit,*) ' Fields will NOT be checked for NaN values when passed from mediator to component'
       endif
    endif

    ! Should target component use all data for first time step?
    do ncomp = 1,ncomps
       if (ncomp /= compmed) then
          call NUOPC_CompAttributeGet(gcomp, name=trim(compname(ncomp))//"_use_data_first_import", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (isPresent .and. isSet) then
             read(cvalue, *) is_local%wrap%med_data_force_first(ncomp)
          else
             is_local%wrap%med_data_force_first(ncomp) = .false.
          endif
          if (maintask) then
             write(logunit,*) trim(compname(ncomp))//'_use_data_first_import is ', is_local%wrap%med_data_force_first(ncomp)
          endif
       end if
    end do

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine AdvertiseFields

  !-----------------------------------------------------------------------------

  subroutine RealizeFieldsWithTransferProvided(gcomp, importState, exportState, clock, rc)

    ! Realize connected Fields with transfer action "provide"

    use ESMF , only : ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_VM, ESMF_SUCCESS
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_TimeInterval
    use ESMF , only : ESMF_VMGet, ESMF_StateIsCreated, ESMF_GridCompGet
    use ESMF , only : ESMF_StateSet, ESMF_StateIntent_Import, ESMF_StateIntent_Export
    use ESMF , only : ESMF_StateIntent_Flag

    ! Input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    type(ESMF_VM)              :: vm
    integer                    :: n
    character(len=*), parameter :: subname = '('//__FILE__//':RealizeFieldsWithTransferProvided)'
    !-----------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

     ! Initialize the internal state mediator vm
     call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     is_local%wrap%vm = vm

    ! Realize States
    do n = 1,ncomps
      if (ESMF_StateIsCreated(is_local%wrap%NStateImp(n), rc=rc)) then
         call ESMF_StateSet(is_local%wrap%NStateImp(n), stateIntent=ESMF_StateIntent_Import, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         call med_fldList_Realize(is_local%wrap%NStateImp(n), med_fldList_GetfldListFr(n), &
              is_local%wrap%flds_scalar_name, is_local%wrap%flds_scalar_num, &
              tag=subname//':Fr_'//trim(compname(n)), rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      if (ESMF_StateIsCreated(is_local%wrap%NStateExp(n), rc=rc)) then
          call ESMF_StateSet(is_local%wrap%NStateExp(n), stateIntent=ESMF_StateIntent_Export, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_fldList_Realize(is_local%wrap%NStateExp(n), med_fldList_getfldListTo(n), &
              is_local%wrap%flds_scalar_name, is_local%wrap%flds_scalar_num, &
              tag=subname//':To_'//trim(compname(n)), rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    enddo

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine RealizeFieldsWithTransferProvided

  !-----------------------------------------------------------------------------

  subroutine ModifyDecompofMesh(gcomp, importState, exportState, clock, rc)

    ! Optionally modify the decomp/distr of transferred Grid/Mesh

    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_GRIDCOMP, ESMF_CLOCK, ESMF_STATE
    use ESMF , only : ESMF_StateIsCreated

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer :: n1
    character(len=*), parameter :: subname = '('//__FILE__//':ModifyDecompofMesh)'
    !-----------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    ! Get the internal state from the mediator gridded component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! Receive Grids
    !------------------

    do n1 = 1,ncomps
       call ESMF_LogWrite(trim(subname)//": calling for component "//trim(compname(n1)), ESMF_LOGMSG_INFO)
       if (ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
          call realizeConnectedGrid(is_local%wrap%NStateImp(n1), trim(compname(n1))//'Imp', rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
       if (ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc)) then
          call realizeConnectedGrid(is_local%wrap%NStateExp(n1), trim(compname(n1))//'Exp', rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
       call ESMF_LogWrite(trim(subname)//": finished for component "//trim(compname(n1)), ESMF_LOGMSG_INFO)
    enddo
    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine realizeConnectedGrid(State,string,rc)

      use ESMF , only : operator(==)
      use ESMF , only : ESMF_STATE, ESMF_Field, ESMF_Grid, ESMF_DistGrid, ESMF_DistGridConnection
      use ESMF , only : ESMF_MAXSTR, ESMF_FieldStatus_Flag, ESMF_GeomType_Flag, ESMF_StateGet
      use ESMF , only : ESMF_FieldGet, ESMF_DistGridGet, ESMF_GridCompGet
      use ESMF , only : ESMF_GeomType_Grid, ESMF_AttributeGet, ESMF_DistGridCreate, ESMF_FieldEmptySet
      use ESMF , only : ESMF_GridCreate, ESMF_LogWrite, ESMF_LogMsg_Info, ESMF_GridGet, ESMF_Failure
      use ESMF , only : ESMF_LogMsg_Warning
      use ESMF , only : ESMF_FieldStatus_Empty, ESMF_FieldStatus_Complete, ESMF_FieldStatus_GridSet
      use ESMF , only : ESMF_GeomType_Mesh, ESMF_MeshGet, ESMF_Mesh, ESMF_MeshEmptyCreate

      ! input/output variables
      type(ESMF_State)   , intent(inout) :: State
      character(len=*)   , intent(in)    :: string
      integer            , intent(out)   :: rc

      ! local variables
      type(ESMF_Field)              :: field
      type(ESMF_Grid)               :: grid, newgrid
      type(ESMF_Mesh)               :: mesh, newmesh
      type(ESMF_DistGrid)           :: distgrid
      type(ESMF_DistGrid)           :: elemdistgrid, newelemdistgrid
      integer                       :: arbDimCount
      integer                       :: dimCount, tileCount
      integer                       :: connectionCount
      integer                       :: fieldCount
      integer                       :: n, n1, i1, i2
      type(ESMF_GeomType_Flag)      :: geomtype
      type(ESMF_FieldStatus_Flag)   :: fieldStatus
      character(len=CX)             :: msgString
      integer                       , allocatable :: minIndexPTile(:,:), maxIndexPTile(:,:)
      character(ESMF_MAXSTR)        , allocatable :: fieldNameList(:)
      type(ESMF_DistGridConnection) , allocatable :: connectionList(:)
      character(len=*),parameter :: subname=' (realizeConnectedGrid) '
      !-----------------------------------------------------------

      ! All of the Fields that set their TransferOfferGeomObject Attribute
      ! to "cannot provide" should now have the accepted Grid available.
      ! Go and pull out this Grid for one of a representative Field and
      ! modify the decomposition and distribution of the Grid to match the Mediator PETs.
      ! On exit from this phase, the connector will transfer the full Grid/Mesh/LocStream
      ! objects (with coordinates) for Field pairs that have a provider and an acceptor side.

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_Success
      if (profile_memory) then
         call ESMF_VMLogMemInfo("Entering "//trim(subname))
      end if
      call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(fieldNameList(fieldCount))
      call ESMF_StateGet(State, itemNameList=fieldNameList, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! do not loop here, assuming that all fields share the
      ! same grid/mesh and because it is more efficient - if
      ! a component has fields on multiple grids/meshes, this
      ! would need to be revisited
      do n=1, min(fieldCount, 1)

         call ESMF_StateGet(State, field=field, itemName=fieldNameList(n), rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         if (fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then

            ! The Mediator is accepting a Grid/Mesh passed to it
            ! through the Connector
            ! While this is still an empty field, it does now hold a Grid/Mesh with DistGrid
            call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            if (geomtype == ESMF_GEOMTYPE_GRID) then

               call ESMF_AttributeGet(field, name="ArbDimCount", value=arbDimCount, &
                    convention="NUOPC", purpose="Instance", rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
               call ESMF_LogWrite(trim(subname)//": geomtype is ESMF_GEOMTYPE_GRID for "//trim(fieldnameList(n)), &
                    ESMF_LOGMSG_INFO)
               write(msgString,'(A,i8)') trim(subname)//':arbdimcount =',arbdimcount
               call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

               ! make decision on whether the incoming Grid is arbDistr or not
               if (arbDimCount>0) then

                  ! The provider defined an arbDistr grid
                  ! - use a regDecomp representation for the grid
                  ! - first get tile min/max, only single tile supported for arbDistr Grid
                  ! - create default regDecomp DistGrid
                  ! - create default regDecomp Grid with just a distgrid
                  call ESMF_LogWrite(trim(subname)//trim(string)//": accept arb2reg grid for "//trim(fieldNameList(n)), &
                       ESMF_LOGMSG_INFO)
                  allocate(minIndexPTile(arbDimCount,1),maxIndexPTile(arbDimCount,1))
                  call ESMF_AttributeGet(field, name="MinIndex", &
                       valueList=minIndexPTile(:,1), &
                       convention="NUOPC", purpose="Instance", rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  call ESMF_AttributeGet(field, name="MaxIndex", &
                       valueList=maxIndexPTile(:,1), convention="NUOPC", purpose="Instance", rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  newgrid = ESMF_GridCreate(distgrid, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  do i1 = 1,arbDimCount
                     write(msgString,'(A,3i8)') trim(subname)//':PTile =',i1,minIndexPTile(i1,1),maxIndexPTile(i1,1)
                     call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
                  enddo

               else   ! arbdimcount <= 0

                  ! The provider sends a non arb grid
                  ! Create a custom DistGrid, based on the minIndex, maxIndex of the accepted DistGrid,
                  ! but with a default regDecomp for the current VM that leads to 1DE/PET.
                  ! - get dimCount and tileCount
                  ! - allocate minIndexPTile and maxIndexPTile according to dimCount and tileCount
                  ! - get minIndex and maxIndex arrays and connectionList
                  ! - create the new DistGrid with the same minIndexPTile and maxIndexPTile
                  ! - create a new Grid on the new DistGrid

                  call ESMF_LogWrite(trim(subname)//trim(string)//": accept reg2reg grid for "//&
                       trim(fieldNameList(n)), ESMF_LOGMSG_INFO)
                  call ESMF_FieldGet(field, grid=grid, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  call ESMF_GridGet(grid, distgrid=distgrid, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  allocate(minIndexPTile(dimCount, tileCount))
                  allocate(maxIndexPTile(dimCount, tileCount))
                  call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
                       maxIndexPTile=maxIndexPTile, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  do i2 = 1,tileCount
                     do i1 = 1,dimCount
                        write(msgString,'(A,4i8)') trim(subname)//':PTile =',i2,i1,minIndexPTile(i1,i2),&
                             maxIndexPTile(i1,i2)
                        call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
                     enddo
                  enddo
                  if (dimcount == 2) then
                     call ESMF_DistGridGet(distgrid, connectionCount=connectionCount, rc=rc)
                     allocate(connectionList(connectionCount))
                     call ESMF_DistGridGet(distgrid, connectionList=connectionList, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
                          maxIndexPTile=maxIndexPTile, connectionList=connectionList, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     call ESMF_LogWrite(trim(subname)//trim(string)//': distgrid with dimcount=2', ESMF_LOGMSG_INFO)
                     newgrid = ESMF_GridCreate(distgrid, gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     deallocate(connectionList)
                  else
                     distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
                          maxIndexPTile=maxIndexPTile, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     call ESMF_LogWrite(trim(subname)//trim(string)//': distgrid with dimcount=1', ESMF_LOGMSG_INFO)
                     newgrid = ESMF_GridCreate(distgrid, gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/), rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  endif

                  ! local clean-up
                  deallocate(minIndexPTile, maxIndexPTile)

               endif  ! arbdimCount

               ! Swap all the Grids in the State
               do n1=1, fieldCount
                  ! access a field in the State and set the Grid
                  call ESMF_StateGet(State, field=field, itemName=fieldNameList(n1), rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  if (fieldStatus==ESMF_FIELDSTATUS_EMPTY .or. fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then
                     call ESMF_FieldEmptySet(field, grid=newgrid, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     call ESMF_LogWrite(trim(subname)//trim(string)//": attach grid for "//trim(fieldNameList(n1)), &
                          ESMF_LOGMSG_INFO)
                     if (dbug_flag > 1) then
                        call Field_GeomPrint(field,trim(fieldNameList(n1))//'_new',rc)
                        if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     end if
                  else
                     call ESMF_LogWrite(trim(subname)//trim(string)//": NOT replacing grid for field: "//&
                          trim(fieldNameList(n1)), ESMF_LOGMSG_WARNING)
                  endif
               enddo

            elseif (geomtype == ESMF_GEOMTYPE_MESH) then

               call ESMF_LogWrite(trim(subname)//": geomtype is ESMF_GEOMTYPE_MESH for "//trim(fieldnameList(n)), &
                    ESMF_LOGMSG_INFO)

               if (dbug_flag > 1) then
                  call Field_GeomPrint(field,trim(fieldNameList(n))//'_orig',rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
               end if
               call ESMF_FieldGet(field, mesh=mesh, rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
               call ESMF_MeshGet(mesh, elementDistGrid=elemDistGrid, rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
               newelemDistGrid = ESMF_DistGridCreate(elemDistGrid, balanceflag=.true., rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
               newmesh = ESMF_MeshEmptyCreate(elementDistGrid=newelemDistGrid, rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return

               ! Swap all the Meshes in the State
               do n1=1, fieldCount
                  ! access a field in the State and set the Mesh
                  call ESMF_StateGet(State, field=field, itemName=fieldNameList(n1), rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  if (fieldStatus==ESMF_FIELDSTATUS_EMPTY .or. fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then
                     call ESMF_FieldEmptySet(field, mesh=newmesh, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     call ESMF_LogWrite(trim(subname)//trim(string)//": attach mesh for "//&
                          trim(fieldNameList(n1)), ESMF_LOGMSG_INFO)
                     if (dbug_flag > 1) then
                        call Field_GeomPrint(field,trim(fieldNameList(n1))//'_new',rc)
                        if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     end if
                  else
                     call ESMF_LogWrite(trim(subname)//trim(string)//": NOT replacing mesh for field: "//&
                          trim(fieldNameList(n1)), ESMF_LOGMSG_WARNING)
                  endif
               enddo

            else  ! geomtype

               call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", ESMF_LOGMSG_INFO)
               rc=ESMF_FAILURE
               return

            endif ! geomtype

         elseif (fieldStatus==ESMF_FIELDSTATUS_EMPTY) then

            call ESMF_LogWrite(trim(subname)//trim(string)//": provide grid for "//trim(fieldNameList(n)), &
                 ESMF_LOGMSG_INFO)

         elseif (fieldStatus==ESMF_FIELDSTATUS_COMPLETE) then

            call ESMF_LogWrite(trim(subname)//trim(string)//": no grid provided for "//trim(fieldNameList(n)), &
                 ESMF_LOGMSG_INFO)

         else

            call ESMF_LogWrite(trim(subname)//": ERROR fieldStatus not supported ", ESMF_LOGMSG_INFO)
            rc=ESMF_FAILURE
            return

         endif   ! fieldStatus

      enddo   ! nflds

      deallocate(fieldNameList)

      if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine realizeConnectedGrid

  end subroutine ModifyDecompofMesh

  !-----------------------------------------------------------------------------

  subroutine RealizeFieldsWithTransferAccept(gcomp, importState, exportState, clock, rc)

    use ESMF , only : ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_LogWrite
    use ESMF , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_StateIsCreated

    !----------------------------------------------------------
    ! realize all Fields with transfer action "accept"
    ! Finish initializing the State Fields
    ! - Fields are partially created when this routine is called.
    ! - Fields contain a geombase object internally created and the geombase object
    !   associates with either a ESMF_Grid, or a ESMF_Mesh, or an or an ESMF_XGrid,
    !   or a ESMF_LocStream.
    ! - Fields containing grids will be transferred! to a Mesh and Realized;
    ! - Fields containg meshes are completed with space allocated internally
    !   for an ESMF_Array based on arrayspec
    !----------------------------------------------------------

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n1
    character(len=*), parameter :: subname = '('//__FILE__//':RealizeFieldsWithTransferAccept)'
    !-----------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do n1 = 1,ncomps
      ! Finish initializing import states and reset state data to spval_init
      if (ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
         call ESMF_LogWrite(trim(subname)//": calling completeFieldInitialize import states from "//trim(compname(n1)), &
              ESMF_LOGMSG_INFO)
        call completeFieldInitialization(is_local%wrap%NStateImp(n1), rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call State_reset(is_local%wrap%NStateImp(n1), value=spval_init, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      ! Finish initializing mediator export states and reset state data to spval_init
      if (ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc)) then
         call ESMF_LogWrite(trim(subname)//": calling completeFieldInitialize export states to "//trim(compname(n1)), &
              ESMF_LOGMSG_INFO)
        call completeFieldInitialization(is_local%wrap%NStateExp(n1), rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call State_reset(is_local%wrap%NStateExp(n1), value=spval_init, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
           call State_GeomPrint(is_local%wrap%NStateExp(n1),'gridExp'//trim(compname(n1)),rc=rc)
           if (ChkErr(rc,__LINE__,u_FILE_u)) return
        end if
      endif
    enddo

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine completeFieldInitialization(State,rc)

      use ESMF  , only : operator(==)
      use ESMF  , only : ESMF_State, ESMF_MAXSTR, ESMF_Grid, ESMF_Mesh, ESMF_Field, ESMF_FieldStatus_Flag
      use ESMF  , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FieldGet
      use ESMF  , only : ESMF_GeomType_Flag, ESMF_FieldCreate, ESMF_MeshCreate, ESMF_GEOMTYPE_GRID
      use ESMF  , only : ESMF_MeshLoc_Element, ESMF_TYPEKIND_R8, ESMF_FIELDSTATUS_GRIDSET
      use ESMF  , only : ESMF_AttributeGet, ESMF_MeshWrite, ESMF_FAILURE
      use NUOPC , only : NUOPC_getStateMemberLists, NUOPC_Realize

      ! input/output variables
      type(ESMF_State)   , intent(inout) :: State
      integer            , intent(out)   :: rc

      ! local varaibles
      integer                     :: n, fieldCount
      character(ESMF_MAXSTR)      :: fieldName
      type(ESMF_Grid)             :: grid
      type(ESMF_Mesh)             :: mesh
      type(ESMF_Field)            :: meshField
      type(ESMF_Field),pointer    :: fieldList(:)
      type(ESMF_FieldStatus_Flag) :: fieldStatus
      type(ESMF_GeomType_Flag)    :: geomtype
      integer                     :: gridToFieldMapCount, ungriddedCount
      integer, allocatable        :: gridToFieldMap(:)
      integer, allocatable        :: ungriddedLBound(:), ungriddedUBound(:)
      logical                     :: isPresent
      logical                     :: meshcreated
      character(len=*), parameter :: subname = '('//__FILE__//':completeFieldInitialization)'
      !-----------------------------------------------------------

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_Success
      if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

      call State_GetNumFields(State, fieldCount, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (fieldCount > 0) then
         nullify(fieldList)
         call NUOPC_getStateMemberLists(State, fieldList=fieldList, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         meshcreated = .false.
         do n=1, fieldCount

            call ESMF_FieldGet(fieldList(n), status=fieldStatus, name=fieldName, geomtype=geomtype, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Input fields contains grid - need to convert to mesh
            if (geomtype == ESMF_GEOMTYPE_GRID .and. fieldName /= is_local%wrap%flds_scalar_name) then

               ! Grab grid
               if (dbug_flag > 1) then
                  call Field_GeomPrint(fieldList(n),trim(fieldName)//'_premesh',rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
               end if
               call ESMF_FieldGet(fieldList(n), grid=grid, rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return

               ! Convert grid to mesh
               if (.not. meshcreated) then
                  if (dbug_flag > 20) then
                     call med_grid_write(grid, trim(fieldName)//'_premesh.nc', rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  end if
                  mesh = ESMF_MeshCreate(grid, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  meshcreated = .true.
                  if (dbug_flag > 20) then
                     call ESMF_MeshWrite(mesh, filename=trim(fieldName)//'_postmesh', rc=rc)
                     if (chkerr(rc,__LINE__,u_FILE_u)) return
                  end if
               end if

               ! Create field on mesh
               meshField = ESMF_FieldCreate(mesh, typekind=ESMF_TYPEKIND_R8, &
                    meshloc=ESMF_MESHLOC_ELEMENT, name=fieldName, rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return

               ! Swap grid for mesh, at this point, only connected fields are in the state
               call NUOPC_Realize(State, field=meshField, rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return

               call ESMF_FieldGet(meshField, status=fieldStatus, rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
               if (fieldStatus == ESMF_FIELDSTATUS_GRIDSET ) then
                 call ESMF_LogWrite(trim(subname)//": ERROR fieldStatus not complete ", ESMF_LOGMSG_INFO)
                 rc = ESMF_FAILURE
                 return
               end if
               call Field_GeomPrint(meshField, trim(subname)//':'//trim(fieldName), rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return

            else

               ! Input fields contain mesh
               if (fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then
                  call ESMF_AttributeGet(fieldList(n), name="GridToFieldMap", convention="NUOPC", &
                       purpose="Instance", itemCount=gridToFieldMapCount, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  allocate(gridToFieldMap(gridToFieldMapCount))
                  call ESMF_AttributeGet(fieldList(n), name="GridToFieldMap", convention="NUOPC", &
                       purpose="Instance", valueList=gridToFieldMap, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  ungriddedCount=0  ! initialize in case it was not set
                  call ESMF_AttributeGet(fieldList(n), name="UngriddedLBound", convention="NUOPC", &
                       purpose="Instance", itemCount=ungriddedCount,  isPresent=isPresent, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  allocate(ungriddedLBound(ungriddedCount), ungriddedUBound(ungriddedCount))
                  if (ungriddedCount > 0) then
                     call ESMF_AttributeGet(fieldList(n), name="UngriddedLBound", convention="NUOPC", &
                          purpose="Instance", valueList=ungriddedLBound, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     call ESMF_AttributeGet(fieldList(n), name="UngriddedUBound", convention="NUOPC", &
                          purpose="Instance", valueList=ungriddedUBound, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  endif
                  call NUOPC_Realize(State, fieldName, typekind=ESMF_TYPEKIND_R8, gridToFieldMap=gridToFieldMap, &
                       ungriddedLbound=ungriddedLbound, ungriddedUbound=ungriddedUbound, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  deallocate(gridToFieldMap, ungriddedLbound, ungriddedUbound)
               end if ! fieldStatus
               call Field_GeomPrint(fieldlist(n), trim(subname)//':'//trim(fieldName), rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return

            end if

         enddo ! end of loop over fields
         deallocate(fieldList)

      endif ! end of fieldcount< 0

      if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine completeFieldInitialization

  end subroutine RealizeFieldsWithTransferAccept

  !-----------------------------------------------------------------------------

  subroutine DataInitialize(gcomp, rc)

    !----------------------------------------------------------
    ! Finish initialization and resolve data dependencies
    ! There will be multiple passes
    ! For first time through:
    !   Do not assume any import fields are connected, just allocate space and such
    !   -- Check present flags
    !   -- Check for active coupling interactions
    !   -- Create FBs: FBImp, FBExp
    !   -- Create mediator specific field bundles (not part of import/export states)
    !   -- Read mediator restarts
    !   -- Initialize route handles field bundles for normalization
    !   -- return!
    ! For second loop:
    !   -- Copy import fields to local FBs
    !   -- Create FBfrac and initialize fractions
    ! Once the ocean is ready:
    !   -- Copy import fields to local FBs
    !   -- Re-initialize fractions
    !   -- Carry out ocnalb_init
    !   -- Carry out aoffluxes_init
    ! Once the atm is ready:
    !   -- Copy import fields to local FBs
    !----------------------------------------------------------

    use ESMF                    , only : ESMF_GridComp, ESMF_Clock, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF                    , only : ESMF_State, ESMF_Time, ESMF_Field, ESMF_StateItem_Flag, ESMF_MAXSTR
    use ESMF                    , only : ESMF_GridCompGet, ESMF_AttributeGet, ESMF_ClockGet, ESMF_Success
    use ESMF                    , only : ESMF_StateIsCreated, ESMF_StateGet, ESMF_FieldBundleIsCreated, ESMF_LogFlush
    use ESMF                    , only : ESMF_FieldBundleGet, ESMF_VM
    use NUOPC                   , only : NUOPC_CompAttributeSet, NUOPC_IsAtTime, NUOPC_SetAttribute
    use NUOPC                   , only : NUOPC_CompAttributeGet
    use med_fraction_mod        , only : med_fraction_init, med_fraction_set
    use med_phases_restart_mod  , only : med_phases_restart_read
    use med_phases_prep_ocn_mod , only : med_phases_prep_ocn_init
    use med_phases_prep_wav_mod , only : med_phases_prep_wav_init
    use med_phases_prep_rof_mod , only : med_phases_prep_rof_init
    use med_phases_prep_glc_mod , only : med_phases_prep_glc_init
    use med_phases_prep_atm_mod , only : med_phases_prep_atm
    use med_phases_post_atm_mod , only : med_phases_post_atm
    use med_phases_post_ice_mod , only : med_phases_post_ice
    use med_phases_post_lnd_mod , only : med_phases_post_lnd
    use med_phases_post_glc_mod , only : med_phases_post_glc
    use med_phases_post_ocn_mod , only : med_phases_post_ocn
    use med_phases_post_rof_mod , only : med_phases_post_rof_init, med_phases_post_rof
    use med_phases_post_wav_mod , only : med_phases_post_wav
    use med_phases_ocnalb_mod   , only : med_phases_ocnalb_run
    use med_phases_aofluxes_mod , only : med_phases_aofluxes_init_fldbuns
    use med_phases_profile_mod  , only : med_phases_profile
    use med_diag_mod            , only : med_diag_zero, med_diag_init
    use med_map_mod             , only : med_map_routehandles_init, med_map_packed_field_create
    use med_io_mod              , only : med_io_init
    use esmFlds                 , only : med_fldList_GetaofluxfldList

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)                :: is_local
    type(ESMF_Clock)                   :: clock
    type(ESMF_State)                   :: importState, exportState
    type(ESMF_Time)                    :: time
    type(ESMF_Field)                   :: field
    type(med_fldList_type), pointer    :: fldListMed_ocnalb
    logical                            :: atCorrectTime
    integer                            :: n1,n2,n
    integer                            :: nsrc,ndst
    integer                            :: fieldCount
    character(ESMF_MAXSTR),allocatable :: fieldNameList(:)
    character(CL), pointer             :: fldnames(:)
    character(CL)                      :: cvalue
    logical                            :: read_restart
    logical                            :: allDone = .false.
    logical,save                       :: first_call = .true.
    real(r8)                           :: real_nx, real_ny, real_ntile
    character(len=CX)                  :: msgString
    character(len=*), parameter :: subname = '('//__FILE__//':DataInitialize)'
    !-----------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="false", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=time, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Beginning  of first_call block
    !---------------------------------------

    if (first_call) then

       ! Allocate module variable
       allocate(compDone(ncomps))

       ! Determine active coupling logical flags
       call med_internalstate_coupling(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

      !----------------------------------------------------------
      ! Create field bundles FBImp, FBExp
      !----------------------------------------------------------

      if (maintask) then
         write(logunit,'(a)') 'Creating mediator field bundles '
      end if

      do n1 = 1,ncomps
         if (is_local%wrap%comp_present(n1) .and. &
              ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc) .and. &
              ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc)) then

            if (maintask) then
               write(logunit,'(a)') trim(subname)//' initializing FBs for '//trim(compname(n1))
            end if

            ! Create FBImp(:) with pointers directly into NStateImp(:)
            call FB_init_pointer(is_local%wrap%NStateImp(n1), is_local%wrap%FBImp(n1,n1), &
                 is_local%wrap%flds_scalar_name, name='FBImp'//trim(compname(n1)), rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Create FBExp(:) with pointers directly into NStateExp(:)
            call FB_init_pointer(is_local%wrap%NStateExp(n1), is_local%wrap%FBExp(n1), &
                 is_local%wrap%flds_scalar_name, name='FBExp'//trim(compname(n1)), rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Create mesh info data
            call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n1), fieldCount=fieldCount, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            if (fieldCount == 0) then
              if (maintask) then
                write(logunit,*) trim(subname)//' '//trim(compname(n1))//' import FB field count is = ', fieldCount
                write(logunit,*) trim(subname)//' '//trim(compname(n1))//' trying to use export FB'
                call ESMF_FieldBundleGet(is_local%wrap%FBExp(n1), fieldCount=fieldCount, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                write(logunit,*) trim(subname)//' '//trim(compname(n1))//' export FB field count is = ', fieldCount
              end if
              call med_meshinfo_create(is_local%wrap%FBExp(n1), &
                   is_local%wrap%mesh_info(n1), is_local%wrap%FBArea(n1), rc=rc)
            else
              call med_meshinfo_create(is_local%wrap%FBImp(n1,n1), &
                   is_local%wrap%mesh_info(n1), is_local%wrap%FBArea(n1), rc=rc)
            end if
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
         end if

         ! The following is FBImp mapped to different grids. FBImp(n1,n1) is handled above
         do n2 = 1,ncomps
            if (n1 /= n2 .and. &
                 is_local%wrap%med_coupling_active(n1,n2) .and. &
                 ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc) .and. &
                 ESMF_StateIsCreated(is_local%wrap%NStateImp(n2),rc=rc)) then

               if (maintask) then
                  write(logunit,'(a)') trim(subname)//' initializing FBs for '//&
                       trim(compname(n1))//'_'//trim(compname(n2))
               end if

               ! Check import FB, if there is no field in it then use export FB
               ! to provide mesh information
               call State_GetNumFields(is_local%wrap%NStateImp(n2), fieldCount, rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
               if (fieldCount == 0) then
                 call FB_init(is_local%wrap%FBImp(n1,n2), is_local%wrap%flds_scalar_name, &
                      STgeom=is_local%wrap%NStateExp(n2), &
                      STflds=is_local%wrap%NStateImp(n1), &
                      name='FBImp'//trim(compname(n1))//'_'//trim(compname(n2)), rc=rc)
               else
                 call FB_init(is_local%wrap%FBImp(n1,n2), is_local%wrap%flds_scalar_name, &
                      STgeom=is_local%wrap%NStateImp(n2), &
                      STflds=is_local%wrap%NStateImp(n1), &
                      name='FBImp'//trim(compname(n1))//'_'//trim(compname(n2)), rc=rc)
               end if
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            endif
         enddo ! loop over n2
      enddo ! loop over n1

      !---------------------------------------
      ! Initialize field bundles needed for ocn albedo calculation
      !---------------------------------------

      ! NOTE: the NStateImp(compocn) or NStateImp(compatm) used below
      ! rather than NStateExp(n2), since the export state might only
      ! contain control data and no grid information if if the target
      ! component (n2) is not prognostic only receives control data back

      ! NOTE: this section must be done BEFORE the second call to esmFldsExchange
      ! Create field bundles for mediator ocean albedo computation

      if ( is_local%wrap%med_coupling_active(compocn,compatm) .or. is_local%wrap%med_coupling_active(compatm,compocn)) then
         ! Create field bundles for mediator ocean albedo computation
         fldListMed_ocnalb => med_fldlist_getocnalbFldList()
         fieldCount = med_fldList_GetNumFlds(fldListMed_ocnalb)
         if (fieldCount > 0) then
            allocate(fldnames(fieldCount))
            call med_fldList_getfldnames(fldListMed_ocnalb%fields, fldnames, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            call FB_init(is_local%wrap%FBMed_ocnalb_a, is_local%wrap%flds_scalar_name, &
                 STgeom=is_local%wrap%NStateImp(compatm), fieldnamelist=fldnames, name='FBMed_ocnalb_a', rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            if (maintask) then
               write(logunit,'(a)') trim(subname)//' initializing FB FBMed_ocnalb_a'
            end if
            call FB_init(is_local%wrap%FBMed_ocnalb_o, is_local%wrap%flds_scalar_name, &
                 STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames, name='FBMed_ocnalb_o', rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            if (maintask) then
               write(logunit,'(a)') trim(subname)//' initializing FB FBMed_ocnalb_o'
            end if
            deallocate(fldnames)
         end if
      end if

      !---------------------------------------
      ! Initialize field bundles needed for atm/ocn flux computation:
      ! is_local%wrap%FBMed_aoflux_a and is_local%wrap%FBMed_aoflux_o
      !---------------------------------------

      ! NOTE: this section must be done BEFORE the second call to esmFldsExchange
      ! Create field bundles for mediator atm/ocean flux computation
      fieldCount = med_fldList_GetNumFlds(med_fldList_getaofluxfldList())
      if ( fieldCount > 0 ) then
         if ( is_local%wrap%med_coupling_active(compocn,compatm) .or. &
              is_local%wrap%med_coupling_active(compatm,compocn)) then
            if ( is_local%wrap%aoflux_grid == 'ogrid' .and. .not. &
                 is_local%wrap%med_coupling_active(compatm,compocn)) then
               is_local%wrap%med_coupling_active(compatm,compocn) = .true.
            end if
            if ( is_local%wrap%aoflux_grid == 'agrid' .and. .not. &
                 is_local%wrap%med_coupling_active(compocn,compatm)) then
               is_local%wrap%med_coupling_active(compocn,compatm) = .true.
            end if
            call med_phases_aofluxes_init_fldbuns(gcomp, rc=rc)
         end if
      end if

      !---------------------------------------
      ! Second call to esmFldsExchange_xxx
      ! Determine mapping and merging info for field exchanges in mediator
      !---------------------------------------

      ! Initialize memory for fldlistFr(:)%flds(:) and fldlistTo(:)%flds(:) - this is needed for
      ! call below for the initialize phase

      if (trim(coupling_mode) == 'cesm') then
         call esmFldsExchange_cesm(gcomp, phase='initialize', rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      else if (coupling_mode(1:3) == 'ufs') then
         call esmFldsExchange_ufs(gcomp, phase='initialize', rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      else if (coupling_mode(1:4) == 'hafs') then
         call esmFldsExchange_hafs(gcomp, phase='initialize', rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if

      if (maintask) then
         call med_fldList_Document_Mapping(logunit, is_local%wrap%med_coupling_active)
         call med_fldList_Document_Merging(logunit, is_local%wrap%med_coupling_active)
      end if

      !---------------------------------------
      ! Initialize route handles and required normalization field bunds
      !---------------------------------------
      call ESMF_LogWrite("before med_map_RouteHandles_init", ESMF_LOGMSG_INFO)
      call med_map_RouteHandles_init(gcomp, is_local%wrap%flds_scalar_name, logunit, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_LogWrite("after  med_map_RouteHandles_init", ESMF_LOGMSG_INFO)

      !---------------------------------------
      ! Initialized packed field data structures
      !---------------------------------------
      do ndst = 1,ncomps
         do nsrc = 1,ncomps
            if (is_local%wrap%med_coupling_active(nsrc,ndst)) then
               call med_map_packed_field_create(ndst, &
                    is_local%wrap%flds_scalar_name, &
                    fieldsSrc=med_fldList_GetfldListFr(nsrc), &
                    FBSrc=is_local%wrap%FBImp(nsrc,nsrc), &
                    FBDst=is_local%wrap%FBImp(nsrc,ndst), &
                    packed_data=is_local%wrap%packed_data(nsrc,ndst,:), rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            end if
         end do
      end do

      if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o) .and. &
            ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_a)) then
          call med_map_packed_field_create(compatm, &
               is_local%wrap%flds_scalar_name, &
               fieldsSrc=med_fldList_getocnalbfldList(), &
               FBSrc=is_local%wrap%FBMed_ocnalb_o, &
               FBDst=is_local%wrap%FBMed_ocnalb_a, &
               packed_data=is_local%wrap%packed_data_ocnalb_o2a(:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
      !---------------------------------------
      ! Initialize ocn export accumulation field bundle
      !---------------------------------------
      if ( is_local%wrap%comp_present(compocn) .and. &
           ESMF_StateIsCreated(is_local%wrap%NStateImp(compocn),rc=rc) .and. &
           ESMF_StateIsCreated(is_local%wrap%NStateExp(compocn),rc=rc)) then
         call med_phases_prep_ocn_init(gcomp, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if

      !---------------------------------------
      ! Initialize wav export accumulation field bundle
      !---------------------------------------
      if ( is_local%wrap%comp_present(compwav) .and. &
           ESMF_StateIsCreated(is_local%wrap%NStateImp(compwav),rc=rc) .and. &
           ESMF_StateIsCreated(is_local%wrap%NStateExp(compwav),rc=rc)) then
         call med_phases_prep_wav_init(gcomp, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if

      !---------------------------------------
      ! Initialize glc module field bundles here if appropriate
      !---------------------------------------
      if (is_local%wrap%lnd2glc_coupling .or. is_local%wrap%ocn2glc_coupling .or. is_local%wrap%accum_lnd2glc) then
         call med_phases_prep_glc_init(gcomp, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if

      !---------------------------------------
      ! Initialize rof module field bundles here if appropriate
      !---------------------------------------
      if (is_local%wrap%med_coupling_active(comprof,complnd)) then
         call med_phases_prep_rof_init(gcomp, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if
      if (is_local%wrap%comp_present(comprof)) then
         call med_phases_post_rof_init(gcomp, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if
      !---------------------------------------
      ! Set the data initialize flag to false
      !---------------------------------------
      call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="false", rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      !---------------------------------------
      ! Set the first call flag to false
      !---------------------------------------

      first_call = .false.

      !---------------------------------------
      ! *** Now return ****
      !---------------------------------------

      ! The Connectors are being "called" for the transfer of Meshes
      ! (or Grids).  However, being "called" can mean different
      ! things! It can mean calling Initialization() phases, or Run()
      ! phases. For most of the initialization hand-shake, only
      ! Initialization() phases are called. This includes the entire
      ! GeomTransfer protocol. However, ONLY the Run phase of a
      ! Connector (full) transfers data AND timestamps!

      ! Once the first time DataInitialize() of CMEPS returns (below),
      ! and NUOPC sees that its InitializeDataComplete is not yet
      ! true, the NUOPC Driver will finally (for the first time!)
      ! execute the Run() phase of all of the Connectors that fit the
      ! *-TO-MED pattern. After that it will call CMEPS
      ! DataInitialize() again. Note that the time stamps are only set
      ! when the Run() phase of all the connectors are run.

      ! The Connectors Run() phase is called before the second call of
      ! the CMEPS DataInitialize phase.  As a result, CMEPS will see
      ! the correct timestamps, which also indicates that the actual
      ! data has been transferred reliably, and CMEPS can safely use it.

      RETURN

    endif  ! end first_call if-block

    !----------------------------------------------------------
    ! Create FBfrac field bundles and initialize fractions
    ! This has some complex dependencies on fractions from import States
    ! and appropriate checks are not implemented. We might need to split
    ! out the fraction FB allocation and the fraction initialization
    !----------------------------------------------------------

    call med_fraction_init(gcomp,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_fraction_set(gcomp,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------
    ! Initialize ocean albedos
    !----------------------------------------------------------

    if (is_local%wrap%comp_present(compocn) .or. is_local%wrap%comp_present(compatm)) then
       call med_phases_ocnalb_run(gcomp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    ! Loop over components and determine if they are at correct time
    !---------------------------------------

    do n1 = 1,ncomps
       compDone(n1) = .true. ! even if component is not present
       if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
          call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemNameList=fieldNameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,fieldCount
             call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemName=fieldNameList(n), field=field, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             atCorrectTime = NUOPC_IsAtTime(field, time, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if (.not. atCorrectTime) then
                compDone(n1) = .false.
             endif
          enddo
          deallocate(fieldNameList)
       endif
    enddo

    !---------------------------------------
    ! Carry out data dependency for atm initialization if needed
    !---------------------------------------

    if (is_local%wrap%comp_present(compatm)) then
       if (.not. compDone(compatm) .and. compDone(compocn)) then
          compDone(compatm) = .true.  ! reset if an item is found that is not done
          call ESMF_StateGet(is_local%wrap%NStateImp(compatm), itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(is_local%wrap%NStateImp(compatm), itemNameList=fieldNameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n=1, fieldCount
             call ESMF_StateGet(is_local%wrap%NStateImp(compatm), itemName=fieldNameList(n), field=field, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             atCorrectTime = NUOPC_IsAtTime(field, time, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if (.not. atCorrectTime) then
                ! If any atm import fields are not time stamped correctly,
                ! then dependency is not satisified - must return to atm
                call ESMF_LogWrite("MED - Initialize-Data-Dependency from ATM NOT YET SATISFIED!!!", &
                     ESMF_LOGMSG_INFO)
                if (maintask) then
                   write(logunit,'(A)') trim(subname)//"MED - Initialize-Data-Dependency from ATM NOT YET SATISFIED!!!"
                end if
                compDone(compatm) = .false.
                exit  ! break out of the loop when first not satisfied found
             endif
          enddo
          deallocate(fieldNameList)

          if (.not. compDone(compatm)) then  ! atmdone is not true
             if (is_local%wrap%comp_present(complnd)) then
                ! map initial lnd->atm
                call med_phases_post_lnd(gcomp, rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             ! do the merge to the atmospheric component
             call med_phases_prep_atm(gcomp, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! change 'Updated' attribute to true for ALL exportState fields
             call ESMF_StateGet(is_local%wrap%NStateExp(compatm), itemCount=fieldCount, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             allocate(fieldNameList(fieldCount))
             call ESMF_StateGet(is_local%wrap%NStateExp(compatm), itemNameList=fieldNameList, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             do n=1, fieldCount
                call ESMF_StateGet(is_local%wrap%NStateExp(compatm), itemName=fieldNameList(n), field=field, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end do
             deallocate(fieldNameList)

             ! Connectors will be automatically called between the mediator and atm until allDone is true
             call ESMF_LogWrite("MED - Initialize-Data-Dependency Sending Data to ATM", ESMF_LOGMSG_INFO)
          endif
       endif
    end if

    !---------------------------------------
    ! Loop over components again and determine if all are at the correct time
    !---------------------------------------

    allDone = .true.
    do n1 = 1,ncomps
       if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
          call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemNameList=fieldNameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n=1, fieldCount
             call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemName=fieldNameList(n), field=field, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             atCorrectTime = NUOPC_IsAtTime(field, time, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if (.not. atCorrectTime) then
                allDone=.false.
                if (dbug_flag > 0) then
                   if (maintask) then
                      write(logunit,'(A)') trim(subname)//" MED - Initialize-Data-Dependency check not yet satisfied for "//&
                           trim(compname(n1))
                   end if
                end if
             endif
          enddo
          deallocate(fieldNameList)
       endif
    enddo
    if (allDone) then
       ! set InitializeDataComplete Component Attribute to "true", indicating
       ! to the driver that this Component has fully initialized its data
       call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite("MED - Initialize-Data-Dependency allDone check Passed", ESMF_LOGMSG_INFO)
    end if

    !---------------------------------------
    ! Create component dimensions in mediator internal state
    !---------------------------------------

    if (allDone) then
       if (maintask) then
          write(logunit,*)
          write(logunit,'(a)') trim(subname)//"Initialize-Data-Dependency allDone check Passed"
       end if
       do n1 = 1,ncomps
          if (maintask) then
             write(logunit,*)
             write(logunit,'(a,2L2)') trim(subname)//" "//trim(compname(n1)), is_local%wrap%comp_present(n1), ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)
          end if
          if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
             call State_GetScalar(scalar_value=real_nx, &
                  scalar_id=is_local%wrap%flds_scalar_index_nx, &
                  state=is_local%wrap%NstateImp(n1), &
                  flds_scalar_name=is_local%wrap%flds_scalar_name, &
                  flds_scalar_num=is_local%wrap%flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call State_GetScalar(scalar_value=real_ny, &
                  scalar_id=is_local%wrap%flds_scalar_index_ny, &
                  state=is_local%wrap%NstateImp(n1), &
                  flds_scalar_name=is_local%wrap%flds_scalar_name, &
                  flds_scalar_num=is_local%wrap%flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if (is_local%wrap%flds_scalar_index_ntile > 0) then
                call State_GetScalar(scalar_value=real_ntile, &
                     scalar_id=is_local%wrap%flds_scalar_index_ntile, &
                     state=is_local%wrap%NstateImp(n1), &
                     flds_scalar_name=is_local%wrap%flds_scalar_name, &
                     flds_scalar_num=is_local%wrap%flds_scalar_num, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                is_local%wrap%ntile(n1) = nint(real_ntile)
             else
                is_local%wrap%ntile(n1) = 0
             end if
             is_local%wrap%nx(n1) = nint(real_nx)
             is_local%wrap%ny(n1) = nint(real_ny)
          endif
          if (is_local%wrap%comp_present(n1)) then
             write(msgString,'(3i8)') is_local%wrap%nx(n1), is_local%wrap%ny(n1), is_local%wrap%ntile(n1)
             if (maintask) then
                write(logunit,'(a)') 'global nx,ny,ntile sizes for '//trim(compname(n1))//":"//trim(msgString)
             end if
             call ESMF_LogWrite(trim(subname)//":"//trim(compname(n1))//":"//trim(msgString), ESMF_LOGMSG_INFO)
          endif
       end do
       if (maintask) write(logunit,*)

       !---------------------------------------
       ! Initialize mediator IO
       !---------------------------------------
       call med_io_init(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       ! Initialize mediator water/heat budget diags
       !---------------------------------------
       call med_diag_init(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_diag_zero(mode='all', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       ! read mediator restarts
       !---------------------------------------
       call NUOPC_CompAttributeGet(gcomp, name="read_restart", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (maintask) then
          write(logunit,*)
          write(logunit,'(a)') trim(subname)//' read_restart = '//trim(cvalue)
       end if
       read(cvalue,*) read_restart
       if (read_restart) then
         call med_phases_restart_read(gcomp, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       !---------------------------------------
       ! Call post routines as part of initialization
       !---------------------------------------
       if (is_local%wrap%comp_present(compatm)) then
          ! map atm->ocn, atm->ice, atm->lnd, atm->wav
          call med_phases_post_atm(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (is_local%wrap%comp_present(compice)) then
          ! call set ice_frac and map ice->ocn and ice->wav
          call med_phases_post_ice(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (allocated(compglc)) then
          ! map initial glc->lnd, glc->ocn and glc->ice
          call med_phases_post_glc(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (is_local%wrap%comp_present(complnd)) then
          ! map initial lnd->atm
          call med_phases_post_lnd(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (is_local%wrap%comp_present(compocn)) then
          ! map initial ocn->ice, ocn->wav
          call med_phases_post_ocn(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (is_local%wrap%comp_present(comprof)) then
          ! map initial rof->lnd, rof->ocn and rof->ice
          call med_phases_post_rof(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (is_local%wrap%comp_present(compwav)) then
          ! map initial wav->ocn, wav->ice, wav->atm
          call med_phases_post_wav(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       call med_phases_profile(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    else ! Not all done
       call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="false", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite("MED - Initialize-Data-Dependency allDone check not yet satisfied, another loop is required", &
            ESMF_LOGMSG_INFO)

    end if

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine DataInitialize

  !-----------------------------------------------------------------------------
  subroutine SetRunClock(gcomp, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_CLOCK, ESMF_Time, ESMF_TimeInterval
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_ClockGet, ESMF_ClockSet
    use ESMF                  , only : ESMF_Success, ESMF_Failure
    use ESMF                  , only : ESMF_Alarm, ESMF_ALARMLIST_ALL, ESMF_ClockGetAlarmList
    use ESMF                  , only : ESMF_AlarmCreate, ESMF_AlarmSet, ESMF_ClockAdvance
    use ESMF                  , only : ESMF_ClockGetAlarmList
    use NUOPC                 , only : NUOPC_CompCheckSetClock, NUOPC_CompAttributeGet
    use NUOPC_Mediator        , only : NUOPC_MediatorGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: mClock  ! mediator clock
    type(ESMF_CLock)        :: dClock  ! driver clock
    type(ESMF_Time)         :: currTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Alarm)        :: stop_alarm
    character(len=CL)       :: cvalue
    character(len=CL)       :: stop_option
    integer                 :: stop_n, stop_ymd
    logical, save           :: stopalarmcreated=.false.
    character(len=*), parameter :: subname = '('//__FILE__//':SetRunClock)'
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    ! query the Mediator for clocks
    call NUOPC_MediatorGet(gcomp, mediatorClock=mClock, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call Clock_TimePrint(dClock, trim(subname)//'driver clock1',rc)
       call Clock_TimePrint(mClock, trim(subname)//'mediat clock1',rc)
    endif

    ! set the mediatorClock to have the current start time as the driverClock
    call ESMF_ClockGet(dClock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockSet(mClock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call Clock_TimePrint(dClock, trim(subname)//'driver clock2',rc)
       call Clock_TimePrint(mClock, trim(subname)//'mediat clock2',rc)
    endif

    ! check and set the component clock against the driver clock
    call NUOPC_CompCheckSetClock(gcomp, dClock, checkTimeStep=.false., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.not. stopalarmcreated) then
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n
       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd
       call med_time_alarmInit(mclock, stop_alarm, stop_option, opt_n=stop_n, opt_ymd=stop_ymd, &
            alarmname='alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       stopalarmcreated = .true.
    end if

    ! Advance med clock to trigger alarms then reset model clock back to currtime
    call ESMF_ClockAdvance(mClock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockSet(mClock, currTime=currtime, timeStep=timestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    endif

  end subroutine SetRunClock

  !-----------------------------------------------------------------------------

  subroutine med_meshinfo_create(FB, mesh_info, FBArea, rc)

    use ESMF , only : ESMF_Array, ESMF_ArrayCreate, ESMF_ArrayDestroy, ESMF_Field, ESMF_FieldGet
    use ESMF , only : ESMF_DistGrid, ESMF_FieldBundle, ESMF_FieldRegridGetArea, ESMF_FieldBundleGet
    use ESMF , only : ESMF_Mesh, ESMF_MeshGet, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
    use ESMF , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF , only : ESMF_FieldCreate, ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
    use med_internalstate_mod , only : mesh_info_type

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    type(mesh_info_type)   , intent(inout) :: mesh_info
    type(ESMF_FieldBundle) , intent(inout) :: FBArea
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)      :: lfield
    type(ESMF_Mesh)       :: lmesh
    integer               :: numOwnedElements
    integer               :: spatialDim
    real(r8), allocatable :: ownedElemCoords(:)
    real(r8), pointer     :: dataptr(:)
    integer               :: n, dimcount, fieldcount
    character(len=*), parameter :: subname = '('//__FILE__//':med_meshinfo_create)'
    !-------------------------------------------------------------------------------

    rc= ESMF_SUCCESS

    ! Loop over number of fields and get mesh from the first field in FB with dimcount==1
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,fieldCount
       call FB_getFieldN(FB, fieldnum=n, field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, mesh=lmesh, dimcount=dimCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dimCount==1) exit
    enddo

    ! Determine dimensions in mesh
    call ESMF_MeshGet(lmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Obtain mesh info areas
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(mesh_info%areas(numOwnedElements))
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    mesh_info%areas(:) = dataptr

    ! Obtain mesh longitudes and latitudes
    allocate(mesh_info%lats(numOwnedElements))
    allocate(mesh_info%lons(numOwnedElements))
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    call ESMF_MeshGet(lmesh, ownedElemCoords=ownedElemCoords)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       mesh_info%lons(n) = ownedElemCoords(2*n-1)
       mesh_info%lats(n) = ownedElemCoords(2*n)
    end do
    deallocate(ownedElemCoords)

    ! Create field bundle with areas so that this can be output to mediator history file
    lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_r8, name='area', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    FBArea = ESMF_FieldBundleCreate(rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBArea, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = mesh_info%areas(:)

  end subroutine med_meshinfo_create

  !-----------------------------------------------------------------------------
  subroutine med_grid_write(grid, fileName, rc)

    use ESMF, only : ESMF_Grid, ESMF_Array, ESMF_ArrayBundle
    use ESMF, only : ESMF_ArrayBundleCreate, ESMF_GridGet
    use ESMF, only : ESMF_GridGetCoord, ESMF_ArraySet, ESMF_ArrayBundleAdd
    use ESMF, only : ESMF_GridGetItem, ESMF_ArrayBundleWrite, ESMF_ArrayBundleDestroy
    use ESMF, only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
    use ESMF, only : ESMF_SUCCESS, ESMF_GRIDITEM_MASK, ESMF_GRIDITEM_AREA

    ! input/output variables
    type(ESMF_Grid) , intent(in)  :: grid
    character(len=*), intent(in)  :: fileName
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_Array)       :: array
    type(ESMF_ArrayBundle) :: arrayBundle
    integer                :: tileCount
    logical                :: isPresent
    character(len=*), parameter :: subname = '('//__FILE__//':med_grid_write)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Check that the grid has tiles or not
    ! Currently only supports single tile
    call ESMF_GridGet(grid, tileCount=tileCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (tileCount .eq. 1) then
      ! Create arraybundle to store grid information
      arrayBundle = ESMF_ArrayBundleCreate(rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Query grid for center stagger
      ! Coordinates
      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
            isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
             staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArraySet(array, name="lon_center", rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayBundleAdd(arrayBundle, (/array/), rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_GridGetCoord(grid, coordDim=2, &
             staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArraySet(array, name="lat_center", rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayBundleAdd(arrayBundle, (/array/), rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      endif

      ! Mask
      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
           staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
             itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArraySet(array, name="mask_center", rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayBundleAdd(arrayBundle, (/array/), rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      endif

      ! Area
      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
           staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
             itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArraySet(array, name="area_center", rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayBundleAdd(arrayBundle, (/array/), rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      endif

      ! Query grid for corner stagger
      ! Coordinates
      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
            isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
             staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArraySet(array, name="lon_corner", rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayBundleAdd(arrayBundle, (/array/), rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_GridGetCoord(grid, coordDim=2, &
             staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArraySet(array, name="lat_corner", rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayBundleAdd(arrayBundle, (/array/), rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      endif

      ! Mask
      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
           staggerLoc=ESMF_STAGGERLOC_CORNER, isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
             itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArraySet(array, name="mask_corner", rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayBundleAdd(arrayBundle, (/array/), rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      endif

      ! Area
      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
           staggerLoc=ESMF_STAGGERLOC_CORNER, isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
             itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArraySet(array, name="area_corner", rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayBundleAdd(arrayBundle, (/array/), rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      endif

      ! Write arraybundle to file
      call ESMF_ArrayBundleWrite(arrayBundle, &
           fileName=trim(fileName), overwrite=.true., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Destroy arraybundle
      call ESMF_ArrayBundleDestroy(arrayBundle, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine med_grid_write

  !-----------------------------------------------------------------------------

  subroutine med_finalize(gcomp, rc)

    use ESMF, only : ESMF_GridComp, ESMF_SUCCESS

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    call memcheck("med_finalize", 0, maintask)
    if (maintask) then
       write(logunit,*)' SUCCESSFUL TERMINATION OF CMEPS'
       call med_phases_profile_finalize()
    end if

  end subroutine med_finalize

end module MED
