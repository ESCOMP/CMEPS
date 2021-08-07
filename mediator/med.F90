module MED

  !-----------------------------------------------------------------------------
  ! Mediator Component.
  !-----------------------------------------------------------------------------
  use ESMF                     , only : ESMF_VMLogMemInfo
  use NUOPC_Model              , only : SetVM
  use med_kind_mod             , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod        , only : dbug_flag          => med_constants_dbug_flag
  use med_constants_mod        , only : spval_init         => med_constants_spval_init
  use med_constants_mod        , only : spval              => med_constants_spval
  use med_constants_mod        , only : czero              => med_constants_czero
  use med_constants_mod        , only : ispval_mask        => med_constants_ispval_mask
  use med_utils_mod            , only : chkerr             => med_utils_ChkErr
  use med_methods_mod          , only : Field_GeomPrint    => med_methods_Field_GeomPrint
  use med_methods_mod          , only : State_GeomPrint    => med_methods_State_GeomPrint
  use med_methods_mod          , only : State_reset        => med_methods_State_reset
  use med_methods_mod          , only : State_getNumFields => med_methods_State_getNumFields
  use med_methods_mod          , only : State_GetScalar    => med_methods_State_GetScalar
  use med_methods_mod          , only : FB_Init            => med_methods_FB_init
  use med_methods_mod          , only : FB_Init_pointer    => med_methods_FB_Init_pointer
  use med_methods_mod          , only : FB_Reset           => med_methods_FB_Reset
  use med_methods_mod          , only : FB_FldChk          => med_methods_FB_FldChk
  use med_methods_mod          , only : FB_diagnose        => med_methods_FB_diagnose
  use med_methods_mod          , only : FB_getFieldN       => med_methods_FB_getFieldN
  use med_methods_mod          , only : clock_timeprint    => med_methods_clock_timeprint
  use med_time_mod             , only : alarmInit          => med_time_alarmInit
  use med_utils_mod            , only : memcheck           => med_memcheck
  use med_internalstate_mod    , only : InternalState
  use med_internalstate_mod    , only : med_coupling_allowed, logunit, mastertask
  use med_phases_profile_mod   , only : med_phases_profile_finalize
  use esmFlds                  , only : ncomps, compname
  use esmFlds                  , only : fldListFr, fldListTo, med_fldList_Realize
  use esmFlds                  , only : ncomps, compname, ncomps
  use esmFlds                  , only : compmed, compatm, compocn, compice, complnd, comprof, compwav ! not arrays
  use esmFlds                  , only : num_icesheets, max_icesheets, compglc, ocn2glc_coupling, lnd2glc_coupling ! compglc is an array
  use esmFlds                  , only : fldListMed_ocnalb
  use esmFlds                  , only : med_fldList_GetNumFlds, med_fldList_GetFldNames, med_fldList_GetFldInfo
  use esmFlds                  , only : med_fldList_Document_Mapping, med_fldList_Document_Merging
  use esmFlds                  , only : coupling_mode
  use esmFlds                  , only : med_name, atm_name, lnd_name, ocn_name
  use esmFlds                  , only : ice_name, rof_name, wav_name, glc_name
  use esmFldsExchange_nems_mod , only : esmFldsExchange_nems
  use esmFldsExchange_cesm_mod , only : esmFldsExchange_cesm
  use esmFldsExchange_hafs_mod , only : esmFldsExchange_hafs

  implicit none
  private

  public  SetServices
  public  SetVM
  private InitializeP0
  private InitializeIPDv03p1 ! advertise fields
  private InitializeIPDv03p3 ! realize connected Fields with transfer action "provide"
  private InitializeIPDv03p4 ! optionally modify the decomp/distr of transferred Grid/Mesh
  private InitializeIPDv03p5 ! realize all Fields with transfer action "accept"
  private DataInitialize     ! finish initialization and resolve data dependencies
  private SetRunClock
  private med_meshinfo_create
  private med_grid_write
  private med_finalize

  character(len=*), parameter :: grid_arbopt = "grid_reg"   ! grid_reg or grid_arb
  character(len=*), parameter :: u_FILE_u  = &
       __FILE__
  logical :: profile_memory = .false.

  character(len=8) :: atm_present, lnd_present
  character(len=8) :: ice_present, rof_present
  character(len=8) :: glc_present, med_present
  character(len=8) :: ocn_present, wav_present

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
    use med_phases_prep_wav_mod , only: med_phases_prep_wav
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

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname=' (module_MED:SetServices) '
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
         phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeIPDv03p1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p3: realize connected Fields with transfer action "provide"
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeIPDv03p3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p4: optionally modify the decomp/distr of transferred Grid/Mesh
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p4"/), userRoutine=InitializeIPDv03p4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p5: realize all Fields with transfer action "accept"
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p5"/), userRoutine=InitializeIPDv03p5, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! attach specializing method for DataInitialize
    !------------------

    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! setup mediator history phase
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_history_write"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_history_write", specRoutine=med_phases_history_write, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_history_write", specRoutine=NUOPC_NoOp, rc=rc)
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
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_post_ocn", specRoutine=NUOPC_NoOp, rc=rc)
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
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_post_ice", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep routines for lnd
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
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_post_lnd", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep and post routines for rof
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
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_post_rof", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep and post routines for wav
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_wav"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_wav", specRoutine=med_phases_prep_wav, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_post_wav"/), userRoutine=mediator_routine_Run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_post_wav", specRoutine=med_phases_post_wav, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_post_wav", specRoutine=NUOPC_NoOp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep and post routines for glc
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
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_post_glc", specRoutine=NUOPC_NoOp, rc=rc)
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
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_ocnalb_run", specRoutine=NUOPC_NoOp, rc=rc)
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
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_TimestampExport, &
         specPhaselabel="med_phases_aofluxes_run", specRoutine=NUOPC_NoOp, rc=rc)
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
    use med_internalstate_mod, only : mastertask, logunit, diagunit
    use esmFlds, only : dststatus_print

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
    character(len=CX) :: do_budgets
    character(len=*),parameter :: subname=' (module_MED:InitializeP0) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    mastertask = .false.
    if (localPet == 0) mastertask=.true.

    ! Determine mediator logunit
    if (mastertask) then
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
       open(newunit=logunit, file=trim(diro)//"/"//trim(logfile))

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
    if (mastertask) then
       write(logunit,'(a)')trim(subname)//": Mediator verbosity is set to "//trim(cvalue)
    end if

    ! Obtain profiling level
    call NUOPC_CompAttributeGet(gcomp, name="Profiling", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (mastertask) then
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

    ! Obtain dststatus_print setting if present
    call NUOPC_CompAttributeGet(gcomp, name='dststatus_print', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) dststatus_print=(trim(cvalue)=="true")
    write(msgString,*) trim(subname)//': Mediator dststatus_print is ',dststatus_print
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine InitializeP0

  !-----------------------------------------------------------------------

  subroutine InitializeIPDv03p1(gcomp, importState, exportState, clock, rc)

    ! Mediator advertises its import and export Fields and sets the
    ! TransferOfferGeomObject Attribute.

    use ESMF  , only : ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_SUCCESS, ESMF_LogFoundAllocError
    use ESMF  , only : ESMF_StateIsCreated
    use ESMF  , only : ESMF_LogMsg_Info, ESMF_LogWrite
    use ESMF  , only : ESMF_END_ABORT, ESMF_Finalize
    use NUOPC , only : NUOPC_AddNamespace, NUOPC_Advertise, NUOPC_AddNestedState
    use NUOPC , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=CS)   :: stdname, shortname
    integer             :: n, n1, n2, ncomp, nflds, ns
    logical             :: isPresent, isSet
    character(len=CS)   :: transferOffer
    character(len=CS)   :: cvalue
    character(len=8)    :: cnum
    type(InternalState) :: is_local
    integer             :: stat
    character(len=CS)   :: attrList(8)
    character(len=*),parameter :: subname=' (module_MED:InitializeIPDv03p1) '
    !-----------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    !------------------
    ! Allocate memory for the internal state and set it in the Component.
    !------------------

    allocate(is_local%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Allocation of the internal state memory failed.", line=__LINE__, file=u_FILE_u)) then
       return  ! bail out
    end if

    call ESMF_GridCompSetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

    ! Only create nested states for active ice sheets
    call NUOPC_CompAttributeGet(gcomp, name='num_icesheets', value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) num_icesheets
    else
       num_icesheets = 0
    end if
    do ns = 1,num_icesheets
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

    !------------------
    ! Initialize mediator flds
    !------------------

    call NUOPC_CompAttributeGet(gcomp, name='coupling_mode', value=coupling_mode, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('coupling_mode = '// trim(coupling_mode), ESMF_LOGMSG_INFO)
    if (mastertask) then
       write(logunit,*) '========================================================'
       write(logunit,'(a)')trim(subname)//' Mediator Coupling Mode is '//trim(coupling_mode)
       write(logunit,*) '========================================================'
       write(logunit,*)
    end if

    if (trim(coupling_mode) == 'cesm') then
       call esmFldsExchange_cesm(gcomp, phase='advertise', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (trim(coupling_mode) == 'nems_orig' .or. trim(coupling_mode) == 'nems_frac' &
       .or. trim(coupling_mode) == 'nems_orig_data') then
       call esmFldsExchange_nems(gcomp, phase='advertise', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (trim(coupling_mode(1:4)) == 'hafs') then
       call esmFldsExchange_hafs(gcomp, phase='advertise', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
        call ESMF_LogWrite(trim(coupling_mode)//' is not a valid coupling_mode', ESMF_LOGMSG_INFO)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    !------------------
    ! Determine component present indices
    !------------------

    call NUOPC_CompAttributeAdd(gcomp, &
         attrList=(/'atm_present','lnd_present','ocn_present','ice_present',&
                    'rof_present','wav_present','glc_present','med_present'/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    med_present = "false"
    atm_present = "false"
    lnd_present = "false"
    ocn_present = "false"
    ice_present = "false"
    rof_present = "false"
    wav_present = "false"
    glc_present = "false"

    ! Note that the present flag is set to true if the component is not stub
    call NUOPC_CompAttributeGet(gcomp, name='ATM_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'satm') atm_present = "true"
       atm_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='LND_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'slnd') lnd_present = "true"
       lnd_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='OCN_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'socn') ocn_present = "true"
       ocn_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ICE_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'sice') ice_present = "true"
       ice_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='ROF_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'srof') rof_present = "true"
       rof_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='WAV_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'swav') wav_present = "true"
       wav_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='GLC_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) /= 'sglc') glc_present = "true"
       glc_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='MED_model', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       med_name = trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='mediator_present', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       med_present = trim(cvalue)
    end if

    call NUOPC_CompAttributeSet(gcomp, name="atm_present", value=atm_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="lnd_present", value=lnd_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="ocn_present", value=ocn_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="ice_present", value=ice_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="rof_present", value=rof_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="wav_present", value=trim(wav_present), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="glc_present", value=trim(glc_present), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name="med_present", value=med_present, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mastertask) then
       write(logunit,*)
       if (trim(atm_present).eq."true") write(logunit,*) "atm_name="//trim(atm_name)
       if (trim(lnd_present).eq."true") write(logunit,*) "lnd_name="//trim(lnd_name)
       if (trim(ocn_present).eq."true") write(logunit,*) "ocn_name="//trim(ocn_name)
       if (trim(ice_present).eq."true") write(logunit,*) "ice_name="//trim(ice_name)
       if (trim(rof_present).eq."true") write(logunit,*) "rof_name="//trim(rof_name)
       if (trim(wav_present).eq."true") write(logunit,*) "wav_name="//trim(wav_name)
       if (trim(glc_present).eq."true") write(logunit,*) "glc_name="//trim(glc_name)
       if (trim(med_present).eq."true") write(logunit,*) "med_name="//trim(med_name)
       write(logunit,*)
    end if

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

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) is_local%wrap%flds_scalar_index_nextsw_cday
    else
       is_local%wrap%flds_scalar_index_nextsw_cday = spval
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
          if (mastertask) write(logunit,*)
          if (ESMF_StateIsCreated(is_local%wrap%NStateImp(ncomp))) then
             nflds = med_fldList_GetNumFlds(fldListFr(ncomp))
             do n = 1,nflds
                call med_fldList_GetFldInfo(fldListFr(ncomp), n, stdname, shortname)
                if (mastertask) then
                   write(logunit,'(a)') trim(subname)//':Fr_'//trim(compname(ncomp))//': '//trim(shortname)
                end if
                if (trim(shortname) == is_local%wrap%flds_scalar_name) then
                   transferOffer = 'will provide'
                else
                   transferOffer = 'cannot provide'
                end if
                call NUOPC_Advertise(is_local%wrap%NStateImp(ncomp), standardName=stdname, shortname=shortname, name=shortname, &
                     TransferOfferGeomObject=transferOffer, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call ESMF_LogWrite(subname//':Fr_'//trim(compname(ncomp))//': '//trim(shortname), ESMF_LOGMSG_INFO)
             end do
          end if
          if (ESMF_StateIsCreated(is_local%wrap%NStateExp(ncomp))) then
             nflds = med_fldList_GetNumFlds(fldListTo(ncomp))
             do n = 1,nflds
                call med_fldList_GetFldInfo(fldListTo(ncomp), n, stdname, shortname)
                if (mastertask) then
                   write(logunit,'(a)') trim(subname)//':To_'//trim(compname(ncomp))//': '//trim(shortname)
                end if
                if (trim(shortname) == is_local%wrap%flds_scalar_name) then
                   transferOffer = 'will provide'
                else
                   transferOffer = 'cannot provide'
                end if
                call NUOPC_Advertise(is_local%wrap%NStateExp(ncomp), standardName=stdname, shortname=shortname, name=shortname, &
                     TransferOfferGeomObject=transferOffer, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call ESMF_LogWrite(subname//':To_'//trim(compname(ncomp))//': '//trim(shortname), ESMF_LOGMSG_INFO)
             end do
          end if
       end if
    end do ! end of ncomps loop

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine InitializeIPDv03p1

  !-----------------------------------------------------------------------------

  subroutine InitializeIPDv03p3(gcomp, importState, exportState, clock, rc)

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
    character(len=*),parameter :: subname=' (module_MED:InitializeIPDv03p3) '
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
         call med_fldList_Realize(is_local%wrap%NStateImp(n), fldListFr(n), &
              is_local%wrap%flds_scalar_name, is_local%wrap%flds_scalar_num, &
              tag=subname//':Fr_'//trim(compname(n)), rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      if (ESMF_StateIsCreated(is_local%wrap%NStateExp(n), rc=rc)) then
          call ESMF_StateSet(is_local%wrap%NStateExp(n), stateIntent=ESMF_StateIntent_Export, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_fldList_Realize(is_local%wrap%NStateExp(n), fldListTo(n), &
              is_local%wrap%flds_scalar_name, is_local%wrap%flds_scalar_num, &
              tag=subname//':To_'//trim(compname(n)), rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    enddo

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine InitializeIPDv03p3

  !-----------------------------------------------------------------------------

  subroutine InitializeIPDv03p4(gcomp, importState, exportState, clock, rc)

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
    integer :: n1,n2
    character(len=*),parameter :: subname=' (module_MED:InitalizeIPDv03p4) '
    !-----------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    ! Get the internal state from the mediator gridded component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! Recieve Grids
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

      type(ESMF_State)   , intent(inout) :: State
      character(len=*)   , intent(in)    :: string
      integer            , intent(out)   :: rc

      ! local variables
      type(ESMF_Field)              :: field
      type(ESMF_Grid)               :: grid
      type(ESMF_Mesh)               :: mesh, newmesh
      integer                       :: localDeCount

      type(ESMF_DistGrid)           :: distgrid
      type(ESMF_DistGrid)           :: nodaldistgrid, newnodaldistgrid
      type(ESMF_DistGrid)           :: elemdistgrid, newelemdistgrid
      type(ESMF_DistGridConnection), allocatable :: connectionList(:)
      integer                       :: arbDimCount
      integer                       :: dimCount, tileCount, petCount
      integer                       :: connectionCount
      integer                       :: deCountPTile, extraDEs
      integer, allocatable          :: minIndexPTile(:,:), maxIndexPTile(:,:)
      integer, allocatable          :: regDecompPTile(:,:)
      integer                       :: i, j, n, n1, fieldCount, nxg, i1, i2
      type(ESMF_GeomType_Flag)      :: geomtype
      character(ESMF_MAXSTR),allocatable :: fieldNameList(:)
      type(ESMF_FieldStatus_Flag)   :: fieldStatus
      character(len=CX)             :: msgString
      character(len=*),parameter :: subname=' (module_MED:realizeConnectedGrid) '
      !-----------------------------------------------------------

      !NOTE: All of the Fields that set their TransferOfferGeomObject Attribute
      !NOTE: to "cannot provide" should now have the accepted Grid available.
      !NOTE: Go and pull out this Grid for one of a representative Field and
      !NOTE: modify the decomposition and distribution of the Grid to match the
      !NOTE: Mediator PETs.

      !TODO: quick implementation, do it for each field one by one
      !TODO: commented out below are application to other fields

      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
      rc = ESMF_Success
      if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

      call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      allocate(fieldNameList(fieldCount))
      call ESMF_StateGet(State, itemNameList=fieldNameList, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_GridCompGet(gcomp, petCount=petCount, rc=rc)
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

         !call NUOPC_GetAttribute(field, name="TransferActionGeomObject", &
         !     value=transferAction, rc=rc)
         !if (ChkErr(rc,__LINE__,u_FILE_u)) return

         if (fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then

            ! The Mediator is accepting a Grid/Mesh passed to it
            ! through the Connector

            ! While this is still an empty field, it does now hold a Grid/Mesh with DistGrid
            call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            if (geomtype == ESMF_GEOMTYPE_GRID) then

               !if (dbug_flag > 1) then
               !   call Field_GeomPrint(field,trim(fieldNameList(n))//'_orig',rc)
               !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
               !end if

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
                  !
                  ! Need to make a choice here to either represent the grid as a
                  ! regDecomp grid on the acceptor side, or to stay with arbDistr grid:
                  !
                  ! Setting the PRECIP_REGDECOMP macro will set up a regDecomp grid on the
                  ! acceptor side.
                  !
                  ! Not setting the PRECIP_REGDECOMP macro will default into keeping the
                  ! original arbDistr Grid.

                  if (grid_arbopt == "grid_reg") then

                     call ESMF_LogWrite(trim(subname)//trim(string)//": accept arb2reg grid for "//trim(fieldNameList(n)), &
                          ESMF_LOGMSG_INFO)

                     ! Use a regDecomp representation for the grid
                     ! first get tile min/max, only single tile supported for arbDistr Grid
                     allocate(minIndexPTile(arbDimCount,1),maxIndexPTile(arbDimCount,1))
                     call ESMF_AttributeGet(field, name="MinIndex", &
                          valueList=minIndexPTile(:,1), &
                          convention="NUOPC", purpose="Instance", rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return

                     call ESMF_AttributeGet(field, name="MaxIndex", &
                          valueList=maxIndexPTile(:,1), &
                          convention="NUOPC", purpose="Instance", rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return

                     ! create default regDecomp DistGrid
                     distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
                          maxIndexPTile=maxIndexPTile, connectionList=connectionList, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return

                     ! Create default regDecomp Grid
                     grid = ESMF_GridCreate(distgrid, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return

                     ! swap out the transferred grid for the newly created one
                     call ESMF_FieldEmptySet(field, grid=grid, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                     do i1 = 1,arbDimCount
                        write(msgString,'(A,3i8)') trim(subname)//':PTile =',i1,minIndexPTile(i1,1),maxIndexPTile(i1,1)
                        call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
                     enddo
                     deallocate(minIndexPTile,maxIndexPTile)

                  elseif (grid_arbopt == "grid_arb") then

                     ! Stick with the arbDistr representation of the grid:
                     ! There is nothing to do here if the same number of DEs is kept on the
                     ! acceptor side. Alternatively, the acceptor side could set up a more
                     ! natural number of DEs (maybe same number as acceptor PETs), and then
                     ! redistribute the arbSeqIndexList. Here simply keep the DEs of the
                     ! provider Grid.
                     call ESMF_LogWrite(trim(subname)//trim(string)//": accept arb2arb grid for "//trim(fieldNameList(n)), &
                          ESMF_LOGMSG_INFO)

                  else   ! grid_arbopt

                     call ESMF_LogWrite(trim(subname)//trim(string)//": ERROR grid_arbopt setting = "//trim(grid_arbopt), &
                          ESMF_LOGMSG_INFO)
                     rc = ESMF_FAILURE
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  endif  ! grid_arbopt


               else   ! arbdimcount <= 0

                  ! The provider defined as non arb grid

                  ! access localDeCount to show this is a real Grid
                  call ESMF_LogWrite(trim(subname)//trim(string)//": accept reg2reg grid for "//&
                       trim(fieldNameList(n)), ESMF_LOGMSG_INFO)

                  call ESMF_FieldGet(field, grid=grid, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  call ESMF_GridGet(grid, localDeCount=localDeCount, distgrid=distgrid, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  ! Create a custom DistGrid, based on the minIndex, maxIndex of the
                  ! accepted DistGrid, but with a default regDecomp for the current VM
                  ! that leads to 1DE/PET.

                  ! get dimCount and tileCount
                  call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, &
                       connectionCount=connectionCount, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
                  allocate(minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount))
                  allocate(connectionList(connectionCount))

                  ! get minIndex and maxIndex arrays, and connectionList
                  call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
                       maxIndexPTile=maxIndexPTile, connectionList=connectionList, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  ! construct a default regDecompPTile -> TODO: move this into ESMF as default

                  allocate(regDecompPTile(dimCount, tileCount))
                  deCountPTile = petCount/tileCount
                  extraDEs = max(0, petCount-deCountPTile)
                  do i=1, tileCount
                     if (i<=extraDEs) then
                        regDecompPTile(1, i) = deCountPTile + 1
                     else
                        regDecompPTile(1, i) = deCountPTile
                     endif
                     do j=2, dimCount
                        regDecompPTile(j, i) = 1
                     enddo
                  enddo

                  do i2 = 1,tileCount
                     do i1 = 1,dimCount
                        write(msgString,'(A,5i8)') trim(subname)//':PTile =',i2,i1,minIndexPTile(i1,i2),&
                             maxIndexPTile(i1,i2),regDecompPTile(i1,i2)
                        call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
                     enddo
                  enddo

                  !--- tcraig, hardwire i direction wraparound, temporary
                  !--- tcraig, now getting info from model distgrid, see above
                  !              allocate(connectionList(1))
                  !              nxg = maxIndexPTile(1,1) - minIndexPTile(1,1) + 1
                  !              write(msgstring,*) trim(subname)//trim(string),': connlist nxg = ',nxg
                  !              call ESMF_LogWrite(trim(msgstring), ESMF_LOGMSG_INFO)
                  !              if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  !              call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
                  !                tileIndexB=1, positionVector=(/nxg, 0/), rc=rc)
                  !              if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  ! create the new DistGrid with the same minIndexPTile and maxIndexPTile,
                  ! but with a default regDecompPTile
                  ! tcraig, force connectionlist and gridEdge arguments to fix wraparound
                  ! need ESMF fixes to implement properly.
                  if (dimcount == 2) then
                     distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
                          maxIndexPTile=maxIndexPTile, regDecompPTile=regDecompPTile, &
                          connectionList=connectionList, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return

                     call ESMF_LogWrite(trim(subname)//trim(string)//': distgrid with dimcount=2', ESMF_LOGMSG_INFO)

                     ! Create a new Grid on the new DistGrid and swap it in the Field
                     grid = ESMF_GridCreate(distgrid, gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  else
                     distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
                          maxIndexPTile=maxIndexPTile, regDecompPTile=regDecompPTile, rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return

                     call ESMF_LogWrite(trim(subname)//trim(string)//': distgrid with dimcount=1', ESMF_LOGMSG_INFO)

                     ! Create a new Grid on the new DistGrid and swap it in the Field
                     grid = ESMF_GridCreate(distgrid, gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/), rc=rc)
                     if (ChkErr(rc,__LINE__,u_FILE_u)) return
                  endif

                  ! local clean-up
                  deallocate(connectionList)
                  deallocate(minIndexPTile, maxIndexPTile, regDecompPTile)

               endif  ! arbdimCount

               ! Swap all the Grids in the State
               do n1=1, fieldCount
                  ! access a field in the State and set the Grid
                  call ESMF_StateGet(State, field=field, itemName=fieldNameList(n1), rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
                  if (ChkErr(rc,__LINE__,u_FILE_u)) return

                  if (fieldStatus==ESMF_FIELDSTATUS_EMPTY .or. fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then
                     call ESMF_FieldEmptySet(field, grid=grid, rc=rc)
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

               ! call ESMF_MeshGet(mesh, nodalDistGrid=nodalDistGrid, rc=rc)
               ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

               ! newnodalDistGrid = ESMF_DistGridCreate(nodalDistGrid, balanceflag=.true., rc=rc)
               ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

               ! Create a new Grid on the new DistGrid and swap it in the Field
               ! newmesh = ESMF_MeshEmptyCreate(elementDistGrid=newelemDistGrid, nodalDistGrid=newnodalDistGrid, rc=rc)
               ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

  end subroutine InitializeIPDv03p4

  !-----------------------------------------------------------------------------

  subroutine InitializeIPDv03p5(gcomp, importState, exportState, clock, rc)

    use ESMF , only : ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_LogWrite
    use ESMF , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_StateIsCreated

    !----------------------------------------------------------
    ! realize all Fields with transfer action "accept"
    !----------------------------------------------------------

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n1,n2
    character(len=*),parameter  :: subname=' (module_MED:InitializeIPDv03p5) '
    !-----------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    rc = ESMF_SUCCESS
    if (profile_memory) call ESMF_VMLogMemInfo("Entering "//trim(subname))

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--- Finish initializing the State Fields
    !--- Write out grid information

    do n1 = 1,ncomps

      if (ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
         call ESMF_LogWrite(trim(subname)//": calling completeFieldInitialize import states from "//trim(compname(n1)), &
              ESMF_LOGMSG_INFO)
        call completeFieldInitialization(is_local%wrap%NStateImp(n1), rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        call State_reset(is_local%wrap%NStateImp(n1), value=spval_init, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

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
      use ESMF  , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FieldGet, ESMF_FieldEmptyComplete
      use ESMF  , only : ESMF_GeomType_Flag, ESMF_FieldCreate, ESMF_MeshCreate, ESMF_GEOMTYPE_GRID
      use ESMF  , only : ESMF_MeshLoc_Element, ESMF_TYPEKIND_R8, ESMF_FIELDSTATUS_GRIDSET
      use ESMF  , only : ESMF_AttributeGet, ESMF_MeshWrite
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
      type(ESMF_Field),pointer    :: fieldList(:) => null()
      type(ESMF_FieldStatus_Flag) :: fieldStatus
      type(ESMF_GeomType_Flag)    :: geomtype
      integer                     :: gridToFieldMapCount, ungriddedCount
      integer, allocatable        :: gridToFieldMap(:)
      integer, allocatable        :: ungriddedLBound(:), ungriddedUBound(:)
      logical                     :: isPresent
      logical                     :: meshcreated
      character(len=*),parameter  :: subname=' (module_MED:completeFieldInitialization) '
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

          call ESMF_FieldGet(fieldList(n), status=fieldStatus, name=fieldName, &
            geomtype=geomtype, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

            meshField = ESMF_FieldCreate(mesh, typekind=ESMF_TYPEKIND_R8, &
                 meshloc=ESMF_MESHLOC_ELEMENT, name=fieldName, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Swap grid for mesh, at this point, only connected fields are in the state
            call NUOPC_Realize(State, field=meshField, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          if (fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then
            call ESMF_LogWrite(subname//" is allocating field memory for field "//trim(fieldName), &
                 ESMF_LOGMSG_INFO)

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
               call ESMF_AttributeGet(fieldList(n), name="UngriddedUBound", convention="NUOPC", &
                    purpose="Instance", valueList=ungriddedUBound, rc=rc)
            endif

            call ESMF_FieldEmptyComplete(fieldList(n), typekind=ESMF_TYPEKIND_R8, gridToFieldMap=gridToFieldMap, &
                 ungriddedLbound=ungriddedLbound, ungriddedUbound=ungriddedUbound, rc=rc)

            deallocate(gridToFieldMap, ungriddedLbound, ungriddedUbound)
          endif   ! fieldStatus

          call Field_GeomPrint(fieldList(n), trim(subname)//':'//trim(fieldName), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

        enddo
        deallocate(fieldList)
      endif

      if (profile_memory) call ESMF_VMLogMemInfo("Leaving "//trim(subname))
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    end subroutine completeFieldInitialization

  end subroutine InitializeIPDv03p5

  !-----------------------------------------------------------------------------

  subroutine DataInitialize(gcomp, rc)

    !----------------------------------------------------------
    ! Finish initialization and resolve data dependencies
    ! There will be multiple passes
    ! For first time through:
    !   Do not assume any import fields are connected, just allocate space and such
    !   -- Check present flags
    !   -- Check for active coupling interactions
    !   -- Create FBs: FBImp, FBExp, FBExpAccum
    !   -- Create mediator specific field bundles (not part of import/export states)
    !   -- Initialize FBExpAccums (to zero), and FBImp (from NStateImp)
    !   -- Read mediator restarts
    !   -- Initialize route handles
    !   -- Initialize field bundles for normalization
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
    use med_phases_prep_glc_mod , only : med_phases_prep_glc_init
    use med_phases_prep_atm_mod , only : med_phases_prep_atm
    use med_phases_post_atm_mod , only : med_phases_post_atm
    use med_phases_post_ice_mod , only : med_phases_post_ice
    use med_phases_post_lnd_mod , only : med_phases_post_lnd_init
    use med_phases_post_glc_mod , only : med_phases_post_glc
    use med_phases_post_ocn_mod , only : med_phases_post_ocn
    use med_phases_post_rof_mod , only : med_phases_post_rof
    use med_phases_post_wav_mod , only : med_phases_post_wav
    use med_phases_ocnalb_mod   , only : med_phases_ocnalb_run
    use med_phases_aofluxes_mod , only : med_phases_aofluxes_run, med_phases_aofluxes_init_fldbuns
    use med_phases_profile_mod  , only : med_phases_profile
    use med_diag_mod            , only : med_diag_zero, med_diag_init
    use med_map_mod             , only : med_map_mapnorm_init, med_map_routehandles_init, med_map_packed_field_create
    use med_io_mod              , only : med_io_init

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)                :: is_local
    type(ESMF_VM)                      :: vm
    type(ESMF_Clock)                   :: clock
    type(ESMF_State)                   :: importState, exportState
    type(ESMF_Time)                    :: time
    type(ESMF_Field)                   :: field
    type(ESMF_StateItem_Flag)          :: itemType
    logical                            :: atCorrectTime, connected
    integer                            :: n1,n2,n,ns
    integer                            :: nsrc,ndst
    integer                            :: cntn1, cntn2
    integer                            :: fieldCount
    character(ESMF_MAXSTR),allocatable :: fieldNameList(:)
    character(CL), pointer             :: fldnames(:) => null()
    character(CL)                      :: cvalue
    character(CL)                      :: cname
    character(CL)                      :: start_type
    logical                            :: read_restart
    logical                            :: allDone = .false.
    logical,save                       :: compDone(ncomps)
    logical,save                       :: first_call = .true.
    real(r8)                           :: real_nx, real_ny
    character(len=CX)                  :: msgString
    character(len=*), parameter        :: subname=' (module_MED:DataInitialize) '
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

      !----------------------------------------------------------
      ! Initialize mediator present flags
      !----------------------------------------------------------

      if (mastertask) then
         write(logunit,'(a)') trim(subname) // "Initializing present flags"
      end if

      do n1 = 1,ncomps
         cname = trim(compname(n1))
         if (cname(1:3) == 'glc') then
            ! Special logic for glc since there can be multiple ice sheets
            call ESMF_AttributeGet(gcomp, name="glc_present", value=cvalue, &
                 convention="NUOPC", purpose="Instance", rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            do ns = 1,max_icesheets
               if (ns <= num_icesheets) then
                  if (trim(cvalue) == 'true') then
                     is_local%wrap%comp_present(compglc(ns)) = .true.
                  else
                     is_local%wrap%comp_present(compglc(ns)) = .false.
                  end if
               end if
            end do
         else
            call ESMF_AttributeGet(gcomp, name=trim(compname(n1))//"_present", value=cvalue, &
                 convention="NUOPC", purpose="Instance", rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            if (trim(cvalue) == "true") then
               is_local%wrap%comp_present(n1) = .true.
            else
               is_local%wrap%comp_present(n1) = .false.
            end if
         end if
         if (mastertask) then
            write(msgString,'(A,L4)') trim(subname)//' comp_present(comp'//trim(compname(n1))//') = ',&
                 is_local%wrap%comp_present(n1)
            write(logunit,'(a)') trim(subname) // trim(msgString)
         end if
      end do

      !----------------------------------------------------------
      ! Check for active coupling interactions
      ! must be allowed, bundles created, and both sides have some fields
      !----------------------------------------------------------

      ! This defines the med_coupling_allowed is a starting point for what is
      ! allowed in this coupled system.  It will be revised further after the system
      ! starts, but any coupling set to false will never be allowed.
      ! are allowed, just update the table below.

      if (mastertask) then
         write(logunit,'(a)') trim(subname) // "Initializing active coupling flags"
      end if

      ! Initialize med_coupling_allowed
      med_coupling_allowed(:,:) = .false.

      ! to atmosphere
      med_coupling_allowed(complnd,compatm) = .true.
      med_coupling_allowed(compice,compatm) = .true.
      med_coupling_allowed(compocn,compatm) = .true.
      med_coupling_allowed(compwav,compatm) = .true.

      ! to land
      med_coupling_allowed(compatm,complnd) = .true.
      med_coupling_allowed(comprof,complnd) = .true.
      do ns = 1,num_icesheets
         med_coupling_allowed(compglc(ns),complnd) = .true.
      end do

      ! to ocean
      med_coupling_allowed(compatm,compocn) = .true.
      med_coupling_allowed(compice,compocn) = .true.
      med_coupling_allowed(comprof,compocn) = .true.
      med_coupling_allowed(compwav,compocn) = .true.
      do ns = 1,num_icesheets
         med_coupling_allowed(compglc(ns),compocn) = .true.
      end do

      ! to ice
      med_coupling_allowed(compatm,compice) = .true.
      med_coupling_allowed(compocn,compice) = .true.
      med_coupling_allowed(comprof,compice) = .true.
      med_coupling_allowed(compwav,compice) = .true.
      do ns = 1,num_icesheets
         med_coupling_allowed(compglc(ns),compice) = .true.
      end do

      ! to river
      med_coupling_allowed(complnd,comprof) = .true.

      ! to wave
      med_coupling_allowed(compatm,compwav) = .true.
      med_coupling_allowed(compocn,compwav) = .true.
      med_coupling_allowed(compice,compwav) = .true.

      ! to land-ice
      do ns = 1,num_icesheets
         med_coupling_allowed(complnd,compglc(ns)) = .true.
         med_coupling_allowed(compocn,compglc(ns)) = .true.
      end do

      ! initialize med_coupling_active table
      is_local%wrap%med_coupling_active(:,:) = .false.
      do n1 = 1,ncomps
        if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
          call State_GetNumFields(is_local%wrap%NStateImp(n1), cntn1, rc=rc) ! Import Field Count
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (cntn1 > 0) then
             do n2 = 1,ncomps
                if (is_local%wrap%comp_present(n2) .and. ESMF_StateIsCreated(is_local%wrap%NStateExp(n2),rc=rc) .and. &
                    med_coupling_allowed(n1,n2)) then
                   call State_GetNumFields(is_local%wrap%NStateExp(n2), cntn2, rc=rc) ! Import Field Count
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   if (cntn2 > 0) then
                      is_local%wrap%med_coupling_active(n1,n2) = .true.
                   endif
                endif
             enddo
          end if
        endif
      enddo

      ! Reset ocn2glc coupling based in input attribute
      if (.not. ocn2glc_coupling) then
         do ns = 1,num_icesheets
            is_local%wrap%med_coupling_active(compocn,compglc(ns)) = .false.
         end do
      end if

      ! create tables of allowed and active coupling flags
      ! - the rows are the destination of coupling
      ! - the columns are the source of coupling
      ! - So, the second column indicates which models the atm is coupled to.
      ! - And the second row indicates which models are coupled to the atm.
      if (mastertask) then
         write(logunit,*) ' '
         write(logunit,'(A)') trim(subname)//' Allowed coupling flags'
         write(logunit,'(2x,A10,20(A5))') '|from to->',(compname(n2),n2=1,ncomps)
         do n1 = 1,ncomps
            write(msgString,'(2x,a1,A,5x,20(L5))') '|',trim(compname(n1)), &
                 (med_coupling_allowed(n1,n2),n2=1,ncomps)
            do n2 = 1,len_trim(msgString)
               if (msgString(n2:n2) == 'F') msgString(n2:n2)='-'
            enddo
            write(logunit,'(A)') trim(msgString)
         enddo

         write(logunit,*) ' '
         write(logunit,'(A)') subname//' Active coupling flags'
         write(logunit,'(2x,A10,20(A5))') '|from to->',(compname(n2),n2=1,ncomps)
         do n1 = 1,ncomps
            write(msgString,'(2x,a1,A,5x,20(L5))') '|',trim(compname(n1)), &
                 (is_local%wrap%med_coupling_active(n1,n2),n2=1,ncomps)
            do n2 = 1,len_trim(msgString)
               if (msgString(n2:n2) == 'F') msgString(n2:n2)='-'
            enddo
            write(logunit,'(A)') trim(msgString)
         enddo
         write(logunit,*) ' '
      endif

      ! Reset the active coupling to permit atm/ocn flux computation on ogrid if needed
      if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a) .and. &
           ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o)) then
         if (  is_local%wrap%aoflux_grid == 'ogrid') then
            if ( is_local%wrap%med_coupling_active(compocn,compatm) .and. .not. &
                 is_local%wrap%med_coupling_active(compatm,compocn)) then
               is_local%wrap%med_coupling_active(compatm,compocn) = .true.
            end if
            if ( is_local%wrap%med_coupling_active(compatm,compocn) .and. .not. &
                 is_local%wrap%med_coupling_active(compocn,compatm)) then
               is_local%wrap%med_coupling_active(compocn,compatm) = .true.
            end if
         end if
      end if

      !----------------------------------------------------------
      ! Create field bundles FBImp, FBExp, FBImpAccum, FBExpAccum
      !----------------------------------------------------------

      if (mastertask) then
         write(logunit,'(a)') 'Creating mediator field bundles '
      end if

      do n1 = 1,ncomps
         if (is_local%wrap%comp_present(n1) .and. &
              ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc) .and. &
              ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc)) then

            if (mastertask) then
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

            ! Create import accumulation field bundles
            call FB_init(is_local%wrap%FBImpAccum(n1,n1), is_local%wrap%flds_scalar_name, &
                 STgeom=is_local%wrap%NStateImp(n1), STflds=is_local%wrap%NStateImp(n1), &
                 name='FBImpAccum'//trim(compname(n1)), rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            call FB_reset(is_local%wrap%FBImpAccum(n1,n1), value=czero, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Create export accumulation field bundles
            call FB_init(is_local%wrap%FBExpAccum(n1), is_local%wrap%flds_scalar_name, &
                 STgeom=is_local%wrap%NStateExp(n1), STflds=is_local%wrap%NStateExp(n1), &
                 name='FBExpAccum'//trim(compname(n1)), rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            call FB_reset(is_local%wrap%FBExpAccum(n1), value=czero, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Create mesh info data
            call ESMF_FieldBundleGet(is_local%wrap%FBImp(n1,n1), fieldCount=fieldCount, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            if (fieldCount == 0) then
              if (mastertask) then
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

         ! The following are FBImp and FBImpAccum mapped to different grids.
         ! FBImp(n1,n1) and FBImpAccum(n1,n1) are handled above
         do n2 = 1,ncomps
            if (n1 /= n2 .and. &
                 is_local%wrap%med_coupling_active(n1,n2) .and. &
                 ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc) .and. &
                 ESMF_StateIsCreated(is_local%wrap%NStateImp(n2),rc=rc)) then

               if (mastertask) then
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

               call FB_init(is_local%wrap%FBImpAccum(n1,n2), is_local%wrap%flds_scalar_name, &
                    STgeom=is_local%wrap%NStateImp(n2), &
                    STflds=is_local%wrap%NStateImp(n1), &
                    name='FBImpAccum'//trim(compname(n1))//'_'//trim(compname(n2)), rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return

               call FB_reset(is_local%wrap%FBImpAccum(n1,n2), value=czero, rc=rc)
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

      ! NOTE: this section must be done BEFORE the call to esmFldsExchange
      ! Create field bundles for mediator ocean albedo computation

      if ( is_local%wrap%med_coupling_active(compocn,compatm) .or. is_local%wrap%med_coupling_active(compatm,compocn)) then
         ! Create field bundles for mediator ocean albedo computation
         fieldCount = med_fldList_GetNumFlds(fldListMed_ocnalb)
         if (fieldCount > 0) then
            allocate(fldnames(fieldCount))
            call med_fldList_getfldnames(fldListMed_ocnalb%flds, fldnames, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return

            call FB_init(is_local%wrap%FBMed_ocnalb_a, is_local%wrap%flds_scalar_name, &
                 STgeom=is_local%wrap%NStateImp(compatm), fieldnamelist=fldnames, name='FBMed_ocnalb_a', rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            if (mastertask) then
               write(logunit,'(a)') trim(subname)//' initializing FB FBMed_ocnalb_a'
            end if

            call FB_init(is_local%wrap%FBMed_ocnalb_o, is_local%wrap%flds_scalar_name, &
                 STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames, name='FBMed_ocnalb_o', rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            if (mastertask) then
               write(logunit,'(a)') trim(subname)//' initializing FB FBMed_ocnalb_o'
            end if
            deallocate(fldnames)
         end if

         ! Create field bundles for mediator atm/ocn flux computation
         call med_phases_aofluxes_init_fldbuns(gcomp, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

      end if

      !---------------------------------------
      ! Determine mapping and merging info for field exchanges in mediator
      !---------------------------------------

      if (trim(coupling_mode) == 'cesm') then
         call esmFldsExchange_cesm(gcomp, phase='initialize', rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      else if (trim(coupling_mode) == 'hafs') then
         call esmFldsExchange_hafs(gcomp, phase='initialize', rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if

      if (mastertask) then
         call med_fldList_Document_Mapping(logunit, is_local%wrap%med_coupling_active)
         call med_fldList_Document_Merging(logunit, is_local%wrap%med_coupling_active)
      end if

      !---------------------------------------
      ! Initialize route handles and required normalization field bunds
      ! Initialized packed field data structures
      !---------------------------------------

      call ESMF_LogWrite("before med_map_RouteHandles_init", ESMF_LOGMSG_INFO)
      call med_map_RouteHandles_init(gcomp, is_local%wrap%flds_scalar_name, logunit, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_LogWrite("after  med_map_RouteHandles_init", ESMF_LOGMSG_INFO)

      call ESMF_LogWrite("before med_map_mapnorm_init", ESMF_LOGMSG_INFO)
      call med_map_mapnorm_init(gcomp, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_LogWrite("after  med_map_mapnorm_init", ESMF_LOGMSG_INFO)

      do ndst = 1,ncomps
         do nsrc = 1,ncomps
            if (is_local%wrap%med_coupling_active(nsrc,ndst)) then
                call med_map_packed_field_create(ndst, &
                     is_local%wrap%flds_scalar_name, &
                     fldsSrc=fldListFr(nsrc)%flds, &
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
               fldsSrc=fldListMed_ocnalb%flds, &
               FBSrc=is_local%wrap%FBMed_ocnalb_o, &
               FBDst=is_local%wrap%FBMed_ocnalb_a, &
               packed_data=is_local%wrap%packed_data_ocnalb_o2a(:), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

      !---------------------------------------
      ! Initialize glc module field bundles here if appropriate
      !---------------------------------------
      do ns = 1,num_icesheets
         if (is_local%wrap%med_coupling_active(complnd,compglc(ns))) then
            lnd2glc_coupling = .true.
            exit
         end if
      end do
      if (lnd2glc_coupling .or. ocn2glc_coupling) then
         call med_phases_prep_glc_init(gcomp, rc=rc)
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
    ! Initialize ocean albedos (this is needed for cesm and hafs)
    !----------------------------------------------------------

    if (trim(coupling_mode(1:5)) /= 'nems_') then
       if (is_local%wrap%comp_present(compocn) .or. is_local%wrap%comp_present(compatm)) then
          call med_phases_ocnalb_run(gcomp, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
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
                if (mastertask) then
                   write(logunit,'(A)') trim(subname)//"MED - Initialize-Data-Dependency from ATM NOT YET SATISFIED!!!"
                end if
                compDone(compatm) = .false.
                exit  ! break out of the loop when first not satisfied found
             endif
          enddo
          deallocate(fieldNameList)

          if (.not. compDone(compatm)) then  ! atmdone is not true
             if (trim(lnd_present) == 'true') then
                ! map initial lnd->atm
                call med_phases_post_lnd_init(gcomp, rc)
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
                   if (mastertask) then
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
       if (mastertask) then
          write(logunit,*)
          write(logunit,'(a)') trim(subname)//"Initialize-Data-Dependency allDone check Passed"
       end if
       do n1 = 1,ncomps
          if (mastertask) then
          write(logunit,*)
          write(logunit,'(a)') trim(subname)//" "//trim(compname(n1))
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
             is_local%wrap%nx(n1) = nint(real_nx)
             is_local%wrap%ny(n1) = nint(real_ny)
             write(msgString,'(2i8,2l4)') is_local%wrap%nx(n1), is_local%wrap%ny(n1)
             if (mastertask) then
                write(logunit,*) 'global nx,ny sizes for '//trim(compname(n1))//":"//trim(msgString)
             end if
             call ESMF_LogWrite(trim(subname)//":"//trim(compname(n1))//":"//trim(msgString), ESMF_LOGMSG_INFO)
          end if
       end do
       if (mastertask) write(logunit,*)

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
       if (mastertask) then
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
       if (trim(atm_present) == 'true') then
          ! map atm->ocn, atm->ice, atm->lnd
          call med_phases_post_atm(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (trim(ice_present) == 'true') then
          ! call set ice_frac and map ice->atm and ice->ocn
          call med_phases_post_ice(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (trim(glc_present) == 'true') then
          ! map initial glc->lnd, glc->ocn and glc->ice
          call med_phases_post_glc(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (trim(lnd_present) == 'true') then
          ! map initial lnd->atm
          call med_phases_post_lnd_init(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (trim(ocn_present) == 'true') then
          ! map initial ocn->ice
          call med_phases_post_ocn(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (trim(rof_present) == 'true') then
          ! map initial rof->lnd, rof->ocn and rof->ice
          call med_phases_post_rof(gcomp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (trim(wav_present) == 'true') then
          ! map initial wav->ocn and wav->ice
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
    type(ESMF_Clock)        :: mediatorClock, driverClock
    type(ESMF_Time)         :: currTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Alarm)        :: stop_alarm
    character(len=CL)       :: cvalue
    character(len=CL)       :: name, stop_option
    integer                 :: stop_n, stop_ymd
    logical                 :: first_time = .true.
    logical, save           :: stopalarmcreated=.false.
    integer                 :: alarmcount

    character(len=*),parameter :: subname=' (module_MED:SetRunClock) '
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    ! query the Mediator for clocks
    call NUOPC_MediatorGet(gcomp, mediatorClock=mediatorClock, driverClock=driverClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call Clock_TimePrint(driverClock  ,trim(subname)//'driver clock1',rc)
       call Clock_TimePrint(mediatorClock,trim(subname)//'mediat clock1',rc)
    endif

    ! set the mediatorClock to have the current start time as the driverClock
    call ESMF_ClockGet(driverClock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockSet(mediatorClock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call Clock_TimePrint(driverClock  ,trim(subname)//'driver clock2',rc)
       call Clock_TimePrint(mediatorClock,trim(subname)//'mediat clock2',rc)
    endif

    ! check and set the component clock against the driver clock
    call NUOPC_CompCheckSetClock(gcomp, driverClock, checkTimeStep=.false., rc=rc)
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
       call alarmInit(mediatorclock, stop_alarm, stop_option, opt_n=stop_n, opt_ymd=stop_ymd, &
            alarmname='alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       stopalarmcreated = .true.
    end if

    !--------------------------------
    ! Advance med clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mediatorClock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mediatorClock, currTime=currtime, timeStep=timestep, rc=rc)
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
    real(r8), pointer     :: dataptr(:) => null()
    integer               :: n, dimcount, fieldcount
    character(len=*),parameter :: subname=' (module_MED:med_meshinfo_create) '
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
    character(len=*), parameter :: subname=' (module_MED_map:med_grid_write) '
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
    call memcheck("med_finalize", 0, mastertask)
    if (mastertask) then
       write(logunit,*)' SUCCESSFUL TERMINATION OF CMEPS'
       call med_phases_profile_finalize()
    end if

  end subroutine med_finalize

end module MED
