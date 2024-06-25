module ESM

  !-----------------------------------------------------------------------------
  ! Code that specializes generic ESM Component code.
  !-----------------------------------------------------------------------------

  use shr_kind_mod , only : r8=>shr_kind_r8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod  , only : shr_sys_abort
  use shr_mpi_mod  , only : shr_mpi_bcast
  use shr_mem_mod  , only : shr_mem_init
  use shr_log_mod  , only : shr_log_setLogunit
  use esm_utils_mod, only : logunit, maintask, dbug_flag, chkerr

  implicit none
  private

  public  :: SetServices
  public  :: ReadAttributes ! used in ensemble_driver

  private :: SetModelServices
  private :: SetRunSequence
  private :: ModifyCplLists
  private :: InitAttributes
  private :: CheckAttributes
  private :: AddAttributes
  private :: InitAdvertize
  private :: esm_init_pelayout
  private :: esm_finalize
  private :: pretty_print_nuopc_freeformat

  real(r8) :: scol_spval = -999._r8

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine SetServices(driver, rc)

    use NUOPC        , only : NUOPC_CompDerive, NUOPC_CompSpecialize, NUOPC_CompSetInternalEntryPoint
    use NUOPC_Driver , only : driver_routine_SS             => SetServices
    use NUOPC_Driver , only : driver_label_SetModelServices => label_SetModelServices
    use NUOPC_Driver , only : driver_label_SetRunSequence   => label_SetRunSequence
    use NUOPC_Driver , only : driver_label_Finalize         => label_Finalize
    use ESMF         , only : ESMF_GridComp, ESMF_Config, ESMF_GridCompSet, ESMF_ConfigLoadFile
    use ESMF         , only : ESMF_METHOD_INITIALIZE
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO

    ! input/output variables
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname = "(esm.F90:SetServices)"
    !---------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    ! NUOPC_Driver registers the generic methods
    call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetRunSequence, &
         specRoutine=SetRunSequence, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! register an internal initialization method
    call NUOPC_CompSetInternalEntryPoint(driver, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p2"/), userRoutine=ModifyCplLists, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !
    ! This prevents the driver trying to "auto" connect to the ensemble_driver
    ! by default the FieldTransferPolicy is "transferall" and we need "transfernone"
    !
    call NUOPC_CompSetInternalEntryPoint(driver, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv05p1"/), userRoutine=InitAdvertize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Set a finalize method
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_Finalize, &
         specRoutine=esm_finalize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create, open and set the config

    call ESMF_GridCompSet(driver, configFile="nuopc.runconfig", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !================================================================================

  subroutine SetModelServices(driver, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_VM, ESMF_Config, ESMF_VMBarrier
    use ESMF         , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_ConfigGetAttribute
    use ESMF         , only : ESMF_ConfigGetLen, ESMF_RC_NOT_VALID, ESMF_LogFoundAllocError
    use ESMF         , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF         , only : ESMF_GridCompSet, ESMF_SUCCESS, ESMF_METHOD_INITIALIZE
    use ESMF         , only : ESMF_VMisCreated, ESMF_GridCompIsPetLocal
    use ESMF         , only : ESMF_RC_FILE_OPEN, ESMF_RC_FILE_READ
    use ESMF         , only : ESMF_AttributeUpdate, ESMF_VMBroadcast
    use ESMF         , only : ESMF_MethodAdd
    use NUOPC        , only : NUOPC_CompSetInternalEntryPoint, NUOPC_CompAttributeGet
    use NUOPC        , only : NUOPC_CompAttributeAdd, NUOPC_CompAttributeSet
    use NUOPC_Driver , only : NUOPC_DriverAddComp, NUOPC_DriverGetComp

    ! input/output variables
    type(ESMF_GridComp)    :: driver
    integer, intent(out)   :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    type(ESMF_Config) :: config
    integer           :: localPet
    character(len=CL) :: meminitStr
    integer           :: global_comm
    integer           :: maxthreads
    character(len=CL) :: msgstr
    integer           :: componentcount
    character(len=*), parameter :: subname = "(esm.F90:SetModelServices)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    !-------------------------------------------
    ! Set the io logunit to the value defined in ensemble_driver
    !-------------------------------------------
    call shr_log_setLogunit(logunit)

    !-------------------------------------------
    ! Get the config and vm objects from the driver
    !-------------------------------------------

    call ESMF_GridCompGet(driver, vm=vm, config=config, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=global_comm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (localPet == 0) then
       maintask=.true.
    else
       maintask = .false.
    end if

    !-------------------------------------------
    ! determine the generic component labels
    !-------------------------------------------

    componentCount = ESMF_ConfigGetLen(config,label="component_list:", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (componentCount == 0) then
      write (msgstr, *) "No models were specified in component_list "
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif

    !-------------------------------------------
    ! Obtain driver attributes
    !-------------------------------------------

    call ReadAttributes(driver, config, "DRIVER_attributes::", formatprint=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "CLOCK_attributes::", formatprint=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "ALLCOMP_attributes::", formatprint=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "PELAYOUT_attributes::", formatprint=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call CheckAttributes(driver, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Initialize other attributes (after initializing driver clock)
    !-------------------------------------------

    call InitAttributes(driver, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Initialize component pe layouts
    !-------------------------------------------

    call esm_init_pelayout(driver, maxthreads, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Memory test
    if (maintask) then
       call shr_mem_init(strbuf=meminitstr)
       write(logunit,*) trim(meminitstr)
    end if

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine SetModelServices

  !================================================================================

  subroutine SetRunSequence(driver, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
    use ESMF         , only : ESMF_Config
    use ESMF         , only : ESMF_GridCompGet, ESMF_ConfigCreate
    use ESMF         , only : ESMF_ConfigLoadFile
    use NUOPC        , only : NUOPC_FreeFormat, NUOPC_FreeFormatDestroy
    use NUOPC        , only : NUOPC_FreeFormatCreate
    use NUOPC_Driver , only : NUOPC_DriverIngestRunSequence, NUOPC_DriverSetRunSequence
    use NUOPC_Driver , only : NUOPC_DriverPrint

    ! input/output variables
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Config)       :: runSeq
    type(NUOPC_FreeFormat)  :: runSeqFF
    character(len=*), parameter :: subname = "(esm.F90:SetRunSequence)"
    !---------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    call shr_log_setLogunit(logunit)

    !--------
    ! Run Sequence and Connectors
    !--------

    ! read free format run sequence

    runSeq = ESMF_ConfigCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ConfigLoadFile(runSeq, "nuopc.runseq", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    runSeqFF = NUOPC_FreeFormatCreate(runSeq, label="runSeq::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_DriverIngestRunSequence(driver, runSeqFF, autoAddConnectors=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
#ifdef DEBUG
    ! Uncomment these to add debugging information for driver
    ! call NUOPC_DriverPrint(driver, orderflag=.true.)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !   line=__LINE__, &
    !   file=__FILE__)) &
    !   return  ! bail out

!    call pretty_print_nuopc_freeformat(runSeqFF, 'run sequence', rc=rc)
!    if (chkerr(rc,__LINE__,u_FILE_u)) return
#endif
    call NUOPC_FreeFormatDestroy(runSeqFF, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine SetRunSequence

  !================================================================================

  subroutine pretty_print_nuopc_freeformat(ffstuff, label, rc)

    use NUOPC, only : NUOPC_FreeFormat, NUOPC_FreeFormatGet, NUOPC_FreeFormatLen
    use ESMF,  only : ESMF_SUCCESS

    ! input/output variables
    type(NUOPC_FreeFormat) , intent(in)  :: ffstuff
    character(len=*)       , intent(in)  :: label
    integer                , intent(out) :: rc

    ! local variables
    integer :: i
    integer :: linecnt
    character(len=NUOPC_FreeFormatLen), pointer :: outstr(:)
    !---------------------------------------

    rc = ESMF_SUCCESS

    if (maintask .or. dbug_flag > 3) then
       write(logunit, *) 'BEGIN: ', trim(label)
       call NUOPC_FreeFormatGet(ffstuff, linecount=linecnt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(outstr(linecnt))
       call NUOPC_FreeFormatGet(ffstuff, stringList=outstr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do i=1,linecnt
          if(len_trim(outstr(i)) > 0) then
             write(logunit, *) trim(outstr(i))
          endif
       enddo
       write(logunit, *) 'END: ', trim(label)
       deallocate(outstr)
    endif

  end subroutine pretty_print_nuopc_freeformat

  !================================================================================

  recursive subroutine ModifyCplLists(driver, importState, exportState, clock, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_LogWrite
    use ESMF         , only : ESMF_LOGMSG_INFO, ESMF_CplComp, ESMF_SUCCESS
    use NUOPC        , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
    use NUOPC_Driver , only : NUOPC_DriverGetComp

    type(ESMF_GridComp)  :: driver
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_CplComp), pointer    :: connectorList(:)
    integer                        :: i, j, cplListSize
    character(len=CL), allocatable :: cplList(:)
    character(len=CL)              :: tempString
    character(len=CL)              :: msgstr
    character(len=*), parameter    :: subname = "(esm.F90:ModifyCplLists)"
    !---------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    call shr_log_setLogunit(logunit)

    call ESMF_LogWrite("Driver is in ModifyCplLists()", ESMF_LOGMSG_INFO, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    nullify(connectorList)
    call NUOPC_DriverGetComp(driver, compList=connectorList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    write (msgstr,*) "Found ", size(connectorList), " Connectors."// " Modifying CplList Attribute...."
    call ESMF_LogWrite(trim(msgstr), ESMF_LOGMSG_INFO, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do i=1, size(connectorList)

      ! query the cplList for connector i
      call NUOPC_CompAttributeGet(connectorList(i), name="CplList", itemCount=cplListSize, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (cplListSize>0) then
        allocate(cplList(cplListSize))

        call NUOPC_CompAttributeGet(connectorList(i), name="CplList", valueList=cplList, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        ! go through all of the entries in the cplList and set the mapping method to "redist"
        do j=1, cplListSize
           !tempString = trim(cplList(j))//":REMAPMETHOD=bilinear"//&
           !":SrcTermProcessing=1:DUMPWEIGHTS=true:TermOrder=SrcSeq"

           tempString = trim(cplList(j))//":remapmethod=redist"
           cplList(j) = trim(tempString)
        enddo

        ! store the modified cplList in CplList attribute of connector i
        call NUOPC_CompAttributeSet(connectorList(i), name="CplList", valueList=cplList, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        deallocate(cplList)
      endif
    enddo

    deallocate(connectorList)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine ModifyCplLists

  !================================================================================

  subroutine InitAttributes(driver, rc)

    use ESMF             , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF             , only : ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
    use ESMF             , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LogSetError, ESMF_LOGMSG_INFO
    use ESMF             , only : ESMF_RC_NOT_VALID
    use ESMF             , only : ESMF_GridCompIsPetLocal, ESMF_VMBroadcast
    use NUOPC            , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use shr_assert_mod   , only : shr_assert_in_domain
    use shr_const_mod    , only : shr_const_tkfrz, shr_const_tktrip
    use shr_const_mod    , only : shr_const_mwwv, shr_const_mwdair
    use shr_frz_mod      , only : shr_frz_freezetemp_init
    use shr_reprosum_mod , only : shr_reprosum_setopts
    use shr_wv_sat_mod   , only : shr_wv_sat_set_default, shr_wv_sat_init
    use shr_wv_sat_mod   , only : shr_wv_sat_make_tables, ShrWVSatTableSpec
    use shr_wv_sat_mod   , only : shr_wv_sat_get_scheme_idx, shr_wv_sat_valid_idx
    use glc_elevclass_mod, only : glc_elevclass_init
   !use shr_scam_mod     , only : shr_scam_checkSurface

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: driver
    integer             , intent(out)   :: rc                 ! return code

    ! local variables
    character(len=CL)            :: errstring
    character(len=CL)            :: cvalue
    logical                      :: reprosum_use_ddpdd    ! setup reprosum, use ddpdd
    real(R8)                     :: reprosum_diffmax      ! setup reprosum, set rel_diff_max
    logical                      :: reprosum_recompute    ! setup reprosum, recompute if tolerance exceeded
    character(LEN=CS)            :: tfreeze_option        ! Freezing point calculation
    integer                      :: glc_nec               ! number of elevation classes in the land component for lnd->glc
    character(LEN=CS)            :: wv_sat_scheme
    real(R8)                     :: wv_sat_transition_start
    logical                      :: wv_sat_use_tables
    real(R8)                     :: wv_sat_table_spacing
    type(ShrWVSatTableSpec)      :: liquid_spec
    type(ShrWVSatTableSpec)      :: ice_spec
    type(ShrWVSatTableSpec)      :: mixed_spec
    integer                      :: localPet, rootpe_med
    integer          , parameter :: ens1=1                ! use first instance of ensemble only
    integer          , parameter :: fix1=1                ! temporary hard-coding to first ensemble, needs to be fixed
    real(R8)         , parameter :: epsilo = shr_const_mwwv/shr_const_mwdair
    character(len=*) , parameter :: subname = '(InitAttributes)'
    !----------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    call shr_log_setLogunit(logunit)

    !----------------------------------------------------------
    ! Initialize options for reproducible sums
    ! TODO: this needs to be moved out of here
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="reprosum_use_ddpdd", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) reprosum_use_ddpdd

    call NUOPC_CompAttributeGet(driver, name="reprosum_diffmax", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) reprosum_diffmax

    call NUOPC_CompAttributeGet(driver, name="reprosum_recompute", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) reprosum_recompute

    call shr_reprosum_setopts(repro_sum_use_ddpdd_in=reprosum_use_ddpdd, &
         repro_sum_rel_diff_max_in=reprosum_diffmax, repro_sum_recompute_in=reprosum_recompute)

    !----------------------------------------------------------
    ! Initialize freezing point calculation for all components
    ! TODO: this needs to be moved out of here
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="tfreeze_option", value=tfreeze_option, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call shr_frz_freezetemp_init(tfreeze_option, maintask)

    call NUOPC_CompAttributeGet(driver, name='cpl_rootpe', value=cvalue, rc=rc)
    read(cvalue, *) rootpe_med
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_GridCompGet(driver, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------
    ! Initialize glc_elevclass_mod module variables
    !----------------------------------------------------------

    ! This must be called on all processors of the driver
    call NUOPC_CompAttributeGet(driver, name='glc_nec', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_nec
    call glc_elevclass_init(glc_nec, logunit=logunit)

    !----------------------------------------------------------
    ! Initialize water vapor info
    !----------------------------------------------------------

    ! TODO: this does not seem to belong here - where should it go?

    call NUOPC_CompAttributeGet(driver, name="wv_sat_scheme", value=wv_sat_scheme, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (.not. shr_wv_sat_valid_idx(shr_wv_sat_get_scheme_idx(trim(wv_sat_scheme)))) then
       call shr_sys_abort(subname//': "'//trim(wv_sat_scheme)//'" is not a recognized saturation vapor pressure scheme name')
    end if
    if (.not. shr_wv_sat_set_default(wv_sat_scheme)) then
       call shr_sys_abort('Invalid wv_sat_scheme.')
    end if

    call NUOPC_CompAttributeGet(driver, name="wv_sat_transition_start", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wv_sat_transition_start

    call shr_assert_in_domain(wv_sat_transition_start, &
         ge=0._R8, le=40._R8, &
         varname="wv_sat_transition_start", msg="Invalid transition temperature range.")

    call NUOPC_CompAttributeGet(driver, name="wv_sat_use_tables", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wv_sat_use_tables

    call NUOPC_CompAttributeGet(driver, name="wv_sat_table_spacing", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wv_sat_table_spacing

    ! A transition range averaging method in CAM is only valid for:
    ! -40 deg C <= T <= 0 deg C
    ! shr_wv_sat_mod itself checks for values with the wrong sign, but we
    ! have to check that the range is no more than 40 deg C here. Even
    ! though this is a CAM-specific restriction, it's not really likely
    ! that any other parameterization will be dealing with mixed-phase
    ! water below 40 deg C anyway.

    call shr_wv_sat_init(shr_const_tkfrz, shr_const_tktrip, wv_sat_transition_start, epsilo, errstring)
    if (errstring /= "") then
       call shr_sys_abort('shr_wv_sat_init: '//trim(errstring))
    end if

    ! The below produces internal lookup tables in the range 175-374K for
    ! liquid water, and 125-274K for ice, with a resolution set by the
    ! option wv_sat_table_spacing.
    ! In theory these ranges could be specified in the namelist, but in
    ! practice users will want to change them *very* rarely if ever, which
    ! is why only the spacing is in the namelist.

    if (wv_sat_use_tables) then
       liquid_spec = ShrWVSatTableSpec(ceiling(200._R8/wv_sat_table_spacing), 175._R8, wv_sat_table_spacing)
       ice_spec    = ShrWVSatTableSpec(ceiling(150._R8/wv_sat_table_spacing), 125._R8, wv_sat_table_spacing)
       mixed_spec  = ShrWVSatTableSpec(ceiling(250._R8/wv_sat_table_spacing), 125._R8, wv_sat_table_spacing)
       call shr_wv_sat_make_tables(liquid_spec, ice_spec, mixed_spec)
    end if

  end subroutine InitAttributes

  !================================================================================

  subroutine CheckAttributes( driver, rc )

    ! !DESCRIPTION: Check that input driver config values have reasonable values

    use ESMF        , only : ESMF_GridComp, ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use NUOPC       , only : NUOPC_CompAttributeGet

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp) , intent(inout) :: driver
    integer             , intent(out)   :: rc

    !----- local -----
    character(len=CS) :: logFilePostFix ! postfix for output log files
    character(len=CL) :: outPathRoot    ! root for output log files
    character(len=CS) :: cime_model
    character(len=*), parameter :: subname = '(driver_attributes_check) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(driver, name="cime_model", value=cime_model, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if ( trim(cime_model) /= 'cesm' .and. trim(cime_model) /= 'ufs') then
       call shr_sys_abort( subname//': cime_model must be set to cesm or ufs, aborting')
    end if

    ! --- LogFile ending name -----
    call NUOPC_CompAttributeGet(driver, name="logFilePostFix", value=logFilePostFix, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if ( len_trim(logFilePostFix) == 0 ) then
       call shr_sys_abort( subname//': logFilePostFix  must be set to something not blank' )
    end if

    ! --- Output path root directory -----
    call NUOPC_CompAttributeGet(driver, name="outPathRoot", value=outPathRoot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if ( len_trim(outPathRoot) == 0 ) then
       call shr_sys_abort( subname//': outPathRoot  must be set' )
    end if
    if ( index(outPathRoot, "/", back=.true.) /= len_trim(outPathRoot) ) then
       call shr_sys_abort( subname//': outPathRoot must end with a slash' )
    end if

  end subroutine CheckAttributes

  !===============================================================================

  subroutine AddAttributes(gcomp, driver, config, compid, compname, inst_suffix, nthrds, rc)

    ! Add specific set of attributes to components from driver attributes

    use ESMF  , only : ESMF_GridComp, ESMF_Config, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogFoundAllocError, ESMF_ConfigGetLen, ESMF_ConfigGetAttribute
    use NUOPC , only : NUOPC_CompAttributeAdd, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet

    ! input/output parameters
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_GridComp) , intent(in)    :: driver
    type(ESMF_Config)   , intent(inout) :: config
    integer             , intent(in)    :: compid
    character(len=*)    , intent(in)    :: compname
    character(len=*)    , intent(in)    :: inst_suffix
    integer             , intent(in)    :: nthrds
    integer             , intent(inout) :: rc
    ! local variables
    integer                        :: inst_index
    logical                        :: computetask
    character(len=CL)              :: cvalue
    character(len=CS)              :: attribute
    character(len=*), parameter    :: subname = "(esm.F90:AddAttributes)"
    !-------------------------------------------
    computetask = .false.
    rc = ESMF_Success
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    call shr_log_setLogunit(logunit)

    !------
    ! Add compid to gcomp attributes
    !------
    write(cvalue,*) compid
    call NUOPC_CompAttributeAdd(gcomp, attrList=(/'MCTID'/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name='MCTID', value=trim(cvalue), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !------
    ! Add driver restart flag to gcomp attributes
    !------
    attribute = 'read_restart'
    call NUOPC_CompAttributeGet(driver, name=trim(attribute), isPresent=computetask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if(.not. computetask) return

    call NUOPC_CompAttributeGet(driver, name=trim(attribute), value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeAdd(gcomp, (/trim(attribute)/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name=trim(attribute), value=trim(cvalue), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !------
    ! Add component specific attributes
    !------
    call ReadAttributes(gcomp, config, trim(compname)//"_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ReadAttributes(gcomp, config, "ALLCOMP_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": call Readattributes for"//trim(compname), ESMF_LOGMSG_INFO)

    call ReadAttributes(gcomp, config, trim(compname)//"_modelio::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) then
       print *,__FILE__,__LINE__,"ERROR reading ",trim(compname)," modelio from runconfig"
       return
    endif
    call ReadAttributes(gcomp, config, "CLOCK_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !------
    ! Add mediator specific attributes - if component is mediator
    !------
    if (compname == 'MED') then
       call ReadAttributes(gcomp, config, "MED_attributes::", rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    !------
    ! Add multi-instance specific attributes
    !------
    call NUOPC_CompAttributeAdd(gcomp, attrList=(/'inst_index'/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! add inst_index attribute (inst_index is not required for cime internal components)
    ! for now hard-wire inst_index to 1
    inst_index = 1
    write(cvalue,*) inst_index
    call NUOPC_CompAttributeSet(gcomp, name='inst_index', value=trim(cvalue), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! add inst_suffix attribute
    if (len_trim(inst_suffix) > 0) then
       call NUOPC_CompAttributeAdd(gcomp, attrList=(/'inst_suffix'/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name='inst_suffix', value=inst_suffix, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    ! Add the nthreads attribute
    call NUOPC_CompAttributeAdd(gcomp, attrList=(/'nthreads'/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(cvalue, *) nthrds
    call NUOPC_CompAttributeSet(gcomp, name='nthreads', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !------
    ! Add single column and single point attributes
    !------
    call esm_set_single_column_attributes(compname, gcomp, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine AddAttributes

  !================================================================================

  subroutine ReadAttributes(gcomp, config, label, relaxedflag, formatprint, rc)

    use ESMF  , only : ESMF_GridComp, ESMF_Config, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use NUOPC , only : NUOPC_FreeFormatCreate, NUOPC_CompAttributeIngest
    use NUOPC , only : NUOPC_FreeFormatDestroy, NUOPC_FreeFormat

    ! input/output arguments
    type(ESMF_GridComp) , intent(inout)        :: gcomp
    type(ESMF_Config)   , intent(in)           :: config
    character(len=*)    , intent(in)           :: label
    logical             , intent(in), optional :: relaxedflag
    logical             , intent(in), optional :: formatprint
    integer             , intent(inout)        :: rc

    ! local variables
    type(NUOPC_FreeFormat)  :: attrFF
    character(len=*), parameter :: subname = "(esm.F90:ReadAttributes)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    call shr_log_setLogunit(logunit)

    if (present(relaxedflag)) then
       attrFF = NUOPC_FreeFormatCreate(config, label=trim(label), relaxedflag=.true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       attrFF = NUOPC_FreeFormatCreate(config, label=trim(label), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    call NUOPC_CompAttributeIngest(gcomp, attrFF, addFlag=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
#ifdef DEBUG
!    if (present (formatprint)) then
!       call pretty_print_nuopc_freeformat(attrFF, trim(label)//' attributes', rc=rc)
!       if (chkerr(rc,__LINE__,u_FILE_u)) return
!    end if
#endif
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine ReadAttributes

  !================================================================================

  subroutine InitAdvertize(driver, importState, exportState, clock, rc)

    ! This empty InitAdvertise is needed because it overrides the behavior
    ! of the default InitAdvertise inside the generic NUOPC_Driver.F90. The
    ! default behavior tries to mirror the fields up the hierarchy (i.e., up
    ! to the ensemble driver). This would be used if we needed to
    ! communicate between the ensemble members. Since we do not need that
    ! right now, we turn it off with this empty subroutine.

    use ESMF, only : ESMF_GridComp, ESMF_State, ESMF_Clock
    use ESMF, only : ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO

    ! input/output variables
    type(ESMF_GridComp)  :: driver
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname = "(esm.F90:InitAdvertize)"
    !---------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

  end subroutine InitAdvertize

  !================================================================================

  subroutine esm_init_pelayout(driver, maxthreads, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_VM, ESMF_VMGet
    use ESMF         , only : ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_Config
    use ESMF         , only : ESMF_ConfigGetLen, ESMF_LogFoundAllocError, ESMF_ConfigGetAttribute
    use ESMF         , only : ESMF_RC_NOT_VALID, ESMF_LogSetError, ESMF_Info, ESMF_InfoSet
    use ESMF         , only : ESMF_GridCompIsPetLocal, ESMF_MethodAdd, ESMF_UtilStringLowerCase
    use ESMF         , only : ESMF_InfoCreate, ESMF_InfoDestroy
    use NUOPC        , only : NUOPC_CompAttributeGet
    use NUOPC_Driver , only : NUOPC_DriverAddComp
#ifndef NO_MPI2
    use mpi          , only : MPI_COMM_NULL, mpi_comm_size
#endif

#ifdef MED_PRESENT
    use med_internalstate_mod , only : med_id
    use med                   , only : MedSetServices => SetServices
#ifdef ESMF_AWARE_THREADING
    use med                   , only : MEDSetVM => SetVM
#endif
#endif
#ifdef ATM_PRESENT
    use atm_comp_nuopc        , only : ATMSetServices => SetServices
#ifdef ESMF_AWARE_THREADING
    use atm_comp_nuopc        , only : ATMSetVM => SetVM
#endif
#endif
#ifdef ICE_PRESENT
    use ice_comp_nuopc        , only : ICESetServices => SetServices
#ifdef ESMF_AWARE_THREADING
    use ice_comp_nuopc        , only : ICESetVM => SetVM
#endif
#endif
#ifdef LND_PRESENT
    use lnd_comp_nuopc        , only : LNDSetServices => SetServices
#ifdef ESMF_AWARE_THREADING
    use lnd_comp_nuopc        , only : LNDSetVM => SetVM
#endif
#endif
#ifdef OCN_PRESENT
    use ocn_comp_nuopc        , only : OCNSetServices => SetServices
#ifdef ESMF_AWARE_THREADING
    use ocn_comp_nuopc        , only : OCNSetVM => SetVM
#endif
#endif
#ifdef WAV_PRESENT
    use wav_comp_nuopc        , only : WAVSetServices => SetServices
#ifdef ESMF_AWARE_THREADING
    use wav_comp_nuopc        , only : WAVSetVM => SetVM
#endif
#endif
#ifdef ROF_PRESENT
    use rof_comp_nuopc        , only : ROFSetServices => SetServices
#ifdef ESMF_AWARE_THREADING
    use rof_comp_nuopc        , only : ROFSetVM => SetVM
#endif
#endif
#ifdef GLC_PRESENT
    use glc_comp_nuopc        , only : GLCSetServices => SetServices
#ifdef ESMF_AWARE_THREADING
    use glc_comp_nuopc        , only : GLCSetVM => SetVM
#endif
#endif
#ifdef NO_MPI2
    include 'mpif.h'
#endif
    ! input/output variables
    type(ESMF_GridComp)            :: driver
    integer, intent(out)           :: maxthreads ! maximum number of threads any component
    integer, intent(out)           :: rc

    ! local variables
    type(ESMF_GridComp)            :: child
    type(ESMF_VM)                  :: vm
    type(ESMF_Config)              :: config
    type(ESMF_Info)                :: info
    integer                        :: PetCount
    integer                        :: ComponentCount
    integer                        :: ntasks, rootpe, nthrds, stride
    integer                        :: ntask
    integer                        :: i
    integer                        :: stat
    character(len=32), allocatable :: compLabels(:)
    character(CS)                  :: namestr
    character(CL)                  :: msgstr
    integer, allocatable           :: petlist(:)
    integer, pointer               :: comms(:), comps(:)
    integer                        :: Global_Comm
    logical                        :: isPresent
    integer, allocatable           :: comp_comm_iam(:)
    logical, allocatable           :: comp_iamin(:)
    character(len=5)               :: inst_suffix
    character(CL)                  :: cvalue
    logical                        :: found_comp
    integer :: rank, nprocs, ierr
    character(len=*), parameter    :: subname = "(esm_pelayout.F90:esm_init_pelayout)"
    !---------------------------------------
    call shr_log_setLogunit(logunit)

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(driver, vm=vm, config=config, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "PELAYOUT_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, petCount=petCount, mpiCommunicator=Global_Comm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    componentCount = ESMF_ConfigGetLen(config,label="component_list:", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(compLabels(componentCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Allocation of compLabels failed.", &
         line=__LINE__, file=u_FILE_u, rcToReturn=rc)) return
    allocate(comp_iamin(componentCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Allocation of compLabels failed.", &
         line=__LINE__, file=u_FILE_u, rcToReturn=rc)) return
    allocate(comp_comm_iam(componentCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Allocation of compLabels failed.", &
         line=__LINE__, file=u_FILE_u, rcToReturn=rc)) return

    call ESMF_ConfigGetAttribute(config, valueList=compLabels, label="component_list:", count=componentCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(driver, name="inst_suffix", value=inst_suffix, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ""
    endif

    allocate(comms(componentCount+1), comps(componentCount+1))
    comps(1) = 1
    comms = MPI_COMM_NULL
    comms(1) = Global_Comm

    maxthreads = 1
    do i=1,componentCount
       namestr = ESMF_UtilStringLowerCase(compLabels(i))
       if (namestr == 'med') namestr = 'cpl'
       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_nthreads', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) nthrds

       if(nthrds > maxthreads) maxthreads = nthrds
    enddo

    do i=1,componentCount
       namestr = ESMF_UtilStringLowerCase(compLabels(i))
       if (namestr == 'med') namestr = 'cpl'
       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_ntasks', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) ntasks

       if (ntasks < 0 .or. ntasks > PetCount) then
          write (msgstr, *) "Invalid NTASKS value specified for component: ",namestr, ' ntasks: ',ntasks, petcount
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif

       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_nthreads', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) nthrds

       info = ESMF_InfoCreate(rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_InfoSet(info, key="/NUOPC/Hint/PePerPet/MaxCount", value=nthrds, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_rootpe', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) rootpe
       if (rootpe < 0 .or. rootpe > PetCount) then
          write (msgstr, *) "Invalid Rootpe value specified for component: ",namestr, ' rootpe: ',rootpe
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif
       if(rootpe+ntasks > PetCount) then
          write (msgstr, *) "Invalid pelayout value specified for component: ",namestr, ' rootpe+ntasks: ',rootpe+ntasks
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif

       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_pestride', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stride
       if (stride < 1 .or. rootpe+(ntasks-1)*stride > PetCount) then
          write (msgstr, *) "Invalid pestride value specified for component: ",namestr,&
               ' rootpe: ',rootpe, ' pestride: ', stride, ' ntasks: ',ntasks, ' PetCount: ', PetCount
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif

       if (allocated(petlist)) then
#ifdef ESMF_AWARE_THREADING
          if(size(petlist) .ne. ntasks*nthrds) then
#else
          if(size(petlist) .ne. ntasks) then
#endif
             deallocate(petlist)
          endif
       endif
       if(.not. allocated(petlist)) then
#ifdef ESMF_AWARE_THREADING
          allocate(petlist(ntasks*nthrds))
#else
          allocate(petlist(ntasks))
#endif
       endif

#ifdef ESMF_AWARE_THREADING
       cnt = 1
       do ntask = rootpe, rootpe+nthrds*ntasks*stride-1, stride
          petlist(cnt) = ntask
          cnt = cnt + 1
       enddo
#else
       do ntask = 1, size(petlist)
          petlist(ntask) = rootpe + (ntask-1)*stride
       enddo
#endif

       comps(i+1) = i+1
       found_comp = .false.
#ifdef MED_PRESENT
       if (trim(compLabels(i)) == 'MED') then
          med_id = i + 1
#ifdef ESMF_AWARE_THREADING
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), MEDSetServices, MEDSetVM, &
               petList=petlist, comp=child, info=info, rc=rc)
#else
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), MEDSetServices,  &
               petList=petlist, comp=child, rc=rc)
#endif
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
#ifdef ATM_PRESENT
       if (trim(compLabels(i)) .eq. 'ATM') then
#ifdef ESMF_AWARE_THREADING
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ATMSetServices, ATMSetVM, &
               petList=petlist, comp=child, info=info, rc=rc)
#else
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ATMSetServices,  &
               petList=petlist, comp=child, rc=rc)
#endif
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
#ifdef LND_PRESENT
       if (trim(compLabels(i)) .eq. 'LND') then
#ifdef ESMF_AWARE_THREADING
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), LNDSetServices, LNDSetVM, &
               PetList=petlist, comp=child, info=info, rc=rc)
#else
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), LNDSetServices, &
               PetList=petlist, comp=child, rc=rc)
#endif
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
#ifdef OCN_PRESENT
       if (trim(compLabels(i)) .eq. 'OCN') then
#ifdef ESMF_AWARE_THREADING
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), OCNSetServices, OCNSetVM, &
               PetList=petlist, comp=child, info=info, rc=rc)
#else
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), OCNSetServices, &
               PetList=petlist, comp=child, rc=rc)
#endif
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
#ifdef ICE_PRESENT
       if (trim(compLabels(i)) .eq. 'ICE') then
#ifdef ESMF_AWARE_THREADING
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ICESetServices, ICESetVM, &
               PetList=petlist, comp=child, info=info, rc=rc)
#else
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ICESetServices, &
               PetList=petlist, comp=child, rc=rc)
#endif
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
#ifdef GLC_PRESENT
       if (trim(compLabels(i)) .eq. 'GLC') then
#ifdef ESMF_AWARE_THREADING
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), GLCSetServices, GLCSetVM, &
               PetList=petlist, comp=child, info=info, rc=rc)
#else
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), GLCSetServices, &
               PetList=petlist, comp=child, rc=rc)
#endif
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
#ifdef ROF_PRESENT
       if (trim(compLabels(i)) .eq. 'ROF') then
#ifdef ESMF_AWARE_THREADING
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ROFSetServices, ROFSetVM, &
               PetList=petlist, comp=child, info=info, rc=rc)
#else
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ROFSetServices,  &
               PetList=petlist, comp=child, rc=rc)
#endif
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
#ifdef WAV_PRESENT
       if (trim(compLabels(i)) .eq. 'WAV') then
#ifdef ESMF_AWARE_THREADING
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), WAVSetServices, WAVSetVM, &
               PetList=petlist, comp=child, info=info, rc=rc)
#else
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), WAVSetServices,  &
               PetList=petlist, comp=child, rc=rc)
#endif
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
#ifdef ESP_PRESENT
       if (trim(compLabels(i)) .eq. 'ESP') then
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ESPSetServices, PetList=petlist, comp=child, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          found_comp = .true.
       end if
#endif
       if (.not. found_comp) then
          write(msgstr,*) 'No component ',trim(compLabels(i)),' found'
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif
       comp_iamin(i) = .false.

       call AddAttributes(child, driver, config, i+1, trim(compLabels(i)), inst_suffix, nthrds, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_GridCompIsPetLocal(child, rc=rc)) then

          call ESMF_GridCompGet(child, vm=vm, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          call ESMF_VMGet(vm, mpiCommunicator=comms(i+1), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (comms(i+1) .ne. MPI_COMM_NULL) then
             call ESMF_VMGet(vm, localPet=comp_comm_iam(i), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             comp_iamin(i) = .true.
             call MPI_Comm_size(comms(i+1), nprocs, ierr)
             call MPI_Comm_rank(comms(i+1), rank, ierr)
             if(nprocs /= ntasks) then
                write(msgstr,*) 'Component ',trim(compLabels(i)),' has mpi task mismatch, do threads align with nodes?'
                call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
                return
             endif
          endif
       endif

       call ESMF_InfoDestroy(info, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    enddo

    deallocate(petlist, comms, comps, comp_iamin, comp_comm_iam)

  end subroutine esm_init_pelayout

  !================================================================================

  subroutine esm_set_single_column_attributes(compname, gcomp, rc)

    ! Generate a mesh for single column

    use netcdf, only : nf90_open, nf90_close, nf90_noerr, nf90_nowrite
    use netcdf, only : nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var
    use NUOPC , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_VM, ESMF_VMGet, ESMF_SUCCESS
    use ESMF  , only : ESMF_Mesh, ESMF_MeshCreate, ESMF_FILEFORMAT_ESMFMESH, ESMF_MeshGet, ESMF_MESHLOC_ELEMENT
    use ESMF  , only : ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet, ESMF_FieldRegridGetArea, ESMF_TYPEKIND_r8

    ! input/output variables
    character(len=*)    , intent(in)    :: compname
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    character(len=CL)      :: single_column_lnd_domainfile
    character(len=CL)      :: single_column_global_meshfile
    real(r8)               :: scol_lon
    real(r8)               :: scol_lat
    real(r8)               :: scol_area
    integer                :: scol_lndmask
    real(r8)               :: scol_lndfrac
    integer                :: scol_ocnmask
    real(r8)               :: scol_ocnfrac
    integer                :: scol_mesh_n
    type(ESMF_Mesh)        :: mesh
    type(ESMF_Field)       :: lfield
    integer                :: lsize
    integer                :: spatialdim
    real(r8), pointer      :: ownedElemCoords(:)
    real(r8), pointer      :: latMesh(:)
    real(r8), pointer      :: lonMesh(:)
    real(r8), pointer      :: dataptr(:)
    integer                :: i,j,ni,nj,n
    integer                :: ncid
    integer                :: dimid
    integer                :: varid_xc
    integer                :: varid_yc
    integer                :: varid_area
    integer                :: varid_mask
    integer                :: varid_frac
    integer                :: start(2)       ! Start index to read in
    integer                :: start3(3)      ! Start index to read in
    integer                :: count3(3)      ! Number of points to read in
    integer                :: status         ! status flag
    real (r8), allocatable :: lats(:)        ! temporary
    real (r8), allocatable :: lons(:)        ! temporary
    real (r8), allocatable :: pos_lons(:)    ! temporary
    real (r8), allocatable :: pos_lats(:)    ! temporary
    real (r8), allocatable :: cols(:)        ! temporary
    real (r8), allocatable :: glob_grid(:,:) ! temporary
    real (r8)              :: pos_scol_lon   ! temporary
    real (r8)              :: pos_scol_lat   ! temporary
    real (r8)              :: scol_data(1)
    integer                :: iscol_data(1)
    integer                :: petcount
    character(len=CL)      :: cvalue
    character(len=*), parameter :: subname= ' (esm_get_single_column_attributes) '
    logical                :: unstructured = .false.
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call shr_log_setLogunit(logunit)
    scol_mesh_n = 0

    ! obtain the single column lon and lat
    call NUOPC_CompAttributeGet(gcomp, name='scol_lon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scol_lon
    call NUOPC_CompAttributeGet(gcomp, name='scol_lat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scol_lat
    call NUOPC_CompAttributeGet(gcomp, name='single_column_lnd_domainfile', value=single_column_lnd_domainfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(gcomp, name='mesh_atm', value=single_column_global_meshfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeAdd(gcomp, attrList=(/'scol_spval'/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if ( (scol_lon < scol_spval .and. scol_lat > scol_spval) .or. &
         (scol_lon > scol_spval .and. scol_lat < scol_spval)) then
       call shr_sys_abort(subname//' ERROR: '//trim(compname)//' both scol_lon and scol_lat must be greater than -999 ')
    end if

    ! Set the special value for single column - if pts_lat or pts_lon are equal to the special value
    ! in the component cap - then single column is not activated
    write(cvalue,*) scol_spval
    call NUOPC_CompAttributeSet(gcomp, name='scol_spval', value=trim(cvalue), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (scol_lon > scol_spval .and.  scol_lat > scol_spval) then

       ! NOTE: currently assume that single column capability is restricted to
       ! ATM, LND, OCN and ICE components only
       ! verify that WAV and LND are not trying to use single column mode
       if (trim(compname) == 'WAV' .or. trim(compname) == 'ROF' .or. trim(compname) == 'GLC') then
          call shr_sys_abort(subname//' ERROR: '//trim(compname)//' does not support single column mode ')
       end if

       ! ensure that single column mode is only run on 1 pet
       call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMGet(vm, petcount=petcount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (petcount > 1) then
          call shr_sys_abort(subname//' ERROR: single column mode must be run on 1 pe')
       endif

       write(logunit,'(a,2(f10.5,2x))')trim(subname)//' single column point for '//trim(compname)//&
            ' has lon and lat = ',scol_lon,scol_lat

       ! This is either a single column or a single point so add attributes
       call NUOPC_CompAttributeAdd(gcomp, &
            attrList=(/'scol_area   ', &
                       'scol_ni     ', &
                       'scol_nj     ', &
                       'scol_lndmask', &
                       'scol_lndfrac', &
                       'scol_ocnmask', &
                       'scol_ocnfrac'/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (trim(single_column_lnd_domainfile) /= 'UNSET') then

          ! In this case the domain file is not a single point file - but normally a
          ! global domain file where a nearest neighbor search will be done to find
          ! the closest point in the domin file to scol_lon and scol_lat

          status = nf90_open(single_column_lnd_domainfile, NF90_NOWRITE, ncid)
          if (status /= nf90_noerr) call shr_sys_abort (trim(subname) //': opening '//&
               trim(single_column_lnd_domainfile))
          status = nf90_inq_dimid (ncid, 'ni', dimid)
          if (status /= nf90_noerr) call shr_sys_abort (trim(subname) //': inq_dimid ni')
          status = nf90_inquire_dimension(ncid, dimid, len=ni)
          if (status /= nf90_noerr) call shr_sys_abort (trim(subname) //': inquire_dimension ni')
          status = nf90_inq_dimid (ncid, 'nj', dimid)
          if (status /= nf90_noerr) call shr_sys_abort (trim(subname) //': inq_dimid nj')
          status = nf90_inquire_dimension(ncid, dimid, len=nj)
          if (status /= nf90_noerr) call shr_sys_abort (trim(subname) //': inquire_dimension nj')

          status = nf90_inq_varid(ncid, 'xc' , varid_xc)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' inq_varid xc')
          status = nf90_inq_varid(ncid, 'yc' , varid_yc)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' inq_varid yc')
          status = nf90_inq_varid(ncid, 'area' , varid_area)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' inq_varid area')
          status = nf90_inq_varid(ncid, 'mask' , varid_mask)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' inq_varid mask')
          status = nf90_inq_varid(ncid, 'frac' , varid_frac)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' inq_varid frac')

          ! Read in domain file for single column
          ! Check for unstructured data ni>1 and nj==1
          if (ni.gt.1 .and. nj == 1) unstructured=.true.

          if (unstructured) then
             allocate(lats(ni))
             allocate(pos_lats(ni))
          else
             allocate(lats(nj))
          end if
          allocate(lons(ni))
          allocate(pos_lons(ni))
          allocate(glob_grid(ni,nj))

          ! The follow assumes that xc and yc are 2 dimensional values
          start3=(/1,1,1/)
          count3=(/ni,nj,1/)
          status = nf90_get_var(ncid, varid_xc, glob_grid, start3, count3)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' get_var xc')
          lons(1:ni) = glob_grid(1:ni,1)
          status = nf90_get_var(ncid, varid_yc, glob_grid, start3, count3)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' get_var yc')
          if (unstructured) then
             lats(1:ni) = glob_grid(1:ni,1)
          else
             lats(1:nj) = glob_grid(1,1:nj)
          end if
          ! find nearest neighbor indices of scol_lon and scol_lat in single_column_lnd_domain file
          ! convert lons array and scol_lon to 0,360 and find index of value closest to 0
          ! and obtain single-column longitude/latitude indices to retrieve
          if (unstructured) then
             allocate(cols(ni))
             pos_lons(:)  = mod(lons(:)  + 360._r8, 360._r8)
             pos_scol_lon = mod(scol_lon + 360._r8, 360._r8)
             pos_lats(:)  = lats(:)  + 90._r8
             pos_scol_lat = scol_lat + 90._r8
             cols=abs(pos_lons - pos_scol_lon)+abs(pos_lats - pos_scol_lat)
             start(1) = MINLOC(cols, dim=1)
             start(2) = 1
             deallocate(cols)
          else
             pos_lons(:)  = mod(lons(:)  + 360._r8, 360._r8)
             pos_scol_lon = mod(scol_lon + 360._r8, 360._r8)
             start(1) = (MINLOC(abs(pos_lons - pos_scol_lon), dim=1))
             start(2) = (MINLOC(abs(lats      -scol_lat    ), dim=1))
          end if
          deallocate(lats)
          deallocate(lons)
          deallocate(pos_lons)
          deallocate(glob_grid)
          ! read in value of nearest neighbor lon and RESET scol_lon and scol_lat
          ! also get area of gridcell, mask and frac
          status = nf90_get_var(ncid, varid_xc, scol_lon, start)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' get_var xc')

          status = nf90_get_var(ncid, varid_yc, scol_lat, start)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' get_var yc')

          status = nf90_get_var(ncid, varid_area, scol_area, start)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' get_var area')

          status = nf90_get_var(ncid, varid_mask, iscol_data, start)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' get_var mask')
          scol_lndmask = iscol_data(1)
          scol_ocnmask = 1 - scol_lndmask

          status = nf90_get_var(ncid, varid_frac, scol_data, start)
          if (status /= nf90_noerr) call shr_sys_abort (subname//' get_var frac')
          scol_lndfrac = scol_data(1)
          scol_ocnfrac = 1._r8 - scol_lndfrac

          if (scol_ocnmask == 0 .and. scol_lndmask == 0) then
             call shr_sys_abort(trim(subname)//' in single column mode '&
                  //' ocean and land mask cannot both be zero')
          end if

          status = nf90_close(ncid)
          if (status /= nf90_noerr) call shr_sys_abort (trim(subname) //': closing '//&
               trim(single_column_lnd_domainfile))

          ! Now read in mesh file to get exact values of scol_lon and scol_lat that will be used
          ! by the models - assume that this occurs only on 1 processor
          mesh = ESMF_MeshCreate(filename=trim(single_column_global_meshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=lsize, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(ownedElemCoords(spatialDim*lsize))
          allocate(lonMesh(lsize), latMesh(lsize))
          call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          scol_mesh_n = 0
          do n = 1,lsize
             lonMesh(n) = ownedElemCoords(2*n-1)
             latMesh(n) = ownedElemCoords(2*n)
             if (abs(lonMesh(n) - scol_lon) < 1.e-4 .and. abs(latMesh(n) - scol_lat) < 1.e-4) then
                scol_mesh_n = n
                exit
             end if
          end do
          scol_lon  = lonMesh(scol_mesh_n)
          scol_lat  = latMesh(scol_mesh_n)

          ! Obtain mesh info areas
          lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_r8, name='area', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldRegridGetArea(lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          scol_area = dataptr(scol_mesh_n)

          ! Set single column attribute values for all components
          write(cvalue,*) scol_lon
          call NUOPC_CompAttributeSet(gcomp, name='scol_lon', value=trim(cvalue), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          write(cvalue,*) scol_lat
          call NUOPC_CompAttributeSet(gcomp, name='scol_lat', value=trim(cvalue), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          write(cvalue,*) scol_area
          call NUOPC_CompAttributeSet(gcomp, name='scol_area', value=trim(cvalue), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Write out diagnostic info
          write(logunit,'(a,2(f13.5,2x))')trim(subname)//' nearest neighbor scol_lon and scol_lat in '&
               //trim(single_column_lnd_domainfile)//' are ',scol_lon,scol_lat
          if (trim(compname) == 'LND') then
             write(logunit,'(a,i4,f13.5)')trim(subname)//' scol_lndmask, scol_lndfrac are ',&
                  scol_lndmask, scol_lndfrac
          else if (trim(compname) == 'OCN') then
             write(logunit,'(a,i4,f13.5)')trim(subname)//' scol_ocnmask, scol_ocnfrac are ',&
                  scol_ocnmask, scol_ocnfrac
          else
             write(logunit,'(a)')trim(subname)//' atm point has unit mask and unit fraction '
          end if
          write(cvalue,*) ni
          call NUOPC_CompAttributeSet(gcomp, name='scol_ni', value=trim(cvalue), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          write(cvalue,*) nj
          call NUOPC_CompAttributeSet(gcomp, name='scol_nj', value=trim(cvalue), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

       else

          ! for single point mode - its assumed that this point is always
          ! either over all ocean or all land and the point is always valid

          scol_lndmask = 1
          scol_lndfrac = 1._r8
          scol_ocnmask = 1
          scol_ocnfrac = 1._r8
          scol_area = 1.e30

          write(cvalue,*) 1
          call NUOPC_CompAttributeSet(gcomp, name='scol_ni', value=trim(cvalue), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeSet(gcomp, name='scol_nj', value=trim(cvalue), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          write(logunit,'(a)')' single point mode is active'
          write(logunit,'(a,f13.5,a,f13.5,a)')' scol_lon is ',scol_lon,' and scol_lat is '

       end if

       write(cvalue,*) scol_area
       call NUOPC_CompAttributeSet(gcomp, name='scol_area', value=trim(cvalue), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write(cvalue,*) scol_lndmask
       call NUOPC_CompAttributeSet(gcomp, name='scol_lndmask', value=trim(cvalue), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write(cvalue,*) scol_lndfrac
       call NUOPC_CompAttributeSet(gcomp, name='scol_lndfrac', value=trim(cvalue), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write(cvalue,*) scol_ocnmask
       call NUOPC_CompAttributeSet(gcomp, name='scol_ocnmask', value=trim(cvalue), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write(cvalue,*) scol_ocnfrac
       call NUOPC_CompAttributeSet(gcomp, name='scol_ocnfrac', value=trim(cvalue), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    end if

  end subroutine esm_set_single_column_attributes

  !================================================================================
  subroutine esm_finalize(driver, rc)

    use ESMF     , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_VM, ESMF_VMGet
    use ESMF     , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_LOGWRITE
    use NUOPC    , only : NUOPC_CompAttributeGet
    use perf_mod , only : t_prf, t_finalizef

    ! input/output variables
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    character(len=CL)    :: timing_dir        ! timing directory
    character(len=5)     :: inst_suffix
    logical              :: isPresent
    type(ESMF_VM)        :: vm
    integer              :: mpicomm
    character(len=*), parameter :: subname = '(esm_finalize) '
    !---------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    call shr_log_setLogunit(logunit)
    call ESMF_GridCompGet(driver, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, mpiCommunicator=mpicomm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="timing_dir",value=timing_dir, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(driver, name="inst_suffix", value=inst_suffix, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ""
    endif
    call t_prf(trim(timing_dir)//'/model_timing'//trim(inst_suffix), mpicom=mpicomm)

    if (maintask) then
       write(logunit,*)' SUCCESSFUL TERMINATION OF CESM'
    end if
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    call t_finalizef()
  end subroutine esm_finalize


end module ESM
