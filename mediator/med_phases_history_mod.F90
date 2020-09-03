module med_phases_history_mod

  !-----------------------------------------------------------------------------
  ! Mediator History control
  !
  ! Each time loop has its own associated clock object. NUOPC manages
  ! these clock objects, i.e. their creation and destruction, as well as
  ! startTime, endTime, timeStep adjustments during the execution. The
  ! outer most time loop of the run sequence is a special case. It uses
  ! the driver clock itself. If a single outer most loop is defined in
  ! the run sequence provided by freeFormat, this loop becomes the driver
  ! loop level directly. Therefore, setting the timeStep or runDuration
  ! for the outer most time loop results modifiying the driver clock
  ! itself. However, for cases with concatenated loops on the upper level
  ! of the run sequence in freeFormat, a single outer loop is added
  ! automatically during ingestion, and the driver clock is used for this
  ! loop instead.
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_VM, ESMF_VMGet
  use ESMF                  , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockGetNextTime, ESMF_ClockGetAlarm
  use ESMF                  , only : ESMF_Calendar
  use ESMF                  , only : ESMF_Time, ESMF_TimeGet
  use ESMF                  , only : ESMF_TimeInterval, ESMF_TimeIntervalGet
  use ESMF                  , only : ESMF_Alarm, ESMF_AlarmIsRinging, ESMF_AlarmRingerOff, ESMF_AlarmGet
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_MAXSTR, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
  use ESMF                  , only : ESMF_Finalize
  use ESMF                  , only : operator(==), operator(-), operator(/=)
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use NUOPC_Model           , only : NUOPC_ModelGet
  use esmFlds               , only : compatm, complnd, compocn, compice, comprof, compglc, ncomps, compname, ncomps
  use esmFlds               , only : fldListFr, fldListTo
  use med_constants_mod     , only : dbug_flag       => med_constants_dbug_flag
  use med_constants_mod     , only : SecPerDay       => med_constants_SecPerDay
  use med_utils_mod         , only : chkerr          => med_utils_ChkErr
  use med_methods_mod       , only : FB_reset        => med_methods_FB_reset
  use med_methods_mod       , only : FB_diagnose     => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_GetFldPtr    => med_methods_FB_GetFldPtr
  use med_methods_mod       , only : FB_accum        => med_methods_FB_accum
  use med_methods_mod       , only : State_GetScalar => med_methods_State_GetScalar
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_time_mod          , only : med_time_alarmInit
  use med_io_mod            , only : med_io_write, med_io_wopen, med_io_enddef
  use med_io_mod            , only : med_io_close, med_io_date2yyyymmdd, med_io_sec2hms
  use med_io_mod            , only : med_io_ymd2date
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public :: med_phases_history_alarms_init
  public :: med_phases_history_write
  public :: med_phases_history_write_atm
  public :: med_phases_history_write_ice
  public :: med_phases_history_write_glc
  public :: med_phases_history_write_lnd
  public :: med_phases_history_write_ocn
  public :: med_phases_history_write_rof
  public :: med_phases_history_write_wav
  public :: med_phases_history_write_aux

  private :: med_phases_history_setauxflds

  ! type(ESMF_FieldBundle) :: FBImpAvg(ncomps)  ! TODO: fill this in
  ! type(ESMF_FieldBundle) :: FBExpAvg(ncomps)  ! TODO: fill this in

  character(*), parameter :: u_FILE_u  = &
       __FILE__

  integer, parameter :: maxfiles = 20
  integer :: nfiles = 0
  type, public :: auxfile_type
     integer                    :: ncomp      ! component index
     character(CS)              :: auxname    ! name for history file creation
     character(CS), allocatable :: flds(:)    ! array of aux field names
     character(CS)              :: alarmname  ! alarm name for output
     integer                    :: deltat     ! interval to write out aux data in seconds
     integer                    :: ntperfile  ! maximum number of time samples per file
     integer                    :: nt = 0     ! time in file
     real(r8)                   :: tbnds(2)   ! CF1.0 time bounds
     character(CS)              :: avgtype    ! instantaneous or time average
     type(ESMF_FieldBundle)     :: FBavg      ! field bundle for time averaging
  end type auxfile_type
  type(auxfile_type) :: auxfiles(maxfiles)

  character(CL) :: case_name  ! case name
  character(CS) :: inst_tag   ! instance tag

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_history_alarms_init(gcomp, rc)

    ! --------------------------------------
    ! Initialize mediator history file alarms
    ! --------------------------------------

    use ESMF        , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF        , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockAdvance, ESMF_ClockSet
    use ESMF        , only : ESMF_Time, ESMF_TimeInterval, ESMF_TimeIntervalGet
    use ESMF        , only : ESMF_Alarm, ESMF_AlarmSet
    use ESMF        , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF        , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF        , only : operator(==), operator(-)
    use NUOPC       , only : NUOPC_CompAttributeGet
    use NUOPC_Model , only : NUOPC_ModelGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(ESMF_VM)           :: vm
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Clock)        :: mclock
    type(ESMF_TimeInterval) :: mtimestep
    type(ESMF_Time)         :: mCurrTime
    type(ESMF_Time)         :: mStartTime
    type(ESMF_TimeInterval) :: timestep
    integer                 :: timestep_length
    character(CS)           :: alarmname       ! alarm name
    character(CL)           :: cvalue          ! attribute string
    character(CL)           :: hist_option     ! freq_option setting (ndays, nsteps, etc)
    integer                 :: hist_n          ! freq_n setting relative to freq_option
    integer                 :: n
    logical                 :: isPresent
    logical                 :: isSet
    character(*),parameter  :: F01 = "(a,2x,i8)"
    character(len=*), parameter :: subname='(med_phases_history_alarms_init)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get model clock, start time, current time and time step
    call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet(mclock, startTime=mStartTime,  currTime=mCurrTime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(mtimestep, s=timestep_length, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(logunit,*)
    write(logunit,F01) trim(subname)//" history clock timestep = ",timestep_length

    ! Determine instantaneous mediator output frequency and type
    call NUOPC_CompAttributeGet(gcomp, name='history_option', isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       call NUOPC_CompAttributeGet(gcomp, name='history_option', value=hist_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name='history_n', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) hist_n
    else
       hist_option = 'none'
       hist_n = -999
    end if

    ! Set alarms for instantaneous mediator history output
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    alarmname = 'alarm_history_inst_all'
    call med_time_alarmInit(mclock, alarm, option=hist_option, opt_n=hist_n, &
         reftime=mStartTime, alarmname=trim(alarmname), rc=rc)
    call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockSet(mclock, currTime=mcurrtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       write(logunit,F01) trim(subname)//" set instantaneous mediator history alarm "//&
            trim(alarmname)//"  with option "//trim(hist_option)//" and frequency ",hist_n
    end if
    do n = 1,ncomps
       if (is_local%wrap%comp_present(n)) then
          alarmname = 'alarm_history_inst_' // trim(compname(n))
          call med_time_alarmInit(mclock, alarm, option=hist_option, opt_n=hist_n, &
               reftime=mStartTime, alarmname=trim(alarmname), rc=rc)
          call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockAdvance(mclock,rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockSet(mclock, currTime=mcurrtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (mastertask) then
             write(logunit,F01) trim(subname)//" set instantaneous mediator history alarm "//&
                  trim(alarmname)//"  with option "//trim(hist_option)//" and frequency ",hist_n
          end if
       end if
    end do

    ! Initialize field bundles for doing time averaged mediator history output
    if (hist_option /= 'none') then
       ! TODO: fill this in
    end if

    ! Determine time average mediator output frequency and type
    call NUOPC_CompAttributeGet(gcomp, name='histavg_option', isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       call NUOPC_CompAttributeGet(gcomp, name='hist_option', value=hist_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name='history_n', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) hist_n
    else
       hist_option = 'none'
       hist_n = -999
    end if

    ! Set alarm for time averaged mediator history output
    alarmname = 'alarm_history_avg_all'
    call med_time_alarmInit(mclock, alarm, option=hist_option, opt_n=hist_n, &
         reftime=mStartTime, alarmname=trim(alarmname), rc=rc)
    call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockSet(mclock, currTime=mcurrtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       write(logunit,F01) trim(subname)//" set average mediator history alarm "//&
            trim(alarmname)//"  with option "//trim(hist_option)//" and frequency ",hist_n
    end if
    do n = 1,ncomps
       if (is_local%wrap%comp_present(n)) then
          alarmname = 'alarm_history_avg_' // trim(compname(n))
          call med_time_alarmInit(mclock, alarm, option=hist_option, opt_n=hist_n, &
               reftime=mStartTime, alarmname=trim(alarmname), rc=rc)
          call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockAdvance(mclock,rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockSet(mclock, currTime=mcurrtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (mastertask) then
             write(logunit,F01) trim(subname)//" set average mediator history alarm "//&
                  trim(alarmname)//"  with option "//trim(hist_option)//" and frequency ",hist_n
          end if
       end if
    end do

    if (mastertask) write(logunit,*)

  end subroutine med_phases_history_alarms_init

  !===============================================================================
  subroutine med_phases_history_write(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for all variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_phases_history_write)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call med_phases_history_write_file_inst(gcomp, 'all', 'alarm_history_inst_all', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)
  end subroutine med_phases_history_write

  !===============================================================================
  subroutine med_phases_history_write_atm(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for atm variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_phases_history_write_atm)'
    !---------------------------------------
    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call med_phases_history_write_file_inst(gcomp, 'atm', 'alarm_history_inst_atm', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)
  end subroutine med_phases_history_write_atm

  !===============================================================================
  subroutine med_phases_history_write_ice(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for ice variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_phases_history_write_ice)'
    !---------------------------------------
    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call med_phases_history_write_file_inst(gcomp, 'ice', 'alarm_history_inst_ice', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)
  end subroutine med_phases_history_write_ice

  !===============================================================================
  subroutine med_phases_history_write_glc(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for glc variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_phases_history_write_glc)'
    !---------------------------------------
    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call med_phases_history_write_file_inst(gcomp, 'glc', 'alarm_history_inst_glc', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)
  end subroutine med_phases_history_write_glc

  !===============================================================================
  subroutine med_phases_history_write_lnd(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for lnd variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_phases_history_write_lnd)'
    !---------------------------------------
    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call med_phases_history_write_file_inst(gcomp, 'lnd', 'alarm_history_inst_lnd', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)
  end subroutine med_phases_history_write_lnd

  !===============================================================================
  subroutine med_phases_history_write_ocn(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for ocn variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_phases_history_write_ocn)'
    !---------------------------------------
    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call med_phases_history_write_file_inst(gcomp, 'ocn', 'alarm_history_inst_ocn', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)
  end subroutine med_phases_history_write_ocn

  !===============================================================================
  subroutine med_phases_history_write_rof(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for rof variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_phases_history_write_rof)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call med_phases_history_write_file_inst(gcomp, 'rof', 'alarm_history_inst_rof', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)
  end subroutine med_phases_history_write_rof

  !===============================================================================
  subroutine med_phases_history_write_wav(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for wav variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(med_phases_history_write_wav)'
    !---------------------------------------
    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    call med_phases_history_write_file_inst(gcomp, 'wav', 'alarm_history_inst_wav', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)
  end subroutine med_phases_history_write_wav

  !===============================================================================
  subroutine med_phases_history_write_file_inst(gcomp, type, alarmname, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    character(len=*)    , intent(in)    :: type
    character(len=*)    , intent(in)    :: alarmname
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Clock)        :: mclock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_VM)           :: vm
    type(ESMF_Calendar)     :: calendar       ! calendar type
    type(InternalState)     :: is_local
    integer                 :: i,j,m,n
    integer                 :: nx,ny          ! global grid size
    character(CL)           :: time_units     ! units of time variable
    character(CL)           :: hist_file
    real(r8)                :: days_since     ! Time interval since reference time
    real(r8)                :: tbnds(2)       ! CF1.0 time bounds
    logical                 :: whead,wdata    ! for writing restart/history cdf files
    integer                 :: iam
    character(len=*), parameter :: subname='(med_phases_history_write_file)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the communicator and localpet
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the history file alarm
    call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetAlarm(mclock, alarmname=trim(alarmname), alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (dbug_flag > 2) then
       call med_phases_history_output_alarminfo(mclock, alarm, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Check if history alarm is ringing - and if so write the mediator history file
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Turn ringer off
       call ESMF_AlarmRingerOff(alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Determine history file name and time units
       call med_phases_history_get_filename(gcomp, type, hist_file, time_units, days_since, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Create history file
       call med_io_wopen(hist_file, vm, iam, clobber=.true.)
       do m = 1,2

          if (m == 1) then
             whead = .true.
             wdata = .false.
          else if (m == 2) then
             whead = .false.
             wdata = .true.
             call med_io_enddef(hist_file)
          end if

          ! write time values (tbnds does not appear in instantaneous output)
          call ESMF_ClockGet(mclock, calendar=calendar, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_write(hist_file, iam, time_units=time_units, calendar=calendar, time_val=days_since, &
               nt=1, whead=whead, wdata=wdata, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! write field bundles
          do n = 1,ncomps
             if (type == 'all' .or. type == trim(compname(n))) then
                if (is_local%wrap%comp_present(n)) then
                   if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
                      nx = is_local%wrap%nx(n)
                      ny = is_local%wrap%ny(n)
                      call med_io_write(hist_file, iam, is_local%wrap%FBimp(n,n), &
                           nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Imp', rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   endif
                   if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(n),rc=rc)) then
                      nx = is_local%wrap%nx(n)
                      ny = is_local%wrap%ny(n)
                      call med_io_write(hist_file, iam, is_local%wrap%FBexp(n), &
                           nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Exp', rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   endif
                   if (ESMF_FieldBundleIsCreated(is_local%wrap%FBFrac(n),rc=rc)) then
                      nx = is_local%wrap%nx(n)
                      ny = is_local%wrap%ny(n)
                      call med_io_write(hist_file, iam, is_local%wrap%FBFrac(n), &
                           nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='Med_frac_'//trim(compname(n)), rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   end if
                endif
             end if
          enddo
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
             nx = is_local%wrap%nx(compocn)
             ny = is_local%wrap%ny(compocn)
             call med_io_write(hist_file, iam, is_local%wrap%FBMed_ocnalb_o, &
                  nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='Med_alb_ocn', rc=rc)
          end if
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o,rc=rc)) then
             nx = is_local%wrap%nx(compocn)
             ny = is_local%wrap%ny(compocn)
             call med_io_write(hist_file, iam, is_local%wrap%FBMed_aoflux_o, &
                  nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='Med_aoflux_ocn', rc=rc)
          end if
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_a,rc=rc)) then
             nx = is_local%wrap%nx(compatm)
             ny = is_local%wrap%ny(compatm)
             call med_io_write(hist_file, iam, is_local%wrap%FBMed_ocnalb_a, &
                  nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='Med_alb_atm', rc=rc)
          end if
          if (type == 'all' .or. type == 'atm') then
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_a,rc=rc)) then
                nx = is_local%wrap%nx(compatm)
                ny = is_local%wrap%ny(compatm)
                call med_io_write(hist_file, iam, is_local%wrap%FBMed_ocnalb_a, &
                     nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='Med_alb_atm', rc=rc)
             end if
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a,rc=rc)) then
                nx = is_local%wrap%nx(compatm)
                ny = is_local%wrap%ny(compatm)
                call med_io_write(hist_file, iam, is_local%wrap%FBMed_aoflux_a, &
                     nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='Med_aoflux_atm', rc=rc)
             end if
          end if
       end do ! end of loop over m

       ! Close file
       call med_io_close(hist_file, iam, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if ! end of if-alarm is ringingblock

  end subroutine med_phases_history_write_file_inst

  !===============================================================================
  subroutine med_phases_history_write_aux(gcomp, rc)
    ! --------------------------------------
    ! Write mediator history file for wav variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(CS) :: alarmname
    character(CL) :: auxflds
    character(CS) :: cvalue
    integer       :: n,nf
    integer       :: ntapes
    logical       :: first_time = .true.
    logical       :: isPresent
    character(CS) :: prefix
    character(len=*), parameter :: subname='(med_phases_history_write_aux)'
    !---------------------------------------
    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)

    if (first_time) then
       nfiles = 0

       ! currently 5 atm2med tapes
       do nf = 1,5
          write(prefix,'(a,i0)') 'histaux_atm2med_file',nf
          call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_flag', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (trim(cvalue) == 'on') then
             nfiles = nfiles + 1
             auxfiles(nfiles)%ncomp = compatm
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_auxname', value=auxfiles(nfiles)%auxname, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_deltat', value=cvalue, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             read(cvalue,*) auxfiles(nfiles)%deltat
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_ntperfile', value=cvalue, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             read(cvalue,*) auxfiles(nfiles)%ntperfile
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_avgtype', value=auxfiles(nfiles)%avgtype, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_flds', value=auxflds, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_phases_history_setauxflds(auxflds, auxfiles(nfiles)%flds, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do

       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', isPresent=isPresent, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if(isPresent) then
          call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', value=inst_tag, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          inst_tag = ""
       endif
    end if

    do n = 1,nfiles
       call med_phases_history_write_aux_file(gcomp, case_name, inst_tag, first_time, n, auxfiles(n), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    first_time = .false.

    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write_aux

  !===============================================================================
  subroutine med_phases_history_write_aux_file(gcomp, case_name, inst_tag, first_time, nfile_index, auxfile, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    character(len=*)    , intent(in)    :: case_name
    character(len=*)    , intent(in)    :: inst_tag
    logical             , intent(in)    :: first_time
    integer             , intent(in)    :: nfile_index
    type(auxfile_type)  , intent(inout) :: auxfile
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: mclock
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: currtime
    type(ESMF_Calendar)     :: calendar          ! calendar type
    type(ESMF_TimeInterval) :: timediff          ! Used to calculate curr_time
    type(ESMF_VM)           :: vm
    character(CS)           :: timestr
    character(CL)           :: time_units        ! units of time variable
    character(CL)           :: hist_file
    real(r8)                :: days_since        ! Time interval since reference time
    integer                 :: nx,ny             ! global grid size
    logical                 :: whead,wdata       ! for writing restart/history cdf files
    integer                 :: iam               ! mpi task 
    integer                 :: start_ymd         ! Starting date YYYYMMDD
    integer                 :: yr,mon,day,sec    ! time units
    integer                 :: diff_day,diff_sec ! time units
    integer                 :: ncomp
    character(len=*), parameter :: subname='(med_phases_history_write_file)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the communicator and localpet
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine time info
    call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet(mclock, currtime=currtime, starttime=starttime, calendar=calendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currtime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    timediff = currtime - starttime
    call ESMF_TimeIntervalGet(timediff, d=diff_day, s=diff_sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! If write history data
    if (currtime /= starttime .and. mod(sec,auxfile%deltat) == 0) then

       ! Increment number of time samples on file
       auxfile%nt = auxfile%nt + 1

       ! Set tbnds(1)
       if (diff_sec == auxfile%deltat) then
          auxfile%tbnds(1) = 0._r8
       else if ( auxfile%nt == 1 ) then
          auxfile%tbnds(1) = day + sec/real(SecPerDay,R8)
       else
          auxfile%tbnds(1) = auxfile%tbnds(2)
       end if

       ! Set tbnds(2)
       days_since = diff_day + diff_sec/real(SecPerDay,R8)
       auxfile%tbnds(2) = days_since

       ! Create history file 
       ncomp = auxfile%ncomp
       nx    = is_local%wrap%nx(ncomp)
       ny    = is_local%wrap%ny(ncomp)

       ! Write  header
       if (auxfile%nt == 1) then

          ! determine history file name
          call ESMF_ClockGet(mclock, currtime=currtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          write(timestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
          write(hist_file, "(8a)") &
               trim(case_name),'.cpl',trim(inst_tag),'.h', trim(auxfile%auxname),'.',trim(timestr), '.nc'

          ! open file
          call med_io_wopen(hist_file, vm, iam, file_ind=nfile_index, clobber=.true.)

          ! define time units
          call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_ymd2date(yr,mon,day,start_ymd)
          time_units = 'days since ' // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(sec, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! define time variables
          call med_io_write(hist_file, iam, time_units, calendar, days_since, &
               nt=auxfile%nt, tbnds=auxfile%tbnds, whead=.true., wdata=.false., file_ind=nfile_index, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! define data variables with a time dimension (include the nt argument below)
          call med_io_write(hist_file, iam, is_local%wrap%FBimp(ncomp,ncomp), &
               nx=nx, ny=ny, nt=auxfile%nt, whead=.true., wdata=.false., pre=trim(compname(ncomp))//'Imp', &
               flds=auxfile%flds, file_ind=nfile_index, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! end definition phase
          call med_io_enddef(hist_file, file_ind=nfile_index)

       end if

       ! Write time variables for time nt
       call med_io_write(hist_file, iam, time_units, calendar, days_since, &
            nt=auxfile%nt, tbnds=auxfile%tbnds, whead=.false., wdata=.true., file_ind=nfile_index, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Write data variables for time nt
       call med_io_write(hist_file, iam, is_local%wrap%FBimp(ncomp,ncomp), &
            nx=nx, ny=ny, nt=auxfile%nt, whead=.false., wdata=.true., pre=trim(compname(ncomp))//'Imp', &
            flds=auxfile%flds, file_ind=nfile_index, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Close file
       if (auxfile%nt == auxfile%ntperfile) then
          call med_io_close(hist_file, iam, file_ind=nfile_index,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          auxfile%nt = 0
       end if

    end if ! end of write_now if-block

  end subroutine med_phases_history_write_aux_file

  !===============================================================================
  subroutine med_phases_history_get_filename(gcomp, type, hist_file, time_units, days_since, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    character(len=*)    , intent(in)    :: type
    character(len=*)    , intent(out)   :: hist_file
    character(len=*)    , intent(out)   :: time_units
    real(r8)            , intent(out)   :: days_since ! Time interval since reference time
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Clock)        :: mclock
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: nexttime
    type(ESMF_TimeInterval) :: timediff       ! Used to calculate curr_time
    type(ESMF_Calendar)     :: calendar       ! calendar type
    character(len=CS)       :: currtimestr
    character(len=CS)       :: nexttimestr
    integer                 :: start_ymd      ! Starting date YYYYMMDD
    integer                 :: yr,mon,day,sec ! time units
    logical                 :: isPresent
    character(CL)           :: case_name      ! case name
    character(CS)           :: inst_tag   ! instance tag
    character(len=*), parameter :: subname='(med_phases_history_get_timeunits)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get case_name and inst_tag
    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', value=inst_tag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_tag = ""
    endif

    ! Get time unit attribute value for variables
    call NUOPC_ModelGet(gcomp, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet(mclock, currtime=currtime, starttime=starttime, calendar=calendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetNextTime(mclock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    call ESMF_TimeGet(nexttime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    timediff = nexttime - starttime
    call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    days_since = day + sec/real(SecPerDay,R8)
    call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_io_ymd2date(yr,mon,day,start_ymd)
    time_units = 'days since ' // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(sec, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine history file name
    ! Use nexttimestr rather than currtimestr here since that is the time at the end of
    ! the timestep and is preferred for history file names
    if (trim(type) == 'all') then
       write(hist_file,"(6a)") trim(case_name), '.cpl',trim(inst_tag),'.hi.', trim(nexttimestr),'.nc'
    else
       write(hist_file,"(6a)") trim(case_name), '.cpl.'//trim(type),trim(inst_tag),'.hi.', trim(nexttimestr),'.nc'
    end if
    if (mastertask) then
       write(logunit,*)
       write(logunit,' (a)') trim(subname)//": writing mediator history file "//trim(hist_file)
       write(logunit,' (a)') trim(subname)//": currtime = "//trim(currtimestr)
       write(logunit,' (a)') trim(subname)//": nexttime = "//trim(nexttimestr)
    end if

  end subroutine med_phases_history_get_filename

  !===============================================================================
  subroutine med_phases_history_output_alarminfo(mclock, alarm, rc)

    ! input/output variables
    type(ESMF_Clock), intent(in)  :: mclock
    type(ESMF_Alarm), intent(in)  :: alarm
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: ringInterval_length
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: nexttime
    integer                 :: yr,mon,day,sec ! time units
    character(len=CS)       :: currtimestr
    character(len=CS)       :: nexttimestr
    character(len=*), parameter :: subname='(med_phases_history_output_alarminfo)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_AlarmGet(alarm, ringInterval=ringInterval, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(ringInterval, s=ringinterval_length, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet(mclock, currtime=currtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    call ESMF_ClockGetNextTime(mclock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
       write(logunit,*)
       write(logunit,*) trim(subname)//": history alarm ringinterval = ", ringInterval_length
       write(logunit,' (a)') trim(subname)//": currtime = "//trim(currtimestr)//" nexttime = "//trim(nexttimestr)
       write(logunit,*) trim(subname) //' history alarm is ringing = ', ESMF_AlarmIsRinging(alarm)
    end if

  end subroutine med_phases_history_output_alarminfo

  !===============================================================================
  subroutine med_phases_history_setauxflds(str, auxflds, rc)

    ! input/output variables
    character(len=*) , intent(in)  :: str        ! colon deliminted string to search
    character(len=*) , allocatable, intent(out) :: auxflds(:) ! memory will be allocate for auxflds
    integer          , intent(out) :: rc

    ! local variables
    integer          :: i,k,n ! generic indecies
    integer          :: nflds ! allocatable size of auxflds
    integer          :: count ! counts occurances of char
    integer          :: kFlds ! number of fields in list
    integer          :: i0,i1 ! name = list(i0:i1)
    integer          :: nChar ! temporary
    logical          :: valid ! check if str is valid
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! check that this is a str is a valid colon dlimited list
    valid = .true.
    nChar = len_trim(str)
    if (nChar < 1) then                     ! list is an empty string
       valid = .false.
    else if (str(1:1) == ':') then          ! first char is delimiter
       valid = .false.
    else if (str(nChar:nChar) == ':') then  ! last  char is delimiter
       valid = .false.
    else if (index(trim(str)," ") > 0) then ! white-space in a field name
       valid = .false.
    end if
    if (.not. valid) then
       write(logunit,*) "ERROR: invalid list = ",trim(str)
       call ESMF_LogWrite("ERROR: invalid list = "//trim(str), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    ! get number of fields in a colon delimited string list
    nflds = 0
    if (len_trim(str) > 0) then
       count = 0
       do n = 1, len_trim(str)
          if (str(n:n) == ':') count = count + 1
       end do
       nflds = count + 1
    endif

    ! allocate memory for auxflds)
    allocate(auxflds(nflds))

    do k = 1,nflds
       ! start with whole list
       i0 = 1
       i1 = len_trim(str)

       ! remove field names before kth field
       do n = 2,k
          i = index(str(i0:i1),':')
          i0 = i0 + i
       end do

       ! remove field names after kth field
       if (k < nFlds) then
          i = index(str(i0:i1),':')
          i1 = i0 + i - 2
       end if

       ! set auxflds(k)
       auxflds(k) = str(i0:i1)//"   "
    end do

  end subroutine med_phases_history_setauxflds

end module med_phases_history_mod
