module esm_time_mod

  use shr_kind_mod        , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, r8=>shr_kind_r8
  use ESMF                , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSet
  use ESMF                , only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockGet, ESMF_ClockSet
  use ESMF                , only : ESMF_ClockAdvance
  use ESMF                , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmGet
  use ESMF                , only : ESMF_Calendar, ESMF_CalKind_Flag, ESMF_CalendarCreate
  use ESMF                , only : ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
  use ESMF                , only : ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF                , only : ESMF_TimeInterval, ESMF_TimeIntervalSet, ESMF_TimeIntervalGet
  use ESMF                , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FAILURE, ESMF_LOGMSG_ERROR
  use ESMF                , only : ESMF_VM, ESMF_VMGet, ESMF_VMBroadcast
  use ESMF                , only : ESMF_VMAllReduce, ESMF_REDUCE_MAX
  use ESMF                , only : ESMF_LOGMSG_INFO, ESMF_FAILURE, ESMF_GridCompIsPetLocal
  use ESMF                , only : operator(<), operator(/=), operator(+)
  use ESMF                , only : operator(-), operator(*) , operator(>=)
  use ESMF                , only : operator(<=), operator(>), operator(==)
  use NUOPC               , only : NUOPC_CompAttributeGet
  use esm_utils_mod       , only : chkerr

  implicit none
  private    ! default private

  public  :: esm_time_clockInit  ! initialize driver clock (assumes default calendar)

!  private :: esm_time_timeInit
  private :: esm_time_alarmInit
  private :: esm_time_date2ymd

  ! Clock and alarm options
  character(len=*), private, parameter :: &
       optNONE           = "none"    , &
       optNever          = "never"   , &
       optNSteps         = "nstep"   , &
       optNSeconds       = "nsecond" , &
       optNMinutes       = "nminute" , &
       optNHours         = "nhour"   , &
       optNDays          = "nday"    , &
       optNMonths        = "nmonth"  , &
       optNYears         = "nyear"   , &
       optMonthly        = "monthly" , &
       optYearly         = "yearly"  , &
       optDate           = "date"    , &
       optGLCCouplingPeriod = "glc_coupling_period"

  ! Module data
  integer, parameter          :: SecPerDay = 86400 ! Seconds per day
  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine esm_time_clockInit(ensemble_driver, instance_driver, logunit, maintask, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: ensemble_driver, instance_driver
    integer, intent(in)  :: logunit
    logical, intent(in)  :: maintask
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_VM)           :: vm, envm
    type(ESMF_Time)         :: StartTime           ! Start time
    type(ESMF_Time)         :: RefTime             ! Reference time
    type(ESMF_Time)         :: CurrTime            ! Current time
    type(ESMF_Time)         :: StopTime            ! Stop time
    type(ESMF_Time)         :: Clocktime           ! Loop time
    type(ESMF_TimeInterval) :: TimeStep            ! Clock time-step
    type(ESMF_Alarm)        :: alarm_stop          ! alarm
    integer                 :: ref_ymd             ! Reference date (YYYYMMDD)
    integer                 :: ref_tod             ! Reference time of day (seconds)
    integer                 :: start_ymd           ! Start date (YYYYMMDD)
    integer                 :: start_tod           ! Start time of day (seconds)
    integer                 :: curr_ymd            ! Current ymd (YYYYMMDD)
    integer                 :: curr_tod            ! Current tod (seconds)
    integer                 :: stop_n              ! Number until stop
    integer                 :: stop_ymd            ! Stop date (YYYYMMDD)
    integer                 :: stop_tod            ! Stop time-of-day
    character(CS)           :: stop_option         ! Stop option units
    integer                 :: atm_cpl_dt          ! Atmosphere coupling interval
    integer                 :: lnd_cpl_dt          ! Land coupling interval
    integer                 :: ice_cpl_dt          ! Sea-Ice coupling interval
    integer                 :: ocn_cpl_dt          ! Ocean coupling interval
    integer                 :: glc_cpl_dt          ! Glc coupling interval
    integer                 :: rof_cpl_dt          ! Runoff coupling interval
    integer                 :: wav_cpl_dt          ! Wav coupling interval
!    integer                 :: esp_cpl_dt          ! Esp coupling interval
    character(CS)           :: glc_avg_period      ! Glc avering coupling period
    logical                 :: read_restart
    character(len=CL)       :: restart_file
    character(len=CL)       :: restart_pfile
    character(len=CL)       :: cvalue
    integer                 :: dtime_drv           ! time-step to use
    integer                 :: yr, mon, day        ! Year, month, day as integers
    integer                 :: unitn               ! unit number
    integer                 :: ierr                ! Return code
    character(CL)           :: tmpstr              ! temporary
    character(CS)           :: inst_suffix
    integer                 :: tmp(4)              ! Array for Broadcast
    integer                 :: myid, bcastID(2)
    logical                 :: isPresent
    logical                 :: inDriver
    logical, save           :: firsttime=.true.
    character(len=*), parameter :: subname = '('//__FILE__//':esm_time_clockInit) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    !---------------------------------------------------------------------------
    ! Determine start time, reference time and current time
    !---------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(ensemble_driver, name="start_ymd", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) start_ymd
    call NUOPC_CompAttributeGet(ensemble_driver, name="start_tod", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) start_tod

    !---------------------------------------------------------------------------
    ! Determine driver clock timestep
    !---------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(ensemble_driver, name="atm_cpl_dt", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) atm_cpl_dt

    call NUOPC_CompAttributeGet(ensemble_driver, name="lnd_cpl_dt", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lnd_cpl_dt

    call NUOPC_CompAttributeGet(ensemble_driver, name="ice_cpl_dt", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ice_cpl_dt

    call NUOPC_CompAttributeGet(ensemble_driver, name="ocn_cpl_dt", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ocn_cpl_dt

    call NUOPC_CompAttributeGet(ensemble_driver, name="glc_cpl_dt", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_cpl_dt

    call NUOPC_CompAttributeGet(ensemble_driver, name="rof_cpl_dt", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) rof_cpl_dt

    call NUOPC_CompAttributeGet(ensemble_driver, name="wav_cpl_dt", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wav_cpl_dt

    dtime_drv = minval((/atm_cpl_dt, lnd_cpl_dt, ocn_cpl_dt, ice_cpl_dt, glc_cpl_dt, rof_cpl_dt, wav_cpl_dt/))
    if(maintask) then
       write(tmpstr,'(i10)') dtime_drv
       call ESMF_LogWrite(trim(subname)//': driver time interval is : '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
       write(logunit,*)   trim(subname)//': driver time interval is : '// trim(tmpstr)
    endif
    call ESMF_GridCompGet(ensemble_driver, vm=envm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    call ESMF_VMGet(envm, localPet=myid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    indriver = ESMF_GridCompIsPetLocal(instance_driver, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if(indriver) then
       call ESMF_GridCompGet(instance_driver, vm=vm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(instance_driver, name='read_restart', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) read_restart

       if (read_restart) then
          
          call NUOPC_CompAttributeGet(instance_driver, name='drv_restart_pointer', value=restart_file, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (trim(restart_file) /= 'none') then

             call NUOPC_CompAttributeGet(instance_driver, name="inst_suffix", isPresent=isPresent, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if(isPresent) then
                call NUOPC_CompAttributeGet(instance_driver, name="inst_suffix", value=inst_suffix, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                inst_suffix = ""
             endif

             restart_pfile = trim(restart_file)//inst_suffix

             if (maintask) then
                call ESMF_LogWrite(trim(subname)//" read rpointer file = "//trim(restart_pfile), &
                     ESMF_LOGMSG_INFO)
                open(newunit=unitn, file=restart_pfile, form='FORMATTED', status='old',iostat=ierr)
                if (ierr < 0) then
                   rc = ESMF_FAILURE
                   call ESMF_LogWrite(trim(subname)//' ERROR rpointer file open returns error', &
                        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
                   return
                end if
                read(unitn,'(a)', iostat=ierr) restart_file
                if (ierr < 0) then
                   rc = ESMF_FAILURE
                   call ESMF_LogWrite(trim(subname)//' ERROR rpointer file read returns error', &
                        ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__)
                   return
                end if
                close(unitn)
                if (maintask) then
                   write(logunit,'(a)') trim(subname)//" reading driver restart from file = "//trim(restart_file)
                end if
                call esm_time_read_restart(restart_file, start_ymd, start_tod, curr_ymd, curr_tod, rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             endif
          else
             if(maintask) then
                write(logunit,*) ' NOTE: the current compset has no mediator - which provides the clock restart information'
                write(logunit,*) '   In this case the restarts are handled solely by the component being used and'
                write(logunit,*) '   and the driver clock will always be starting from the initial date on restart'
             end if
             curr_ymd = start_ymd
             curr_tod = start_tod
          endif
       else
          curr_ymd = start_ymd
          curr_tod = start_tod                                            
       end if ! end if read_restart
    endif
    if(maintask) then
       bcastID(1) = myid
       tmp(1) = start_ymd ; tmp(2) = start_tod
       tmp(3) = curr_ymd  ; tmp(4) = curr_tod
    else
       bcastID(1) = 0
       tmp = 0
    endif
    call ESMF_VMAllReduce(envm, bcastID(1:1), bcastID(2:2), 1, ESMF_REDUCE_MAX,rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMBroadcast(envm, tmp, 4, bcastID(2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    start_ymd = tmp(1) ; start_tod = tmp(2)
    curr_ymd  = tmp(3) ; curr_tod  = tmp(4)

    ! Determine start time (THE FOLLOWING ASSUMES THAT THE DEFAULT CALENDAR IS SET in the driver)

    call esm_time_date2ymd(start_ymd, yr, mon, day)
    call ESMF_TimeSet( StartTime, yy=yr, mm=mon, dd=day, s=start_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if(maintask) then
       write(tmpstr,'(i10)') start_ymd
       call ESMF_LogWrite(trim(subname)//': driver start_ymd: '// trim(tmpstr), ESMF_LOGMSG_INFO)
       write(logunit,*)   trim(subname)//': driver start_ymd: '// trim(tmpstr)
       write(tmpstr,'(i10)') start_tod
       call ESMF_LogWrite(trim(subname)//': driver start_tod: '// trim(tmpstr), ESMF_LOGMSG_INFO)
       write(logunit,*)   trim(subname)//': driver start_tod: '// trim(tmpstr)
    endif

    ! Determine current time
    call esm_time_date2ymd(curr_ymd, yr, mon, day)
    call ESMF_TimeSet( CurrTime, yy=yr, mm=mon, dd=day, s=curr_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if(maintask) then
       write(tmpstr,'(i10)') curr_ymd
       call ESMF_LogWrite(trim(subname)//': driver curr_ymd: '// trim(tmpstr), ESMF_LOGMSG_INFO)
       write(logunit,*)   trim(subname)//': driver curr_ymd: '// trim(tmpstr)
       write(tmpstr,'(i10)') curr_tod
       call ESMF_LogWrite(trim(subname)//': driver curr_tod: '// trim(tmpstr), ESMF_LOGMSG_INFO)
       write(logunit,*)   trim(subname)//': driver curr_tod: '// trim(tmpstr)
    endif
    ! Set reference time - HARD-CODED TO START TIME
    ref_ymd = start_ymd
    ref_tod = start_tod
    call esm_time_date2ymd(ref_ymd, yr, mon, day)
    call ESMF_TimeSet( RefTime, yy=yr, mm=mon, dd=day, s=ref_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    call ESMF_TimeIntervalSet( TimeStep, s=dtime_drv, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Create an instance_driver clock
    !---------------------------------------------------------------------------

    ! Create the clock
    clock = ESMF_ClockCreate(TimeStep, StartTime, refTime=RefTime, name='ESMF Driver Clock', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Advance the clock to the current time (in case of a restart)
    call ESMF_ClockGet(clock, currTime=clocktime, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do while( clocktime < CurrTime)
       call ESMF_ClockAdvance( clock, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet( clock, currTime=clocktime, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! Set the driver gridded component clock to the created clock
    if (indriver) then
       call ESMF_GridCompSet(instance_driver, clock=clock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! Set driver clock stop time
    call NUOPC_CompAttributeGet(ensemble_driver, name="stop_option", value=stop_option, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(ensemble_driver, name="stop_n", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_n
    call NUOPC_CompAttributeGet(ensemble_driver, name="stop_ymd", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_ymd
    call NUOPC_CompAttributeGet(ensemble_driver, name="stop_tod", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_tod
    if ( stop_ymd < 0) then
       stop_ymd = 99990101
       stop_tod = 0
    endif

    if (maintask) then
       write(tmpstr,'(i10)') stop_ymd
       call ESMF_LogWrite(trim(subname)//': driver stop_ymd: '// trim(tmpstr), ESMF_LOGMSG_INFO)
       write(logunit,*)   trim(subname)//': driver stop_ymd: '// trim(tmpstr)
       write(tmpstr,'(i10)') stop_tod
       call ESMF_LogWrite(trim(subname)//': driver stop_tod: '// trim(tmpstr), ESMF_LOGMSG_INFO)
       write(logunit,*)   trim(subname)//': driver stop_tod: '// trim(tmpstr)
    endif

    call esm_time_alarmInit(clock, &
         alarm   = alarm_stop,           &
         option  = stop_option,          &
         opt_n   = stop_n,               &
         opt_ymd = stop_ymd,             &
         opt_tod = stop_tod,             &
         RefTime = CurrTime,             &
         alarmname = 'alarm_stop', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AlarmGet(alarm_stop, RingTime=StopTime, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(clock, StopTime=StopTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Create the ensemble driver clock
    !---------------------------------------------------------------------------
    if(firsttime) then
       TimeStep = StopTime - ClockTime
       clock = ESMF_ClockCreate(TimeStep, ClockTime, StopTime=StopTime, &
            refTime=RefTime, name='ESMF ensemble Driver Clock', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_GridCompSet(ensemble_driver, clock=clock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       firsttime = .false.
    endif
 end subroutine esm_time_clockInit

 !===============================================================================

 subroutine esm_time_alarmInit( clock, alarm, option, &
      opt_n, opt_ymd, opt_tod, RefTime, alarmname, rc)

   ! Setup an alarm in a clock
   ! Notes: The ringtime sent to AlarmCreate MUST be the next alarm
   ! time.  If you send an arbitrary but proper ringtime from the
   ! past and the ring interval, the alarm will always go off on the
   ! next clock advance and this will cause serious problems.  Even
   ! if it makes sense to initialize an alarm with some reference
   ! time and the alarm interval, that reference time has to be
   ! advance forward to be >= the current time.  In the logic below
   ! we set an appropriate "NextAlarm" and then we make sure to
   ! advance it properly based on the ring interval.

   ! input/output variables
   type(ESMF_Clock)            , intent(inout) :: clock     ! clock
   type(ESMF_Alarm)            , intent(inout) :: alarm     ! alarm
   character(len=*)            , intent(in)    :: option    ! alarm option
   integer          , optional , intent(in)    :: opt_n     ! alarm freq
   integer          , optional , intent(in)    :: opt_ymd   ! alarm ymd
   integer          , optional , intent(in)    :: opt_tod   ! alarm tod (sec)
   type(ESMF_Time)  , optional , intent(in)    :: RefTime   ! ref time
   character(len=*) , optional , intent(in)    :: alarmname ! alarm name
   integer                     , intent(inout) :: rc        ! Return code

   ! local variables
   type(ESMF_Calendar)     :: cal              ! calendar
   integer                 :: lymd             ! local ymd
   integer                 :: ltod             ! local tod
   integer                 :: cyy,cmm,cdd,csec ! time info
   character(len=64)       :: lalarmname       ! local alarm name
   logical                 :: update_nextalarm ! update next alarm
   type(ESMF_Time)         :: CurrTime         ! Current Time
   type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
   type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
   character(len=*), parameter :: subname = '(med_time_alarmInit): '
   !-------------------------------------------------------------------------------

   rc = ESMF_SUCCESS

   lalarmname = 'alarm_unknown'
   if (present(alarmname)) lalarmname = trim(alarmname)
   ltod = 0
   if (present(opt_tod)) ltod = opt_tod
   lymd = -1
   if (present(opt_ymd)) lymd = opt_ymd

   call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
   if (ChkErr(rc,__LINE__,u_FILE_u)) return

   call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
   if (ChkErr(rc,__LINE__,u_FILE_u)) return

   ! initial guess of next alarm, this will be updated below
   if (present(RefTime)) then
      NextAlarm = RefTime
   else
      NextAlarm = CurrTime
   endif

   ! Get calendar from clock
   call ESMF_ClockGet(clock, calendar=cal, rc=rc)
   if (ChkErr(rc,__LINE__,u_FILE_u)) return

   ! Error checks
   if (trim(option) == optdate) then
      if (.not. present(opt_ymd)) then
         call ESMF_LogWrite(trim(subname)//trim(option)//' requires opt_ymd', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
      end if
      if (lymd < 0 .or. ltod < 0) then
         call ESMF_LogWrite(subname//trim(option)//'opt_ymd, opt_tod invalid', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
      end if
    else if (&
         trim(option) == optNSteps   .or. trim(option) == trim(optNSteps)//'s'   .or. &
         trim(option) == optNSeconds .or. trim(option) == trim(optNSeconds)//'s' .or. &
         trim(option) == optNMinutes .or. trim(option) == trim(optNMinutes)//'s' .or. &
         trim(option) == optNHours   .or. trim(option) == trim(optNHours)//'s'   .or. &
         trim(option) == optNDays    .or. trim(option) == trim(optNDays)//'s'    .or. &
         trim(option) == optNMonths  .or. trim(option) == trim(optNMonths)//'s'  .or. &
         trim(option) == optNYears   .or. trim(option) == trim(optNYears)//'s' ) then
      if (.not.present(opt_n)) then
         call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
      end if
      if (opt_n <= 0) then
         call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
      end if
    end if

   ! Determine inputs for call to create alarm
   selectcase (trim(option))

   case (optNONE)
      call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      update_nextalarm  = .false.

   case (optDate)
      call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call esm_time_date2ymd(opt_ymd, cyy, cmm, cdd)

      call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=cdd, s=ltod, calendar=cal, rc=rc )
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      update_nextalarm  = .false.

   case (optNever)
      call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      update_nextalarm  = .false.

   case (optNSteps,trim(optNSteps)//'s')
      call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      AlarmInterval = AlarmInterval * opt_n
      update_nextalarm  = .true.

   case (optNSeconds,trim(optNSeconds)//'s')
      call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      AlarmInterval = AlarmInterval * opt_n
      update_nextalarm  = .true.

   case (optNMinutes,trim(optNMinutes)//'s')
      call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
      AlarmInterval = AlarmInterval * opt_n
      update_nextalarm  = .true.

   case (optNHours,trim(optNHours)//'s')
      call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      AlarmInterval = AlarmInterval * opt_n
      update_nextalarm  = .true.

   case (optNDays,trim(optNDays)//'s')
      call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      AlarmInterval = AlarmInterval * opt_n
      update_nextalarm  = .true.

   case (optNMonths,trim(optNMonths)//'s')
      call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      AlarmInterval = AlarmInterval * opt_n
      update_nextalarm  = .true.

   case (optMonthly)
      call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=cal, rc=rc )
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      update_nextalarm  = .true.

   case (optNYears, trim(optNYears)//'s')
      call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      AlarmInterval = AlarmInterval * opt_n
      update_nextalarm  = .true.

   case (optYearly)
      call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=cal, rc=rc )
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      update_nextalarm  = .true.

   case default
      call ESMF_LogWrite(subname//'unknown option '//trim(option), ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return

   end select

   ! --------------------------------------------------------------------------------
   ! --- AlarmInterval and NextAlarm should be set ---
   ! --------------------------------------------------------------------------------

   ! --- advance Next Alarm so it won't ring on first timestep for
   ! --- most options above. go back one alarminterval just to be careful

   if (update_nextalarm) then
      NextAlarm = NextAlarm - AlarmInterval
      do while (NextAlarm <= CurrTime)
         NextAlarm = NextAlarm + AlarmInterval
      enddo
   endif

   alarm = ESMF_AlarmCreate( name=lalarmname, clock=clock, ringTime=NextAlarm, &
        ringInterval=AlarmInterval, rc=rc)
   if (ChkErr(rc,__LINE__,u_FILE_u)) return

 end subroutine esm_time_alarmInit

 !===============================================================================
#ifdef UNUSEDFUNCTION
 subroutine esm_time_timeInit( Time, ymd, cal, tod, desc, logunit )

   !  Create the ESMF_Time object corresponding to the given input time, given in
   !  YMD (Year Month Day) and TOD (Time-of-day) format.
   !  Set the time by an integer as YYYYMMDD and integer seconds in the day

   ! input/output parameters:
   type(ESMF_Time)     , intent(inout)        :: Time ! ESMF time
   integer             , intent(in)           :: ymd  ! year, month, day YYYYMMDD
   type(ESMF_Calendar) , intent(in)           :: cal  ! ESMF calendar
   integer             , intent(in), optional :: tod  ! time of day in seconds
   character(len=*)    , intent(in), optional :: desc ! description of time to set
   integer             , intent(in), optional :: logunit

   ! local variables
   integer                     :: yr, mon, day ! Year, month, day as integers
   integer                     :: ltod         ! local tod
   character(len=256)          :: ldesc        ! local desc
   integer                     :: rc           ! return code
   character(len=*), parameter :: subname = '(esm_time_m_ETimeInit) '
   !-------------------------------------------------------------------------------

   ltod = 0
   if (present(tod)) ltod = tod
   ldesc = ''
   if (present(desc)) ldesc = desc

   if ( (ymd < 0) .or. (ltod < 0) .or. (ltod > SecPerDay) )then
      if (present(logunit)) then
         write(logunit,*) subname//': ERROR yymmdd is a negative number or '// &
              'time-of-day out of bounds', ymd, ltod
      end if
      call ESMF_LogWrite( subname//'ERROR: Bad input' , ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if

   call esm_time_date2ymd (ymd,yr,mon,day)

   call ESMF_TimeSet( Time, yy=yr, mm=mon, dd=day, s=ltod, calendar=cal, rc=rc )
   if (ChkErr(rc,__LINE__,u_FILE_u)) return

 end subroutine esm_time_timeInit
#endif
 !===============================================================================

 subroutine esm_time_date2ymd (date, year, month, day)

   ! input/output variables
   integer, intent(in)  :: date             ! coded-date (yyyymmdd)
   integer, intent(out) :: year,month,day   ! calendar year,month,day

   ! local variables
   integer :: tdate   ! temporary date
   character(*),parameter :: subName = "(esm_time_date2ymd)"
   !-------------------------------------------------------------------------------

   tdate = abs(date)
   year = int(tdate/10000)
   if (date < 0) then
      year = -year
   end if
   month = int( mod(tdate,10000)/  100)
   day = mod(tdate,  100)

 end subroutine esm_time_date2ymd

 !===============================================================================

 subroutine esm_time_read_restart(restart_file, start_ymd, start_tod, curr_ymd, curr_tod, rc)

   use netcdf , only : nf90_open, nf90_nowrite, nf90_noerr
   use netcdf , only : nf90_inq_varid, nf90_get_var, nf90_close
   use ESMF   , only : ESMF_LogWrite, ESMF_LOGMSG_INFO

   ! input/output variables
   character(len=*), intent(in) :: restart_file
   integer, intent(out)         :: start_ymd           ! Start date (YYYYMMDD)
   integer, intent(out)         :: start_tod           ! Start time of day (seconds)
   integer, intent(out)         :: curr_ymd            ! Current ymd (YYYYMMDD)
   integer, intent(out)         :: curr_tod            ! Current tod (seconds)
   integer, intent(out)         :: rc

   ! local variables
   integer                 :: status, ncid, varid ! netcdf stuff
   character(CL)           :: tmpstr              ! temporary
   character(len=*), parameter :: subname = "(esm_time_read_restart)"
   !----------------------------------------------------------------

   ! use netcdf here since it's serial
   rc = ESMF_SUCCESS
   status = nf90_open(restart_file, NF90_NOWRITE, ncid)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_open: '//trim(restart_file), ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   endif

   status = nf90_inq_varid(ncid, 'start_ymd', varid)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid start_ymd', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if
   status = nf90_get_var(ncid, varid, start_ymd)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var start_ymd', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if

   status = nf90_inq_varid(ncid, 'start_tod', varid)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid start_tod', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if
   status = nf90_get_var(ncid, varid, start_tod)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var start_tod', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if

   status = nf90_inq_varid(ncid, 'curr_ymd', varid)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid curr_ymd', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if
   status = nf90_get_var(ncid, varid, curr_ymd)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var curr_ymd', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if

   status = nf90_inq_varid(ncid, 'curr_tod', varid)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid curr_tod', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if
   status = nf90_get_var(ncid, varid, curr_tod)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var curr_tod', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if

   status = nf90_close(ncid)
   if (status /= nf90_NoErr) then
      call ESMF_LogWrite(trim(subname)//' ERROR: nf90_close', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
   end if

   write(tmpstr,*) trim(subname)//" read start_ymd = ",start_ymd
   call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

   write(tmpstr,*) trim(subname)//" read start_tod = ",start_tod
   call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

   write(tmpstr,*) trim(subname)//" read curr_ymd  = ",curr_ymd
   call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

   write(tmpstr,*) trim(subname)//" read curr_tod  = ",curr_tod
   call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

 end subroutine esm_time_read_restart

end module esm_time_mod
