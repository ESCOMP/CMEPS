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
  use ESMF                  , only : ESMF_Alarm
  use med_constants_mod     , only : dbug_flag       => med_constants_dbug_flag
  use med_constants_mod     , only : SecPerDay       => med_constants_SecPerDay
  use med_utils_mod         , only : chkerr          => med_utils_ChkErr
  use med_methods_mod       , only : FB_reset        => med_methods_FB_reset
  use med_methods_mod       , only : FB_diagnose     => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_GetFldPtr    => med_methods_FB_GetFldPtr
  use med_methods_mod       , only : FB_accum        => med_methods_FB_accum
  use med_methods_mod       , only : State_GetScalar => med_methods_State_GetScalar
  use med_map_mod           , only : med_map_FB_Regrid_Norm
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_time_mod          , only : med_time_alarmInit
  use med_io_mod            , only : med_io_write, med_io_wopen, med_io_enddef
  use med_io_mod            , only : med_io_close, med_io_date2yyyymmdd, med_io_sec2hms
  use med_io_mod            , only : med_io_ymd2date
  use perf_mod              , only : t_startf, t_stopf
  use esmFlds               , only : ncomps

  implicit none
  private

  public :: med_phases_history_alarm_init
  public :: med_phases_history_write

  ! type(ESMF_Alarm) :: alarm_hist_inst
  ! type(ESMF_Alarm) :: alarm_hist_avg

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_history_alarm_init(gcomp, rc)

    ! --------------------------------------
    ! Initialize mediator history file alarms (module variables)
    ! --------------------------------------

    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF  , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockAdvance, ESMF_ClockSet
    use ESMF  , only : ESMF_Time
    use ESMF  , only : ESMF_TimeInterval, ESMF_TimeIntervalGet
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF  , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF  , only : operator(==), operator(-)
    use ESMF  , only : ESMF_ALARMLIST_ALL, ESMF_Alarm, ESMF_AlarmSet
    use NUOPC , only : NUOPC_CompAttributeGet
    use NUOPC_Model, only : NUOPC_ModelGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Clock)        :: mclock, dclock
    type(ESMF_TimeInterval) :: mtimestep, dtimestep
    type(ESMF_Time)         :: mCurrTime
    type(ESMF_Time)         :: mStartTime
    type(ESMF_TimeInterval) :: timestep
    integer                 :: alarmcount
    integer                 :: timestep_length
    character(CL)           :: cvalue          ! attribute string
    character(CL)           :: histinst_option ! freq_option setting (ndays, nsteps, etc)
    character(CL)           :: histavg_option  ! freq_option setting (ndays, nsteps, etc)
    integer                 :: histinst_n      ! freq_n setting relative to freq_option
    integer                 :: histavg_n       ! freq_n setting relative to freq_option
    character(len=*), parameter :: subname='(med_phases_history_alarm_init)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    ! -----------------------------
    ! Get model clock
    ! -----------------------------

    call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get start time
    call ESMF_ClockGet(mclock, startTime=mStartTime,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------------
    ! Set alarm for instantaneous mediator history output
    ! -----------------------------

    call NUOPC_CompAttributeGet(gcomp, name='history_option', value=histinst_option, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(gcomp, name='history_n', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) histinst_n

    call med_time_alarmInit(mclock, alarm, option=histinst_option, opt_n=histinst_n, &
         reftime=mStartTime, alarmname='alarm_history_inst', rc=rc)

    call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------------
    ! Set alarm for averaged mediator history output
    ! -----------------------------

    !TODO: add isSet and isPresent flags to reading these and other config attributes
    !call NUOPC_CompAttributeGet(gcomp, name='histavg_option', value=histavg_option, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call NUOPC_CompAttributeGet(gcomp, name='histavg_n', value=cvalue, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) histavg_n

    !call med_time_alarmInit(mclock, alarm, option=histavg_option, opt_n=histavg_n, &
    !     reftime=mStartTime, alarmname='alarm_history_avg', rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockGet(mclock, currTime=mCurrTime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(mtimestep, s=timestep_length, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=mcurrtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------------
    ! Write mediator diagnostic output
    ! -----------------------------

    if (mastertask) then
       write(logunit,*)
       write(logunit,100) trim(subname)//" history clock timestep = ",timestep_length
       write(logunit,100) trim(subname)//" set instantaneous mediator history alarm with option "//&
            trim(histinst_option)//" and frequency ",histinst_n
       !write(logunit,100) trim(subname)//" set averaged mediator history alarm with option "//&
       !     trim(histavg_option)//" and frequency ",histavg_n
100    format(a,2x,i8)
       write(logunit,*)
    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": exited", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_phases_history_alarm_init

  !===============================================================================

  subroutine med_phases_history_write(gcomp, rc)

    ! --------------------------------------
    ! Write mediator history file
    ! --------------------------------------

    use ESMF    , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF    , only : ESMF_VM, ESMF_VMGet
    use ESMF    , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockGetNextTime, ESMF_ClockGetAlarm
    use ESMF    , only : ESMF_Calendar
    use ESMF    , only : ESMF_Time, ESMF_TimeGet
    use ESMF    , only : ESMF_TimeInterval, ESMF_TimeIntervalGet
    use ESMF    , only : ESMF_Alarm, ESMF_AlarmIsRinging, ESMF_AlarmRingerOff, ESMF_AlarmGet
    use ESMF    , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF    , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF    , only : ESMF_FieldBundleIsCreated, ESMF_MAXSTR, ESMF_ClockPrint, ESMF_AlarmIsCreated
    use ESMF    , only : operator(==), operator(-)
    use ESMF    , only : ESMF_ALARMLIST_ALL, ESMF_ClockGetAlarmList
    use NUOPC   , only : NUOPC_CompAttributeGet
    use esmFlds , only : compatm, complnd, compocn, compice, comprof, compglc, ncomps, compname
    use esmFlds , only : fldListFr, fldListTo
    use NUOPC_Model, only : NUOPC_ModelGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: mclock, dclock
    type(ESMF_TimeInterval) :: mtimestep, dtimestep
    integer                 :: timestep_length
    type(ESMF_Alarm)        :: alarm
    integer                 :: alarmCount
    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: reftime
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: nexttime
    type(ESMF_TimeInterval) :: timediff       ! Used to calculate curr_time
    type(ESMF_Calendar)     :: calendar       ! calendar type
    character(len=64)       :: currtimestr
    character(len=64)       :: nexttimestr
    type(InternalState)     :: is_local
    character(CS)           :: histavg_option ! Histavg option units
    integer                 :: i,j,m,n,n1,ncnt
    integer                 :: start_ymd      ! Starting date YYYYMMDD
    integer                 :: start_tod      ! Starting time-of-day (s)
    integer                 :: nx,ny          ! global grid size
    integer                 :: yr,mon,day,sec ! time units
    real(r8)                :: rval           ! real tmp value
    real(r8)                :: dayssince      ! Time interval since reference time
    integer                 :: fk             ! index
    character(CL)           :: time_units     ! units of time variable
    character(CL)           :: case_name      ! case name
    character(CL)           :: hist_file      ! Local path to history filename
    character(CS)           :: cpl_inst_tag   ! instance tag
    character(CL)           :: cvalue         ! attribute string
    real(r8)                :: tbnds(2)       ! CF1.0 time bounds
    logical                 :: whead,wdata    ! for writing restart/history cdf files
    integer                 :: iam
    logical                 :: isPresent
    type(ESMF_TimeInterval) :: RingInterval
    integer                 :: ringInterval_length
    logical                 :: first_time = .true.
    character(len=*), parameter :: subname='(med_phases_history_write)'
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the communicator and localpet
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', value=cpl_inst_tag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       cpl_inst_tag = ""
    endif

    call NUOPC_ModelGet(gcomp, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (first_time) then
       call med_phases_history_alarm_init(gcomp, rc)
    end if

    !---------------------------------------
    ! Check if history alarm is ringing - and if so write the mediator history file
    !---------------------------------------

    ! TODO: Add history averaging functionality and Determine if history average alarm is on
    ! if (ESMF_AlarmIsRinging(AlarmHistAvg, rc=rc)) then
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !    alarmIsOn = .true.
    !    call ESMF_AlarmRingerOff( AlarmHist, rc=rc )
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else
    !    alarmisOn = .false.
    ! endif

    call ESMF_ClockGetAlarm(mclock, alarmname='alarm_history_inst', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 2) then
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
    end if

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Turn ringer off
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Get time info for history file
       call ESMF_GridCompGet(gcomp, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGet(mclock, currtime=currtime, reftime=reftime, starttime=starttime, calendar=calendar, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGetNextTime(mclock, nextTime=nexttime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec

       call ESMF_TimeGet(nexttime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
       timediff = nexttime - reftime
       call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       dayssince = day + sec/real(SecPerDay,R8)

       call ESMF_TimeGet(reftime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_io_ymd2date(yr,mon,day,start_ymd)
       start_tod = sec
       time_units = 'days since ' // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(start_tod, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Use nexttimestr rather than currtimestr here since that is the time at the end of
       ! the timestep and is preferred for history file names
       write(hist_file,"(6a)") trim(case_name), '.cpl',trim(cpl_inst_tag),'.hi.', trim(nexttimestr),'.nc'

       if (mastertask) then
          write(logunit,*)
          write(logunit,' (a)') trim(subname)//": writing mediator history file "//trim(hist_file)
          write(logunit,' (a)') trim(subname)//": currtime = "//trim(currtimestr)
          write(logunit,' (a)') trim(subname)//": nexttime = "//trim(nexttimestr)
       end if

       call med_io_wopen(hist_file, vm, iam, clobber=.true.)
       do m = 1,2
          whead=.false.
          wdata=.false.
          if (m == 1) then
             whead=.true.
          elseif (m == 2) then
             wdata=.true.
             call med_io_enddef(hist_file)
          endif

          tbnds = dayssince

          if (tbnds(1) >= tbnds(2)) then
             call med_io_write(hist_file, iam, time_units=time_units, calendar=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call med_io_write(hist_file, iam, time_units=time_units, calendar=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata, tbnds=tbnds, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          do n = 1,ncomps
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
          enddo
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
             nx = is_local%wrap%nx(compocn)
             ny = is_local%wrap%ny(compocn)
             call med_io_write(hist_file, iam, is_local%wrap%FBMed_ocnalb_o, &
                  nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='Med_alb_ocn', rc=rc)
          end if
          !TODO: don't write aoflux_(oa) when they're not being used
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
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a,rc=rc)) then
             nx = is_local%wrap%nx(compatm)
             ny = is_local%wrap%ny(compatm)
             call med_io_write(hist_file, iam, is_local%wrap%FBMed_aoflux_a, &
                  nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='Med_aoflux_atm', rc=rc)
          end if
       enddo

       call med_io_close(hist_file, iam, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    endif

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

    first_time = .false.

  end subroutine med_phases_history_write

  !===============================================================================

end module med_phases_history_mod
