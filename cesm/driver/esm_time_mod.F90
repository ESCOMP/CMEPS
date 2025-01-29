module esm_time_mod

  use shr_kind_mod        , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, r8=>shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort
  use ESMF                , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSet
  use ESMF                , only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockGet, ESMF_ClockSet
  use ESMF                , only : ESMF_ClockAdvance
  use ESMF                , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmGet
  use ESMF                , only : ESMF_Calendar, ESMF_CalKind_Flag, ESMF_CalendarCreate
  use ESMF                , only : ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
  use ESMF                , only : ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF                , only : ESMF_TimeInterval, ESMF_TimeIntervalSet, ESMF_TimeIntervalGet
  use ESMF                , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FAILURE
  use ESMF                , only : ESMF_VM, ESMF_VMGet, ESMF_VMBroadcast
  use ESMF                , only : ESMF_VMAllReduce, ESMF_REDUCE_MAX, ESMF_ClockGetAlarm
  use ESMF                , only : ESMF_LOGMSG_INFO, ESMF_FAILURE, ESMF_GridCompIsPetLocal
  use ESMF                , only : operator(<), operator(/=), operator(+)
  use ESMF                , only : operator(-), operator(*) , operator(>=)
  use ESMF                , only : operator(<=), operator(>), operator(==)
  use NUOPC               , only : NUOPC_CompAttributeGet
  use esm_utils_mod       , only : chkerr
  use nuopc_shr_methods   , only : AlarmInit
  
  implicit none
  private    ! default private

  public  :: esm_time_clockinit  ! initialize driver clock (assumes default calendar)

  private :: esm_time_date2ymd

  ! Clock and alarm options
  character(len=*), private, parameter :: &
       optNONE           = "none"    , &
       optNever          = "never"   , &
       optNSeconds       = "nsecond" , &
       optNMinutes       = "nminute" , &
       optNHours         = "nhour"   , &
       optNDays          = "nday"    , &
       optNMonths        = "nmonth"  , &
       optNYears         = "nyear"   , &
       optMonthly        = "monthly" , &
       optYearly         = "yearly"  , &
       optDate           = "date"    , &
       optEnd            = "end"     , &
       optGLCCouplingPeriod = "glc_coupling_period"

  ! Module data
  integer, parameter          :: SecPerDay = 86400 ! Seconds per day
  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine esm_time_clockinit(ensemble_driver, instance_driver, logunit, maintask, rc)
    use nuopc_shr_methods, only : get_minimum_timestep, dtime_drv
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
    character(CS)           :: glc_avg_period      ! Glc avering coupling period
    logical                 :: read_restart
    character(len=CL)       :: restart_file
    character(len=CL)       :: restart_pfile
    character(len=CL)       :: cvalue
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
    logical                 :: exists
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
          
          call NUOPC_CompAttributeGet(instance_driver, name='drv_restart_pointer', value=restart_pfile, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (trim(restart_pfile) /= 'none') then

             if (maintask) then
                write(logunit,*) " read rpointer file = "//trim(restart_pfile)
                inquire( file=trim(restart_pfile), exist=exists)
                if (.not. exists) then
                   call shr_sys_abort(trim(subname)//' ERROR rpointer file '//trim(restart_pfile)//' not found',&
                        line=__LINE__, file=__FILE__)
                endif
                call ESMF_LogWrite(trim(subname)//" read rpointer file = "//trim(restart_pfile), &
                     ESMF_LOGMSG_INFO)
                open(newunit=unitn, file=restart_pfile, form='FORMATTED', status='old',iostat=ierr)
                if (ierr < 0) then
                   call shr_sys_abort(trim(subname)//' ERROR rpointer file open returns error', &
                        line=__LINE__, file=__FILE__)
                end if
                read(unitn,'(a)', iostat=ierr) restart_file
                if (ierr < 0) then
                   call shr_sys_abort(trim(subname)//' ERROR rpointer file read returns error', &
                        line=__LINE__, file=__FILE__)
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

    dtime_drv = get_minimum_timestep(ensemble_driver, rc)
    if(maintask) then
       write(tmpstr,'(i10)') dtime_drv
       call ESMF_LogWrite(trim(subname)//': driver time interval is : '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
       write(logunit,*)   trim(subname)//': driver time interval is : '// trim(tmpstr)
    endif
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

    call alarmInit(clock, &
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
  end subroutine esm_time_clockinit

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
      call shr_sys_abort(trim(subname)//' ERROR: nf90_open: '//trim(restart_file),&
           file=__FILE__, line=__LINE__)
   endif

   status = nf90_inq_varid(ncid, 'start_ymd', varid)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid start_ymd', &
           file=__FILE__, line=__LINE__)
   end if
   status = nf90_get_var(ncid, varid, start_ymd)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var start_ymd', &
           file=__FILE__, line=__LINE__)
   end if

   status = nf90_inq_varid(ncid, 'start_tod', varid)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid start_tod', &
           file=__FILE__, line=__LINE__)
   end if
   status = nf90_get_var(ncid, varid, start_tod)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var start_tod', &
           file=__FILE__, line=__LINE__)
   end if

   status = nf90_inq_varid(ncid, 'curr_ymd', varid)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid curr_ymd', &
           file=__FILE__, line=__LINE__)
   end if
   status = nf90_get_var(ncid, varid, curr_ymd)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var curr_ymd', &
           file=__FILE__, line=__LINE__)
   end if

   status = nf90_inq_varid(ncid, 'curr_tod', varid)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid curr_tod', &
           file=__FILE__, line=__LINE__)
   end if
   status = nf90_get_var(ncid, varid, curr_tod)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var curr_tod', &
           file=__FILE__, line=__LINE__)
   end if

   status = nf90_close(ncid)
   if (status /= nf90_NoErr) then
      call shr_sys_abort(trim(subname)//' ERROR: nf90_close', &
           file=__FILE__, line=__LINE__)
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
