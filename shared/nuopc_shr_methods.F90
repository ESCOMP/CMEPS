module nuopc_shr_methods
  use ESMF         , only : operator(<), operator(/=), operator(+)
  use ESMF         , only : operator(-), operator(*) , operator(>=)
  use ESMF         , only : operator(<=), operator(>), operator(==), MOD
  use ESMF         , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError, ESMF_MAXSTR
  use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FAILURE
  use ESMF         , only : ESMF_State, ESMF_StateGet
  use ESMF         , only : ESMF_Field, ESMF_FieldGet
  use ESMF         , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSet
  use ESMF         , only : ESMF_GeomType_Flag, ESMF_FieldStatus_Flag
  use ESMF         , only : ESMF_Mesh, ESMF_MeshGet, ESMF_AlarmSet
  use ESMF         , only : ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_COMPLETE
  use ESMF         , only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockGet, ESMF_ClockSet
  use ESMF         , only : ESMF_ClockPrint, ESMF_ClockAdvance
  use ESMF         , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmGet, ESMF_AlarmSet
  use ESMF         , only : ESMF_Calendar, ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
  use ESMF         , only : ESMF_Time, ESMF_TimeGet, ESMF_TimeSet, ESMF_ClockGetAlarm
  use ESMF         , only : ESMF_TimeInterval, ESMF_TimeIntervalSet, ESMF_TimeIntervalGet
  use ESMF         , only : ESMF_VM, ESMF_VMGet, ESMF_VMBroadcast, ESMF_VMGetCurrent
  use ESMF         , only : ESMF_ClockGetNextTime
  use NUOPC        , only : NUOPC_CompAttributeGet
  use NUOPC_Model  , only : NUOPC_ModelGet
  use shr_kind_mod , only : r8 => shr_kind_r8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod  , only : shr_sys_abort
  use shr_log_mod  , only : shr_log_setLogUnit

  implicit none
  private

  public  :: memcheck
  public  :: get_component_instance
  public  :: set_component_logging
  public  :: log_clock_advance
  public  :: state_getscalar
  public  :: state_setscalar
  public  :: state_diagnose
  public  :: alarmInit
  public  :: get_minimum_timestep
  public  :: chkerr
  public  :: shr_get_rpointer_name
  private :: timeInit
  private :: field_getfldptr
  
  ! Module data
  
  ! Clock and alarm options shared with esm_time_mod along with dtime_driver which is initialized there.
  ! Dtime_driver is needed here for setting alarm options which use the nstep option and is a module variable
  ! to avoid requiring a change in all model caps. 
  character(len=*), public, parameter :: &
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
       optEnd            = "end"       , &
       optGLCCouplingPeriod = "glc_coupling_period"

  integer, public :: dtime_drv ! initialized in esm_time_mod.F90

  integer, parameter :: SecPerDay = 86400 ! Seconds per day
  integer, parameter :: memdebug_level=1
  character(len=1024)                   :: msgString
  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine memcheck(string, level, maintask)

    ! input/output variables
    character(len=*) , intent(in) :: string
    integer          , intent(in) :: level
    logical          , intent(in) :: maintask

    ! local variables
#ifdef CESMCOUPLED
    integer, external :: GPTLprint_memusage
    integer :: ierr = 0
#endif
    !-----------------------------------------------------------------------

#ifdef CESMCOUPLED
    if ((maintask .and. memdebug_level > level) .or. memdebug_level > level+1) then
       ierr = GPTLprint_memusage(string)
    endif
#endif

  end subroutine memcheck

!===============================================================================

  subroutine get_component_instance(gcomp, inst_suffix, inst_index, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(out) :: inst_suffix
    integer          , intent(out) :: inst_index
    integer          , intent(out) :: rc

    ! local variables
    logical          :: isPresent
    character(len=4) :: cvalue
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       cvalue = inst_suffix(2:)
       read(cvalue, *) inst_index
    else
       inst_suffix = ""
       inst_index=1
    endif

  end subroutine get_component_instance

!===============================================================================
  subroutine set_component_logging(gcomp, maintask, logunit, shrlogunit, rc)
    use NUOPC, only: NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd  
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    logical, intent(in)  :: maintask
    integer, intent(out) :: logunit
    integer, intent(out) :: shrlogunit
    integer, intent(out) :: rc

    ! local variables
    character(len=CL) :: diro
    character(len=CL) :: logfile
    character(len=CL) :: inst_suffix
    integer :: inst_index ! Not used here
    integer :: n
    character(len=CL) :: name
    character(len=*), parameter :: subname = "("//__FILE__//": set_component_logging)"   
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (maintask) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call get_component_instance(gcomp, inst_suffix, inst_index, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! Multiinstance logfile name needs a correction
       if(len_trim(inst_suffix) > 0) then
          n = index(logfile, '.')
          logfile = logfile(1:n-1)//trim(inst_suffix)//logfile(n:)
       endif

       open(newunit=logunit,file=trim(diro)//"/"//trim(logfile))

    else
       logUnit = 6
    endif
    
    call ESMF_GridCompGet(gcomp, name=name, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": setting logunit for component: "//trim(name), ESMF_LOGMSG_INFO)
    call NUOPC_CompAttributeAdd(gcomp, (/"logunit"/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, "logunit", logunit, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call shr_log_setLogUnit (logunit)
    ! Still need to set this return value
    shrlogunit = logunit
    call ESMF_LogWrite(trim(subname)//": done for component "//trim(name), ESMF_LOGMSG_INFO)
  end subroutine set_component_logging

!===============================================================================

  subroutine log_clock_advance(clock, component, logunit, rc)

    ! input/output variables
    type(ESMF_Clock)               :: clock
    character(len=*) , intent(in)  :: component
    integer          , intent(in)  :: logunit
    integer          , intent(out) :: rc

    ! local variables
    character(len=CL) :: cvalue, prestring
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    write(prestring, *) "------>Advancing ",trim(component)," from: "
    call ESMF_ClockPrint(clock, options="currTime", unit=cvalue, preString=trim(prestring), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

    call ESMF_ClockPrint(clock, options="stopTime", unit=cvalue, &
         preString="--------------------------------> to: ", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

  end subroutine log_clock_advance

!===============================================================================

  subroutine state_getscalar(state, scalar_id, scalar_value, flds_scalar_name, flds_scalar_num, rc)

    ! ----------------------------------------------
    ! Get scalar data from State for a particular name and broadcast it to all other pets
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)     :: state
    integer,          intent(in)     :: scalar_id
    real(r8),         intent(out)    :: scalar_value
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_VM)     :: vm
    type(ESMF_Field)  :: field
    real(r8), pointer :: farrayptr(:,:)
    real(r8)          :: tmp(1)
    character(len=*), parameter :: subname='(state_getscalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
        call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u)
        rc = ESMF_FAILURE
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
      endif
      tmp(:) = farrayptr(scalar_id,:)
    endif
    call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    scalar_value = tmp(1)

  end subroutine state_getscalar

!================================================================================

  subroutine state_setscalar(scalar_value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)

    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------

    ! input/output arguments
    real(r8),         intent(in)     :: scalar_value
    integer,          intent(in)     :: scalar_id
    type(ESMF_State), intent(inout)  :: State
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_Field)  :: lfield
    type(ESMF_VM)     :: vm
    real(r8), pointer :: farrayptr(:,:)
    character(len=*), parameter :: subname='(state_setscalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
       call ESMF_FieldGet(lfield, farrayPtr = farrayptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       endif
       farrayptr(scalar_id,1) = scalar_value
    endif

  end subroutine state_setscalar

!===============================================================================

  subroutine state_diagnose(State, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    integer                         :: n
    type(ESMf_Field)                :: lfield
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(r8), pointer               :: dataPtr1d(:)
    real(r8), pointer               :: dataPtr2d(:,:)
    character(len=*),parameter      :: subname='(state_diagnose)'
    ! ----------------------------------------------

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    call ESMF_StateGet(state, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount

       call ESMF_StateGet(state, itemName=lfieldnamelist(n), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call field_getfldptr(lfield, fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data
       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,a)') trim(string)//': for 1d field '//trim(lfieldnamelist(n))
             call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
             write(msgString,'(A,3g14.7,i8)') trim(string)//': 1d field '//trim(lfieldnamelist(n)), &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
             call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
             call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
          endif
       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,a)') trim(string)//': for 2d field '//trim(lfieldnamelist(n))
             call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
             write(msgString,'(A,3g14.7,i8)') trim(string)//': 2d field '//trim(lfieldnamelist(n)), &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
             call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
             call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
          endif
       else
          call shr_sys_abort(trim(subname)//": ERROR rank not supported ")
       endif
    enddo

    deallocate(lfieldnamelist)

  end subroutine state_diagnose

!===============================================================================

  subroutine field_getfldptr(field, fldptr1, fldptr2, rank, abort, rc)

    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_Field)  , intent(in)              :: field
    real(r8), pointer , intent(inout), optional :: fldptr1(:)
    real(r8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort
    integer           , intent(out)  , optional :: rc

    ! local variables
    type(ESMF_GeomType_Flag)    :: geomtype
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Mesh)             :: lmesh
    integer                     :: lrank, nnodes, nelements
    logical                     :: labort
    character(len=*), parameter :: subname='(field_getfldptr)'
    ! ----------------------------------------------

    if (.not.present(rc)) then
       call shr_sys_abort(trim(subname)//": ERROR rc not present ", &
            line=__LINE__, file=u_FILE_u)
    endif

    rc = ESMF_SUCCESS

    labort = .true.
    if (present(abort)) then
       labort = abort
    endif
    lrank = -99

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       lrank = 0
       if (labort) then
          call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       else
          call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       endif
    else

       call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (geomtype == ESMF_GEOMTYPE_GRID) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (geomtype == ESMF_GEOMTYPE_MESH) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (nnodes == 0 .and. nelements == 0) lrank = 0
       else
          call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", &
               ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       endif ! geomtype

       if (lrank == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
               ESMF_LOGMSG_INFO)
       elseif (lrank == 1) then
          if (.not.present(fldptr1)) then
             call shr_sys_abort(trim(subname)//": ERROR missing rank=1 array ", &
                  line=__LINE__, file=u_FILE_u)
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (lrank == 2) then
          if (.not.present(fldptr2)) then
             call shr_sys_abort(trim(subname)//": ERROR missing rank=2 array ", &
                  line=__LINE__, file=u_FILE_u)
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call shr_sys_abort(trim(subname)//": ERROR in rank ", &
               line=__LINE__, file=u_FILE_u)
       endif

    endif  ! status

    if (present(rank)) then
       rank = lrank
    endif

  end subroutine field_getfldptr

!===============================================================================

  subroutine alarmInit( clock, alarm, option, &
       opt_n, opt_ymd, opt_tod, RefTime, alarmname, advance_clock, rc)
    use ESMF, only : ESMF_AlarmPrint, ESMF_ClockGetAlarm
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
    logical          , optional , intent(in)    :: advance_clock ! advance clock to trigger alarm
    integer                     , intent(inout) :: rc        ! Return code

    ! local variables
    type(ESMF_Calendar)     :: cal                ! calendar
    integer                 :: lymd             ! local ymd
    integer                 :: ltod             ! local tod
    integer                 :: cyy,cmm,cdd,csec ! time info
    character(len=64)       :: lalarmname       ! local alarm name
    logical                 :: update_nextalarm ! update next alarm
    type(ESMF_Time)         :: CurrTime         ! Current Time
    type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
    type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
    type(ESMF_TimeInterval) :: TimeStepInterval ! Component timestep interval
    character(len=*), parameter :: subname = '(alarmInit): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lalarmname = 'alarm_unknown'
    if (present(alarmname)) lalarmname = trim(alarmname)
    ltod = 0
    if (present(opt_tod)) ltod = opt_tod
    lymd = -1
    if (present(opt_ymd)) lymd = opt_ymd

    call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! initial guess of next alarm, this will be updated below
    if (present(RefTime)) then
       NextAlarm = RefTime
    else
       NextAlarm = CurrTime
    endif
    
    ! Determine calendar
    call ESMF_ClockGet(clock, calendar=cal)

    ! Error checks
    if (trim(option) == optdate) then
       if (.not. present(opt_ymd)) then
          call shr_sys_abort(trim(subname)//trim(option)//' requires opt_ymd')
       end if
       if (lymd < 0 .or. ltod < 0) then
          call shr_sys_abort(subname//trim(option)//'opt_ymd, opt_tod invalid')
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
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
    end if
    call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine inputs for call to create alarm
    selectcase (trim(option))

    case (optNONE)
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optNever)
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optDate)
       call timeInit(NextAlarm, lymd, cal, ltod, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

   case (optNSteps,trim(optNSteps)//'s')
      call ESMF_ClockGet(clock, TimeStep=TimestepInterval, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if(dtime_drv > 0) then
         call ESMF_TimeIntervalSet(AlarmInterval, s=dtime_drv, rc=rc )
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      else
         call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      
      AlarmInterval = AlarmInterval * opt_n
      ! timestepinterval*0 is 0 of kind ESMF_TimeStepInterval
      if (mod(AlarmInterval, TimestepInterval) /= (TimestepInterval*0)) then
         call shr_sys_abort(subname//'illegal Alarm setting for '//trim(alarmname))
      endif
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

    case (optEnd)
       call ESMF_ClockGetAlarm(clock, alarmname="alarm_stop", alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmGet(alarm, ringTime=NextAlarm, rc=rc) 
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
    case default
       call shr_sys_abort(subname//'unknown option '//trim(option))

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
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    
    ! Advance model clock to trigger alarm then reset model clock back to currtime
    if (present(advance_clock)) then
       if (advance_clock) then
          call ESMF_AlarmSet(alarm, clock=clock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockGet(clock, currTime=CurrTime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockAdvance(clock,rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockSet(clock, currTime=currtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

  end subroutine alarmInit

!===============================================================================

  subroutine timeInit( Time, ymd, cal, tod, rc)

    ! Create the ESMF_Time object corresponding to the given input time,
    ! given in YMD (Year Month Day) and TOD (Time-of-day) format.
    ! Set the time by an integer as YYYYMMDD and integer seconds in the day

    ! input/output parameters:
    type(ESMF_Time)     , intent(inout) :: Time ! ESMF time
    integer             , intent(in)    :: ymd  ! year, month, day YYYYMMDD
    type(ESMF_Calendar) , intent(in)    :: cal  ! ESMF calendar
    integer             , intent(in)    :: tod  ! time of day in seconds
    integer             , intent(out)   :: rc

    ! local variables
    integer :: year, mon, day ! year, month, day as integers
    integer :: tdate          ! temporary date
    character(len=*), parameter :: subname='(timeInit)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if ( (ymd < 0) .or. (tod < 0) .or. (tod > SecPerDay) )then
       call shr_sys_abort( subname//'ERROR yymmdd is a negative number or time-of-day out of bounds' )
    end if

    tdate = abs(ymd)
    year = int(tdate/10000)
    if (ymd < 0) year = -year
    mon = int( mod(tdate,10000)/  100)
    day = mod(tdate,  100)

    call ESMF_TimeSet( Time, yy=year, mm=mon, dd=day, s=tod, calendar=cal, rc=rc )
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine timeInit

!===============================================================================

    integer function get_minimum_timestep(gcomp, rc)
    ! Get the minimum timestep interval in seconds based on the nuopc.config variables *_cpl_dt,
    ! if none of these variables are defined this routine will throw an error
    type(ESMF_GridComp), intent(in) :: gcomp
    integer, intent(out)            :: rc

    character(len=CS)          :: cvalue
    integer                    :: comp_dt             ! coupling interval of component
    integer, parameter         :: ncomps = 8
    character(len=3),dimension(ncomps) :: compname
    character(len=10)          :: comp
    logical                    :: is_present, is_set  ! determine if these variables are used
    integer                    :: i
    !---------------------------------------------------------------------------
    ! Determine driver clock timestep
    !---------------------------------------------------------------------------
    compname = (/"atm", "lnd", "ice", "ocn", "rof", "glc", "wav", "esp"/)
    get_minimum_timestep = huge(1)

    do i=1,ncomps
       comp = compname(i)//"_cpl_dt"
    
       call NUOPC_CompAttributeGet(gcomp, name=comp, isPresent=is_present, isSet=is_set, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (is_present .and. is_set) then
          call NUOPC_CompAttributeGet(gcomp, name=comp, value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) comp_dt
          get_minimum_timestep = min(comp_dt, get_minimum_timestep)
       endif
    enddo

    if(get_minimum_timestep == huge(1)) then
       call shr_sys_abort('minimum_timestep_error: this option is not supported ')
    endif
    if(get_minimum_timestep <= 0) then
       call shr_sys_abort('minimum_timestep_error ERROR ')
    endif
  end function get_minimum_timestep

  subroutine shr_get_rpointer_name(gcomp, compname, ymd, time, rpfile, mode, rc)
    type(ESMF_GRIDCOMP), intent(in) :: gcomp
    character(len=3), intent(in) :: compname
    integer, intent(in) :: ymd
    integer, intent(in) :: time
    character(len=*), intent(out) :: rpfile
    character(len=*), intent(in) :: mode
    integer, intent(out) :: rc
    
    ! local vars
    integer :: yr, mon, day
    character(len=16) timestr
    logical :: isPresent
    character(len=ESMF_MAXSTR)  :: inst_suffix
    character(len=*), parameter :: subname='shr_get_rpointer_name'
    
    rc = ESMF_SUCCESS

    inst_suffix = ""
    call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(ispresent) call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', value=inst_suffix, rc=rc)
    
    yr = ymd/10000
    mon = (ymd - yr*10000)/100
    day = (ymd - yr*10000 - mon*100)
    write(timestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',time
    write(rpfile,*) "rpointer."//compname//trim(inst_suffix)//'.'//trim(timestr)
    rpfile = adjustl(rpfile)
    if (mode.eq.'read') then
       inquire(file=trim(rpfile), exist=isPresent)
       if(.not. isPresent) then
          rpfile = "rpointer."//compname//trim(inst_suffix)
          inquire(file=trim(rpfile), exist=isPresent)
          if(.not. isPresent) then
             call shr_sys_abort( subname//'ERROR no rpointer file found in '//rpfile//' or in '//rpfile//'.'//timestr )
          endif
       endif
    endif
  end subroutine shr_get_rpointer_name

  logical function chkerr(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

end module nuopc_shr_methods
