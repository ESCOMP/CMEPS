module med_phases_history_mod

  !-----------------------------------------------------------------------------
  ! Mediator History control
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_VM
  use ESMF                  , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockSet, ESMF_ClockAdvance
  use ESMF                  , only : ESMF_ClockGetNextTime, ESMF_ClockGetAlarm, ESMF_ClockIsCreated
  use ESMF                  , only : ESMF_Calendar, ESMF_Time, ESMF_TimeGet
  use ESMF                  , only : ESMF_TimeInterval, ESMF_TimeIntervalGet, ESMF_TimeIntervalSet
  use ESMF                  , only : ESMF_Alarm, ESMF_AlarmIsRinging, ESMF_AlarmRingerOff, ESMF_AlarmGet
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_MAXSTR, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
  use ESMF                  , only : ESMF_Finalize
  use ESMF                  , only : operator(-), operator(+)
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use NUOPC_Model           , only : NUOPC_ModelGet
  use med_utils_mod         , only : chkerr => med_utils_ChkErr
  use med_internalstate_mod , only : ncomps, compname
  use med_internalstate_mod , only : InternalState, maintask, logunit
  use med_time_mod          , only : med_time_alarmInit
  use med_io_mod            , only : med_io_write, med_io_wopen, med_io_enddef, med_io_close
  use perf_mod              , only : t_startf, t_stopf
  use pio                   , only : file_desc_t

  implicit none
  private

  ! Public routine called from med_internal_state_init
  public :: med_phases_history_init

  ! Public routine called from the run sequence
  public :: med_phases_history_write         ! inst only - for all variables

  ! Public routines called from post phases
  public :: med_phases_history_write_comp    ! inst, avg, aux for component
  public :: med_phases_history_write_med     ! inst only, med aoflux and ocn albedoes
  public :: med_phases_history_write_lnd2glc ! inst only, yearly average of lnd->glc data on lnd grid

  ! Private routines
  private :: med_phases_history_write_comp_inst  ! write instantaneous file for a given component
  private :: med_phases_history_write_comp_avg   ! write averaged file for a given component
  private :: med_phases_history_write_comp_aux   ! write auxiliary file for a given component
  private :: med_phases_history_init_histclock
  private :: med_phases_history_query_ifwrite
  private :: med_phases_history_set_timeinfo
  private :: med_phases_history_fldbun_accum
  private :: med_phases_history_fldbun_average

  ! ----------------------------
  ! Instantaneous history files all components
  ! ----------------------------
  character(CL)  :: hist_option_all_inst  ! freq_option setting (ndays, nsteps, etc)
  integer        :: hist_n_all_inst       ! freq_n setting relative to freq_option

  ! ----------------------------
  ! Instantaneous history files datatypes/variables per component
  ! ----------------------------
  type, public :: instfile_type
     type(file_desc_t):: io_file
     logical          :: write_inst
     character(CS)    :: hist_option
     integer          :: hist_n
     type(ESMF_Clock) :: clock
     type(ESMF_Alarm) :: alarm
     character(CS)    :: alarmname
     logical          :: is_clockset = .false.
     logical          :: is_active = .false.
  end type instfile_type
  type(instfile_type) , allocatable, public :: instfiles(:)

  ! ----------------------------
  ! Time averaging history files
  ! ----------------------------
  type, public :: avgfile_type
     type(file_desc_t)      :: io_file
     logical                :: write_avg
     type(ESMF_FieldBundle) :: FBaccum_import    ! field bundle for time averaging
     integer                :: accumcnt_import   ! field bundle accumulation counter
     type(ESMF_FieldBundle) :: FBaccum_export    ! field bundle for time averaging
     integer                :: accumcnt_export   ! field bundle accumulation counter
     character(CS)          :: hist_option
     integer                :: hist_n
     type(ESMF_Clock)       :: clock
     type(ESMF_Alarm)       :: alarm
     character(CS)          :: alarmname
     logical                :: is_clockset = .false.
     logical                :: is_active = .false.
  end type avgfile_type
  type(avgfile_type), allocatable :: avgfiles(:)

  ! ----------------------------
  ! Auxiliary history files
  ! ----------------------------
  type, public :: auxfile_type
     type(file_desc_t)          :: io_file
     character(CS), allocatable :: flds(:)       ! array of aux field names
     character(CS)              :: auxname       ! name for history file creation
     character(CL)              :: histfile = '' ! current history file name
     integer                    :: ntperfile     ! maximum number of time samples per file
     integer                    :: nt = 0        ! time in file
     logical                    :: doavg         ! if true, time average, otherwise instantaneous
     type(ESMF_FieldBundle)     :: FBaccum       ! field bundle for time averaging
     integer                    :: accumcnt      ! field bundle accumulation counter
     type(ESMF_Clock)           :: clock         ! auxiliary history clock
     type(ESMF_Alarm)           :: alarm         ! auxfile alarm
     character(CS)              :: alarmname     ! name of write alarm
  end type auxfile_type

  integer, parameter :: max_auxfiles = 10
  type, public :: auxcomp_type
     type(auxfile_type) :: files(max_auxfiles)
     integer            :: num_auxfiles  = 0       ! actual number of auxiliary files
     logical            :: init_auxfiles = .false. ! if auxfile initial has occured
  end type auxcomp_type
  type(auxcomp_type), allocatable, public :: auxcomp(:)

  ! ----------------------------
  ! Other private module variables
  ! ----------------------------

  logical :: whead(2) = (/.true. , .false./)
  logical :: wdata(2) = (/.false., .true. /)

  character(CL) :: case_name = 'unset'  ! case name
  character(CS) :: inst_tag = 'unset'   ! instance tag
  logical       :: debug_alarms = .true.
  character(*), parameter :: u_FILE_u  = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_history_init()
    ! allocate module memory
    allocate(instfiles(ncomps))
    allocate(avgfiles(ncomps))
    allocate(auxcomp(ncomps))
  end subroutine med_phases_history_init

  !===============================================================================
  subroutine med_phases_history_write(gcomp, rc)

    ! --------------------------------------
    ! Write instantaneous mediator history file for all variables
    ! --------------------------------------

    use med_io_mod, only : med_io_write_time, med_io_define_time
    use ESMF      , only : ESMF_Alarm, ESMF_AlarmSet
    use ESMF      , only : ESMF_FieldBundleIsCreated
    use med_internalstate_mod, only : compocn, compatm

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(file_desc_t)       :: io_file
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: mclock
    type(ESMF_Alarm)        :: alarm
    character(CS)           :: alarmname
    character(CL)           :: cvalue       ! attribute string
    logical                 :: isPresent
    logical                 :: isSet
    type(ESMF_VM)           :: vm
    type(ESMF_Calendar)     :: calendar     ! calendar type
    integer                 :: m,n        ! indices
    character(CL)           :: time_units   ! units of time variable
    character(CL)           :: hist_file    ! history file name
    real(r8)                :: time_val     ! time coordinate output
    real(r8)                :: time_bnds(2) ! time bounds output
    logical                 :: write_now    ! true => write to history type
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: nexttime
    character(len=CS)       :: currtimestr
    character(len=CS)       :: nexttimestr
    integer                 :: yr,mon,day,sec    ! time units
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: ringInterval_length
    logical                 :: first_time = .true.
    character(len=*), parameter :: subname='(med_phases_history_write)'
    !---------------------------------------

    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    alarmname='alarm_history_inst_all'

    if (first_time) then
       call NUOPC_CompAttributeGet(gcomp, name='history_option', isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          call NUOPC_CompAttributeGet(gcomp, name='history_option', value=hist_option_all_inst, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeGet(gcomp, name='history_n', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) hist_n_all_inst
       else
          ! If attribute is not present - don't write history output
          hist_option_all_inst = 'none'
          hist_n_all_inst = -999
       end if

       ! Set alarm name and initialize clock and alarm for instantaneous history output
       ! The alarm for the full history write is set on the mediator clock not as a separate alarm
       if (hist_option_all_inst /= 'none' .and. hist_option_all_inst /= 'never') then

          ! Initialize alarm on mediator clock for instantaneous mediator history output for all variables
          call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockGet(mclock, startTime=starttime,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_time_alarmInit(mclock, alarm, option=hist_option_all_inst, opt_n=hist_n_all_inst, &
               reftime=starttime, alarmname=alarmname, rc=rc)
          call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Advance model clock to trigger alarms then reset model clock back to currtime
          call ESMF_ClockGet(mclock, currTime=CurrTime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockAdvance(mclock,rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockSet(mclock, currTime=currtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Write diagnostic info
          if (maintask) then
             write(logunit,'(a,2x,i8)') trim(subname) // "  initialized history alarm "//&
                  trim(alarmname)//"  with option "//trim(hist_option_all_inst)//" and frequency ",hist_n_all_inst
          end if
       end if
       first_time = .false.
    end if

    write_now = .false.
    if (hist_option_all_inst /= 'none' .and. hist_option_all_inst /= 'never') then
       call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGetAlarm(mclock, alarmname=trim(alarmname), alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          ! Set write flag to .true. and turn ringer off
          write_now = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Write diagnostic info if appropriate
          if (maintask .and. debug_alarms) then
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
             write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec

             if (maintask) then
                write(logunit,*)
                write(logunit,'(a,i8)') trim(subname)//" : history alarmname "//trim(alarmname)//&
                     ' is ringing, interval length is ', ringInterval_length
                write(logunit,'(a)') trim(subname)//" : mclock currtime = "//trim(currtimestr)//&
                     " mclock nexttime = "//trim(nexttimestr)
             end if
          end if
       end if

       ! If write now flag is true
       if (write_now) then

          ! Determine time_val and tbnds data for history as well as history file name
          call med_phases_history_set_timeinfo(gcomp, mclock, alarmname, &
               time_val, time_bnds, time_units, hist_file, doavg=.false., compname='all', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Create history file
          call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_wopen(hist_file, io_file, vm, rc, clobber=.true.)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Loop over whead/wdata phases
          do m = 1,2
             if (m == 2) then
                call med_io_enddef(io_file)
             end if

             ! Write time values
             if (whead(m)) then
                call ESMF_ClockGet(mclock, calendar=calendar, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call med_io_define_time(io_file, time_units, calendar, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                call med_io_write_time(io_file, time_val, time_bnds, nt=1, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             do n = 2,ncomps ! skip the mediator here
                ! Write import and export field bundles
                if (is_local%wrap%comp_present(n)) then
                   if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
                      call med_io_write(io_file, is_local%wrap%FBimp(n,n), whead(m), wdata(m), &
                           is_local%wrap%nx(n), is_local%wrap%ny(n), nt=1, pre=trim(compname(n))//'Imp', &
                           ntile=is_local%wrap%ntile(n), rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   endif
                   if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(n),rc=rc)) then
                      call med_io_write(io_file, is_local%wrap%FBexp(n), whead(m), wdata(m), &
                           is_local%wrap%nx(n), is_local%wrap%ny(n), nt=1, pre=trim(compname(n))//'Exp', &
                           ntile=is_local%wrap%ntile(n), rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   endif
                end if
                ! Write mediator fraction field bundles
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBFrac(n),rc=rc)) then
                   call med_io_write(io_file, is_local%wrap%FBFrac(n), whead(m), wdata(m), &
                        is_local%wrap%nx(n), is_local%wrap%ny(n), nt=1, pre='Med_frac_'//trim(compname(n)), &
                        ntile=is_local%wrap%ntile(n), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
                ! Write component mediator area field bundles
                call med_io_write(io_file, is_local%wrap%FBArea(n), whead(m), wdata(m), &
                     is_local%wrap%nx(n), is_local%wrap%ny(n), nt=1, pre='MED_'//trim(compname(n)), &
                     ntile=is_local%wrap%ntile(n), rc=rc)
             end do

             ! Write atm/ocn fluxes and ocean albedoes if field bundles are created
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
                call med_io_write(io_file, is_local%wrap%FBMed_ocnalb_o, whead(m), wdata(m), &
                     is_local%wrap%nx(compocn), is_local%wrap%ny(compocn), nt=1, pre='Med_alb_ocn', rc=rc)
             end if
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o,rc=rc)) then
                call med_io_write(io_file, is_local%wrap%FBMed_aoflux_o, whead(m), wdata(m), &
                     is_local%wrap%nx(compocn), is_local%wrap%ny(compocn), nt=1, pre='Med_aoflux_ocn', rc=rc)
             end if
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_a,rc=rc)) then
                call med_io_write(io_file, is_local%wrap%FBMed_ocnalb_a, whead(m), wdata(m), &
                     is_local%wrap%nx(compatm), is_local%wrap%ny(compatm), nt=1, pre='Med_alb_atm', rc=rc)
             end if
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a,rc=rc)) then
                call med_io_write(io_file, is_local%wrap%FBMed_aoflux_a, whead(m), wdata(m), &
                     is_local%wrap%nx(compatm), is_local%wrap%ny(compatm), nt=1, pre='Med_aoflux_atm', rc=rc)
             end if

          end do ! end of loop over whead/wdata m index phases

          ! Close file
          call med_io_close(io_file, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       end if ! end of write_now if-block
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write

  !===============================================================================
  subroutine med_phases_history_write_med(gcomp, rc)

    ! Write mediator history file for med variables - only instantaneous files are written
    ! This writes out ocean albedoes and atm/ocean fluxes computed by the mediator
    ! along with the fractions computed by the mediator

    use ESMF      , only : ESMF_FieldBundleIsCreated
    use med_io_mod, only : med_io_write_time, med_io_define_time
    use med_internalstate_mod, only : compmed, compocn, compatm

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_VM)       :: vm
    type(ESMF_Calendar) :: calendar     ! calendar type
    integer             :: m            ! indices
    character(CL)       :: time_units   ! units of time variable
    character(CL)       :: hist_file    ! history file name
    real(r8)            :: time_val     ! time coordinate output
    real(r8)            :: time_bnds(2) ! time bounds output
    logical             :: write_now    ! true => write to history type
    character(CL)       :: cvalue       ! attribute string
    character(CL)       :: hist_option  ! freq_option setting (ndays, nsteps, etc)
    integer             :: hist_n       ! freq_n setting relative to freq_option
    character(CL)       :: hist_option_in
    character(CL)       :: hist_n_in
    logical             :: isPresent
    logical             :: isSet
    character(len=*), parameter :: subname='(med_phases_history_write_med)'
    !---------------------------------------
    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! alarm is not set determine hist_option and hist_n
    if (.not. instfiles(compmed)%is_clockset) then
       ! Determine attribute prefix
       write(hist_option_in,'(a)') 'history_option_med_inst'
       write(hist_n_in,'(a)') 'history_n_med_inst'

       ! Determine instantaneous mediator output frequency and type
       call NUOPC_CompAttributeGet(gcomp, name=trim(hist_option_in), isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          call NUOPC_CompAttributeGet(gcomp, name=trim(hist_option_in), value=hist_option, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeGet(gcomp, name=trim(hist_n_in), value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) hist_n
       else
          ! If attribute is not present - don't write history output
          hist_option = 'none'
          hist_n = -999
       end if

       ! Set alarm name and initialize clock and alarm for instantaneous history output
       if (hist_option /= 'none' .and. hist_option /= 'never') then
          instfiles(compmed)%alarmname =  'alarm_history_inst_med'
          call med_phases_history_init_histclock(gcomp, instfiles(compmed)%clock, &
               instfiles(compmed)%alarm, instfiles(compmed)%alarmname, hist_option, hist_n, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          instfiles(compmed)%is_active = .true.
          instfiles(compmed)%is_clockset = .true.
       else
          instfiles(compmed)%is_active = .false.
          ! this is set to true here even if history file is not active
          instfiles(compmed)%is_clockset = .true.
       end if
    end if

    ! if history file is active and history clock is initialized - process history file
    if (instfiles(compmed)%is_active .and. instfiles(compmed)%is_clockset) then

       ! Determine if will write to history file
       call med_phases_history_query_ifwrite(gcomp, instfiles(compmed)%clock, instfiles(compmed)%alarmname, &
            write_now, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! If write now flag is true
       if (write_now) then

          ! Determine time_val and tbnds data for history as well as history file name
          call med_phases_history_set_timeinfo(gcomp, instfiles(compmed)%clock, instfiles(compmed)%alarmname, &
               time_val, time_bnds, time_units, hist_file, doavg=.false., compname='med', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Create history file
          call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_wopen(hist_file, instfiles(compmed)%io_file, vm, rc, clobber=.true.)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do m = 1,2
             ! Write time values
             if (whead(m)) then
                call ESMF_ClockGet(instfiles(compmed)%clock, calendar=calendar, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call med_io_define_time(instfiles(compmed)%io_file, time_units, calendar, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                call med_io_enddef(instfiles(compmed)%io_file)
                call med_io_write_time(instfiles(compmed)%io_file, time_val, time_bnds, nt=1, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             ! Write aoflux fields computed in mediator
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o,rc=rc)) then
                call med_io_write(instfiles(compmed)%io_file, is_local%wrap%FBMed_aoflux_o, whead(m), wdata(m), &
                     is_local%wrap%nx(compocn), is_local%wrap%ny(compocn), nt=1, pre='Med_aoflux_ocn', rc=rc)
             end if
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a,rc=rc)) then
                call med_io_write(instfiles(compmed)%io_file, is_local%wrap%FBMed_aoflux_a, whead(m), wdata(m), &
                     is_local%wrap%nx(compatm), is_local%wrap%ny(compatm), nt=1, pre='Med_aoflux_atm', rc=rc)
             end if

             ! If appropriate - write ocn albedos computed in mediator
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
                call med_io_write(instfiles(compmed)%io_file, is_local%wrap%FBMed_ocnalb_o, whead(m), wdata(m), &
                     is_local%wrap%nx(compocn), is_local%wrap%ny(compocn), nt=1, pre='Med_alb_ocn', rc=rc)
             end if
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_a,rc=rc)) then
                call med_io_write(instfiles(compmed)%io_file, is_local%wrap%FBMed_ocnalb_a, whead(m), wdata(m), &
                     is_local%wrap%nx(compatm), is_local%wrap%ny(compatm), nt=1, pre='Med_alb_atm', rc=rc)
             end if
          end do ! end of loop over m

          ! Close file
          call med_io_close(instfiles(compmed)%io_file, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       end if ! end of if-write_now block
    end if  ! end of if-active block

  end subroutine med_phases_history_write_med

  !===============================================================================
  subroutine med_phases_history_write_lnd2glc(gcomp, fldbun_lnd, rc, fldbun_glc)

    ! Write yearly average of lnd -> glc fields on both land and glc grids

    use med_internalstate_mod, only : complnd, compglc
    use med_constants_mod , only : SecPerDay => med_constants_SecPerDay
    use med_io_mod        , only : med_io_write_time, med_io_define_time
    use med_io_mod        , only : med_io_date2yyyymmdd, med_io_sec2hms, med_io_ymd2date

    ! input/output variables
    type(ESMF_GridComp)    , intent(in)  :: gcomp
    type(ESMF_FieldBundle) , intent(in)  :: fldbun_lnd
    integer                , intent(out) :: rc
    type(ESMF_FieldBundle) , intent(in), optional :: fldbun_glc(:)

    ! local variables
    type(file_desc_t)       :: io_file
    type(InternalState)     :: is_local
    type(ESMF_VM)           :: vm
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: nexttime
    type(ESMF_Calendar)     :: calendar     ! calendar type
    type(ESMF_TimeInterval) :: timediff(2)  ! time bounds upper and lower relative to start
    character(len=CS)       :: nexttime_str
    integer                 :: yr,mon,day,sec
    integer                 :: start_ymd    ! starting date YYYYMMDD
    character(CL)           :: time_units   ! units of time variable
    real(r8)                :: time_val     ! time coordinate output
    real(r8)                :: time_bnds(2) ! time bounds output
    character(len=CL)       :: hist_file
    integer                 :: m,n
    logical                 :: isPresent
    character(len=*), parameter :: subname='(med_phases_history_write_lnd2glc)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the model clock
    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine starttime, currtime and nexttime
    call ESMF_ClockGet(clock, currtime=currtime, starttime=starttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine time units
    call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_io_ymd2date(yr,mon,day,start_ymd)
    time_units = 'days since ' // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(sec, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set time bounds and time coord
    timediff(1) = nexttime - starttime
    call ESMF_TimeIntervalGet(timediff(1), d=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    time_val = day + sec/real(SecPerDay,R8)
    time_bnds(1) = time_val
    time_bnds(2) = time_val

    ! Determine history file name
    if (trim(case_name) == 'unset') then
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
    call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(nexttime_str,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    write(hist_file, "(6a)") trim(case_name),'.cpl',trim(inst_tag),'.hx.1yr2glc.',trim(nexttime_str),'.nc'

    ! Create history file
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_io_wopen(hist_file, io_file, vm, rc, clobber=.true.)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Write data to history file
    do m = 1,2
       if (whead(m)) then
          call ESMF_ClockGet(clock, calendar=calendar, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_define_time(io_file, time_units, calendar, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call med_io_enddef(io_file)
          call med_io_write_time(io_file, time_val, time_bnds, nt=1, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       call med_io_write(io_file, fldbun_lnd, whead(m), wdata(m), &
            is_local%wrap%nx(complnd), is_local%wrap%ny(complnd), &
            nt=1, pre=trim(compname(complnd))//'Imp', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (present(fldbun_glc)) then
          do n = 1,size(fldbun_glc)
             call med_io_write(io_file, fldbun_glc(n), whead(m), wdata(m), &
                  is_local%wrap%nx(compglc(n)), is_local%wrap%ny(compglc(n)), &
                  nt=1, pre=trim(compname(compglc(n)))//'Exp', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end do
       end if

    end do ! end of loop over m

    ! Close history file
    call med_io_close(io_file, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_history_write_lnd2glc

  !===============================================================================
  subroutine med_phases_history_write_comp(gcomp, compid, rc)

    ! Write mediator history file for compid variables

    ! input/output variables
    type(ESMF_GridComp), intent(inout) :: gcomp
    integer            , intent(in)    :: compid
    integer            , intent(out)   :: rc
    !---------------------------------------
    rc = ESMF_SUCCESS

    call med_phases_history_write_comp_inst(gcomp, compid, instfiles(compid), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_comp_avg(gcomp, compid, avgfiles(compid), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_comp_aux(gcomp, compid, auxcomp(compid), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_history_write_comp

  !===============================================================================
  subroutine med_phases_history_write_comp_inst(gcomp, compid, instfile, rc)

    ! Write instantaneous mediator history file for component compid

    use med_io_mod, only : med_io_write_time, med_io_define_time
    use ESMF      , only : ESMF_FieldBundleIsCreated

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(in)    :: compid
    type(instfile_type) , intent(inout) :: instfile
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState) :: is_local
    character(CL)       :: cvalue       ! attribute string
    character(CL)       :: hist_option  ! freq_option setting (ndays, nsteps, etc)
    integer             :: hist_n       ! freq_n setting relative to freq_option
    character(CL)       :: hist_option_in
    character(CL)       :: hist_n_in
    logical             :: isPresent
    logical             :: isSet
    type(ESMF_VM)       :: vm
    type(ESMF_Calendar) :: calendar     ! calendar type
    integer             :: m            ! indices
    integer             :: nx,ny        ! global grid size
    integer             :: ntile        ! number of tiles for tiled domain eg CSG
    character(CL)       :: time_units   ! units of time variable
    character(CL)       :: hist_file    ! history file name
    real(r8)            :: time_val     ! time coordinate output
    real(r8)            :: time_bnds(2) ! time bounds output
    logical             :: write_now    ! true => write to history type
    character(len=*), parameter :: subname='(med_phases_history_write_inst_comp)'
    !---------------------------------------

    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! alarm is not set determine hist_option and hist_n
    if (.not. instfile%is_clockset) then

       ! Determine attribute name
       write(hist_option_in,'(a)') 'history_option_'//trim(compname(compid))//'_inst'
       write(hist_n_in,'(a)') 'history_n_'//trim(compname(compid))//'_inst'

       ! Determine instantaneous mediator output frequency and type
       call NUOPC_CompAttributeGet(gcomp, name=trim(hist_option_in), isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          call NUOPC_CompAttributeGet(gcomp, name=trim(hist_option_in), value=hist_option, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeGet(gcomp, name=trim(hist_n_in), value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) hist_n
       else
          ! If attribute is not present - don't write history output
          hist_option = 'none'
          hist_n = -999
       end if

       ! Set alarm name and initialize clock and alarm for instantaneous history output
       if (hist_option /= 'none' .and. hist_option /= 'never') then
          instfile%alarmname =  'alarm_history_inst_'//trim(compname(compid))
          call med_phases_history_init_histclock(gcomp, instfile%clock, &
               instfile%alarm, instfile%alarmname, hist_option, hist_n, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          instfile%is_active = .true.
          instfile%is_clockset = .true.
       else
          instfile%is_active = .false.
          ! this is set to true here even if history file is not active
          instfile%is_clockset = .true.
       end if
    end if ! end of if-clock set if block

    ! if history file is active and history clock is initialized - process history file
    if (instfile%is_active .and. instfile%is_clockset) then

       ! Determine if should write to history file
       call med_phases_history_query_ifwrite(gcomp, instfile%clock, instfile%alarmname, write_now, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! If write now flag is true
       if (write_now) then
          ! Determine time_val and tbnds data for history as well as history file name
          call med_phases_history_set_timeinfo(gcomp, instfile%clock, instfile%alarmname, &
               time_val, time_bnds, time_units, hist_file, doavg=.false., compname=compname(compid), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Create history file
          call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_wopen(hist_file, instfile%io_file, vm, rc, clobber=.true.)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do m = 1,2
             ! Write time values
             if (whead(m)) then
                call ESMF_ClockGet(instfile%clock, calendar=calendar, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call med_io_define_time(instfile%io_file, time_units, calendar, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                call med_io_enddef(instfile%io_file)
                call med_io_write_time(instfile%io_file, time_val, time_bnds, nt=1, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             nx = is_local%wrap%nx(compid)
             ny = is_local%wrap%ny(compid)
             ntile = is_local%wrap%ntile(compid)
             ! Define/write import field bundle
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(compid,compid),rc=rc)) then
                call med_io_write(instfile%io_file, is_local%wrap%FBimp(compid,compid), whead(m), wdata(m), nx, ny, &
                     nt=1, pre=trim(compname(compid))//'Imp', ntile=ntile, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             endif
             ! Define/write import export bundle
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(compid),rc=rc)) then
                call med_io_write(instfile%io_file, is_local%wrap%FBexp(compid), whead(m), wdata(m), nx, ny, &
                     nt=1, pre=trim(compname(compid))//'Exp', ntile=ntile, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             endif
             ! Define/Write mediator fractions
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBFrac(compid),rc=rc)) then
                call med_io_write(instfile%io_file, is_local%wrap%FBFrac(compid), whead(m), wdata(m), nx, ny, &
                     nt=1, pre='Med_frac_'//trim(compname(compid)), ntile=ntile, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

          end do ! end of loop over m

          ! Close file
          call med_io_close(instfile%io_file, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       end if
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write_comp_inst

  !===============================================================================
  subroutine med_phases_history_write_comp_avg(gcomp, compid, avgfile, rc)

    ! Write mediator average history file variables for component compid

    use ESMF              , only : ESMF_FieldBundleIsCreated
    use med_constants_mod , only : czero => med_constants_czero
    use med_methods_mod   , only : med_methods_FB_init, med_methods_FB_reset
    use med_io_mod        , only : med_io_write_time, med_io_define_time

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(in)    :: compid
    type(avgfile_type)  , intent(inout) :: avgfile
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState)     :: is_local
    character(CL)           :: cvalue        ! attribute string
    character(CL)           :: hist_option   ! freq_option setting (ndays, nsteps, etc)
    integer                 :: hist_n        ! freq_n setting relative to freq_option
    character(CL)           :: hist_option_in
    character(CL)           :: hist_n_in
    logical                 :: isPresent
    logical                 :: isSet
    type(ESMF_VM)           :: vm
    type(ESMF_Calendar)     :: calendar          ! calendar type
    integer                 :: m                 ! indices
    integer                 :: nx,ny             ! global grid size
    integer                 :: ntile             ! number of tiles for tiled domain eg CSG
    character(CL)           :: time_units        ! units of time variable
    character(CL)           :: hist_file         ! history file name
    real(r8)                :: time_val          ! time coordinate output
    real(r8)                :: time_bnds(2)      ! time bounds output
    logical                 :: write_now         ! true => write to history type
    character(CS)           :: scalar_name
    character(len=*), parameter :: subname='(med_phases_history_write_comp_avg)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! alarm is not set determine hist_option and hist_n
    if (.not. avgfile%is_clockset) then

       ! Determine attribute name
       write(hist_option_in,'(a)') 'history_option_'//trim(compname(compid))//'_avg'
       write(hist_n_in,'(a)') 'history_n_'//trim(compname(compid))//'_avg'

       ! Determine time average mediator output frequency and type
       call NUOPC_CompAttributeGet(gcomp, name=trim(hist_option_in), isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          call NUOPC_CompAttributeGet(gcomp, name=trim(hist_option_in), value=hist_option, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeGet(gcomp, name=trim(hist_n_in), value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) hist_n
       else
          hist_option = 'none'
          hist_n = -999
       end if
       if (hist_option /= 'never' .and. hist_option /= 'none') then

          ! Set alarm name, initialize clock and alarm for average history output and
          avgfile%alarmname =  'alarm_history_avg_'//trim(compname(compid))
          call med_phases_history_init_histclock(gcomp, avgfile%clock, &
               avgfile%alarm, avgfile%alarmname, hist_option, hist_n, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          avgfile%is_active = .true.
          avgfile%is_clockset = .true.

          ! Initialize accumulation import/export field bundles
          scalar_name = trim(is_local%wrap%flds_scalar_name)
          if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(compid,compid)) .and. .not. &
               ESMF_FieldBundleIsCreated(avgfile%FBaccum_import)) then
             call med_methods_FB_init(avgfile%FBaccum_import, scalar_name, &
                  STgeom=is_local%wrap%NStateImp(compid), STflds=is_local%wrap%NStateImp(compid), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call med_methods_FB_reset(avgfile%FBaccum_import, czero, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             avgfile%accumcnt_import = 0
          end if
          if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(compid)) .and. .not. &
               ESMF_FieldBundleIsCreated(avgfile%FBaccum_export)) then
             call med_methods_FB_init(avgfile%FBaccum_export, scalar_name, &
                  STgeom=is_local%wrap%NStateExp(compid), STflds=is_local%wrap%NStateExp(compid), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call med_methods_FB_reset(avgfile%FBaccum_export, czero, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             avgfile%accumcnt_export = 0
          end if

       else

          avgfile%is_active = .false.
          ! this is set to true here even if history file is not active
          avgfile%is_clockset = .true.

       end if
    end if ! end of if-clock set if block

    ! if history file is active and history clock is initialized - process history file
    if (avgfile%is_active .and. avgfile%is_clockset) then

       ! Determine if will write to history file
       call med_phases_history_query_ifwrite(gcomp, avgfile%clock, avgfile%alarmname, write_now, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Accumulate and then average if write_now flag is true
       if (ESMF_FieldBundleIsCreated(avgfile%FBaccum_import)) then
          call med_phases_history_fldbun_accum(is_local%wrap%FBImp(compid,compid), &
               avgfile%FBaccum_import, avgfile%accumcnt_import,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (write_now) then
             call med_phases_history_fldbun_average(avgfile%FBaccum_import, avgfile%accumcnt_import, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
       if (ESMF_FieldBundleIsCreated(avgfile%FBaccum_export)) then
          call med_phases_history_fldbun_accum(is_local%wrap%FBExp(compid), &
               avgfile%FBaccum_export, avgfile%accumcnt_export, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (write_now) then
             call med_phases_history_fldbun_average(avgfile%FBaccum_export, avgfile%accumcnt_export, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if

       ! If write now flag is true
       if (write_now) then

          ! Determine time_val and tbnds data for history as well as history file name
          call med_phases_history_set_timeinfo(gcomp, avgfile%clock, avgfile%alarmname, &
               time_val, time_bnds, time_units, hist_file, doavg=.true., compname=trim(compname(compid)), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Create history file
          call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_wopen(hist_file, avgfile%io_file, vm, rc, clobber=.true.)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do m = 1,2
             ! Write time values
             if (whead(m)) then
                call ESMF_ClockGet(avgfile%clock, calendar=calendar, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call med_io_define_time(avgfile%io_file, time_units, calendar, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                call med_io_enddef(avgfile%io_file)
                call med_io_write_time(avgfile%io_file, time_val, time_bnds, nt=1, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             ! Write import and export field bundles
             if (is_local%wrap%comp_present(compid)) then
                nx = is_local%wrap%nx(compid)
                ny = is_local%wrap%ny(compid)
                ntile = is_local%wrap%ntile(compid)
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(compid,compid),rc=rc)) then
                   call med_io_write(avgfile%io_file, avgfile%FBaccum_import, whead(m), wdata(m), nx, ny, &
                        nt=1, pre=trim(compname(compid))//'Imp', ntile=ntile, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   if (wdata(m)) then
                      call med_methods_FB_reset(avgfile%FBAccum_import, czero, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                   end if
                endif
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(compid),rc=rc)) then
                   call med_io_write(avgfile%io_file, avgfile%FBaccum_export, whead(m), wdata(m), nx, ny, &
                        nt=1, pre=trim(compname(compid))//'Exp', ntile=ntile, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   if (wdata(m)) then
                      call med_methods_FB_reset(avgfile%FBAccum_export, czero, rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                   end if
                endif
             end if
          end do ! end of loop over m

          ! Close file
          call med_io_close(avgfile%io_file, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       end if ! end of write_now if-block
    end if ! end of clock created if-block

    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write_comp_avg

  !===============================================================================
  subroutine med_phases_history_write_comp_aux(gcomp, compid, auxcomp, rc)

    ! -----------------------------
    ! Write mediator auxiliary history file for auxcomp component
    ! Initialize auxiliary history file
    ! Each time this routine is called the routine SetRunClock in med.F90 is called
    ! at the beginning and the mediator clock current time and time step is set to the
    ! driver current time and time step
    ! -----------------------------

    use ESMF             , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundleRemove
    use ESMF             , only : ESMF_Field, ESMF_FieldGet !DEBUG
    use med_constants_mod, only : czero => med_constants_czero
    use med_io_mod       , only : med_io_write_time, med_io_define_time
    use med_methods_mod  , only : med_methods_FB_init
    use med_methods_mod  , only : med_methods_FB_reset
    use med_methods_mod  , only : med_methods_FB_fldchk

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(in)    :: compid
    type(auxcomp_type)  , intent(inout) :: auxcomp
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(ESMF_VM)           :: vm
    type(ESMF_Calendar)     :: calendar    ! calendar type
    logical                 :: isPresent   ! is attribute present
    logical                 :: isSet       ! is attribute set
    character(CL)           :: hist_option ! freq_option setting (ndays, nsteps, etc)
    integer                 :: hist_n      ! freq_n setting relative to freq_option
    integer                 :: nfcnt
    integer                 :: nfile
    integer                 :: nfld
    integer                 :: n,n1,nf
    character(CL)           :: prefix
    character(CL)           :: cvalue
    character(CL)           :: auxflds
    integer                 :: fieldCount
    logical                 :: found
    logical                 :: enable_auxfile
    character(CL)           :: time_units        ! units of time variable
    integer                 :: nx,ny             ! global grid size
    logical                 :: write_now         ! if true, write time sample to file
    real(r8)                :: time_val          ! time coordinate output
    real(r8)                :: time_bnds(2)      ! time bounds output
    character(CS), allocatable  :: fieldNameList(:)
    character(len=*), parameter :: subname='(med_phases_history_write_comp_aux)'
    !---------------------------------------

    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.not. auxcomp%init_auxfiles) then

       ! Initialize number of aux files for this component to zero
       nfcnt = 0
       do nfile = 1,max_auxfiles
          ! Determine attribute prefix
          write(prefix,'(a,i0)') 'histaux_'//trim(compname(compid))//'2med_file',nfile

          ! Determine if will write the file
          call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_enabled', isPresent=isPresent, isSet=isSet, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (isPresent .and. isSet) then
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_enabled', value=cvalue, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             read(cvalue,'(l7)') enable_auxfile
          else
             enable_auxfile = .false.
          end if

          ! If file will be written - then initialize auxcomp%files(nfcnt)
          if (enable_auxfile) then
             ! Increment nfcnt
             nfcnt = nfcnt + 1

             ! Determine number of time samples per file
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_ntperfile', value=cvalue, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             read(cvalue,*) auxcomp%files(nfcnt)%ntperfile
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! Determine if will do time average for aux file
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_doavg', value=cvalue, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             read(cvalue,*) auxcomp%files(nfcnt)%doavg

             ! Determine the colon delimited field names for this file
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_flds', value=auxflds, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! Determine fields that will be output to auxhist files
             if (trim(auxflds) == 'all') then

                ! Output all fields sent to the mediator from ncomp to the auxhist files
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(compid,compid), fieldCount=fieldCount, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                allocate(auxcomp%files(nfcnt)%flds(fieldcount))
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(compid,compid), fieldNameList=auxcomp%files(nfcnt)%flds, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

             else

                ! Translate the colon deliminted string (auxflds) into a character array (fieldnamelist)
                ! Note that the following call allocates the memory for fieldnamelist
                call get_auxflds(auxflds, fieldnamelist, rc)

                ! TODO: print warning statement if remove field
                ! TODO: if request field that is NOT in the field definition file - then quit
                ! Remove all fields from fieldnamelist that are not in FBImp(compid,compid)
                fieldCount = size(fieldnamelist)
                do n = 1,fieldcount
                   if (.not. med_methods_FB_fldchk(is_local%wrap%FBImp(compid,compid), trim(fieldnamelist(n)), rc)) then
                      do n1 = n, fieldCount-1
                         fieldnamelist(n1) = fieldnamelist(n1+1)
                      end do
                      fieldCount = fieldCount - 1
                   end if
                end do

                ! Create auxcomp%files(nfcnt)%flds array
                allocate(auxcomp%files(nfcnt)%flds(fieldcount))
                do n = 1,fieldcount
                   auxcomp%files(nfcnt)%flds(n) = trim(fieldnamelist(n))
                end do

                ! Deallocate memory from fieldnamelist
                deallocate(fieldnamelist) ! this was allocated in med_phases_history_get_auxflds

             end if ! end of if auxflds is set to 'all'

             if (maintask) then
                write(logunit,*)
                write(logunit,'(a,i4,a)') trim(subname) // '   Writing the following fields to auxfile ',nfcnt,&
                     ' for component '//trim(compname(compid))
                do nfld = 1,size(auxcomp%files(nfcnt)%flds)
                   write(logunit,'(8x,a)') trim(auxcomp%files(nfcnt)%flds(nfld))
                end do
             end if

             ! Create FBaccum if averaging is on
             if (auxcomp%files(nfcnt)%doavg) then

                ! First duplicate all fields in FBImp(compid,compid)
                call ESMF_LogWrite(trim(subname)// ": initializing FBaccum(compid)", ESMF_LOGMSG_INFO)
                if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compid,compid)) .and. .not. &
                     ESMF_FieldBundleIsCreated(auxcomp%files(nfcnt)%FBaccum)) then
                   call med_methods_FB_init(auxcomp%files(nfcnt)%FBaccum, is_local%wrap%flds_scalar_name, &
                        STgeom=is_local%wrap%NStateImp(compid), STflds=is_local%wrap%NStateImp(compid), &
                        rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   call med_methods_FB_reset(auxcomp%files(nfcnt)%FBaccum, czero, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   auxcomp%files(nfcnt)%accumcnt = 0
                end if

                ! Now remove all fields from FBAccum that are not in the input flds list
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(compid,compid), fieldCount=fieldCount, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                allocate(fieldNameList(fieldCount))
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(compid,compid), fieldNameList=fieldNameList, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                do n = 1,size(fieldnamelist)
                   found = .false.
                   do n1 = 1,size(auxcomp%files(nfcnt)%flds)
                      if (trim(fieldnamelist(n)) == trim(auxcomp%files(nfcnt)%flds(n1))) then
                         found = .true.
                         exit
                      end if
                   end do
                   if (.not. found) then
                      call ESMF_FieldBundleRemove(auxcomp%files(nfcnt)%FBaccum, fieldnamelist(n:n), rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                   end if
                end do
                deallocate(fieldnameList)

                ! Check that FBAccum has at least one field left - if not exit
                call ESMF_FieldBundleGet(auxcomp%files(nfcnt)%FBAccum, fieldCount=nfld, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                if (nfld == 0) then
                   call ESMF_LogWrite(subname//'FBAccum is zero for '//trim(auxcomp%files(nfcnt)%auxname), &
                        ESMF_LOGMSG_ERROR)
                   rc = ESMF_FAILURE
                   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) then
                      call ESMF_Finalize(endflag=ESMF_END_ABORT)
                   end if
                end if

             end if

             ! Determine auxiliary file output frequency and type
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_history_option', value=hist_option, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_history_n', value=cvalue, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             read(cvalue,*) hist_n

             ! Determine alarmname
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_auxname', value=auxcomp%files(nfcnt)%auxname, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             write(auxcomp%files(nfcnt)%alarmname,'(a,i0)') 'alarm_'//trim(prefix)

             ! Initialize clock and alarm for instantaneous history output
             call med_phases_history_init_histclock(gcomp, auxcomp%files(nfcnt)%clock, &
                  auxcomp%files(nfcnt)%alarm, auxcomp%files(nfcnt)%alarmname, hist_option, hist_n, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

          end if ! end of isPresent and isSet and  if flag is on for file n
       end do ! end of loop over nfile

       ! Set number of aux files for this component - this is a module variable
       auxcomp%num_auxfiles = nfcnt

       ! Set initialization flags to .true.
       auxcomp%init_auxfiles = .true.

    end if ! end of initialization if-block

    ! Write auxiliary history files for component compid
    do nf = 1,auxcomp%num_auxfiles

       ! Determine if will write to history file
       call med_phases_history_query_ifwrite(gcomp, auxcomp%files(nf)%clock, auxcomp%files(nf)%alarmname, write_now, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Do accumulation and average if required
       if (auxcomp%files(nf)%doavg) then
          call med_phases_history_fldbun_accum(is_local%wrap%FBImp(compid,compid), &
               auxcomp%files(nf)%FBaccum, auxcomp%files(nf)%accumcnt, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (write_now) then
             call med_phases_history_fldbun_average(auxcomp%files(nf)%FBaccum, auxcomp%files(nf)%accumcnt, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       end if

       ! Write time sample to file
       if ( write_now ) then

          ! Set shorthand variables
          nx = is_local%wrap%nx(compid)
          ny = is_local%wrap%ny(compid)

          ! Increment number of time samples on file
          auxcomp%files(nf)%nt = auxcomp%files(nf)%nt + 1

          ! Determine time_val and tbnds data for history as well as history file name
          if (auxcomp%files(nf)%nt == 1) then
             call med_phases_history_set_timeinfo(gcomp, auxcomp%files(nf)%clock, auxcomp%files(nf)%alarmname, &
                  time_val, time_bnds, time_units, auxcomp%files(nf)%histfile, auxcomp%files(nf)%doavg, &
                  auxname=auxcomp%files(nf)%auxname, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call med_phases_history_set_timeinfo(gcomp, auxcomp%files(nf)%clock, auxcomp%files(nf)%alarmname, &
                  time_val, time_bnds, time_units, auxcomp%files(nf)%histfile, auxcomp%files(nf)%doavg, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Write  header
          if (auxcomp%files(nf)%nt == 1) then

             ! open file
             call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_wopen(auxcomp%files(nf)%histfile, auxcomp%files(nf)%io_file, vm, rc, file_ind=nf, clobber=.true.)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! define time variables
             call ESMF_ClockGet(auxcomp%files(nf)%clock, calendar=calendar, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_define_time(auxcomp%files(nf)%io_file, time_units, calendar, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! define data variables with a time dimension (include the nt argument below)
             call med_io_write(auxcomp%files(nf)%io_file, is_local%wrap%FBimp(compid,compid), &
                  whead(1), wdata(1), nx, ny, nt=auxcomp%files(nf)%nt, &
                  pre=trim(compname(compid))//'Imp', flds=auxcomp%files(nf)%flds, &
                  use_float=.true., rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! end definition phase
             call med_io_enddef(auxcomp%files(nf)%io_file)
          end if

          ! Write time variables for time nt
          call med_io_write_time(auxcomp%files(nf)%io_file, time_val, time_bnds, nt=auxcomp%files(nf)%nt, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Write data variables for time nt
          if (auxcomp%files(nf)%doavg) then
             call med_io_write(auxcomp%files(nf)%io_file, auxcomp%files(nf)%FBaccum, whead(2), wdata(2), nx, ny, &
                  nt=auxcomp%files(nf)%nt, pre=trim(compname(compid))//'Imp', flds=auxcomp%files(nf)%flds, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_methods_FB_reset(auxcomp%files(nf)%FBaccum, value=czero, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call med_io_write(auxcomp%files(nf)%io_file, is_local%wrap%FBimp(compid,compid), whead(2), wdata(2), nx, ny, &
                  nt=auxcomp%files(nf)%nt, pre=trim(compname(compid))//'Imp', flds=auxcomp%files(nf)%flds, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Close file
          if (auxcomp%files(nf)%nt == auxcomp%files(nf)%ntperfile) then
             call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_close(auxcomp%files(nf)%io_file, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             auxcomp%files(nf)%nt = 0
          end if

       end if ! end of write_now if-block

    end do
    call t_stopf('MED:'//subname)

  contains

    subroutine get_auxflds(str, flds, rc)
      ! input/output variables
      character(len=*)               , intent(in)  :: str     ! colon deliminted string to search
      character(len=*) , allocatable , intent(out) :: flds(:) ! memory will be allocate for flds
      integer                        , intent(out) :: rc
      ! local variables
      integer          :: i,k,n ! generic indecies
      integer          :: nflds ! allocatable size of flds
      integer          :: count ! counts occurances of char
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
         if (maintask) write(logunit,*) "ERROR: invalid list = ",trim(str)
         call ESMF_LogWrite("ERROR: invalid list = "//trim(str), ESMF_LOGMSG_ERROR)
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
      ! allocate memory for flds)
      allocate(flds(nflds))
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
         ! set flds(k)
         flds(k) = str(i0:i1)//"   "
      end do
    end subroutine get_auxflds

  end subroutine med_phases_history_write_comp_aux

  !===============================================================================
  subroutine med_phases_history_fldbun_accum(fldbun, fldbun_accum, count, rc)

    use ESMF, only : ESMF_Field, ESMF_FieldGet

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: fldbun
    type(ESMF_FieldBundle) , intent(inout) :: fldbun_accum
    integer                , intent(out)   :: count
    integer                , intent(out)   :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: lfield
    type(ESMF_Field)       :: lfield_accum
    integer                :: fieldCount_accum
    character(CL), pointer :: fieldnames_accum(:)
    real(r8), pointer      :: dataptr1d(:)
    real(r8), pointer      :: dataptr2d(:,:)
    real(r8), pointer      :: dataptr1d_accum(:)
    real(r8), pointer      :: dataptr2d_accum(:,:)
    integer                :: ungriddedUBound_accum(1)
    integer                :: ungriddedUBound(1)
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Loop over field names in fldbun_accum

    call ESMF_FieldBundleGet(fldbun_accum, fieldCount=fieldCount_accum, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldnames_accum(fieldCount_accum))
    call ESMF_FieldBundleGet(fldbun_accum, fieldCount=fieldCount_accum, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(fldbun_accum, fieldNameList=fieldnames_accum, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldcount_accum
       call ESMF_FieldBundleGet(fldbun, fieldName=trim(fieldnames_accum(n)), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldBundleGet(fldbun_accum, fieldName=trim(fieldnames_accum(n)), field=lfield_accum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield_accum, ungriddedUBound=ungriddedUBound_accum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (ungriddedUBound(1) /= ungriddedUBound_accum(1)) then
          call ESMF_LogWrite(" upper bounds for field and field_accum do not match", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
       end if

       if (ungriddedUBound(1) > 0) then
          call ESMF_FieldGet(lfield, farrayptr=dataptr2d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield_accum, farrayptr=dataptr2d_accum, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr2d_accum(:,:) = dataptr2d_accum(:,:) + dataptr2d(:,:)
       else
          call ESMF_FieldGet(lfield, farrayptr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield_accum, farrayptr=dataptr1d_accum, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          dataptr1d_accum(:) = dataptr1d_accum(:) + dataptr1d(:)
       end if
    end do
    deallocate(fieldnames_accum)

    ! Accumulate counter
    count = count + 1

  end subroutine med_phases_history_fldbun_accum

  !===============================================================================
  subroutine med_phases_history_fldbun_average(fldbun_accum, count, rc)

    use ESMF              , only : ESMF_Field, ESMF_FieldGet
    use med_constants_mod , only : czero => med_constants_czero

    ! input/output variables
    type(ESMF_FieldBundle) , intent(inout) :: fldbun_accum
    integer                , intent(inout) :: count
    integer                , intent(out)   :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: lfield_accum
    integer                :: fieldCount
    character(CL), pointer :: fieldnames(:)
    real(r8), pointer      :: dataptr1d_accum(:)
    real(r8), pointer      :: dataptr2d_accum(:,:)
    integer                :: ungriddedUBound(1)
    !---------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fldbun_accum, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldnames(fieldCount))
    call ESMF_FieldBundleGet(fldbun_accum, fieldNameList=fieldnames, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1, fieldcount
       call ESMF_FieldBundleGet(fldbun_accum, fieldName=trim(fieldnames(n)), field=lfield_accum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield_accum, ungriddedUBound=ungriddedUBound, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (ungriddedUBound(1) > 0) then
          call ESMF_FieldGet(lfield_accum, farrayptr=dataptr2d_accum, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (count == 0) then
             dataptr2d_accum(:,:) = czero
          else
             dataptr2d_accum(:,:) = dataptr2d_accum(:,:) / real(count, r8)
          end if
       else
          call ESMF_FieldGet(lfield_accum, farrayptr=dataptr1d_accum, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (count == 0) then
             dataptr1d_accum(:) = czero
          else
             dataptr1d_accum(:) = dataptr1d_accum(:) / real(count, r8)
          end if
       end if
    end do
    deallocate(fieldnames)

    ! Reset counter
    count = 0

  end subroutine med_phases_history_fldbun_average

  !===============================================================================
  subroutine med_phases_history_init_histclock(gcomp, hclock, alarm, alarmname, hist_option, hist_n, rc)

    use NUOPC_Mediator, only : NUOPC_MediatorGet
    use ESMF          , only : ESMF_ClockCreate, ESMF_ClockGet, ESMF_ClockSet
    use med_time_mod  , only : med_time_alarmInit

    ! input/output variables
    type(ESMF_GridComp) , intent(in)    :: gcomp
    type(ESMF_Clock)    , intent(inout) :: hclock
    type(ESMF_Alarm)    , intent(inout) :: alarm
    character(len=*)    , intent(in)    :: alarmname
    character(len=*)    , intent(in)    :: hist_option  ! freq_option setting (ndays, nsteps, etc)
    integer             , intent(in)    :: hist_n       ! freq_n setting relative to freq_option
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Clock)        :: mclock, dclock
    type(ESMF_Time)         :: StartTime
    type(ESMF_TimeInterval) :: mtimestep, dtimestep
    integer                 :: msec, dsec
    character(len=*), parameter :: subname='(med_phases_history_init_histclock) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_MediatorGet(gcomp, mediatorClock=mClock, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, timeStep=mtimestep, starttime=starttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(mtimestep, s=msec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(dtimestep, s=dsec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (maintask) then
       write(logunit,'(a,2x,i8,2x,i8)') trim(subname) // "  mediator, driver timesteps for " &
            //trim(alarmname),msec,dsec
    end if

    ! Create history clock from mediator clock - THIS CALL DOES NOT COPY ALARMS
    hclock = ESMF_ClockCreate(mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize history alarm and advance history clock to trigger
    ! alarms then reset history clock back to mcurrtime
    call med_time_alarmInit(hclock, alarm, option=hist_option, opt_n=hist_n, &
         reftime=StartTime, alarmname=trim(alarmname), advance_clock=.true., rc=rc)

    ! Write diagnostic info
    if (maintask) then
       write(logunit,'(a,2x,i8)') trim(subname) // "  initialized history alarm "//&
            trim(alarmname)//"  with option "//trim(hist_option)//" and frequency ",hist_n
    end if

  end subroutine med_phases_history_init_histclock

  !===============================================================================
  subroutine med_phases_history_query_ifwrite(gcomp, hclock, alarmname, write_now, rc)

    use NUOPC_Mediator, only : NUOPC_MediatorGet

    ! input/output variables
    type(ESMF_GridComp) , intent(in)    :: gcomp
    type(ESMF_Clock)    , intent(inout) :: hclock    ! write clock
    character(len=*)    , intent(in)    :: alarmname ! write alarmname
    logical             , intent(out)   :: write_now ! if true => write now
    integer             , intent(out)   :: rc        ! error code

    ! local variables
    type(ESMF_Clock)        :: mclock         ! mediator clock
    type(ESMF_Alarm)        :: alarm          ! write alarm
    type(ESMF_Time)         :: currtime       ! current time
    character(len=CS)       :: currtimestr    ! current time string
    type(ESMF_Time)         :: nexttime       ! next time
    character(len=CS)       :: nexttimestr    ! next time string
    integer                 :: yr,mon,day,sec ! time units
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: ringInterval_length
    character(len=*), parameter :: subname='(med_phases_history_query_ifwrite) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Update hclock to trigger alarm
    call ESMF_ClockAdvance(hclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the history file alarm and determine if alarm is ringing
    call ESMF_ClockGetAlarm(hclock, alarmname=trim(alarmname), alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set write_now flag and turn ringer off if appropriate
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write_now = .true.
       call ESMF_AlarmRingerOff(alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write_now = .false.
    end if

    ! Write diagnostic output
    if (write_now) then
       if (maintask .and. debug_alarms) then
          ! output alarm info
          call ESMF_AlarmGet(alarm, ringInterval=ringInterval, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeIntervalGet(ringInterval, s=ringinterval_length, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_ClockGet(hclock, currtime=currtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
          call ESMF_ClockGetNextTime(hclock, nextTime=nexttime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec

          if (maintask) then
             write(logunit,*)
             write(logunit,'(a,i8)') trim(subname)//" : history alarmname "//trim(alarmname)//&
                  ' is ringing, interval length is ', ringInterval_length
             write(logunit,'(a)') trim(subname)//" : hclock currtime = "//trim(currtimestr)//&
                  " hclock nexttime = "//trim(nexttimestr)
          end if

          call NUOPC_MediatorGet(gcomp, mediatorClock=mClock, rc=rc)
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
          write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec

          if (maintask) then
             write(logunit,'(a)') trim(subname)//" : mclock currtime = "//trim(currtimestr)//&
                  " mclock nexttime = "//trim(nexttimestr)
          end if

       end if
    end if

  end subroutine med_phases_history_query_ifwrite

  !===============================================================================
  subroutine med_phases_history_set_timeinfo(gcomp, hclock, alarmname, &
       time_val, time_bnds, time_units, histfile, doavg, auxname, compname, rc)

    use NUOPC_Mediator    , only : NUOPC_MediatorGet
    use ESMF              , only : ESMF_GridComp, ESMF_Clock, ESMF_Alarm, ESMF_Time, ESMF_TimeInterval
    use ESMF              , only : ESMF_ClockGet, ESMF_ClockGetNextTime, ESMF_ClockGetAlarm
    use ESMF              , only : ESMF_AlarmGet, ESMF_TimeIntervalGet, ESMF_TimeGet
    use med_constants_mod , only : SecPerDay => med_constants_SecPerDay
    use med_io_mod        , only : med_io_ymd2date, med_io_date2yyyymmdd, med_io_sec2hms

    ! input/output variables
    type(ESMF_GridComp)         , intent(in)  :: gcomp
    type(ESMF_Clock)            , intent(in)  :: hclock
    character(len=*)            , intent(in)  :: alarmname
    real(r8)                    , intent(out) :: time_val
    real(r8)                    , intent(out) :: time_bnds(2)
    character(len=*)            , intent(out) :: time_units
    character(len=*)            , intent(out) :: histfile
    logical                     , intent(in)  :: doavg
    character(len=*) , optional , intent(in)  :: auxname
    character(len=*) , optional , intent(in)  :: compname
    integer                     , intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: mclock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: nexttime
    type(ESMF_TimeInterval) :: ringInterval   ! alarm interval
    type(ESMF_TimeInterval) :: timediff(2)    ! time bounds upper and lower relative to start
    character(len=CL)       :: currtime_str
    character(len=CL)       :: nexttime_str
    character(len=CL)       :: hist_str
    integer                 :: yr,mon,day,sec ! time units
    integer                 :: start_ymd      ! Starting date YYYYMMDD
    logical                 :: isPresent
    character(len=*), parameter :: subname='(med_phases_history_set_timeinfo) '
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Determine starttime, currtime and nexttime from the mediator clock rather than the input history clock
    call NUOPC_MediatorGet(gcomp, mediatorClock=mClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet(mclock, currtime=currtime, starttime=starttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetNextTime(mclock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine time units
    call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_io_ymd2date(yr,mon,day,start_ymd)
    time_units = 'days since ' // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(sec, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set time bounds and time coord
    if (doavg) then
       call ESMF_ClockGetAlarm(hclock, alarmname=trim(alarmname), alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmGet(alarm, ringInterval=ringInterval, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       timediff(2) = nexttime - starttime
       timediff(1) = nexttime - starttime - ringinterval
       call ESMF_TimeIntervalGet(timediff(2), d_r8=time_bnds(2), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeIntervalGet(timediff(1), d_r8=time_bnds(1), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       time_val = 0.5_r8 * (time_bnds(1) + time_bnds(2))
    else
       timediff(1) = nexttime - starttime
       call ESMF_TimeIntervalGet(timediff(1), d=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       time_val = day + sec/real(SecPerDay,R8)
       time_bnds(1) = time_val
       time_bnds(2) = time_val
    end if

    ! Determine history file name
    ! Use nexttime_str rather than currtime_str here since that is the time at the end of
    ! the timestep and is preferred for history file names

    call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(nexttime_str,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec

    if (trim(case_name) == 'unset') then
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

    if (present(auxname)) then
       write(histfile, "(8a)") trim(case_name),'.cpl' ,trim(inst_tag),'.hx.',trim(auxname),'.',&
            trim(nexttime_str),'.nc'
    else if (present(compname)) then
       if (doavg) then
          hist_str = '.ha.'
       else
          hist_str = '.hi.'
       end if
       if (trim(compname) /= 'all') then
          hist_str = trim(hist_str) // trim(compname) // '.'
       end if
       write(histfile, "(6a)") trim(case_name),'.cpl',trim(inst_tag),trim(hist_str),trim(nexttime_str),'.nc'
    end if

    if (maintask) then
       call ESMF_TimeGet(currtime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(currtime_str,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
       write(logunit,*)
       write(logunit,' (a)') trim(subname) // " writing mediator history file "//trim(histfile)
       write(logunit,' (a)') trim(subname) // "   currtime = "//trim(currtime_str)//" nexttime = "//trim(nexttime_str)
    end if

  end subroutine med_phases_history_set_timeinfo

end module med_phases_history_mod
