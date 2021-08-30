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
  ! itself. However, for cases with cocnatenated loops on the upper level
  ! of the run sequence in freeFormat, a single outer loop is added
  ! automatically during ingestion, and the driver clock is used for this
  ! loop instead.
  !-----------------------------------------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_VM, ESMF_VMGet
  use ESMF                  , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockSet, ESMF_ClockAdvance, ESMF_ClockCreate
  use ESMF                  , only : ESMF_ClockGetNextTime, ESMF_ClockGetAlarm, ESMF_ClockIsCreated
  use ESMF                  , only : ESMF_Calendar
  use ESMF                  , only : ESMF_Time, ESMF_TimeGet
  use ESMF                  , only : ESMF_TimeInterval, ESMF_TimeIntervalGet, ESMF_TimeIntervalSet
  use ESMF                  , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmSet
  use ESMF                  , only : ESMF_AlarmIsRinging, ESMF_AlarmRingerOff, ESMF_AlarmGet
  use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet
  use ESMF                  , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundleRemove
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_MAXSTR, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
  use ESMF                  , only : ESMF_Finalize
  use ESMF                  , only : operator(==), operator(-), operator(+), operator(/=), operator(<=)
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use NUOPC_Model           , only : NUOPC_ModelGet
  use esmFlds               , only : compmed, compatm, complnd, compocn, compice, comprof, compglc, compwav
  use esmFlds               , only : ncomps, compname, num_icesheets
  use esmFlds               , only : fldListFr, fldListTo
  use med_constants_mod     , only : SecPerDay => med_constants_SecPerDay
  use med_utils_mod         , only : chkerr    => med_utils_ChkErr
  use med_methods_mod       , only : med_methods_FB_reset
  use med_methods_mod       , only : med_methods_FB_fldchk
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_time_mod          , only : med_time_alarmInit
  use med_io_mod            , only : med_io_write, med_io_wopen, med_io_enddef
  use med_io_mod            , only : med_io_close, med_io_date2yyyymmdd, med_io_sec2hms
  use med_io_mod            , only : med_io_ymd2date
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  ! Public routines called from the run sequence
  public :: med_phases_history_write     ! inst only
  public :: med_phases_history_write_med ! inst only
  public :: med_phases_history_write_atm ! inst, avg, aux
  public :: med_phases_history_write_ice ! inst, avg, aux
  public :: med_phases_history_write_glc ! inst, avg, aux
  public :: med_phases_history_write_lnd ! inst, avg, aux
  public :: med_phases_history_write_ocn ! inst, avg, aux
  public :: med_phases_history_write_rof ! inst, avg, aux
  public :: med_phases_history_write_wav ! inst, avg, aux

  ! Private routines
  private :: med_phases_history_write_inst_comp ! write instantaneous file for a given component
  private :: med_phases_history_write_avg_comp  ! write averaged file for a given component
  private :: med_phases_history_write_aux_comp  ! write auxiliary file for a given component
  private :: med_phases_history_write_hfile
  private :: med_phases_history_write_hfileaux
  private :: med_phases_history_get_filename
  private :: med_phases_history_get_auxflds
  private :: med_phases_history_output_alarminfo
  private :: med_phases_history_ymds2rday_offset

  character(CL) :: case_name = 'unset'  ! case name
  character(CS) :: inst_tag = 'unset'   ! instance tag

  ! ----------------------------
  ! Time averaging history files
  ! ----------------------------
  type, public :: avgfile_type
     type(ESMF_FieldBundle) :: FBaccum    ! field bundle for time averaging
     integer                :: accumcnt   ! field bundle accumulation counter
  end type avgfile_type
  type(avgfile_type) :: avgfiles_import(ncomps)
  type(avgfile_type) :: avgfiles_export(ncomps)
  type(avgfile_type) :: avgfiles_aoflux_ocn
  type(avgfile_type) :: avgfiles_ocnalb_ocn
  type(avgfile_type) :: avgfiles_aoflux_atm
  type(avgfile_type) :: avgfiles_ocnalb_atm
  type(ESMF_Clock)   :: hclock_avg_comp(ncomps)

  ! ----------------------------
  ! Auxiliary history files
  ! ----------------------------
  integer, parameter :: max_auxfiles = 10
  type, public :: auxfile_type
     character(CS), allocatable :: flds(:)       ! array of aux field names
     character(CS)              :: auxname       ! name for history file creation
     character(CL)              :: histfile = '' ! current history file name
     character(CS)              :: alarmname     ! name of write alarm
     integer                    :: ntperfile     ! maximum number of time samples per file
     integer                    :: nt = 0        ! time in file
     logical                    :: useavg        ! if true, time average, otherwise instantaneous
     type(ESMF_FieldBundle)     :: FBaccum       ! field bundle for time averaging
     integer                    :: accumcnt      ! field bundle accumulation counter
     type(ESMF_Clock)           :: hclock        ! auxiliary history clock
  end type auxfile_type
  integer            , public :: num_auxfiles(ncomps) = 0
  type(auxfile_type) , public :: auxfiles(max_auxfiles,ncomps)

  ! ----------------------------
  ! Instantaneous history files
  ! ----------------------------
  type(ESMF_Clock) :: hclock_inst_all
  type(ESMF_Clock) :: hclock_inst_comp(ncomps)

  logical       :: debug_alarms = .true.
  character(*), parameter :: u_FILE_u  = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_history_write(gcomp, rc)
    ! --------------------------------------
    ! Write instantaneous mediator history file for all variables
    ! --------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(CL)               :: alarmname
    type(ESMF_Clock)            :: mclock
    type(ESMF_Alarm)            :: alarm
    type(ESMF_Time)             :: CurrTime
    type(ESMF_Time)             :: StartTime
    type(ESMF_TimeInterval)     :: timestep
    integer                     :: timestep_length
    character(CL)               :: hist_option   ! freq_option setting (ndays, nsteps, etc)
    integer                     :: hist_n        ! freq_n setting relative to freq_option
    character(CL)               :: cvalue        ! attribute string
    logical                     :: isPresent
    logical                     :: isSet
    logical                     :: first_time = .true.
    character(len=*), parameter :: subname='(med_phases_history_write)'
    !---------------------------------------
    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    alarmname =  'alarm_history_inst_all'
    if (first_time) then
       call NUOPC_CompAttributeGet(gcomp, name='history_option', isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          call NUOPC_CompAttributeGet(gcomp, name='history_option', value=hist_option, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeGet(gcomp, name='history_n', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) hist_n
       else
          ! If attribute is not present - don't write history output
          hist_option = 'none'
          hist_n = -999
       end if

       if (hist_option /= 'none' .and. hist_option /= 'never') then
          ! First create hclock from mclock - THIS CALL DOES NOT COPY ALARMS
          call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          hclock_inst_all = ESMF_ClockCreate(mclock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Set alarm for instantaneous history output
          ! Advance history clock to trigger alarms then reset history clock back to mcurrtime
          call ESMF_ClockGet(hclock_inst_all, startTime=StartTime,  currTime=CurrTime, timeStep=timestep, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_time_alarmInit(hclock_inst_all, alarm, option=hist_option, opt_n=hist_n, &
               reftime=StartTime, alarmname=trim(alarmname), rc=rc)
          call ESMF_AlarmSet(alarm, clock=hclock_inst_all, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockAdvance(hclock_inst_all,rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockSet(hclock_inst_all, currTime=currtime)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Write diagnostic info
          if (mastertask) then
             call ESMF_TimeIntervalGet(timestep, s=timestep_length, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             write(logunit,'(a,2x,i8)') "  initialized instantaneous history alarm "//&
                  trim(alarmname)//"  with option "//trim(hist_option)//" and frequency ",hist_n
             write(logunit,'(a,2x,i8)') "  history clock timestep = ",timestep_length
          end if
       end if

       first_time = .false.
    end if

    if (ESMF_ClockIsCreated(hclock_inst_all)) then
       ! Advance the clock
       call ESMF_ClockAdvance(hclock_inst_all, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Write the instantaneous history file for all relevant components
       call med_phases_history_write_hfile(gcomp, 'all', hclock_inst_all, 'alarm_history_inst_all', .false., rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write

  !===============================================================================
  subroutine med_phases_history_write_med(gcomp, rc)
    ! Write mediator history file for med variables - only instantaneous files are written
    ! This writes out ocean albedoes and atm/ocean fluxes computed by the mediator
    ! along with the fractions computed by the mediator
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    logical              :: first_time = .true.
    !---------------------------------------
    rc = ESMF_SUCCESS
    call med_phases_history_write_inst_comp(gcomp, compmed, &
         first_time, 'med_phases_history_write_inst_med', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (first_time) first_time = .false.
  end subroutine med_phases_history_write_med

  !===============================================================================
  subroutine med_phases_history_write_atm(gcomp, rc)
    ! Write mediator history file for atm variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    logical              :: first_time = .true.
    !---------------------------------------
    rc = ESMF_SUCCESS
    call med_phases_history_write_inst_comp(gcomp, compatm, &
         first_time, 'med_phases_history_write_inst_atm', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_avg_comp(gcomp, compatm, &
         first_time, 'med_phases_history_write_avg_atm', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_aux_comp(gcomp, compatm, &
         first_time, 'med_phases_history_write_aux_atm', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (first_time) first_time = .false.
  end subroutine med_phases_history_write_atm

  !===============================================================================
  ! Write mediator history file for ice variables
  subroutine med_phases_history_write_ice(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    logical              :: first_time = .true.
    !---------------------------------------
    rc = ESMF_SUCCESS
    call med_phases_history_write_inst_comp(gcomp, compice, &
         first_time, 'med_phases_history_write_inst_ice', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_avg_comp(gcomp, compice, &
         first_time, 'med_phases_history_write_avg_ice', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_aux_comp(gcomp, compice, &
         first_time, 'med_phases_history_write_aux_ice', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (first_time) first_time = .false.
  end subroutine med_phases_history_write_ice

  !===============================================================================
  ! Write mediator history file for glc variables
  subroutine med_phases_history_write_glc(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    integer              :: ns
    character(len=CS)    :: cns
    logical              :: first_time = .true.
    !---------------------------------------
    rc = ESMF_SUCCESS
    do ns = 1,num_icesheets
       write(cns,*) ns
       call med_phases_history_write_inst_comp(gcomp, compglc(ns), &
            first_time, 'med_phases_history_write_inst_glc'//trim(cns), rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_phases_history_write_avg_comp(gcomp, compglc(ns), &
            first_time, 'med_phases_history_write_avg_glc'//trim(cns), rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_phases_history_write_aux_comp(gcomp, compglc(ns), &
            first_time, 'med_phases_history_write_aux_glc'//trim(cns), rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
    if (first_time) first_time = .false.
  end subroutine med_phases_history_write_glc

  !===============================================================================
  ! Write mediator history file for lnd variables
  subroutine med_phases_history_write_lnd(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    logical              :: first_time = .true.
    !---------------------------------------
    rc = ESMF_SUCCESS
    call med_phases_history_write_inst_comp(gcomp, complnd, &
         first_time, 'med_phases_history_write_inst_lnd', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_avg_comp(gcomp, complnd, &
         first_time, 'med_phases_history_write_avg_lnd', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_aux_comp(gcomp, complnd, &
         first_time, 'med_phases_history_write_aux_lnd', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (first_time) first_time = .false.
  end subroutine med_phases_history_write_lnd

  !===============================================================================
  ! Write mediator history file for ocn variables
  subroutine med_phases_history_write_ocn(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    logical              :: first_time = .true.
    !---------------------------------------
    rc = ESMF_SUCCESS
    call med_phases_history_write_inst_comp(gcomp, compocn, &
         first_time, 'med_phases_history_write_inst_ocn', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_avg_comp(gcomp, compocn, &
         first_time, 'med_phases_history_write_avg_ocn', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_aux_comp(gcomp, compocn, &
         first_time, 'med_phases_history_write_aux_ocn', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (first_time) first_time = .false.
  end subroutine med_phases_history_write_ocn

  !===============================================================================
  ! Write mediator history file for rof variables
  subroutine med_phases_history_write_rof(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    logical              :: first_time = .true.
    !---------------------------------------
    rc = ESMF_SUCCESS
    call med_phases_history_write_inst_comp(gcomp, comprof, &
         first_time, 'med_phases_history_write_inst_rof', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_avg_comp(gcomp, comprof, &
         first_time, 'med_phases_history_write_avg_rof', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_aux_comp(gcomp, comprof, &
         first_time, 'med_phases_history_write_aux_rof', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (first_time) first_time = .false.
  end subroutine med_phases_history_write_rof

  !===============================================================================
  ! Write mediator history file for wav variables
  subroutine med_phases_history_write_wav(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    logical              :: first_time = .true.
    !---------------------------------------
    rc = ESMF_SUCCESS
    call med_phases_history_write_inst_comp(gcomp, compwav, &
         first_time, 'med_phases_history_write_inst_wav', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_avg_comp(gcomp, compwav, &
         first_time, 'med_phases_history_write_avg_wav', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_phases_history_write_aux_comp(gcomp, compwav, &
         first_time, 'med_phases_history_write_aux_wav', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (first_time) first_time = .false.
  end subroutine med_phases_history_write_wav

  !===============================================================================
  subroutine med_phases_history_write_inst_comp(gcomp, compid, first_time, subname, rc)
    ! Write instantaneous mediator history file for component compid

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(in)    :: compid
    logical             , intent(in)    :: first_time
    character(len=*)    , intent(in)    :: subname
    integer             , intent(out)   :: rc

    ! local variables
    character(CL)               :: alarmname
    type(ESMF_Clock)            :: mclock
    type(ESMF_Alarm)            :: alarm
    type(ESMF_Time)             :: CurrTime
    type(ESMF_Time)             :: StartTime
    type(ESMF_TimeInterval)     :: timestep
    integer                     :: timestep_length
    character(CL)               :: hist_option   ! freq_option setting (ndays, nsteps, etc)
    integer                     :: hist_n        ! freq_n setting relative to freq_option
    character(CL)               :: hist_option_in
    character(CL)               :: hist_n_in
    character(CL)               :: cvalue        ! attribute string
    logical                     :: isPresent
    logical                     :: isSet
    !---------------------------------------
    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    alarmname =  'alarm_history_inst_'//trim(compname(compid))

    if (first_time) then

       ! Determine attribute prefix
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

       if (hist_option /= 'none' .and. hist_option /= 'never') then
          ! First create hclock from mclock - THIS CALL DOES NOT COPY ALARMS
          call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          hclock_inst_comp(compid) = ESMF_ClockCreate(mclock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Set alarm for instantaneous history output
          ! Advance history clock to trigger alarms then reset history clock back to mcurrtime
          call ESMF_ClockGet(hclock_inst_comp(compid), startTime=StartTime,  currTime=CurrTime, timeStep=timestep, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_time_alarmInit(hclock_inst_comp(compid), alarm, option=hist_option, opt_n=hist_n, &
               reftime=StartTime, alarmname=trim(alarmname), rc=rc)
          call ESMF_AlarmSet(alarm, clock=hclock_inst_comp(compid), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockAdvance(hclock_inst_comp(compid),rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockSet(hclock_inst_comp(compid), currTime=currtime)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Write diagnostic info
          if (mastertask) then
             call ESMF_TimeIntervalGet(timestep, s=timestep_length, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             write(logunit,'(a,2x,i8)') "    initialized instantaneous history alarm "//&
                  trim(alarmname)//"  with option "//trim(hist_option)//" and frequency ",hist_n
             write(logunit,'(a,2x,i8)') "    history clock timestep = ",timestep_length
          end if
       end if
    end if

    if (ESMF_ClockIsCreated(hclock_inst_comp(compid))) then
       ! Advance the clock
       call ESMF_ClockAdvance(hclock_inst_comp(compid), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Write the instantaneous history file
       call med_phases_history_write_hfile(gcomp, trim(compname(compid)), hclock_inst_comp(compid), &
            trim(alarmname), .false., rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write_inst_comp

  !===============================================================================
  subroutine med_phases_history_write_avg_comp(gcomp, compid, first_time, subname, rc)

    ! Write mediator average history file variables for component compid

    use med_constants_mod, only : czero => med_constants_czero

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(in)    :: compid
    logical             , intent(in)    :: first_time
    character(len=*)    , intent(in)    :: subname
    integer             , intent(out)   :: rc
    ! local variables
    integer                 :: n
    character(CL)           :: alarmname
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: mclock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: CurrTime
    type(ESMF_Time)         :: StartTime
    type(ESMF_TimeInterval) :: timestep
    integer                 :: timestep_length
    character(CL)           :: cvalue        ! attribute string
    character(CL)           :: hist_option   ! freq_option setting (ndays, nsteps, etc)
    integer                 :: hist_n        ! freq_n setting relative to freq_option
    character(CL)           :: hist_option_in
    character(CL)           :: hist_n_in
    logical                 :: isPresent
    logical                 :: isSet
    !---------------------------------------

    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    ! Set alarm name
    alarmname =  'alarm_history_avg_'//trim(compname(compid))

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (first_time) then

       ! Determine attribute prefix
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
       ! Create time average field bundles (module variables)
       if (hist_option /= 'never' .and. hist_option /= 'none') then

          ! First create hclock from mclock - THIS CALL DOES NOT COPY ALARMS
          call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          hclock_avg_comp(compid) = ESMF_ClockCreate(mclock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Set alarm for time averaged history output
          ! Advance history clock to trigger alarms then reset history clock back to mcurrtime
          call ESMF_ClockGet(hclock_avg_comp(compid), startTime=StartTime,  currTime=CurrTime, timeStep=timestep, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_time_alarmInit(hclock_avg_comp(compid), alarm, option=hist_option, opt_n=hist_n, &
               reftime=StartTime, alarmname=trim(alarmname), rc=rc)
          call ESMF_AlarmSet(alarm, clock=hclock_avg_comp(compid), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockAdvance(hclock_avg_comp(compid),rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockSet(hclock_avg_comp(compid), currTime=currtime)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (compid /= compmed) then ! component is not mediator
             ! create accumulated import and export field bundles
             call med_phases_history_init_fldbun_accum(is_local%wrap%FBimp(compid,compid), &
                  is_local%wrap%flds_scalar_name, avgfiles_import(compid)%FBaccum, avgfiles_import(compid)%accumcnt, rc=rc)
             call med_phases_history_init_fldbun_accum(is_local%wrap%FBExp(compid), &
                  is_local%wrap%flds_scalar_name, avgfiles_export(compid)%FBaccum, avgfiles_export(compid)%accumcnt, rc=rc)
          else ! component is mediator
             ! create accumulated atm/ocn and ocnalb field bundles
             call med_phases_history_init_fldbun_accum(is_local%wrap%FBMed_aoflux_o, &
                  is_local%wrap%flds_scalar_name, avgfiles_aoflux_ocn%FBaccum, avgfiles_aoflux_ocn%accumcnt, rc=rc)
             call med_phases_history_init_fldbun_accum(is_local%wrap%FBMed_aoflux_a, &
                  is_local%wrap%flds_scalar_name, avgfiles_aoflux_atm%FBaccum, avgfiles_aoflux_atm%accumcnt, rc=rc)
             call med_phases_history_init_fldbun_accum(is_local%wrap%FBMed_ocnalb_o, &
                  is_local%wrap%flds_scalar_name, avgfiles_ocnalb_ocn%FBaccum, avgfiles_ocnalb_ocn%accumcnt, rc=rc)
             call med_phases_history_init_fldbun_accum(is_local%wrap%FBMed_ocnalb_a, &
                  is_local%wrap%flds_scalar_name, avgfiles_ocnalb_atm%FBaccum, avgfiles_ocnalb_atm%accumcnt, rc=rc)
          end if
       end if

    end if ! end of initialization (first_time) if block

    if (ESMF_ClockIsCreated(hclock_avg_comp(compid))) then
       ! Update clock
       call ESMF_ClockAdvance(hclock_avg_comp(compid), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Write history file
       call med_phases_history_write_hfile(gcomp, trim(compname(compid)), hclock_avg_comp(compid), &
            trim(alarmname), .true., rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write_avg_comp

  !===============================================================================
  subroutine med_phases_history_write_aux_comp(gcomp, compid, first_time, subname, rc)

    ! -----------------------------
    ! Write mediator history file for component compid
    ! Initialize auxiliary history file
    ! Each time this routine is called the routine SetRunClock in med.F90 is called
    ! at the beginning and the mediator clock current time and time step is set to the
    ! driver current time and time step
    ! -----------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(in)    :: compid
    logical             , intent(in)    :: first_time
    character(len=*)    , intent(in)    :: subname
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: mclock      ! mediator clock
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: currtime
    type(ESMF_TimeInterval) :: timestep
    type(ESMF_Alarm)        :: alarm
    logical                 :: isPresent   ! is attribute present
    logical                 :: isSet       ! is attribute set
    character(CL)           :: hist_option ! freq_option setting (ndays, nsteps, etc)
    integer                 :: hist_n      ! freq_n setting relative to freq_option
    integer                 :: nfcnt
    integer                 :: nfile
    integer                 :: nfld
    integer                 :: n,n1
    character(CL)           :: prefix
    character(CL)           :: cvalue
    character(CL)           :: auxflds
    integer                 :: fieldCount
    logical                 :: found
    logical                 :: enable_auxfile
    character(CS), allocatable  :: fieldNameList(:)
    !---------------------------------------
    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    if (first_time) then

       ! Get the internal state
       nullify(is_local%wrap)
       call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
             read(cvalue,'(l)') enable_auxfile
          else
             enable_auxfile = .false.
          end if

          ! If file will be written - then initialize auxfiles(nfcnt,compid)
          if (enable_auxfile) then

             ! Increment nfcnt
             nfcnt = nfcnt + 1

             ! Determine number of time samples per file
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_ntperfile', value=cvalue, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             read(cvalue,*) auxfiles(nfcnt,compid)%ntperfile
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! Determine if will do time average for aux file
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_useavg', value=cvalue, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             read(cvalue,*) auxfiles(nfcnt,compid)%useavg

             ! Determine the colon delimited field names for this file
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_flds', value=auxflds, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! Determine fields that will be output to auxhist files
             if (trim(auxflds) == 'all') then

                ! Output all fields sent to the mediator from ncomp to the auxhist files
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(compid,compid), &
                     fieldCount=fieldCount, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                allocate(auxfiles(nfcnt,compid)%flds(fieldcount))
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(compid,compid), &
                     fieldNameList=auxfiles(nfcnt,compid)%flds, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

             else

                ! Translate the colon deliminted string (auxflds) into a character array (fieldnamelist)
                ! Note that the following call allocates the memory for fieldnamelist
                call med_phases_history_get_auxflds(auxflds, fieldnamelist, rc)

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

                ! Create auxfiles(nfcnt,compid)%flds array
                allocate(auxfiles(nfcnt,compid)%flds(fieldcount))
                do n = 1,fieldcount
                   auxfiles(nfcnt,compid)%flds(n) = trim(fieldnamelist(n))
                end do

                ! Deallocate memory from fieldnamelist
                deallocate(fieldnamelist) ! this was allocated in med_phases_history_get_auxflds

             end if ! end of if auxflds is set to 'all'

             if (mastertask) then
                write(logunit,*)
                write(logunit,'(a,i4,a)') '   Writing the following fields to auxfile ',nfcnt,&
                     ' for component '//trim(compname(compid))
                do nfld = 1,size(auxfiles(nfcnt,compid)%flds)
                   write(logunit,'(8x,a)') trim(auxfiles(nfcnt,compid)%flds(nfld))
                end do
             end if

             ! Create FBaccum if averaging is on
             if (auxfiles(nfcnt,compid)%useavg) then

                ! First duplicate all fields in FBImp(compid,compid)
                call ESMF_LogWrite(trim(subname)// ": initializing FBaccum(compid)", ESMF_LOGMSG_INFO)
                call  med_phases_history_init_fldbun_accum(is_local%wrap%FBImp(compid,compid), &
                     is_local%wrap%flds_scalar_name, auxfiles(nfcnt,compid)%FBaccum, auxfiles(nfcnt,compid)%accumcnt, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! Now remove all fields from FBAccum that are not in the input flds list
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(compid,compid), fieldCount=fieldCount, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                allocate(fieldNameList(fieldCount))
                call ESMF_FieldBundleGet(is_local%wrap%FBImp(compid,compid), fieldNameList=fieldNameList, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                do n = 1,size(fieldnamelist)
                   found = .false.
                   do n1 = 1,size(auxfiles(nfcnt,compid)%flds)
                      if (trim(fieldnamelist(n)) == trim(auxfiles(nfcnt,compid)%flds(n1))) then
                         found = .true.
                         exit
                      end if
                   end do
                   if (.not. found) then
                      call ESMF_FieldBundleRemove(auxfiles(nfcnt,compid)%FBaccum, fieldnamelist(n:n), rc=rc)
                      if (chkerr(rc,__LINE__,u_FILE_u)) return
                   end if
                end do
                deallocate(fieldnameList)

                ! Check that FBAccum has at least one field left - if not exit
                call ESMF_FieldBundleGet(auxfiles(nfcnt,compid)%FBAccum, fieldCount=nfld, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                if (nfld == 0) then
                   call ESMF_LogWrite(subname//'FBAccum is zero for '//trim(auxfiles(nfcnt,compid)%auxname), &
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

             ! First create hclock from mclock - THIS CALL DOES NOT COPY ALARMS
             call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             auxfiles(nfcnt,compid)%hclock = ESMF_ClockCreate(mclock, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! Set alarm for auxiliary history output
             ! Advance history clock to trigger alarms then reset history clock back to mcurrtime
             call NUOPC_CompAttributeGet(gcomp, name=trim(prefix)//'_auxname', &
                  value=auxfiles(nfcnt,compid)%auxname, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             write(auxfiles(nfcnt,compid)%alarmname,'(a,i0)') 'alarm_'//trim(prefix)
             call ESMF_ClockGet(auxfiles(nfcnt,compid)%hclock, &
                  startTime=starttime,  currTime=currtime, timeStep=timestep, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_time_alarmInit(auxfiles(nfcnt,compid)%hclock, alarm, option=hist_option, opt_n=hist_n, &
                  reftime=starttime, alarmname=trim(auxfiles(nfcnt,compid)%alarmname), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_AlarmSet(alarm, clock=auxfiles(nfcnt,compid)%hclock, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_ClockAdvance(auxfiles(nfcnt,compid)%hclock, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_ClockSet(auxfiles(nfcnt,compid)%hclock, currtime=currtime)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if (mastertask) then
                write(logunit,'(a,2x,i8)') "    created auxiliary history alarm "//&
                     trim(auxfiles(nfcnt,compid)%alarmname)//"  with option "//trim(hist_option)//" and frequency ",hist_n
             end if
          end if ! end of isPresent and isSet and  if flag is on for file n
       end do ! end of loop over nfile

       ! Set number of aux files for this component
       num_auxfiles(compid) = nfcnt

    end if ! end of initialization (first time) block

    ! Write auxiliary history files for component compid
    do n = 1,num_auxfiles(compid)
       call ESMF_ClockAdvance(auxfiles(n,compid)%hclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_phases_history_write_hfileaux(gcomp, n, compid, auxfiles(n,compid), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write_aux_comp

  !===============================================================================
  subroutine med_phases_history_write_hfile(gcomp, comptype, hclock, alarmname, doavg, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    character(len=*)    , intent(in)    :: comptype
    type(ESMF_Clock)    , intent(in)    :: hclock
    character(len=*)    , intent(in)    :: alarmname
    logical             , intent(in)    :: doavg
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(ESMF_VM)           :: vm
    type(ESMF_Clock)        :: mclock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: nexttime
    type(ESMF_Calendar)     :: calendar          ! calendar type
    type(ESMF_TimeInterval) :: timediff(2)       ! time bounds upper and lower relative to start
    type(ESMF_TimeInterval) :: ringInterval      ! alarm interval
    real(r8)                :: tbnds(2)          ! CF1.0 time bounds
    integer                 :: i,j,m,n
    integer                 :: nx,ny             ! global grid size
    character(CL)           :: time_units        ! units of time variable
    character(CL)           :: hist_file
    real(r8)                :: days_since        ! Time interval since reference time
    real(r8)                :: avg_time          ! Time coordinate output
    logical                 :: whead,wdata       ! for writing restart/history cdf files
    integer                 :: iam
    logical                 :: write_now
    integer                 :: yr,mon,day,sec    ! time units
    character(len=*), parameter :: subname='(med_phases_history_write_hfile)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the history file alarm and determine if alarm is ringing
    ! call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetAlarm(hclock, alarmname=trim(alarmname), alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) then
          write(logunit,*)'DEBUG: alarm ',trim(alarmname),' is ringing'
       end if
       ! Set write_now flag
       write_now = .true.
       ! Turn ringer off
       call ESMF_AlarmRingerOff(alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! Write diagnostic output
       if (debug_alarms) then
          call med_phases_history_output_alarminfo(hclock, alarm, alarmname, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    else
       write_now = .false.
    end if

    ! Accumulate if alarm is not on and then average if write_now flag is true
    if (doavg) then
       do n = 1,ncomps
          if (comptype == 'all' .or. comptype == trim(compname(n))) then
             if (mastertask) then
                write(logunit,*)'DEBUG: write_now ',write_now,' for comp ' ,trim(compname(n))
             end if
             if (ESMF_FieldBundleIsCreated(avgfiles_import(n)%FBaccum)) then
                call med_phases_history_fldbun_accum(is_local%wrap%FBImp(n,n), avgfiles_import(n)%FBaccum, &
                     avgfiles_import(n)%accumcnt,  rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                if (write_now) then
                   call med_phases_history_fldbun_average(avgfiles_import(n)%FBaccum, &
                        avgfiles_import(n)%accumcnt, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             end if
             if (ESMF_FieldBundleIsCreated(avgfiles_export(n)%FBaccum)) then
                call med_phases_history_fldbun_accum(is_local%wrap%FBExp(n), avgfiles_export(n)%FBaccum, &
                     avgfiles_export(n)%accumcnt,  rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                if (write_now) then
                   call med_phases_history_fldbun_average(avgfiles_export(n)%FBaccum, &
                        avgfiles_export(n)%accumcnt, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             end if
          end if
       end do
    end if

    ! Check if history alarm is ringing - and if so write the mediator history file
    if (write_now) then

       ! Determine history file name and time units
       call med_phases_history_get_filename(gcomp, doavg, comptype, hist_file, time_units, days_since, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Set tbnds and avg_time if doing averaging
       if (doavg) then
          call ESMF_ClockGet(hclock, currtime=currtime, starttime=starttime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ClockGetNextTime(hclock, nextTime=nexttime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AlarmGet(alarm, ringInterval=ringInterval, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          timediff(2) = nexttime - starttime
          timediff(1) = nexttime - ringinterval - starttime
          call ESMF_TimeIntervalGet(timediff(2), d_r8=tbnds(2), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeIntervalGet(timediff(1), d_r8=tbnds(1), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          avg_time = 0.5_r8 * (tbnds(1) + tbnds(2))
       end if

       ! Create history file
       call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMGet(vm, localPet=iam, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
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

          ! Write time values (tbnds does not appear in instantaneous output)
          call ESMF_ClockGet(hclock, calendar=calendar, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (doavg) then
             call med_io_write(hist_file, iam, time_units=time_units, calendar=calendar, time_val=avg_time, &
                  nt=1, tbnds=tbnds, whead=whead, wdata=wdata, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call med_io_write(hist_file, iam, time_units=time_units, calendar=calendar, time_val=days_since, &
                  nt=1, whead=whead, wdata=wdata, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Write import and export field bundles
          do n = 2,ncomps ! skip the mediator here
             if (comptype == 'all' .or. comptype == trim(compname(n))) then
                if (is_local%wrap%comp_present(n)) then
                   nx = is_local%wrap%nx(n)
                   ny = is_local%wrap%ny(n)
                   if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
                      if (doavg) then
                         call med_io_write(hist_file, iam, avgfiles_import(n)%FBaccum, &
                              nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, &
                              pre=trim(compname(n))//'Imp', rc=rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      else
                         call med_io_write(hist_file, iam, is_local%wrap%FBimp(n,n), &
                              nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, &
                              pre=trim(compname(n))//'Imp', rc=rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      end if
                   endif
                   if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(n),rc=rc)) then
                      if (doavg) then
                         call med_io_write(hist_file, iam, avgfiles_export(n)%FBaccum, &
                              nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, &
                              pre=trim(compname(n))//'Exp', rc=rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      else
                         call med_io_write(hist_file, iam, is_local%wrap%FBexp(n), &
                              nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, &
                              pre=trim(compname(n))//'Exp', rc=rc)
                         if (ChkErr(rc,__LINE__,u_FILE_u)) return
                      end if
                   endif
                endif
             end if
          enddo

          ! Write mediator fractions
          ! Also write atm/ocn fluxes and ocean albedoes if field bundles are created
          if (.not. doavg) then
             if (comptype == 'all' .or. comptype == 'med') then
                do n = 2,ncomps ! skip the mediator here
                   if (ESMF_FieldBundleIsCreated(is_local%wrap%FBFrac(n),rc=rc)) then
                      call med_io_write(hist_file, iam, is_local%wrap%FBFrac(n), &
                           nx=is_local%wrap%nx(n), ny=is_local%wrap%ny(n), nt=1, whead=whead, wdata=wdata, &
                           pre='Med_frac_'//trim(compname(n)), rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   end if
                end do
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
                   call med_io_write(hist_file, iam, is_local%wrap%FBMed_ocnalb_o, &
                        nx=is_local%wrap%nx(compocn), ny=is_local%wrap%ny(compocn), nt=1, whead=whead, wdata=wdata, &
                        pre='Med_alb_ocn', rc=rc)
                end if
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o,rc=rc)) then
                   call med_io_write(hist_file, iam, is_local%wrap%FBMed_aoflux_o, &
                        nx=is_local%wrap%nx(compocn), ny=is_local%wrap%ny(compocn), nt=1, whead=whead, wdata=wdata, &
                        pre='Med_aoflux_ocn', rc=rc)
                end if
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_a,rc=rc)) then
                   call med_io_write(hist_file, iam, is_local%wrap%FBMed_ocnalb_a, &
                        nx=is_local%wrap%nx(compatm), ny=is_local%wrap%ny(compatm), nt=1, whead=whead, wdata=wdata, &
                        pre='Med_alb_atm', rc=rc)
                end if
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a,rc=rc)) then
                   call med_io_write(hist_file, iam, is_local%wrap%FBMed_aoflux_a, &
                        nx=is_local%wrap%nx(compatm), ny=is_local%wrap%ny(compatm), nt=1, whead=whead, wdata=wdata, &
                        pre='Med_aoflux_atm', rc=rc)
                end if
             end if
          end if

       end do ! end of loop over m

       ! Close file
       call med_io_close(hist_file, iam, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if ! end of if-alarm is ringingblock

  end subroutine med_phases_history_write_hfile

  !===============================================================================
  subroutine med_phases_history_write_hfileaux(gcomp, nfile_index, comp_index, auxfile, rc)

    use med_constants_mod, only : czero => med_constants_czero

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    integer             , intent(in)    :: nfile_index
    integer             , intent(in)    :: comp_index
    type(auxfile_type)  , intent(inout) :: auxfile
    integer             , intent(out)   :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(ESMF_VM)           :: vm
    type(ESMF_Clock)        :: mclock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: nexttime
    type(ESMF_Calendar)     :: calendar          ! calendar type
    type(ESMF_TimeInterval) :: timediff(2)       ! time bounds upper and lower relative to start
    type(ESMF_TimeInterval) :: ringInterval      ! alarm interval
    character(CS)           :: timestr           ! yr-mon-day-sec string
    character(CL)           :: time_units        ! units of time variable
    real(r8)                :: avg_time          ! Time coordinate output
    integer                 :: nx,ny             ! global grid size
    logical                 :: whead,wdata       ! for writing restart/history cdf files
    logical                 :: write_now         ! if true, write time sample to file
    integer                 :: iam               ! mpi task
    integer                 :: start_ymd         ! Starting date YYYYMMDD
    integer                 :: yr,mon,day,sec    ! time units
    real(r8)                :: tbnds(2)          ! CF1.0 time bounds
    character(len=*), parameter :: subname='(med_phases_history_write_hfileaux)'
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
    ! Use nexttime rather than currtime for the time difference form
    ! start since that is the time at the end of the time step
    call ESMF_ClockGet(auxfile%hclock, currtime=currtime, starttime=starttime, calendar=calendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetNextTime(auxfile%hclock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetAlarm(auxfile%hclock, alarmname=trim(auxfile%alarmname), alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    write_now = .false.
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       write_now = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (debug_alarms) then
          call med_phases_history_output_alarminfo(auxfile%hclock, alarm, auxfile%alarmname, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_AlarmGet(alarm, ringInterval=ringInterval, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       timediff(2) = currtime - starttime
       timediff(1) = currtime - starttime - ringinterval
       call ESMF_TimeIntervalGet(timediff(2), d_r8=tbnds(2), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeIntervalGet(timediff(1), d_r8=tbnds(1), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       avg_time = 0.5_r8 * (tbnds(1) + tbnds(2))
    end if

    if (mastertask .and. debug_alarms) then
       if (write_now) then
          write(logunit,'(a)')' alarmname = '//trim(auxfile%alarmname)//' is ringing'
          write(logunit,'(a,f13.5,a,f13.5)')' tbnds(1) = ',tbnds(1),' tbnds(2) = ',tbnds(2)
       else
          write(logunit,'(a)')' alarmname = '//trim(auxfile%alarmname)//' is not ringing'
       end if
       call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,'(a,4(i6,2x))')' nexttime is ',yr,mon,day,sec
       call ESMF_TimeGet(currtime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,'(a,4(i6,2x))')' currtime is ',yr,mon,day,sec
       call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,'(a,4(i6,2x))')' starttime is ',yr,mon,day,sec
    end if

    ! Do accumulation and average if required
    if (auxfile%useavg) then
       call med_phases_history_fldbun_accum(is_local%wrap%FBImp(comp_index,comp_index), auxfile%FBaccum, &
            auxfile%accumcnt, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (write_now) then
          call med_phases_history_fldbun_average(auxfile%FBaccum, auxfile%accumcnt, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    end if

    ! Write time sample to file
    if ( write_now ) then

       ! Increment number of time samples on file
       auxfile%nt = auxfile%nt + 1

       ! Set shorthand variables
       nx = is_local%wrap%nx(comp_index)
       ny = is_local%wrap%ny(comp_index)

       ! Write  header
       if (auxfile%nt == 1) then

          ! determine history file name
          call ESMF_ClockGet(auxfile%hclock, currtime=currtime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          write(timestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
          if (trim(case_name) == 'unset') then
             call med_phases_history_set_casename(gcomp, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          write(auxfile%histfile, "(8a)") &
               trim(case_name),'.cpl',trim(inst_tag),'.hx.', trim(auxfile%auxname),'.',trim(timestr), '.nc'

          ! open file
          call med_io_wopen(auxfile%histfile, vm, iam, file_ind=nfile_index, clobber=.true.)

          ! define time units
          call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_ymd2date(yr,mon,day,start_ymd)
          time_units = 'days since ' // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(sec, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! define time variables
          call med_io_write(auxfile%histfile, iam, time_units, calendar, avg_time, &
               nt=auxfile%nt, tbnds=tbnds, whead=.true., wdata=.false., file_ind=nfile_index, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! define data variables with a time dimension (include the nt argument below)
          call med_io_write(auxfile%histfile, iam, is_local%wrap%FBimp(comp_index,comp_index), &
               nx=nx, ny=ny, nt=auxfile%nt, whead=.true., wdata=.false., pre=trim(compname(comp_index))//'Imp', &
               flds=auxfile%flds, file_ind=nfile_index, use_float=.true., rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! end definition phase
          call med_io_enddef(auxfile%histfile, file_ind=nfile_index)

       end if

       ! Write time variables for time nt
       call med_io_write(auxfile%histfile, iam, time_units, calendar, avg_time, &
            nt=auxfile%nt, tbnds=tbnds, whead=.false., wdata=.true., file_ind=nfile_index, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Write data variables for time nt
       if (auxfile%useavg) then
          call med_io_write(auxfile%histfile, iam, auxfile%FBaccum, &
               nx=nx, ny=ny, nt=auxfile%nt, whead=.false., wdata=.true., pre=trim(compname(comp_index))//'Imp', &
               flds=auxfile%flds, file_ind=nfile_index, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_methods_FB_reset(auxfile%FBaccum, value=czero, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call med_io_write(auxfile%histfile, iam, is_local%wrap%FBimp(comp_index,comp_index), &
               nx=nx, ny=ny, nt=auxfile%nt, whead=.false., wdata=.true., pre=trim(compname(comp_index))//'Imp', &
               flds=auxfile%flds, file_ind=nfile_index, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Close file
       if (auxfile%nt == auxfile%ntperfile) then
          call med_io_close(auxfile%histfile, iam, file_ind=nfile_index,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          auxfile%nt = 0
       end if

    end if ! end of write_now if-block

  end subroutine med_phases_history_write_hfileaux

  !===============================================================================
  subroutine med_phases_history_get_filename(gcomp, doavg, comptype, hist_file, time_units, days_since, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    logical             , intent(in)    :: doavg
    character(len=*)    , intent(in)    :: comptype
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
    character(len=CS)       :: histstr
    character(len=*), parameter :: subname='(med_phases_history_get_filename)'
    !---------------------------------------

    rc = ESMF_SUCCESS

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
    if (doavg) then
       histstr = 'ha.'
    else
       histstr = 'hi.'
    end if
    if (trim(comptype) /= 'all') then
       histstr = trim(histstr) // trim(comptype) // '.'
    end if
    if (trim(case_name) == 'unset') then
       call med_phases_history_set_casename(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    write(hist_file,"(6a)") trim(case_name),'.cpl.',trim(inst_tag),trim(histstr),trim(nexttimestr),'.nc'
    if (mastertask) then
       write(logunit,*)
       write(logunit,' (a)') "writing mediator history file "//trim(hist_file)
       write(logunit,' (a)') "currtime = "//trim(currtimestr)//" nexttime = "//trim(nexttimestr)
    end if

  end subroutine med_phases_history_get_filename

  !===============================================================================
  subroutine med_phases_history_get_auxflds(str, flds, rc)

    ! input/output variables
    character(len=*)               , intent(in)  :: str     ! colon deliminted string to search
    character(len=*) , allocatable , intent(out) :: flds(:) ! memory will be allocate for flds
    integer                        , intent(out) :: rc

    ! local variables
    integer          :: i,k,n ! generic indecies
    integer          :: nflds ! allocatable size of flds
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
       if (mastertask) write(logunit,*) "ERROR: invalid list = ",trim(str)
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

  end subroutine med_phases_history_get_auxflds

  !===============================================================================
  subroutine med_phases_history_output_alarminfo(mclock, alarm, alarmname, rc)

    ! input/output variables
    type(ESMF_Clock), intent(in)  :: mclock
    type(ESMF_Alarm), intent(in)  :: alarm
    character(len=*), intent(in)  :: alarmname
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: ringInterval_length
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: nexttime
    character(len=CS)       :: currtimestr
    character(len=CS)       :: nexttimestr
    integer                 :: yr,mon,day,sec ! time units
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
       write(logunit,'(a,i8)') trim(subname)//": history alarmname "//trim(alarmname)//&
            ' is ringing, interval length is ', ringInterval_length
       write(logunit,'(a)') trim(subname)//": currtime = "//trim(currtimestr)//" nexttime = "//trim(nexttimestr)
    end if

  end subroutine med_phases_history_output_alarminfo

  !===============================================================================
  subroutine med_phases_history_ymds2rday_offset(currtime, rdays_offset, &
       years_offset, months_offset, days_offset, seconds_offset, rc)

    ! Given the current time and optional year, month, day and seconds offsets
    ! from the current time: Return an offset from the current time given in fractional days.
    ! For example, if day_offset = -2 and seconds_offset = -21600, rday_offset will be -2.25.
    ! One or more of the following optional arguments should be provided:

    ! input/output variables
    type(ESMF_Time) , intent(in)           :: currtime       ! current time
    real(r8)        , intent(out)          :: rdays_offset   ! offset from current time in fractional days
    integer         , intent(in), optional :: years_offset   ! number of years offset from current time
    integer         , intent(in), optional :: months_offset  ! number of months offset from current time
    integer         , intent(in), optional :: days_offset    ! number of days offset from current time
    integer         , intent(in), optional :: seconds_offset ! number of seconds offset from current time
    integer         , intent(out)          :: rc

    ! local variables
    type(ESMF_TimeInterval) :: timeinterval
    !---------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_TimeIntervalSet(timeinterval=timeinterval, startTime=currtime, &
         YY=years_offset, MM=months_offset, D=days_offset, S=seconds_offset, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet(timeinterval=timeinterval, d_r8=rdays_offset, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_history_ymds2rday_offset

  !===============================================================================
  subroutine med_phases_history_set_casename(gcomp, rc)

    ! Set module variables case_name and inst_tag

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    logical :: isPresent
    logical :: isSet
    !---------------------------------------
    rc = ESMF_SUCCESS

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

  end subroutine med_phases_history_set_casename

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
    integer                :: fieldCount
    character(CL), pointer :: fieldnames(:) => null()
    real(r8), pointer      :: dataptr1d(:) => null()
    real(r8), pointer      :: dataptr2d(:,:) => null()
    real(r8), pointer      :: dataptr1d_accum(:) => null()
    real(r8), pointer      :: dataptr2d_accum(:,:) => null()
    integer                :: ungriddedUBound(1)
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Accumulate field
    call ESMF_FieldBundleGet(fldbun_accum, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldnames(fieldCount))
    call ESMF_FieldBundleGet(fldbun_accum, fieldNameList=fieldnames, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1, fieldcount
       call ESMF_FieldBundleGet(fldbun, fieldName=trim(fieldnames(n)), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(fldbun_accum, fieldName=trim(fieldnames(n)), field=lfield_accum, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
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
    deallocate(fieldnames)

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
    character(CL), pointer :: fieldnames(:) => null()
    real(r8), pointer      :: dataptr1d_accum(:) => null()
    real(r8), pointer      :: dataptr2d_accum(:,:) => null()
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
  subroutine med_phases_history_init_fldbun_accum(fldbun, scalar_name, fldbun_accum, count, rc)

    use ESMF              , only : ESMF_FieldBundleIsCreated
    use med_constants_mod , only : czero => med_constants_czero
    use med_methods_mod   , only : med_methods_FB_init
    use med_methods_mod   , only : med_methods_FB_reset

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: fldbun
    character(len=*)       , intent(in)    :: scalar_name
    type(ESMF_FieldBundle) , intent(inout) :: fldbun_accum
    integer                , intent(out)   :: count
    integer                , intent(out)   :: rc
    !---------------------------------------

    rc = ESMF_SUCCESS

    if (ESMF_FieldBundleIsCreated(fldbun) .and. .not. ESMF_FieldBundleIsCreated(fldbun_accum)) then
       call med_methods_FB_init(fldbun_accum, scalar_name, FBgeom=fldbun, FBflds=fldbun, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call med_methods_FB_reset(fldbun_accum, czero, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       count = 0
    end if

  end subroutine med_phases_history_init_fldbun_accum

end module med_phases_history_mod
