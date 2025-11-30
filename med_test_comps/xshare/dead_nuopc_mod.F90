module dead_nuopc_mod

  use ESMF              , only : ESMF_Gridcomp, ESMF_State, ESMF_StateGet
  use ESMF              , only : ESMF_Field, ESMF_FieldGet
  use ESMF              , only : ESMF_Clock, ESMF_Time, ESMF_TimeInterval, ESMF_Alarm
  use ESMF              , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_ClockSet, ESMF_ClockAdvance, ESMF_AlarmSet
  use ESMF              , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_METHOD_INITIALIZE
  use ESMF              , only : ESMF_FAILURE, ESMF_LOGMSG_ERROR
  use ESMF              , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMBroadcast, ESMF_VMGet
  use ESMF              , only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VmGet
  use ESMF              , only : operator(/=), operator(==), operator(+)
  use shr_kind_mod      , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod       , only : shr_sys_abort
  use shr_string_mod    , only : shr_string_withoutSuffix
  use shr_wtracers_mod  , only : WTRACERS_SUFFIX, shr_wtracers_get_initial_ratio
  use dead_methods_mod  , only : chkerr, alarmInit

  implicit none
  private

  public :: dead_read_inparms
  public :: ModelInitPhase
  public :: ModelSetRunClock
  public :: fld_list_add
  public :: fld_list_realize
  public :: set_all_export_fields

  private :: set_wtracer_field

  ! !PUBLIC DATA MEMBERS:
  integer, parameter, public :: fldname_maxlen = 128

  type fld_list_type
     character(len=fldname_maxlen) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0

     ! Water tracer fields are handled via ungridded dimensions, but we track the size of
     ! this dimension separately to better distinguish between the ungridded dimension
     ! used for water tracers vs. the ungridded dimension used for other purposes -
     ! particularly for the case of fields that have both. For fields that are not water
     ! tracer fields, num_wtracers will be 0; for fields that are water tracer fields,
     ! num_wtracers will be the number of water tracers in this simulation.
     integer :: num_wtracers = 0
  end type fld_list_type
  public :: fld_list_type

  integer, parameter, public :: fldsMax = 100
  integer                    :: dbug_flag = 0
  character(*), parameter    :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dead_read_inparms(model, inst_suffix, logunit, nxg, nyg)

    ! input/output variables
    character(len=*) , intent(in)    :: model
    character(len=*) , intent(in)    :: inst_suffix ! char string associated with instance
    integer          , intent(in)    :: logunit     ! logging unit number
    integer          , intent(out)   :: nxg         ! global dim i-direction
    integer          , intent(out)   :: nyg         ! global dim j-direction

    ! local variables
    type(ESMF_VM)           :: vm
    character(CL)           :: fileName ! generic file name
    integer                 :: nunit    ! unit number
    integer                 :: unitn    ! Unit for namelist file
    integer                 :: tmp(2)   ! array for broadcast
    integer                 :: localPet ! mpi id of current task in current context
    integer                 :: rc       ! return code
    character(*), parameter :: F00   = "('(dead_read_inparms) ',8a)"
    character(*), parameter :: F01   = "('(dead_read_inparms) ',a,a,4i8)"
    character(*), parameter :: F03   = "('(dead_read_inparms) ',a,a,i8,a)"
    character(*), parameter :: subName = "(dead_read_inpamrs) "
    !-------------------------------------------------------------------------------

    ! read the input parms (used to configure model)
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    nxg = -9999
    nyg = -9999

    if (localPet==0) then
       open(newunit=unitn, file='x'//model//'_in'//trim(inst_suffix), status='old' )
       read(unitn,*) nxg
       read(unitn,*) nyg
       close (unitn)
    endif

    tmp(1) = nxg
    tmp(2) = nyg
    call ESMF_VMBroadcast(vm, tmp, 3, 0, rc=rc)
    nxg = tmp(1)
    nyg = tmp(2)

    if (localPet==0) then
       write(logunit,*)' Read in X'//model//' input from file= x'//model//'_in'
       write(logunit,F00) model
       write(logunit,F00) model,'         Model  :  ',model
       write(logunit,F01) model,'           NGX  :  ',nxg
       write(logunit,F01) model,'           NGY  :  ',nyg
       write(logunit,F00) model,'    inst_suffix :  ',trim(inst_suffix)
       write(logunit,F00) model
    end if

  end subroutine dead_read_inparms

  !===============================================================================
  subroutine fld_list_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer                    , intent(inout) :: num
    type(fld_list_type)        , intent(inout) :: fldlist(:)
    character(len=*)           , intent(in)    :: stdname
    integer,          optional , intent(in)    :: ungridded_lbound
    integer,          optional , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(dead_nuopc_mod:fld_list_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information
    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fld_list_add

  !===============================================================================
  subroutine fld_list_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer           :: n
    type(ESMF_Field)  :: field
    character(len=80) :: stdname
    character(len=*),parameter  :: subname='(dead_nuopc_mod:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. &
                 fldlist(n)%ungridded_ubound > 0 .and. &
                 fldlist(n)%num_wtracers > 0) then
                ! This field has two ungridded dimensions: one for water tracers and one
                ! for some other purpose. The first ungridded dimension will be for water
                ! tracers.
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/1, fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%num_wtracers, fldlist(n)%ungridded_ubound/), &
                     gridToFieldMap=(/3/), rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             else if (fldlist(n)%ungridded_lbound > 0 .and. &
                      fldlist(n)%ungridded_ubound > 0) then
                ! This field has one ungridded dimension
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                     gridToFieldMap=(/2/), rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             else if (fldlist(n)%num_wtracers > 0) then
                ! This field has one ungridded dimension, for water tracers
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/1/), &
                     ungriddedUbound=(/fldlist(n)%num_wtracers/), &
                     gridToFieldMap=(/2/), rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             else
                ! This field has no ungridded dimensions
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             end if
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------

      use ESMF, only : ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(dead_nuopc_mod:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fld_list_realize

  !================================================================================
  subroutine set_all_export_fields(exportState, flds, fld_min, fld_max, lon, lat, field_setexport, rc, fld_num_save)

    ! ----------------------------------------------
    ! Set all export fields for a given component's state
    !
    ! This accepts a procedure argument for the subroutine that does the actual setting
    ! for each field, since this procedure can differ between different xcomps.
    !
    ! Water tracer fields are handled specially: these are set equal to the corresponding
    ! bulk field times the initial ratio for this tracer.
    ! ----------------------------------------------

    ! input/output arguments
    type(ESMF_State), intent(inout)  :: exportState
    type(fld_list_type), intent(in)  :: flds(:)
    integer, intent(in)              :: fld_min  ! first index in flds to set
    integer, intent(in)              :: fld_max  ! last index in flds to set
    real(r8), intent(in)             :: lon(:)
    real(r8), intent(in)             :: lat(:)
    integer, intent(out)             :: rc

    ! fld_num_save can be provided to continue where we left off from the last call.
    ! This is useful for multiple ice sheets, for example, where we want different field
    ! values for each ice sheet. It should generally be set to 1 for the initial call from
    ! a component (but could be set to some other value if desired).
    integer, optional, intent(inout) :: fld_num_save

    interface
       subroutine field_setexport(exportState, fldname, lon, lat, nf, ungridded_index, rc)
          import :: ESMF_State
          import :: r8

          type(ESMF_State), intent(inout) :: exportState
          character(len=*), intent(in)    :: fldname
          real(r8), intent(in)            :: lon(:)
          real(r8), intent(in)            :: lat(:)
          integer, intent(in)             :: nf
          integer, optional, intent(in)   :: ungridded_index
          integer, intent(out)            :: rc
       end subroutine field_setexport
    end interface

    ! local variables
    integer :: nf, nind, fld_num
    character(len=*), parameter :: subname='(dead_nuopc_mod:set_all_export_fields)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (present(fld_num_save)) then
       fld_num = fld_num_save
    else
       fld_num = 1
    end if

    do nf = fld_min,fld_max
       if (flds(nf)%num_wtracers > 0) then
          ! We'll handle water tracers specially, below. A few notes about this:
          ! - We handle water tracers after we are done setting all non-water tracer
          !   fields, because the setting of water tracer fields depends on the
          !   corresponding non-tracer fields.
          ! - We do *not* increment fld_num for the water tracer fields. This ensures that
          !   values put in the non-tracer fields remain the same even when introducing
          !   water tracers.
          cycle
       end if

       if (flds(nf)%ungridded_ubound == 0) then
          call field_setexport(exportState, trim(flds(nf)%stdname), lon, lat, nf=fld_num, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          fld_num = fld_num + 1
       else
          do nind = 1,flds(nf)%ungridded_ubound
             call field_setexport(exportState, trim(flds(nf)%stdname), lon, lat, nf=fld_num, &
                  ungridded_index=nind, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             fld_num = fld_num + 1
          end do
       end if
    end do

    ! Now handle water tracers.
    do nf = fld_min,fld_max
       if (flds(nf)%num_wtracers > 0) then
          call set_wtracer_field(exportState, flds(nf), rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    if (present(fld_num_save)) then
       fld_num_save = fld_num
    end if

  end subroutine set_all_export_fields

  !================================================================================
  subroutine set_wtracer_field(exportState, fld, rc)

    ! ----------------------------------------------
    ! Sets a single water tracer field (for all tracers), based on the corresponding bulk
    ! field and the initial ratio of this tracer.
    ! ----------------------------------------------

    ! input/output arguments
    type(ESMF_State), intent(inout) :: exportState
    type(fld_list_type), intent(in) :: fld
    integer, intent(out)            :: rc

    ! local variables
    logical :: has_suffix
    character(len=fldname_maxlen) :: wtracer_bulk_fldname

    type(ESMF_Field) :: field_wtracers
    type(ESMF_Field) :: field_bulk

    ! If there is no ungridded dimension other than the water tracer dimension, we'll use
    ! these variables:
    real(r8), pointer :: data_bulk_1d(:)
    real(r8), pointer :: data_wtracers_2d(:,:)

    ! If there is an additional ungridded dimension in addition to the water tracer
    ! dimension, we'll use these variables:
    real(r8), pointer :: data_bulk_2d(:,:)
    real(r8), pointer :: data_wtracers_3d(:,:,:)

    integer :: n

    character(len=*), parameter :: subname='(dead_nuopc_mod:set_wtracer_field)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call shr_string_withoutSuffix( &
         in_str = fld%stdname, &
         suffix = WTRACERS_SUFFIX, &
         has_suffix = has_suffix, &
         out_str = wtracer_bulk_fldname)
    if (.not. has_suffix) then
       call ESMF_LogWrite(subname//": ERROR: "//trim(fld%stdname)// &
            " does not end with the expected suffix for a water tracer field", &
            ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    call ESMF_StateGet(exportState, itemName=trim(fld%stdname), field=field_wtracers, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(exportState, itemName=trim(wtracer_bulk_fldname), field=field_bulk, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (fld%ungridded_lbound > 0 .and. fld%ungridded_ubound > 0) then
       ! There is an additional ungridded dimension in addition to the water tracer
       ! dimension. Note that we assume that the bulk field matches the tracer field in
       ! terms of the size of this ungridded dimension.
       call ESMF_FieldGet(field_wtracers, farrayPtr=data_wtracers_3d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_bulk, farrayPtr=data_bulk_2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       do n = 1, fld%num_wtracers
          data_wtracers_3d(n,:,:) = data_bulk_2d(:,:) * shr_wtracers_get_initial_ratio(n)
       end do
    else
       ! No additional ungridded dimension
       call ESMF_FieldGet(field_wtracers, farrayPtr=data_wtracers_2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_bulk, farrayPtr=data_bulk_1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       do n = 1, fld%num_wtracers
          data_wtracers_2d(n,:) = data_bulk_1d(:) * shr_wtracers_get_initial_ratio(n)
       end do
    end if

  end subroutine set_wtracer_field

  !===============================================================================
  subroutine ModelInitPhase(gcomp, importState, exportState, clock, rc)

    use NUOPC, only : NUOPC_CompFilterPhaseMap

    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelInitPhase

  !===============================================================================
  subroutine ModelSetRunClock(gcomp, rc)

    use ESMF        , only : ESMF_ClockGetAlarmList, ESMF_ALARMLIST_ALL
    use NUOPC_Model , only : NUOPC_ModelGet
    use NUOPC       , only : NUOPC_CompAttributeGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option       ! Restart option units
    integer                  :: restart_n            ! Number until restart interval
    integer                  :: restart_ymd          ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname='dead_nuopc_mod:(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart alarm
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelSetRunClock

end module dead_nuopc_mod
