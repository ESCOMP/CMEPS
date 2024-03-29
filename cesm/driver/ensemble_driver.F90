module Ensemble_driver

  !-----------------------------------------------------------------------------
  ! Code that creates the ensemble driver layer above the esm driver instance.
  ! The ensmeble driver is configured to run a single clock cycle in nuopc with time step
  ! length of stop_time - start_time.  It's purpose is to instantiate NINST copies of the
  ! esm driver and its components layed out concurently across mpi tasks.
  !-----------------------------------------------------------------------------

  use shr_kind_mod  , only : cl=>shr_kind_cl, cs=>shr_kind_cs, cx=>shr_kind_cx
  use shr_log_mod   , only : shr_log_setLogUnit
  use esm_utils_mod , only : maintask, logunit, chkerr

  implicit none
  private

  public  :: SetServices
  private :: SetModelServices
  private :: ensemble_finalize

  integer, allocatable :: asyncio_petlist(:)
  logical :: asyncio_task=.false.
  logical :: asyncIO_available=.false.
  integer :: number_of_members
  integer :: inst  ! ensemble instance containing this task
  character(*),parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine SetServices(ensemble_driver, rc)

    use NUOPC        , only : NUOPC_CompDerive, NUOPC_CompSpecialize, NUOPC_CompAttributeSet
    use NUOPC        , only : NUOPC_CompAttributeGet
    use NUOPC_Driver , only : driver_routine_SS             => SetServices
    use NUOPC_Driver , only : ensemble_label_SetModelServices => label_SetModelServices
    use NUOPC_Driver , only : ensemble_label_PostChildrenAdvertise => label_PostChildrenAdvertise
    use NUOPC_Driver , only : label_Finalize
    use ESMF         , only : ESMF_GridComp, ESMF_GridCompSet
    use ESMF         , only : ESMF_Config, ESMF_ConfigCreate, ESMF_ConfigLoadFile
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO

    type(ESMF_GridComp)  :: ensemble_driver
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Config) :: config
    logical           :: isPresent
    character(len=*), parameter :: subname = "(ensemble_driver.F90:SetServices)"
    !---------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    ! NUOPC_Driver registers the generic methods
    call NUOPC_CompDerive(ensemble_driver, driver_routine_SS, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(ensemble_driver, specLabel=ensemble_label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! PostChildrenAdvertise is a NUOPC specialization which happens after Advertize but before Realize
    ! We have overloaded this specialization location to initilize IO.
    ! So after all components have called Advertise but before any component calls Realize
    ! IO will be initialized and any async IO tasks will be split off to the PIO async IO driver.
    call NUOPC_CompSpecialize(ensemble_driver, specLabel=ensemble_label_PostChildrenAdvertise, &
         specRoutine=InitializeIO, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create, open and set the config
    config = ESMF_ConfigCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ConfigLoadFile(config, "nuopc.runconfig", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSet(ensemble_driver, config=config, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! NUOPC component drivers end the initialization process with an internal call to InitializeDataResolution.
    ! The ensemble_driver does not need to InitializeDataResolution and doing so will cause a hang
    ! if asyncronous IO is used.  This attribute is available after ESMF8.4.0b03 to toggle that control.
    ! Cannot use asyncIO with older ESMF versions.
    call NUOPC_CompAttributeGet(ensemble_driver, name="InitializeDataResolution", &
         isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if(isPresent) then
       call ESMF_LogWrite(trim(subname)//": setting InitializeDataResolution false", ESMF_LOGMSG_INFO)
       call NUOPC_CompAttributeSet(ensemble_driver, name="InitializeDataResolution", value="false", rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       asyncIO_available = .true.
       call ESMF_LogWrite(trim(subname)//": asyncio is available", ESMF_LOGMSG_INFO)
    endif
    ! Set a finalize method, it calls pio_finalize
    call NUOPC_CompSpecialize(ensemble_driver, specLabel=label_Finalize, &
         specRoutine=ensemble_finalize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !================================================================================

  subroutine SetModelServices(ensemble_driver, rc)

    use ESMF          , only : ESMF_GridComp, ESMF_VM, ESMF_Config, ESMF_Clock, ESMF_VMGet
    use ESMF          , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_ConfigGetAttribute
    use ESMF          , only : ESMF_ConfigGetLen, ESMF_RC_NOT_VALID, ESMF_LogFoundAllocError
    use ESMF          , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF          , only : ESMF_GridCompSet, ESMF_SUCCESS, ESMF_METHOD_INITIALIZE, ESMF_RC_ARG_BAD
    use ESMF          , only : ESMF_CalendarSetDefault
    use ESMF          , only : ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
    use NUOPC         , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use NUOPC_Driver  , only : NUOPC_DriverAddComp, NUOPC_DriverGetComp
    use esm           , only : ESMSetServices => SetServices, ReadAttributes
    use esm_time_mod  , only : esm_time_clockInit
    use perf_mod      , only : t_startf, t_stopf, t_initf
    ! input/output variables
    type(ESMF_GridComp)    :: ensemble_driver
    integer, intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_GridComp)    :: driver
    type(ESMF_Config)      :: config
    integer                :: n
    integer, pointer       :: petList(:)
    integer                :: petCount
    integer                :: localPet
    character(len=512)     :: diro
    character(len=512)     :: logfile
    logical                :: read_restart
    character(len=CS)      :: read_restart_string
    integer                :: ntasks_per_member
    integer                :: iopetcnt
    integer                :: petcnt
    logical                :: comp_task
    integer                :: pio_asyncio_ntasks
    integer                :: pio_asyncio_stride
    integer                :: pio_asyncio_rootpe
    integer                :: Global_Comm
    character(len=CL)      :: start_type     ! Type of startup
    character(len=7)       :: drvrinst
    character(len=5)       :: inst_suffix
    character(len=CX)      :: msgstr
    character(len=CL)      :: cvalue
    character(len=CL)      :: calendar
    character(len=*) , parameter :: start_type_start = "startup"
    character(len=*) , parameter :: start_type_cont  = "continue"
    character(len=*) , parameter :: start_type_brnch = "branch"
    character(len=*) , parameter :: subname = "(ensemble_driver.F90:SetModelServices)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    call ESMF_GridCompGet(ensemble_driver, config=config, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=global_comm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (localPet == 0) then
       maintask=.true.
    else
       maintask = .false.
    end if

    call t_initf('drv_in', LogPrint=.true., LogUnit=logunit, mpicom=global_comm, mastertask=maintask)
    call t_startf(subname)

    !-------------------------------------------
    ! Initialize clocks
    !-------------------------------------------

    call ReadAttributes(ensemble_driver, config, "ALLCOMP_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(ensemble_driver, config, "CLOCK_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(ensemble_driver, 'calendar', calendar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (calendar == 'NO_LEAP') then
       call ESMF_CalendarSetDefault(ESMF_CALKIND_NOLEAP, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (calendar == 'GREGORIAN') then
       call ESMF_CalendarSetDefault(ESMF_CALKIND_GREGORIAN, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       write (msgstr, *) "Only NO_LEAP and GREGORIAN calendars currently supported"
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    end if

    ! Check valid values of start type
    call NUOPC_CompAttributeGet(ensemble_driver, name="start_type", value=start_type, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if ((trim(start_type) /= start_type_start) .and.  &
        (trim(start_type) /= start_type_cont ) .and.  &
        (trim(start_type) /= start_type_brnch)) then
       write (msgstr, *) subname//': start_type invalid = '//trim(start_type)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    end if

    if (trim(start_type) == trim(start_type_cont) .or. trim(start_type) == trim(start_type_brnch)) then
       read_restart = .true.
    else
       read_restart = .false.
    endif
    write(read_restart_string,*) read_restart

    ! Add read_restart to ensemble_driver attributes
    call NUOPC_CompAttributeAdd(ensemble_driver, attrList=(/'read_restart'/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(ensemble_driver, name='read_restart', value=trim(read_restart_string), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Extract the config object from the ensemble_driver
    !-------------------------------------------

    call ReadAttributes(ensemble_driver, config, "PELAYOUT_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Determine number of ensemble members and the number of tasks per member
    !-------------------------------------------

    call NUOPC_CompAttributeGet(ensemble_driver, name="ninst", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) number_of_members

    call NUOPC_CompAttributeGet(ensemble_driver, name="pio_asyncio_ntasks", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) pio_asyncio_ntasks

    call NUOPC_CompAttributeGet(ensemble_driver, name="pio_asyncio_stride", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) pio_asyncio_stride

    call NUOPC_CompAttributeGet(ensemble_driver, name="pio_asyncio_rootpe", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) pio_asyncio_rootpe

    call ESMF_VMGet(vm, localPet=localPet, PetCount=PetCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ntasks_per_member = PetCount/number_of_members - pio_asyncio_ntasks
    if(modulo(PetCount-pio_asyncio_ntasks*number_of_members, number_of_members) .ne. 0) then
       write (msgstr,'(a,i5,a,i3,a,i3,a)') &
            "PetCount (",PetCount,") - Async IOtasks (",pio_asyncio_ntasks*number_of_members,") must be evenly divisable by number of members (",number_of_members,")"
       call ESMF_LogSetError(ESMF_RC_ARG_BAD, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    endif

    !-------------------------------------------
    ! Loop over number of ensemblel members
    !-------------------------------------------

    allocate(petList(ntasks_per_member))
    allocate(asyncio_petlist(pio_asyncio_ntasks))
    !
    ! Logic for asyncio variables is handled in cmeps buildnml.
    ! here we assume that pio_asyncio_stride and pio_asyncio_ntasks are only set
    ! if asyncio is enabled.
    !
    logunit = 6
    do inst=1,number_of_members
       petcnt=1
       iopetcnt = 1
       comp_task = .false.
       asyncio_task = .false.
       ! Determine pet list for driver instance
       if(pio_asyncio_ntasks > 0) then
          do n=pio_asyncio_rootpe,pio_asyncio_rootpe+pio_asyncio_stride*(pio_asyncio_ntasks-1),pio_asyncio_stride
             asyncio_petlist(iopetcnt) = (inst-1)*(ntasks_per_member+pio_asyncio_ntasks) + n
             if(asyncio_petlist(iopetcnt) == localPet) asyncio_task = .true.
             iopetcnt = iopetcnt+1
          enddo
          iopetcnt = 1
       endif
       do n=0,ntasks_per_member+pio_asyncio_ntasks-1
          if(pio_asyncio_ntasks > 0) then
             if( asyncio_petlist(iopetcnt)==(inst-1)*(ntasks_per_member+pio_asyncio_ntasks) + n) then
                ! Here if asyncio is true and this is an io task
                iopetcnt = iopetcnt+1
             else if(petcnt <= ntasks_per_member) then
                ! Here if this is a compute task
                petList(petcnt) = n + (inst-1)*(ntasks_per_member + pio_asyncio_ntasks)
                if (petList(petcnt) == localPet) then
                   comp_task=.true.
                endif
                petcnt = petcnt+1
             else
                msgstr = "ERROR task cannot be neither a compute task nor an asyncio task"
                call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
                return  ! bail out
             endif
          else
             ! Here if asyncio is false
             petList(petcnt) = (inst-1)*ntasks_per_member + n
             if (petList(petcnt) == localPet) comp_task=.true.
             petcnt = petcnt+1
          endif
       enddo
       if(inst == localPet/(ntasks_per_member+pio_asyncio_ntasks) + 1) then
          if(comp_task .and. asyncio_task) then
             write(msgstr,*) "ERROR task cannot be both a compute task and an asyncio task", inst, petlist
             call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
             return  ! bail out
          elseif (.not. comp_task .and. .not. asyncio_task) then
             write(msgstr,*) "ERROR task is nether a compute task nor an asyncio task", inst, petlist
             call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
             return  ! bail out
          endif
       endif
       ! Add driver instance to ensemble driver
       write(drvrinst,'(a,i4.4)') "ESM",inst
    
       call NUOPC_DriverAddComp(ensemble_driver, drvrinst, ESMSetServices, petList=petList, comp=driver, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       write(msgstr, *) ": driver added on PETS ",petlist(1),' to ',petlist(petcnt-1), comp_task, asyncio_task
       call ESMF_LogWrite(trim(subname)//msgstr)
       ! Set the driver log to the driver task 0
       if (comp_task) then
          if(number_of_members > 1) then
             call NUOPC_CompAttributeAdd(driver, attrList=(/'inst_suffix'/), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             write(inst_suffix,'(a,i4.4)') '_',inst
             call NUOPC_CompAttributeSet(driver, name='inst_suffix', value=inst_suffix, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          else
             inst_suffix = ''
          endif
          ! CESM does not use this ESMF feature and at large processor counts it can be expensive to have it on.
          call NUOPC_CompAttributeSet(driver, name="HierarchyProtocol", value="off", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          ! Set the driver instance attributes
          call NUOPC_CompAttributeAdd(driver, attrList=(/'read_restart'/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeSet(driver, name='read_restart', value=trim(read_restart_string), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          call ReadAttributes(driver, config, "CLOCK_attributes::", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          call ReadAttributes(driver, config, "DRIVER_attributes::", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          call ReadAttributes(driver, config, "DRV_modelio::", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          write(msgStr, *) trim(subname), ' instance = ',inst, 'attributes read'
          call ESMF_LogWrite(msgStr)
          if (localPet == petList(1)) then
             call NUOPC_CompAttributeGet(driver, name="diro", value=diro, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call NUOPC_CompAttributeGet(driver, name="logfile", value=logfile, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             ! Multiinstance logfile name needs a correction
             if(len_trim(inst_suffix) > 0) then
                n = index(logfile, '.')
                logfile = logfile(1:n-1)//trim(inst_suffix)//logfile(n:)
             endif
             open (newunit=logunit,file=trim(diro)//"/"//trim(logfile))
             maintask = .true.
          endif

       endif
       call shr_log_setLogUnit (logunit)
       ! Create a clock for each driver instance
       call esm_time_clockInit(ensemble_driver, driver, logunit, localpet==petList(1), rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    enddo
    inst = localPet/(ntasks_per_member+pio_asyncio_ntasks) + 1

    deallocate(petList)
    call t_stopf(subname)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine SetModelServices

  subroutine InitializeIO(ensemble_driver, rc)
    use ESMF, only: ESMF_GridComp, ESMF_LOGMSG_INFO, ESMF_LogWrite
    use ESMF, only: ESMF_SUCCESS, ESMF_VM, ESMF_GridCompGet, ESMF_VMGet
    use ESMF, only: ESMF_CONFIG, ESMF_GridCompIsPetLocal, ESMF_State, ESMF_Clock
    use NUOPC, only: NUOPC_CompAttributeGet, NUOPC_CompGet
    use NUOPC_DRIVER, only: NUOPC_DriverGetComp
    use driver_pio_mod   , only: driver_pio_init, driver_pio_component_init
#ifndef NO_MPI2
    use MPI,  only : MPI_Comm_split, MPI_UNDEFINED
#endif
    type(ESMF_GridComp) :: ensemble_driver
    type(ESMF_VM) :: ensemble_vm
    integer, intent(out) :: rc
    character(len=*), parameter :: subname = '('//__FILE__//':InitializeIO)'
    type(ESMF_GridComp), pointer :: dcomp(:)
    integer :: iam
    integer :: Global_Comm, Instance_Comm
    integer :: drv
    integer :: PetCount
    integer :: key, color, i
    type(ESMF_GridComp) :: driver
    character(len=7)       :: drvrinst
    character(len=8) :: compname

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    call shr_log_setLogUnit (logunit)

    call ESMF_GridCompGet(ensemble_driver, vm=ensemble_vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(ensemble_vm, localpet=iam, mpiCommunicator=Global_Comm, PetCount=PetCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if(number_of_members > 1) then
       color = inst
       key = modulo(iam, PetCount/number_of_members)
#ifndef NO_MPI2
       call MPI_Comm_split(Global_Comm, color, key, Instance_Comm, rc)
#endif
       do i=1,size(asyncio_petlist)
          asyncio_petList(i) = modulo(asyncio_petList(i), PetCount/number_of_members)
       enddo
    else
       Instance_Comm = Global_Comm
    endif
    write(drvrinst,'(a,i4.4)') "ESM",inst
    call NUOPC_DriverGetComp(ensemble_driver, drvrinst, comp=driver, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": call driver_pio_init "//compname, ESMF_LOGMSG_INFO)
    call driver_pio_init(driver, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": call driver_pio_component_init "//compname, ESMF_LOGMSG_INFO)
    call driver_pio_component_init(driver, Instance_Comm, asyncio_petlist, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(subname)//": driver_pio_component_init done "//compname, ESMF_LOGMSG_INFO)

    deallocate(asyncio_petlist)
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
  end subroutine InitializeIO

  subroutine ensemble_finalize(ensemble_driver, rc)
    use ESMF, only : ESMF_GridComp, ESMF_SUCCESS
    use driver_pio_mod, only: driver_pio_finalize
    type(ESMF_GridComp) :: Ensemble_driver
    integer, intent(out) :: rc
    rc = ESMF_SUCCESS
    call shr_log_setLogUnit (logunit)
    call driver_pio_finalize()

  end subroutine ensemble_finalize
end module Ensemble_driver
