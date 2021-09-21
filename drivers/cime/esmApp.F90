program esmApp

  !-----------------------------------------------------------------------------
  ! Generic ESM application driver
  !-----------------------------------------------------------------------------

  use ESMF,            only : ESMF_Initialize, ESMF_CALKIND_GREGORIAN, ESMF_LOGKIND_MULTI
  use ESMF,            only : ESMF_END_ABORT, ESMF_LogFoundError, ESMF_Finalize, ESMF_LOGERR_PASSTHRU
  use ESMF,            only : ESMF_GridCompSetServices, ESMF_GridCompFinalize, ESMF_LogSet, ESMF_LogWrite
  use ESMF,            only : ESMF_GridCompDestroy, ESMF_LOGMSG_INFO, ESMF_GridComp, ESMF_GridCompRun
  use ESMF,            only : ESMF_GridCompFinalize, ESMF_GridCompCreate, ESMF_GridCompInitialize
  use ESMF,            only : ESMF_LOGKIND_MULTI_ON_ERROR, ESMF_LogKind_Flag
  use ESMF,            only : ESMF_VMGet, ESMF_VM, ESMF_InitializePreMPI

  use mpi
  use NUOPC,           only : NUOPC_FieldDictionarySetup
  use ensemble_driver, only : SetServices
  use shr_pio_mod,     only : shr_pio_init1
  use shr_sys_mod,     only : shr_sys_abort
  !
  ! The CrayLabs SmartSim interface is provided in directory share  https://github.com/ESCOMP/CESM_share
  ! Please see file cime/tools/smartsim/README.md for a complete explanation of the CESM interface to smartsim
  ! create_smartsim_cluster is set to true if the database is using 3 or more nodes, false if its using 1 and 
  ! (in the pbs interface at least) 2 is not allowed.
  !
  use nuopc_shr_methods, only : sr_client, use_smartredis

  implicit none

  ! local variables
  integer                 :: COMP_COMM
  integer                 :: rc, urc
  type(ESMF_LogKind_Flag) :: logkindflag
  type(ESMF_GridComp)     :: ensemble_driver_comp
  logical                 :: create_esmf_pet_files = .false.
  integer                 :: iam, ier
  integer                 :: fileunit
  integer                 :: provided
  type(ESMF_VM)           :: vm
  logical                 :: create_smartsim_cluster = .false.

  namelist /debug_inparm / create_esmf_pet_files
  !
  ! The CrayLabs SmartSim interface is provided in directory share  https://github.com/ESCOMP/CESM_share
  ! Please see file cime/tools/smartsim/README.md for a complete explanation of the CESM interface to smartsim
  ! The use_smartredis variable is set in file drv_in and if true the variable sr_client is initialized in esmApp.F90
  !
  namelist /smartsim_inparm/ use_smartredis, create_smartsim_cluster
  !-----------------------------------------------------------------------------
  ! Initiallize MPI
  !-----------------------------------------------------------------------------
#ifndef NO_MPI2
  call ESMF_InitializePreMPI()
  call MPI_init_thread(MPI_THREAD_SERIALIZED, provided, rc)
#else
  call MPI_init(rc)
#endif
  COMP_COMM = MPI_COMM_WORLD

  !-----------------------------------------------------------------------------
  ! Initialize PIO
  !-----------------------------------------------------------------------------

  ! For planned future use of async io using pio2.  The IO tasks are seperated from the compute tasks here
  ! and COMP_COMM will be MPI_COMM_NULL on the IO tasks which then call shr_pio_init2 and do not return until
  ! the model completes.  All other tasks call ESMF_Initialize.  8 is the maximum number of component models
  ! supported

  call shr_pio_init1(8, "drv_in", COMP_COMM)

  !-----------------------------------------------------------------------------
  ! Initialize ESMF
  !-----------------------------------------------------------------------------

  ! by default, ESMF_LOGKIND_MULTI_ON_ERROR does not create files PET[N*].ESMF_LogFile unless there is an error
  ! if want those files, comment out the following line and uncomment the line logkindflag = ESMF_LOGKIND_MULTI
  call mpi_comm_rank(COMP_COMM, iam, ier)
  if (iam==0) then
     open(newunit=fileunit, status="old", file="drv_in")
     read(fileunit, debug_inparm, iostat=ier)
     if (ier > 0) then
	call shr_sys_abort('esmApp: error reading in debug_inparm namelist from drv_in')
     end if
     close(fileunit)
  end if
  call mpi_bcast (create_esmf_pet_files, 1, MPI_LOGICAL, 0, COMP_COMM, ier)

  if (create_esmf_pet_files) then
     logkindflag = ESMF_LOGKIND_MULTI
  else
     logkindflag = ESMF_LOGKIND_MULTI_ON_ERROR
  end if
  call ESMF_Initialize(mpiCommunicator=COMP_COMM, logkindflag=logkindflag, logappendflag=.false., &
       defaultCalkind=ESMF_CALKIND_GREGORIAN, ioUnitLBound=5001, ioUnitUBound=5101, vm=vm, rc=rc)

  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_VMGet(vm, mpiCommunicator=COMP_COMM, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_LogSet(flush=.true.)

  call ESMF_LogWrite("esmApp STARTING", ESMF_LOGMSG_INFO, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Initialize the CrayLabs SmartRedis client, a stub is provided in share if 
  ! smartsim is not used.  This client shall be used by all component models.
  !-----------------------------------------------------------------------------

  if (iam==0) then
     open(newunit=fileunit, status="old", file="drv_in")
     read(fileunit, smartsim_inparm, iostat=ier)
     if (ier > 0) then
	call shr_sys_abort('esmApp: error reading in smartsim_inparm namelist from drv_in')
     end if
     close(fileunit)
  end if
  call mpi_bcast (use_smartredis, 1, MPI_LOGICAL, 0, COMP_COMM, ier)
  call mpi_bcast (create_smartsim_cluster, 1, MPI_LOGICAL, 0, COMP_COMM, ier)
  if (use_smartredis) then
     call ESMF_Logwrite("Using SmartSim interface", ESMF_LOGMSG_INFO, rc=rc)
     call sr_client%initialize(create_smartsim_cluster)
  endif

  !-----------------------------------------------------------------------------
  ! Operate on the NUOPC Field dictionary
  !-----------------------------------------------------------------------------

  call NUOPC_FieldDictionarySetup("fd.yaml", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Create the earth system ensemble driver Component
  !-----------------------------------------------------------------------------

  ensemble_driver_comp = ESMF_GridCompCreate(name="ensemble", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! SetServices for the ensemble driver Component
  !-----------------------------------------------------------------------------

  call ESMF_GridCompSetServices(ensemble_driver_comp, SetServices, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Call Initialize for the earth system ensemble Component
  !-----------------------------------------------------------------------------

  call ESMF_GridCompInitialize(ensemble_driver_comp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Call Run  for the ensemble driver
  !-----------------------------------------------------------------------------
  call ESMF_GridCompRun(ensemble_driver_comp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Call Finalize for the ensemble driver
  ! Destroy the ensemble driver
  !-----------------------------------------------------------------------------

  call ESMF_GridCompFinalize(ensemble_driver_comp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_LogWrite("ESMF_GridCompDestroy called", ESMF_LOGMSG_INFO, rc=rc)

  ! call ESMF_LogSet(flush=.true., trace=.true., rc=rc)
  call ESMF_GridCompDestroy(ensemble_driver_comp, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_LogWrite("ESMF_GridCompDestroy finished", ESMF_LOGMSG_INFO, rc=rc)

  !-----------------------------------------------------------------------------
  ! Finalize ESMF
  !-----------------------------------------------------------------------------

  call ESMF_LogWrite("esmApp FINISHED", ESMF_LOGMSG_INFO, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_Finalize()

end program
