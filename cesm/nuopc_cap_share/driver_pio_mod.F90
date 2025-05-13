module driver_pio_mod
  use pio         , only : pio_offset_kind, pio_rearr_opt_t, PIO_REARR_COMM_UNLIMITED_PEND_REQ
  use pio         , only : pio_iotype_netcdf, pio_iotype_pnetcdf, pio_iotype_netcdf4c, pio_iotype_netcdf4p
  use pio         , only : iosystem_desc_t, PIO_64BIT_DATA, PIO_64BIT_OFFSET, PIO_REARR_COMM_COLL
  use pio         , only : PIO_REARR_COMM_P2P, pio_init, pio_set_log_level
  use pio         , only : pio_set_blocksize, pio_set_buffer_size_limit, pio_finalize
  use shr_pio_mod,  only : io_compname, pio_comp_settings, iosystems, io_compid, shr_pio_getindex
  use shr_kind_mod, only : CS=>shr_kind_CS, shr_kind_cl, shr_kind_in
  use shr_log_mod,  only : shr_log_getLogUnit
  use shr_sys_mod,  only : shr_sys_abort
#ifndef NO_MPI2
  use mpi, only : mpi_comm_null, mpi_comm_world, mpi_finalize
#endif
  use esm_utils_mod, only : chkerr
  implicit none
#ifdef NO_MPI2
#include <mpif.h>
#endif
  private
  public :: driver_pio_init
  public :: driver_pio_component_init
  public :: driver_pio_finalize
  private :: driver_pio_log_comp_settings

  integer :: io_comm
  integer :: pio_debug_level=0, pio_blocksize=0
  integer(kind=pio_offset_kind) :: pio_buffer_size_limit=-1

  type(pio_rearr_opt_t) :: pio_rearr_opts
  logical               :: pio_async_interface

  integer :: total_comps
  logical :: maintask
#define DEBUGI 1

#ifdef DEBUGI
  integer :: drank
#endif

  character(*), parameter :: u_FILE_u = &
       __FILE__

contains

!>
!! @public
!! @brief if pio_async_interface is true, tasks in io_comm do not return from this subroutine.
!!
!! if pio_async_interface is false each component namelist pio_inparm is read from compname_modelio.nml
!! Then a subset of each components compute tasks are Identified as IO tasks using the root, stride and count
!! variables to select the tasks.
!!
!<

  subroutine driver_pio_init(driver, rc)
    use ESMF, only : ESMF_GridComp, ESMF_VM, ESMF_Config, ESMF_GridCompGet
    use ESMF, only : ESMF_VMGet, ESMF_RC_NOT_VALID, ESMF_LogSetError
    use NUOPC, only: NUOPC_CompAttributeGet
    use shr_string_mod, only : shr_string_toLower
    type(ESMF_GridComp)            :: driver
    integer, intent(out) :: rc

    type(ESMF_VM)     :: vm
    integer :: i
    character(len=shr_kind_cl) :: nlfilename, cname
    integer :: ret
    integer :: localPet
    integer :: logunit
    character(len=CS) :: pio_rearr_comm_type, pio_rearr_comm_fcd
    character(CS) :: msgstr

    character(*), parameter :: subName = '(driver_pio_init) '
    
    call shr_log_getLogUnit(logunit)
    call ESMF_GridCompGet(driver, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    maintask = (localPet == 0)

    call NUOPC_CompAttributeGet(driver, name="pio_buffer_size_limit", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname,*) pio_buffer_size_limit

    ! 0 is a valid value of pio_buffer_size_limit
    if(pio_buffer_size_limit>=0) then
       if(maintask) write(logunit,*) 'Setting pio_buffer_size_limit : ',pio_buffer_size_limit
       call pio_set_buffer_size_limit(pio_buffer_size_limit)
    endif

    call NUOPC_CompAttributeGet(driver, name="pio_blocksize", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname, *) pio_blocksize
    
    if(pio_blocksize>0) then
       if(maintask) write(logunit,*) 'Setting pio_blocksize : ',pio_blocksize
       call pio_set_blocksize(pio_blocksize)
    endif

    call NUOPC_CompAttributeGet(driver, name="pio_debug_level", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname, *) pio_debug_level

    if(pio_debug_level > 0) then
       if(maintask) write(logunit,*) 'Setting pio_debug_level : ',pio_debug_level
       ret = pio_set_log_level(pio_debug_level)
    endif
       
    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_type", value=pio_rearr_comm_type, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if(trim(pio_rearr_comm_type) .eq. 'p2p') then
       pio_rearr_opts%comm_type = PIO_REARR_COMM_P2P
    else
       pio_rearr_opts%comm_type = PIO_REARR_COMM_COLL
    endif
    
    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_fcd", value=pio_rearr_comm_fcd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    
    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_enable_hs_comp2io", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_enable_hs_comp2io", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    pio_rearr_opts%comm_fc_opts_comp2io%enable_hs = logical((trim(cname) .eq. '.true.'), kind=1)

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_enable_hs_io2comp", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    pio_rearr_opts%comm_fc_opts_io2comp%enable_hs = logical((trim(cname) .eq. '.true.'), kind=1)

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_enable_isend_comp2io", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    pio_rearr_opts%comm_fc_opts_comp2io%enable_isend = logical((trim(cname) .eq. '.true.'), kind=1)

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_enable_isend_io2comp", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    pio_rearr_opts%comm_fc_opts_io2comp%enable_isend = logical((trim(cname) .eq. '.true.'), kind=1)

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_max_pend_req_comp2io", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname, *) pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_max_pend_req_io2comp", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname, *) pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req

    if(maintask) then
       ! Log the rearranger options
       write(logunit, *) "PIO rearranger options:"
       write(logunit, *) "  comm type     = ", pio_rearr_opts%comm_type, " (",trim(pio_rearr_comm_type),")"
       write(logunit, *) "  comm fcd      = ", pio_rearr_opts%fcd, " (",trim(pio_rearr_comm_fcd),")"
       if(pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req == PIO_REARR_COMM_UNLIMITED_PEND_REQ) then
          write(logunit, *) "  max pend req (comp2io)  = PIO_REARR_COMM_UNLIMITED_PEND_REQ (-1)"
       else
          write(logunit, *) "  max pend req (comp2io)  = ", pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req
       end if
       write(logunit, *) "  enable_hs (comp2io)     = ", pio_rearr_opts%comm_fc_opts_comp2io%enable_hs
       write(logunit, *) "  enable_isend (comp2io)  = ", pio_rearr_opts%comm_fc_opts_comp2io%enable_isend
       if(pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req == PIO_REARR_COMM_UNLIMITED_PEND_REQ) then
          write(logunit, *) "  max pend req (io2comp)  = PIO_REARR_COMM_UNLIMITED_PEND_REQ (-1)"
       else
          write(logunit, *) "  max pend req (io2comp)  = ", pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req
       end if
       write(logunit, *) "  enable_hs (io2comp)    = ", pio_rearr_opts%comm_fc_opts_io2comp%enable_hs
       write(logunit, *) "  enable_isend (io2comp)  = ", pio_rearr_opts%comm_fc_opts_io2comp%enable_isend
    end if

  end subroutine driver_pio_init

  subroutine driver_pio_component_init(driver, inst_comm, asyncio_petlist, rc)
    use ESMF, only : ESMF_GridComp, ESMF_LogSetError, ESMF_RC_NOT_VALID, ESMF_GridCompIsCreated, ESMF_VM, ESMF_VMGet
    use ESMF, only : ESMF_GridCompGet, ESMF_GridCompIsPetLocal, ESMF_VMIsCreated, ESMF_Finalize, ESMF_PtrInt1D
    use ESMF, only : ESMF_LOGMSG_INFO, ESMF_LOGWRITE
    use NUOPC, only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use NUOPC_Driver, only : NUOPC_DriverGetComp
    use mpi, only :  MPI_INTEGER, MPI_MAX, MPI_IN_PLACE, MPI_LOR, MPI_LOGICAL

    type(ESMF_GridComp) :: driver
    integer, intent(in) :: asyncio_petlist(:) 
    integer, intent(in) :: Inst_comm ! The communicator associated with the driver
    integer, intent(out) :: rc

    type(ESMF_VM) :: vm
    integer :: i, npets, default_stride
    integer :: j, myid
    integer :: k
    integer :: comp_comm, comp_rank
    integer, allocatable :: procs_per_comp(:), async_procs_per_comp(:)
    integer, allocatable :: io_proc_list(:), asyncio_tasks(:), comp_proc_list(:,:)

    type(ESMF_GridComp), pointer :: gcomp(:)

    character(CS) :: cval
    character(CS) :: msgstr
    integer :: do_async_init
    integer :: totalpes
    integer :: asyncio_ntasks
    integer :: asyncio_stride
    integer :: pecnt
    integer :: ierr
    integer :: iocomm
    integer :: pp
    integer :: async_rearr
    integer :: maxprocspercomp, driver_myid
    integer, allocatable :: driverpetlist(:)
    integer, allocatable :: asyncio_comp_comm(:)
    integer :: logunit
    integer :: ioproc
    integer :: n
    logical :: asyncio_task
    logical, allocatable :: petlocal(:)
    type(ESMF_PtrInt1D), pointer :: petLists(:)
    type(iosystem_desc_t), allocatable :: async_iosystems(:)
    character(len=*), parameter :: subname = '('//__FILE__//':shr_pio_component_init)'

    asyncio_ntasks = size(asyncio_petlist)

    call shr_log_getLogUnit(logunit)
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call MPI_Comm_rank(Inst_comm, myid, rc)
    call MPI_Comm_size(Inst_comm, totalpes, rc)

    asyncio_task=.false.

    do i=1,asyncio_ntasks
       ! asyncio_petlist is in 
       if(asyncio_petlist(i) == myid) then
          asyncio_task = .true.
          exit
       endif
    enddo
    write(msgstr,*) 'asyncio_task = ', asyncio_task, myid, asyncio_petlist
    call ESMF_LogWrite(trim(subname)//msgstr, ESMF_LOGMSG_INFO, rc=rc)
    nullify(gcomp)
    nullify(petLists)
    if (.not. asyncio_task) then
       call ESMF_GridCompGet(gridcomp=driver, vm=vm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMGet(vm, localPet=driver_myid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_DriverGetComp(driver, compList=gcomp, petLists=petLists, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif
    if(associated(gcomp)) then
       total_comps = size(gcomp)       
    else
       total_comps = 0
    endif

    call ESMF_LogWrite(trim(subname)//": share total_comps and driverpecount", ESMF_LOGMSG_INFO)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if(totalpes > 1) then
       call MPI_AllReduce(MPI_IN_PLACE, total_comps, 1, MPI_INTEGER, &
            MPI_MAX, Inst_comm, rc)
    endif

    allocate(pio_comp_settings(total_comps))
    allocate(procs_per_comp(total_comps))
    allocate(io_compid(total_comps))
    allocate(io_compname(total_comps))
    allocate(iosystems(total_comps))
    allocate(petlocal(total_comps))
    do_async_init = 0
    procs_per_comp = 0

    do i=1,total_comps
       if(associated(gcomp)) then
          petlocal(i) = ESMF_GridCompIsPetLocal(gcomp(i), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          call NUOPC_CompAttributeGet(gcomp(i), name="pio_async_interface", value=cval, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          pio_comp_settings(i)%pio_async_interface = (trim(cval) == '.true.')
          
          call NUOPC_CompAttributeGet(gcomp(i), name="pio_rearranger", value=cval, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          read(cval, *) pio_comp_settings(i)%pio_rearranger          
       else
          petlocal(i) = .false.
       endif
       pio_comp_settings(i)%pio_async_interface = .false.
       io_compid(i) = i+1

       if (petlocal(i)) then
          call NUOPC_CompAttributeAdd(gcomp(i), attrList=(/'MCTID'/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          write(cval, *) io_compid(i)
          call NUOPC_CompAttributeSet(gcomp(i), name="MCTID", value=trim(cval), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          call ESMF_GridCompGet(gcomp(i), vm=vm, name=cval, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite(trim(subname)//": initialize component: "//trim(cval), ESMF_LOGMSG_INFO)
          io_compname(i) = trim(cval)

          call ESMF_VMGet(vm, mpiCommunicator=comp_comm, localPet=comp_rank, petCount=npets, &
               ssiLocalPetCount=default_stride, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          procs_per_comp(i) = npets

          if(.not. pio_comp_settings(i)%pio_async_interface) then
             call NUOPC_CompAttributeGet(gcomp(i), name="pio_stride", value=cval, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             read(cval, *) pio_comp_settings(i)%pio_stride
             if(pio_comp_settings(i)%pio_stride <= 0 .or. pio_comp_settings(i)%pio_stride > npets) then
                pio_comp_settings(i)%pio_stride = min(npets, default_stride)
             endif
          
             call NUOPC_CompAttributeGet(gcomp(i), name="pio_numiotasks", value=cval, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             read(cval, *) pio_comp_settings(i)%pio_numiotasks
          
             if(pio_comp_settings(i)%pio_numiotasks < 0 .or. pio_comp_settings(i)%pio_numiotasks > npets) then
                pio_comp_settings(i)%pio_numiotasks = max(1,npets/pio_comp_settings(i)%pio_stride)
             endif


             call NUOPC_CompAttributeGet(gcomp(i), name="pio_root", value=cval, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             read(cval, *) pio_comp_settings(i)%pio_root
          
             if(pio_comp_settings(i)%pio_root < 0 .or. pio_comp_settings(i)%pio_root > npets) then
                pio_comp_settings(i)%pio_root = 0
             endif
          endif
         
          call NUOPC_CompAttributeGet(gcomp(i), name="pio_typename", value=cval, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          select case (trim(cval))
          case ('pnetcdf')
             pio_comp_settings(i)%pio_iotype = PIO_IOTYPE_PNETCDF
          case ('netcdf')
             pio_comp_settings(i)%pio_iotype = PIO_IOTYPE_NETCDF
          case ('netcdf4p')
             pio_comp_settings(i)%pio_iotype = PIO_IOTYPE_NETCDF4P
          case ('netcdf4c')
             pio_comp_settings(i)%pio_iotype = PIO_IOTYPE_NETCDF4C
          case DEFAULT
             write (msgstr, *) "Invalid PIO_TYPENAME Setting for component ", trim(cval)
             call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
             return
          end select
                       
          call NUOPC_CompAttributeGet(gcomp(i), name="pio_netcdf_format", value=cval, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call driver_pio_getioformatfromname(cval, pio_comp_settings(i)%pio_netcdf_ioformat, PIO_64BIT_DATA)
          
          if (.not. pio_comp_settings(i)%pio_async_interface) then
             if(pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req < PIO_REARR_COMM_UNLIMITED_PEND_REQ) then
                pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req = pio_comp_settings(i)%pio_numiotasks
             endif
             if(pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req < PIO_REARR_COMM_UNLIMITED_PEND_REQ) then
                pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req = pio_comp_settings(i)%pio_numiotasks
             endif
             
             call pio_init(comp_rank ,comp_comm ,pio_comp_settings(i)%pio_numiotasks, 0, pio_comp_settings(i)%pio_stride, &
                  pio_comp_settings(i)%pio_rearranger, iosystems(i), pio_comp_settings(i)%pio_root, &
                  pio_rearr_opts)
          endif
          ! Write the PIO settings to the beggining of each component log
          if(comp_rank == 0) call driver_pio_log_comp_settings(gcomp(i), rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

       endif
    enddo

    call ESMF_LogWrite(trim(subname)//": check for async", ESMF_LOGMSG_INFO)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do i=1,total_comps
       call MPI_AllReduce(MPI_IN_PLACE, pio_comp_settings(i)%pio_async_interface, 1, MPI_LOGICAL, &
            MPI_LOR, Inst_comm, rc)
       if(pio_comp_settings(i)%pio_async_interface) then
          do_async_init = do_async_init + 1
       endif
    enddo

!
!   Get the PET list for each component using async IO
!

    call MPI_Allreduce(MPI_IN_PLACE, do_async_init, 1, MPI_INTEGER, MPI_MAX, Inst_comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, procs_per_comp, total_comps, MPI_INTEGER, MPI_MAX, Inst_comm, ierr)
    if (do_async_init > 0) then
       maxprocspercomp = 0
       do i=1,total_comps
          if(procs_per_comp(i) > maxprocspercomp) maxprocspercomp = procs_per_comp(i)
       enddo
       call MPI_AllReduce(MPI_IN_PLACE, maxprocspercomp, 1, MPI_INTEGER, &
            MPI_MAX, Inst_comm, rc)

       allocate(asyncio_comp_comm(do_async_init))
       allocate(comp_proc_list(maxprocspercomp, do_async_init))
       j = 1
       k = 1
       comp_proc_list = -1
       if(.not. asyncio_task) then
          do i=1,total_comps
             if(pio_comp_settings(i)%pio_async_interface) then
                comp_proc_list(1:procs_per_comp(i), j) = petLists(i)%ptr
                ! IO tasks are not in the driver comp so we need to correct the comp_proc_list
                do k=1,size(asyncio_petlist)
                  ioproc = asyncio_petlist(k)
                  do n=1,procs_per_comp(i)
                     if(petLists(i)%ptr(n) >= (ioproc-k+1)) comp_proc_list(n,j) = comp_proc_list(n,j) + 1
                  enddo
                enddo
                j = j+1
             endif
!             deallocate(petLists(i)%ptr)
          enddo
       endif
       ! Copy comp_proc_list to io tasks
       do i=1,do_async_init
          call MPI_AllReduce(MPI_IN_PLACE, comp_proc_list(:,i), maxprocspercomp, MPI_INTEGER, MPI_MAX, Inst_comm, ierr)
       enddo
       if(asyncio_ntasks == 0) then
          call shr_sys_abort(subname//' ERROR: ASYNC IO Requested but no IO PES assigned')
       endif

       allocate(async_iosystems(do_async_init))
       allocate(async_procs_per_comp(do_async_init))
       j=1
       async_rearr = 0
       do i=1,total_comps
          if(pio_comp_settings(i)%pio_async_interface) then
             async_procs_per_comp(j) = procs_per_comp(i)
             j = j+1
             if(.not.asyncio_task) then
                if(async_rearr == 0) then
                   async_rearr = pio_comp_settings(i)%pio_rearranger
                elseif(async_rearr .ne. pio_comp_settings(i)%pio_rearranger .and. pio_comp_settings(i)%pio_rearranger > 0) then
                   write(msgstr,*) i,async_rearr,pio_comp_settings(i)%pio_rearranger
                   call shr_sys_abort(subname//' ERROR: all async component rearrangers must match '//msgstr)
                endif
             endif
          endif
       enddo
       
       ! IO tasks should not return until the run is completed
       !ierr = pio_set_log_level(1)
       call ESMF_LogWrite(trim(subname)//": call async pio_init", ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call MPI_AllReduce(MPI_IN_PLACE, async_rearr, 1, MPI_INTEGER, &
            MPI_MAX, Inst_comm, rc)
       call pio_init(async_iosystems, Inst_comm, async_procs_per_comp, &
            comp_proc_list, asyncio_petlist, &
            async_rearr, asyncio_comp_comm, io_comm)
       if(.not. asyncio_task) then
          j=1
          do i=1,total_comps
             if(pio_comp_settings(i)%pio_async_interface) then
                iosystems(i) = async_iosystems(j)
                j = j+1
             endif
          enddo
       endif
    endif
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    if(associated(petLists)) deallocate(petLists)
    if(associated(gcomp)) deallocate(gcomp)
  end subroutine driver_pio_component_init

  subroutine driver_pio_log_comp_settings(gcomp, rc)
    use ESMF, only : ESMF_GridComp, ESMF_GridCompGet, ESMF_SUCCESS
    use NUOPC, only: NUOPC_CompAttributeGet
    use, intrinsic :: iso_fortran_env, only: output_unit
 
    type(ESMF_GridComp) :: gcomp
    integer, intent(out) :: rc
    integer :: compid
    character(len=CS) :: name, cval
    integer :: i
    integer :: logunit
    logical :: isPresent

    rc = ESMF_SUCCESS
    call ESMF_GridCompGet(gcomp, name=name, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="MCTID", value=cval, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if(isPresent) then
       read(cval, *) compid
       i = shr_pio_getindex(compid)
    endif

    logunit = 6
    call NUOPC_CompAttributeGet(gcomp, name="logunit", value=logunit, isPresent=ispresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if(.not. isPresent) then
       logunit = output_unit
       if(maintask) write(logunit,*) 'Attribute logunit not set for ',trim(name)
    endif
    if(pio_comp_settings(i)%pio_async_interface) then
       write(logunit,*) trim(name),': using ASYNC IO interface'
    else
       write(logunit,*) trim(name),': PIO numiotasks=',   pio_comp_settings(i)%pio_numiotasks
       write(logunit, *) trim(name), ': PIO stride=',pio_comp_settings(i)%pio_stride
       write(logunit, *) trim(name),': PIO rearranger=',pio_comp_settings(i)%pio_rearranger
       write(logunit, *) trim(name),': PIO root=',pio_comp_settings(i)%pio_root
    endif
  end subroutine driver_pio_log_comp_settings

!===============================================================================
  subroutine driver_pio_finalize(  )
    integer :: ierr
    integer :: i
    do i=1,size(iosystems)
       call pio_finalize(iosystems(i), ierr)
    end do

  end subroutine driver_pio_finalize

!===============================================================================

  subroutine driver_pio_getioformatfromname(pio_netcdf_format, pio_netcdf_ioformat, pio_default_netcdf_ioformat)
    use shr_string_mod, only : shr_string_toupper
    character(len=*), intent(inout) :: pio_netcdf_format
    integer, intent(out) :: pio_netcdf_ioformat
    integer, intent(in) :: pio_default_netcdf_ioformat

    pio_netcdf_format = shr_string_toupper(pio_netcdf_format)
    if ( pio_netcdf_format .eq. 'CLASSIC' ) then
       pio_netcdf_ioformat = 0
    elseif ( pio_netcdf_format .eq. '64BIT_OFFSET' ) then
       pio_netcdf_ioformat = PIO_64BIT_OFFSET
    elseif ( pio_netcdf_format .eq. '64BIT_DATA' ) then
       pio_netcdf_ioformat = PIO_64BIT_DATA
    else
       pio_netcdf_ioformat = pio_default_netcdf_ioformat
    endif

  end subroutine driver_pio_getioformatfromname


  subroutine driver_pio_getiotypefromname(typename, iotype, defaulttype)
    use shr_string_mod, only : shr_string_toupper
    character(len=*), intent(inout) :: typename
    integer, intent(out) :: iotype
    integer, intent(in) :: defaulttype

    integer :: logunit
    
    call shr_log_getLogUnit(logunit)

    typename = shr_string_toupper(typename)
    if      ( typename .eq. 'NETCDF' ) then
       iotype = pio_iotype_netcdf
    else if ( typename .eq. 'PNETCDF') then
       iotype = pio_iotype_pnetcdf
    else if ( typename .eq. 'NETCDF4P') then
       iotype = pio_iotype_netcdf4p
    else if ( typename .eq. 'NETCDF4C') then
       iotype = pio_iotype_netcdf4c
    else if ( typename .eq. 'NOTHING') then
       iotype = defaulttype
    else if ( typename .eq. 'DEFAULT') then
       iotype = defaulttype
    else
       write(logunit,*) 'driver_pio_mod: WARNING Bad io_type argument - using iotype_netcdf'
       iotype=pio_iotype_netcdf
    end if

  end subroutine driver_pio_getiotypefromname

!===============================================================================

end module driver_pio_mod
