module driver_pio_mod
  use pio
  use shr_pio_mod,  only : io_compname, pio_comp_settings, iosystems, io_compid, shr_pio_getindex
  use shr_kind_mod, only : CS=>shr_kind_CS, shr_kind_cl, shr_kind_in
  use shr_file_mod, only : shr_file_getunit, shr_file_freeunit
  use shr_log_mod,  only : shr_log_unit
  use shr_mpi_mod,  only : shr_mpi_bcast, shr_mpi_chkerr
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
  public :: driver_pio_log_comp_settings

  integer :: io_comm
  integer :: pio_debug_level=0, pio_blocksize=0
  integer(kind=pio_offset_kind) :: pio_buffer_size_limit=-1

  type(pio_rearr_opt_t) :: pio_rearr_opts
  logical, allocatable :: pio_async_interface(:)

  integer :: total_comps
  logical :: mastertask
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
    character(len=CS) :: pio_rearr_comm_type, pio_rearr_comm_fcd
    character(CS) :: msgstr

    character(*), parameter :: subName = '(driver_pio_init) '

    call ESMF_GridCompGet(driver, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    mastertask = (localPet == 0)

    call NUOPC_CompAttributeGet(driver, name="pio_buffer_size_limit", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname,*) pio_buffer_size_limit

    ! 0 is a valid value of pio_buffer_size_limit
    if(pio_buffer_size_limit>=0) then
       if(mastertask) write(shr_log_unit,*) 'Setting pio_buffer_size_limit : ',pio_buffer_size_limit
       call pio_set_buffer_size_limit(pio_buffer_size_limit)
    endif

    call NUOPC_CompAttributeGet(driver, name="pio_blocksize", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname, *) pio_blocksize
    
    if(pio_blocksize>0) then
       if(mastertask) write(shr_log_unit,*) 'Setting pio_blocksize : ',pio_blocksize
       call pio_set_blocksize(pio_blocksize)
    endif

    call NUOPC_CompAttributeGet(driver, name="pio_debug_level", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname, *) pio_debug_level

    if(pio_debug_level > 0) then
       if(mastertask) write(shr_log_unit,*) 'Setting pio_debug_level : ',pio_debug_level
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

    pio_rearr_opts%comm_fc_opts_comp2io%enable_hs = (trim(cname) .eq. '.true.')

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_enable_hs_io2comp", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    pio_rearr_opts%comm_fc_opts_io2comp%enable_hs = (trim(cname) .eq. '.true.')

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_enable_isend_comp2io", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    pio_rearr_opts%comm_fc_opts_comp2io%enable_isend = (trim(cname) .eq. '.true.')

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_enable_isend_io2comp", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    pio_rearr_opts%comm_fc_opts_io2comp%enable_isend = (trim(cname) .eq. '.true.')

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_max_pend_req_comp2io", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname, *) pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req

    call NUOPC_CompAttributeGet(driver, name="pio_rearr_comm_max_pend_req_io2comp", value=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cname, *) pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req

    if(mastertask) then
       ! Log the rearranger options
       write(shr_log_unit, *) "PIO rearranger options:"
       write(shr_log_unit, *) "  comm type     = ", pio_rearr_opts%comm_type, " (",trim(pio_rearr_comm_type),")"
       write(shr_log_unit, *) "  comm fcd      = ", pio_rearr_opts%fcd, " (",trim(pio_rearr_comm_fcd),")"
       if(pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req == PIO_REARR_COMM_UNLIMITED_PEND_REQ) then
          write(shr_log_unit, *) "  max pend req (comp2io)  = PIO_REARR_COMM_UNLIMITED_PEND_REQ (-1)"
       else
          write(shr_log_unit, *) "  max pend req (comp2io)  = ", pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req
       end if
       write(shr_log_unit, *) "  enable_hs (comp2io)     = ", pio_rearr_opts%comm_fc_opts_comp2io%enable_hs
       write(shr_log_unit, *) "  enable_isend (comp2io)  = ", pio_rearr_opts%comm_fc_opts_comp2io%enable_isend
       if(pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req == PIO_REARR_COMM_UNLIMITED_PEND_REQ) then
          write(shr_log_unit, *) "  max pend req (io2comp)  = PIO_REARR_COMM_UNLIMITED_PEND_REQ (-1)"
       else
          write(shr_log_unit, *) "  max pend req (io2comp)  = ", pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req
       end if
       write(shr_log_unit, *) "  enable_hs (io2comp)    = ", pio_rearr_opts%comm_fc_opts_io2comp%enable_hs
       write(shr_log_unit, *) "  enable_isend (io2comp)  = ", pio_rearr_opts%comm_fc_opts_io2comp%enable_isend
    end if

  end subroutine driver_pio_init

  subroutine driver_pio_component_init(driver, ncomps, rc)
    use ESMF, only : ESMF_GridComp, ESMF_LogSetError, ESMF_RC_NOT_VALID, ESMF_GridCompIsCreated, ESMF_VM, ESMF_VMGet
    use ESMF, only : ESMF_GridCompGet, ESMF_GridCompIsPetLocal, ESMF_VMIsCreated
    use NUOPC, only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use NUOPC_Driver, only : NUOPC_DriverGetComp

    type(ESMF_GridComp) :: driver
    type(ESMF_VM) :: vm
    integer, intent(in) :: ncomps
    integer, intent(out) :: rc

    integer :: i, npets, default_stride
    integer :: j
    integer :: comp_comm, comp_rank
    type(ESMF_GridComp), pointer :: gcomp(:)
    character(CS) :: cval
    character(CS) :: msgstr
    integer :: do_async_init
    type(iosystem_desc_t), allocatable :: async_iosystems(:)

    allocate(pio_comp_settings(ncomps))
    allocate(gcomp(ncomps))

    allocate(io_compid(ncomps))
    allocate(io_compname(ncomps))
    allocate(iosystems(ncomps))

    allocate(pio_async_interface(ncomps))

    nullify(gcomp)
    do_async_init = 0

    call NUOPC_DriverGetComp(driver, compList=gcomp, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    total_comps = size(gcomp)
    
    do i=1,total_comps
       io_compid(i) = i+1

       if (ESMF_GridCompIsPetLocal(gcomp(i), rc=rc)) then
          call ESMF_GridCompGet(gcomp(i), vm=vm, name=cval, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          io_compname(i) = trim(cval)
          
          call NUOPC_CompAttributeAdd(gcomp(i), attrList=(/'MCTID'/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          
          write(cval, *) io_compid(i)
          call NUOPC_CompAttributeSet(gcomp(i), name="MCTID", value=trim(cval), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          call ESMF_VMGet(vm, mpiCommunicator=comp_comm, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if(comp_comm .ne. MPI_COMM_NULL) then
             call ESMF_VMGet(vm, petCount=npets, localPet=comp_rank, ssiLocalPetCount=default_stride, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          
             call NUOPC_CompAttributeGet(gcomp(i), name="pio_stride", value=cval, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             read(cval, *) pio_comp_settings(i)%pio_stride
             if(pio_comp_settings(i)%pio_stride <= 0 .or. pio_comp_settings(i)%pio_stride > npets) then
                pio_comp_settings(i)%pio_stride = min(npets, default_stride)
             endif
             
             call NUOPC_CompAttributeGet(gcomp(i), name="pio_rearranger", value=cval, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             read(cval, *) pio_comp_settings(i)%pio_rearranger
             
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
             
             call NUOPC_CompAttributeGet(gcomp(i), name="pio_async_interface", value=cval, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             pio_async_interface(i) = (trim(cval) == '.true.')
             
             call NUOPC_CompAttributeGet(gcomp(i), name="pio_netcdf_format", value=cval, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call driver_pio_getioformatfromname(cval, pio_comp_settings(i)%pio_netcdf_ioformat, PIO_64BIT_DATA)
             
             if (pio_async_interface(i)) then
                do_async_init = do_async_init + 1
             else
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
          endif
       endif
    enddo
    if (do_async_init > 0) then
       allocate(async_iosystems(do_async_init))
       j=1
       do i=1,total_comps
          if(pio_async_interface(i)) then
             iosystems(i) = async_iosystems(j)
             j = j+1
          endif
       enddo
       
    endif

    deallocate(gcomp)
  end subroutine driver_pio_component_init

  subroutine driver_pio_log_comp_settings(gcomp, logunit)
    use ESMF, only : ESMF_GridComp, ESMF_GridCompGet
    use NUOPC, only: NUOPC_CompAttributeGet

    type(ESMF_GridComp) :: gcomp
    integer, intent(in) :: logunit

    integer :: compid
    character(len=CS) :: name, cval
    integer :: i
    integer :: rc
    logical :: isPresent

    call ESMF_GridCompGet(gcomp, name=name, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="MCTID", value=cval, isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if(isPresent) then
       read(cval, *) compid
       i = shr_pio_getindex(compid)
    endif
    write(logunit,*) trim(name),': PIO numiotasks=',   pio_comp_settings(i)%pio_numiotasks
    
    write(logunit, *) trim(name), ': PIO stride=',pio_comp_settings(i)%pio_stride
    
    write(logunit, *) trim(name),': PIO rearranger=',pio_comp_settings(i)%pio_rearranger

    write(logunit, *) trim(name),': PIO root=',pio_comp_settings(i)%pio_root
        
  end subroutine driver_pio_log_comp_settings

!===============================================================================
  subroutine driver_pio_finalize(  )
    integer :: ierr
    integer :: i
    do i=1,total_comps
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
       write(shr_log_unit,*) 'driver_pio_mod: WARNING Bad io_type argument - using iotype_netcdf'
       iotype=pio_iotype_netcdf
    end if

  end subroutine driver_pio_getiotypefromname

!===============================================================================

end module driver_pio_mod
