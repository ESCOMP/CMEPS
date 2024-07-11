module med_io_mod

  !------------------------------------------
  ! Create mediator history files
  !------------------------------------------

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, I8=>SHR_KIND_I8, R8=>SHR_KIND_R8
  use med_kind_mod          , only : R4=>SHR_KIND_R4
  use med_constants_mod     , only : fillvalue => SHR_CONST_SPVAL
  use ESMF                  , only : ESMF_VM, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LogFoundError, ESMF_LOGMSG_ERROR
  use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_END_ABORT, ESMF_LOGERR_PASSTHRU
  use ESMF                  , only : ESMF_VMGetCurrent, ESMF_VMGet, ESMF_VMBroadCast, ESMF_Finalize
  use NUOPC                 , only : NUOPC_FieldDictionaryGetEntry
  use NUOPC                 , only : NUOPC_FieldDictionaryHasEntry
  use pio                   , only : file_desc_t, iosystem_desc_t
  use med_internalstate_mod , only : logunit, med_id, maintask
  use med_constants_mod     , only : dbug_flag    => med_constants_dbug_flag
  use med_methods_mod       , only : FB_getFieldN => med_methods_FB_getFieldN
  use med_methods_mod       , only : FB_getFldPtr => med_methods_FB_getFldPtr
  use med_methods_mod       , only : FB_getNameN  => med_methods_FB_getNameN
  use med_utils_mod         , only : chkerr       => med_utils_ChkErr

  implicit none
  private

  ! public member functions:
  public :: med_io_wopen
  public :: med_io_close
  public :: med_io_redef
  public :: med_io_enddef
  public :: med_io_sec2hms
  public :: med_io_read
  public :: med_io_define_time
  public :: med_io_write_time
  public :: med_io_write
  public :: med_io_init
  public :: med_io_date2yyyymmdd
  public :: med_io_datetod2string
  public :: med_io_ymd2date

  ! private member functions
  private :: med_io_file_exists

  ! public data members:
  interface med_io_read
     module procedure med_io_read_FB
     module procedure med_io_read_int
     module procedure med_io_read_int1d
     module procedure med_io_read_r8
     module procedure med_io_read_r81d
     module procedure med_io_read_char
  end interface med_io_read
  interface med_io_write
     module procedure med_io_write_FB
     module procedure med_io_write_int
     module procedure med_io_write_int1d
     module procedure med_io_write_r8
     module procedure med_io_write_r81d
     module procedure med_io_write_char
  end interface med_io_write
  interface med_io_date2ymd
     module procedure med_io_date2ymd_int
     module procedure med_io_date2ymd_long
  end interface med_io_date2ymd
  interface  med_io_datetod2string
     module procedure med_io_datetod2string_int
     module procedure med_io_datetod2string_long
  end interface med_io_datetod2string
  interface med_io_ymd2date
     module procedure med_io_ymd2date_int
     module procedure med_io_ymd2date_long
  end interface med_io_ymd2date

  ! module data
  character(*),parameter         :: prefix    = "med_io_"
  character(*),parameter         :: modName   = "(med_io_mod) "
  character(*),parameter         :: version   = "cmeps0"

  integer                        :: pio_iotype
  integer                        :: pio_ioformat
  type(iosystem_desc_t), pointer :: io_subsystem
  character(*),parameter         :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  logical function med_io_file_exists(vm, filename)

    !---------------
    ! inquire if i/o file exists
    !---------------

    ! input/output variables
    type(ESMF_VM)                :: vm
    character(len=*), intent(in) :: filename

    ! local variables
    integer :: tmp(1)
    integer :: iam
    integer :: rc
    !-------------------------------------------------------------------------------

    tmp(1) = 0

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    med_io_file_exists = .false.
    if (iam==0) then
       inquire(file=trim(filename),exist=med_io_file_exists)
       if (med_io_file_exists) tmp(1) = 1
    end if
    call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if(tmp(1) == 1) med_io_file_exists = .true.

  end function med_io_file_exists

  !===============================================================================
  subroutine med_io_init(gcomp, rc)

    !---------------
    ! initialize pio
    !---------------

    use ESMF , only : ESMF_GridComp, ESMF_UtilStringUpperCase
#ifdef CESMCOUPLED
    use shr_pio_mod , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
#else
    use pio  , only : pio_init, pio_setdebuglevel, pio_set_rearr_opts
    use pio  , only : PIO_64BIT_OFFSET, PIO_64BIT_DATA
    use pio  , only : PIO_IOTYPE_NETCDF, PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF4C, PIO_IOTYPE_NETCDF4P
    use pio  , only : PIO_REARR_BOX, PIO_REARR_SUBSET
    use pio  , only : PIO_REARR_COMM_P2P, PIO_REARR_COMM_COLL
    use pio  , only : PIO_REARR_COMM_FC_2D_ENABLE, PIO_REARR_COMM_FC_2D_DISABLE
    use pio  , only : PIO_REARR_COMM_FC_1D_COMP2IO, PIO_REARR_COMM_FC_1D_IO2COMP
    use NUOPC, only : NUOPC_CompAttributeGet
#endif

    ! input/output arguments
    type(ESMF_GridComp), intent(in)    :: gcomp
    integer            , intent(out)   :: rc

#ifndef CESMCOUPLED
    ! local variables
    type(ESMF_VM)           :: vm
    integer                 :: ret
    integer                 :: comm
    integer                 :: localPet, petCount
    integer                 :: pio_numiotasks
    integer                 :: pio_stride
    integer                 :: pio_rearranger
    integer                 :: pio_root
    integer                 :: pio_debug_level
    integer                 :: pio_rearr_comm_type
    integer                 :: pio_rearr_comm_fcd
    logical                 :: pio_rearr_comm_enable_hs_comp2io
    logical                 :: pio_rearr_comm_enable_isend_comp2io
    integer                 :: pio_rearr_comm_max_pend_req_comp2io
    logical                 :: pio_rearr_comm_enable_hs_io2comp
    logical                 :: pio_rearr_comm_enable_isend_io2comp
    integer                 :: pio_rearr_comm_max_pend_req_io2comp
    logical                 :: isPresent, isSet
    character(len=CS)       :: cvalue
    character(*), parameter :: subname = '(med_io_init)'
    !-------------------------------------------------------------------------------
#endif

#ifdef CESMCOUPLED
    io_subsystem => shr_pio_getiosys(med_id)
    pio_iotype   =  shr_pio_getiotype(med_id)
    pio_ioformat =  shr_pio_getioformat(med_id)
#else
    ! query VM
    call ESMF_VMGetCurrent(vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=comm, localPet=localPet, petCount=petCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! query component specific PIO attributes
    ! pio_netcdf_format
    call NUOPC_CompAttributeGet(gcomp, name='pio_netcdf_format', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'CLASSIC') then
         pio_ioformat = 0
       else if (trim(cvalue) .eq. '64BIT_OFFSET') then
         pio_ioformat = PIO_64BIT_OFFSET
       else if (trim(cvalue) .eq. '64BIT_DATA') then
         pio_ioformat = PIO_64BIT_DATA
       else
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_ioformat (CLASSIC|64BIT_OFFSET|64BIT_DATA)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       cvalue = '64BIT_OFFSET'
       pio_ioformat = PIO_64BIT_OFFSET
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_netcdf_format = ', trim(cvalue), pio_ioformat

    ! pio_typename
    call NUOPC_CompAttributeGet(gcomp, name='pio_typename', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'NETCDF') then
         pio_iotype = PIO_IOTYPE_NETCDF
       else if (trim(cvalue) .eq. 'PNETCDF') then
         pio_iotype = PIO_IOTYPE_PNETCDF
       else if (trim(cvalue) .eq. 'NETCDF4C') then
         pio_iotype = PIO_IOTYPE_NETCDF4C
       else if (trim(cvalue) .eq. 'NETCDF4P') then
         pio_iotype = PIO_IOTYPE_NETCDF4P
       else
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_typename (NETCDF|PNETCDF|NETCDF4C|NETCDF4P)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       cvalue = 'NETCDF'
       pio_iotype = PIO_IOTYPE_NETCDF
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_typename = ', trim(cvalue), pio_iotype

    ! pio_root
    call NUOPC_CompAttributeGet(gcomp, name='pio_root', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_root
       if (pio_root < 0) then
          pio_root = 1
       endif
       pio_root = min(pio_root, petCount-1)
    else
       pio_root = 1
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_root = ', pio_root

    ! pio_stride
    call NUOPC_CompAttributeGet(gcomp, name='pio_stride', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_stride
    else
       pio_stride = -99
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_stride = ', pio_stride

    ! pio_numiotasks
    call NUOPC_CompAttributeGet(gcomp, name='pio_numiotasks', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_numiotasks
    else
       pio_numiotasks = -99
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks

    ! check for parallel IO, it requires at least two io pes
    if (petCount > 1 .and. pio_numiotasks == 1 .and. &
       (pio_iotype .eq. PIO_IOTYPE_PNETCDF .or. pio_iotype .eq. PIO_IOTYPE_NETCDF4P)) then
       pio_numiotasks = 2
       pio_stride = min(pio_stride, petCount/2)
       if (localPet == 0) then
          write(logunit,*) ' parallel io requires at least two io pes - following parameters are updated:'
          write(logunit,*) trim(subname), ' : pio_stride = ', pio_stride
          write(logunit,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks
       end if
    endif

    ! check/set/correct io pio parameters
    if (pio_stride > 0 .and. pio_numiotasks < 0) then
       pio_numiotasks = max(1, petCount/pio_stride)
       if (localPet == 0) write(logunit,*) trim(subname), ' : update pio_numiotasks = ', pio_numiotasks
    else if(pio_numiotasks > 0 .and. pio_stride < 0) then
       pio_stride = max(1, petCount/pio_numiotasks)
       if (localPet == 0) write(logunit,*) trim(subname), ' : update pio_stride = ', pio_stride
    else if(pio_numiotasks < 0 .and. pio_stride < 0) then
       pio_stride = max(1,petCount/4)
       pio_numiotasks = max(1,petCount/pio_stride)
       if (localPet == 0) write(logunit,*) trim(subname), ' : update pio_numiotasks = ', pio_numiotasks
       if (localPet == 0) write(logunit,*) trim(subname), ' : update pio_stride = ', pio_stride
    end if
    if (pio_stride == 1) then
       pio_root = 0
    endif

    if (pio_root + (pio_stride)*(pio_numiotasks-1) >= petCount .or. &
       pio_stride <= 0 .or. pio_numiotasks <= 0 .or. pio_root < 0 .or. pio_root > petCount-1) then
       if (petCount < 100) then
          pio_stride = max(1, petCount/4)
       else if(petCount < 1000) then
          pio_stride = max(1, petCount/8)
       else
          pio_stride = max(1, petCount/16)
       end if
       if(pio_stride > 1) then
          pio_numiotasks = petCount/pio_stride
          pio_root = min(1, petCount-1)
       else
          pio_numiotasks = petCount
          pio_root = 0
       end if
       if (localPet == 0) then
          write(logunit,*) 'pio_stride, iotasks or root out of bounds - resetting to defaults:'
          write(logunit,*) trim(subname), ' : pio_root = ', pio_root
          write(logunit,*) trim(subname), ' : pio_stride = ', pio_stride
          write(logunit,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks
       end if
    end if

    ! pio_rearranger
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearranger', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'BOX') then
         pio_rearranger = PIO_REARR_BOX
       else if (trim(cvalue) .eq. 'SUBSET') then
         pio_rearranger = PIO_REARR_SUBSET
       else
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_rearranger (BOX|SUBSET)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       cvalue = 'SUBSET'
       pio_rearranger = PIO_REARR_SUBSET
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_rearranger = ', trim(cvalue), pio_rearranger

    ! init PIO
    allocate(io_subsystem)
    if (localPet == 0) write(logunit,*) trim(subname),' calling pio init'
    call pio_init(localPet, comm, pio_numiotasks, 0, pio_stride, pio_rearranger, io_subsystem, base=pio_root)

    ! PIO debug related options
    ! pio_debug_level
    call NUOPC_CompAttributeGet(gcomp, name='pio_debug_level', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_debug_level
       if (pio_debug_level < 0 .or. pio_debug_level > 6) then
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_debug_level (0-6)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       pio_debug_level = 0
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_debug_level = ', pio_debug_level

    ! set PIO debug level
    call pio_setdebuglevel(pio_debug_level)

    ! query shared PIO rearranger attributes
    ! pio_rearr_comm_type
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_type', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'P2P') then
          pio_rearr_comm_type = PIO_REARR_COMM_P2P
       else if (trim(cvalue) .eq. 'COLL') then
          pio_rearr_comm_type = PIO_REARR_COMM_COLL
       else
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_rearr_comm_type (P2P|COLL)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       cvalue = 'P2P'
       pio_rearr_comm_type = PIO_REARR_COMM_P2P
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_rearr_comm_type = ', trim(cvalue), pio_rearr_comm_type

    ! pio_rearr_comm_fcd
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_fcd', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. '2DENABLE') then
          pio_rearr_comm_fcd = PIO_REARR_COMM_FC_2D_ENABLE
       else if (trim(cvalue) .eq. 'IO2COMP') then
          pio_rearr_comm_fcd = PIO_REARR_COMM_FC_1D_IO2COMP
       else if (trim(cvalue) .eq. 'COMP2IO') then
          pio_rearr_comm_fcd = PIO_REARR_COMM_FC_1D_COMP2IO
       else if (trim(cvalue) .eq. '2DDISABLE') then
          pio_rearr_comm_fcd = PIO_REARR_COMM_FC_2D_DISABLE
       else
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_rearr_comm_fcd (2DENABLE|IO2COMP|COMP2IO|2DDISABLE)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       cvalue = '2DENABLE'
       pio_rearr_comm_fcd = PIO_REARR_COMM_FC_2D_ENABLE
    end if
    if (localPet == 0) write(logunit,*) trim(subname), ' : pio_rearr_comm_fcd = ', trim(cvalue), pio_rearr_comm_fcd

    ! pio_rearr_comm_enable_hs_comp2io
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_enable_hs_comp2io', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_enable_hs_comp2io
    else
       pio_rearr_comm_enable_hs_comp2io = .true.
    end if

    ! pio_rearr_comm_enable_isend_comp2io
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_enable_isend_comp2io', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_enable_isend_comp2io
    else
       pio_rearr_comm_enable_isend_comp2io = .false.
    end if

    ! pio_rearr_comm_max_pend_req_comp2io
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_max_pend_req_comp2io', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_max_pend_req_comp2io
    else
       pio_rearr_comm_max_pend_req_comp2io = 0
    end if

    ! pio_rearr_comm_enable_hs_io2comp
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_enable_hs_io2comp', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_enable_hs_io2comp
    else
       pio_rearr_comm_enable_hs_io2comp = .false.
    end if

    ! pio_rearr_comm_enable_isend_io2comp
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_enable_isend_io2comp', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_enable_isend_io2comp
    else
       pio_rearr_comm_enable_isend_io2comp = .true.
    end if

    ! pio_rearr_comm_max_pend_req_io2comp
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_max_pend_req_io2comp', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_max_pend_req_io2comp
    else
       pio_rearr_comm_max_pend_req_io2comp = 64
    end if

    ! print out PIO rearranger parameters
    if (localPet == 0) then
       write(logunit,*) trim(subname), ' : pio_rearr_comm_enable_hs_comp2io = ', pio_rearr_comm_enable_hs_comp2io
       write(logunit,*) trim(subname), ' : pio_rearr_comm_enable_isend_comp2io = ', pio_rearr_comm_enable_isend_comp2io
       write(logunit,*) trim(subname), ' : pio_rearr_comm_max_pend_req_comp2io = ', pio_rearr_comm_max_pend_req_comp2io
       write(logunit,*) trim(subname), ' : pio_rearr_comm_enable_hs_io2comp = ', pio_rearr_comm_enable_hs_io2comp
       write(logunit,*) trim(subname), ' : pio_rearr_comm_enable_isend_io2comp = ', pio_rearr_comm_enable_isend_io2comp
       write(logunit,*) trim(subname), ' : pio_rearr_comm_max_pend_req_io2comp = ', pio_rearr_comm_max_pend_req_io2comp
    end if

    ! set PIO rearranger options
    if (localPet == 0) write(logunit,*) trim(subname),' calling pio_set_rearr_opts'
    ret = pio_set_rearr_opts(io_subsystem, pio_rearr_comm_type, pio_rearr_comm_fcd, &
                             pio_rearr_comm_enable_hs_comp2io, pio_rearr_comm_enable_isend_comp2io, &
                             pio_rearr_comm_max_pend_req_comp2io, &
                             pio_rearr_comm_enable_hs_io2comp, pio_rearr_comm_enable_isend_io2comp, &
                             pio_rearr_comm_max_pend_req_io2comp)
#endif

  end subroutine med_io_init

  !===============================================================================
  subroutine med_io_wopen(filename, io_file, vm, rc, clobber, file_ind, model_doi_url)

    !---------------
    ! open netcdf file
    !---------------

    use pio , only : PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR
    use pio , only : pio_openfile, pio_createfile, PIO_GLOBAL, pio_enddef
    use pio , only : pio_put_att, pio_redef, pio_get_att
    use pio , only : pio_seterrorhandling, pio_file_is_open, pio_clobber, pio_write, pio_noclobber

    ! input/output arguments
    character(*),            intent(in) :: filename
    type(file_desc_t),       intent(inout) :: io_file
    type(ESMF_VM)                       :: vm
    integer,                 intent(out) :: rc
    logical,       optional, intent(in) :: clobber
    integer,       optional, intent(in) :: file_ind
    character(CL), optional, intent(in) :: model_doi_url
    ! local variables
    logical       :: lclobber
    integer       :: rcode
    integer       :: nmode
    integer       :: lfile_ind
    integer       :: iam
    character(CL) :: lversion
    character(CL) :: lmodel_doi_url
    character(*),parameter :: subName = '(med_io_wopen) '
    !-------------------------------------------------------------------------------

    lversion=trim(version)

    lclobber = .false.
    if (present(clobber)) lclobber=clobber

    lmodel_doi_url = 'unset'
    if (present(model_doi_url)) lmodel_doi_url = model_doi_url

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    if (.not. pio_file_is_open(io_file)) then

       if (med_io_file_exists(vm, filename)) then
          if (lclobber) then
             nmode = pio_clobber
             ! only applies to classic NETCDF files.
             if(pio_iotype == PIO_IOTYPE_NETCDF .or. pio_iotype == PIO_IOTYPE_PNETCDF) then
                nmode = ior(nmode,pio_ioformat)
             endif
             rcode = pio_createfile(io_subsystem, io_file, pio_iotype, trim(filename), nmode)
             if(iam==0) write(logunit,'(a)') trim(subname)//' creating file '//trim(filename)
             rcode = pio_put_att(io_file,pio_global,"file_version",version)
             rcode = pio_put_att(io_file,pio_global,"model_doi_url",lmodel_doi_url)
          else
             rcode = pio_openfile(io_subsystem, io_file, pio_iotype, trim(filename), pio_write)
             if (iam==0) write(logunit,'(a)') trim(subname)//' opening file '//trim(filename)
             call pio_seterrorhandling(io_file,PIO_BCAST_ERROR)
             rcode = pio_get_att(io_file,pio_global,"file_version",lversion)
             call pio_seterrorhandling(io_file,PIO_INTERNAL_ERROR)
             if (trim(lversion) /= trim(version)) then
                rcode = pio_redef(io_file)
                rcode = pio_put_att(io_file,pio_global,"file_version",version)
                rcode = pio_enddef(io_file)
             endif
          endif
       else
          nmode = pio_noclobber
          ! only applies to classic NETCDF files.
          if(pio_iotype == PIO_IOTYPE_NETCDF .or. pio_iotype == PIO_IOTYPE_PNETCDF) then
             nmode = ior(nmode,pio_ioformat)
          endif

          rcode = pio_createfile(io_subsystem, io_file, pio_iotype, trim(filename), nmode)
          if (iam==0) write(logunit,'(a)') trim(subname) //' creating file '// trim(filename)
          rcode = pio_put_att(io_file,pio_global,"file_version",version)
          rcode = pio_put_att(io_file,pio_global,"model_doi_url",lmodel_doi_url)
       endif

    else
       ! filename is already open, just return
    endif

  end subroutine med_io_wopen

  !===============================================================================
  subroutine med_io_close(io_file, rc)

    !---------------
    ! close netcdf file
    !---------------

    use pio, only: pio_file_is_open, pio_closefile

    ! input/output variables
    type(file_desc_t) :: io_file
    integer          , intent(out) :: rc

    ! local variables

    character(*),parameter :: subName = '(med_io_close) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (pio_file_is_open(io_file)) then
       call pio_closefile(io_file)
    endif

  end subroutine med_io_close

  !===============================================================================
  subroutine med_io_redef(io_file)

    use pio, only : pio_redef

    ! input/output variables
    type(file_desc_t) :: io_file
    ! local variables
    integer :: rcode
    !-------------------------------------------------------------------------------

    rcode = pio_redef(io_file)

  end subroutine med_io_redef

  !===============================================================================
  subroutine med_io_enddef(io_file)

    use pio, only : pio_enddef

    ! input/output variables
    type(file_desc_t) :: io_file

    ! local variables

    integer :: rcode
    !-------------------------------------------------------------------------------

    rcode = pio_enddef(io_file)

  end subroutine med_io_enddef

  !===============================================================================
  character(len=24) function med_io_date2yyyymmdd (date)

    integer, intent(in) :: date  ! date expressed as an integer: yyyymmdd

    call med_io_datetod2string(date_str = med_io_date2yyyymmdd, ymd = date)

  end function med_io_date2yyyymmdd

  !===============================================================================
  character(len=8) function med_io_sec2hms (seconds, rc)

    ! input arguments
    integer, intent(in)  :: seconds
    integer, intent(out) :: rc

    ! local variables
    integer :: hours     ! hours of hh:mm:ss
    integer :: minutes   ! minutes of hh:mm:ss
    integer :: secs      ! seconds of hh:mm:ss
    !----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (seconds < 0 .or. seconds > 86400) then
       write(logunit,*)'med_io_sec2hms: bad input seconds:', seconds
       call ESMF_LogWrite('med_io_sec2hms: bad input seconds', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    hours   = seconds / 3600
    minutes = (seconds - hours*3600) / 60
    secs    = (seconds - hours*3600 - minutes*60)

    if (minutes < 0 .or. minutes > 60) then
       write(logunit,*)'med_io_sec2hms: bad minutes = ',minutes
       call ESMF_LogWrite('med_io_sec2hms: bad minutes', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    if (secs < 0 .or. secs > 60) then
       write(logunit,*)'med_io_sec2hms: bad secs = ',secs
       call ESMF_LogWrite('med_io_sec2hms: bad secs', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    write(med_io_sec2hms,80) hours, minutes, secs
80  format(i2.2,':',i2.2,':',i2.2)

  end function med_io_sec2hms

  !===============================================================================
  subroutine med_io_write_FB(io_file, FB, whead, wdata, nx, ny, nt, &
       fillval, pre, flds, tavg, use_float, ntile, rc)

    !---------------
    ! Write FB to netcdf file
    !---------------

    use ESMF,  only : operator(==)
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE, ESMF_END_ABORT
    use ESMF , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundle, ESMF_Mesh, ESMF_DistGrid
    use ESMF , only : ESMF_FieldBundleGet, ESMF_FieldGet, ESMF_MeshGet, ESMF_DistGridGet
    use ESMF , only : ESMF_Field, ESMF_FieldGet, ESMF_AttributeGet
    use ESMF , only : ESMF_CoordSys_Flag, ESMF_COORDSYS_SPH_DEG, ESMF_COORDSYS_SPH_RAD, ESMF_COORDSYS_CART
    use pio  , only : var_desc_t, io_desc_t, pio_offset_kind
    use pio  , only : pio_def_dim, pio_inq_dimid, pio_real, pio_def_var, pio_put_att, pio_double
    use pio  , only : pio_inq_varid, pio_setframe, pio_write_darray, pio_initdecomp, pio_freedecomp
    use pio  , only : pio_syncfile

    ! input/output variables
    type(file_desc_t)                       :: io_file
    type(ESMF_FieldBundle)     , intent(in) :: FB        ! data to be written
    logical                    , intent(in) :: whead     ! write header
    logical                    , intent(in) :: wdata     ! write data
    integer                    , intent(in) :: nx        ! 2d grid size if available
    integer                    , intent(in) :: ny        ! 2d grid size if available
    integer ,         optional , intent(in) :: nt        ! time sample
    real(r8),         optional , intent(in) :: fillval   ! fill value
    character(len=*), optional , intent(in) :: pre       ! prefix to variable name
    character(len=*), optional , intent(in) :: flds(:)   ! specific fields to write out
    logical,          optional , intent(in) :: tavg      ! is this a tavg
    logical,          optional , intent(in) :: use_float ! write output as float rather than double
    integer,          optional , intent(in) :: ntile     ! number of nx * ny tiles
    integer                    , intent(out):: rc

    ! local variables
    type(ESMF_Field)              :: field
    type(ESMF_Mesh)               :: mesh
    type(ESMF_Distgrid)           :: distgrid
    type(ESMF_CoordSys_Flag)      :: coordsys
    integer                       :: rcode
    integer                       :: nf,ns,ng
    integer                       :: k,n
    integer                       :: ndims, nelements
    integer    ,target            :: dimid2(2)
    integer    ,target            :: dimid3(3)
    integer    ,target            :: dimid4(4)
    integer    ,pointer           :: dimid(:)
    type(var_desc_t)              :: varid
    type(io_desc_t)               :: iodesc
    integer(kind=Pio_Offset_Kind) :: frame
    character(CL)                 :: itemc       ! string converted to char
    character(CL)                 :: name1       ! var name
    character(CL)                 :: cunit       ! var units
    character(CL)                 :: lpre        ! local prefix
    character(CS)                 :: coordvarnames(2)   ! coordinate variable names
    character(CS)                 :: coordnames(2)      ! coordinate long names
    character(CS)                 :: coordunits(2)      ! coordinate units
    integer                       :: lnx,lny,lntile
    logical                       :: luse_float
    real(r8)                      :: lfillvalue
    integer, pointer              :: minIndexPTile(:,:)
    integer, pointer              :: maxIndexPTile(:,:)
    integer                       :: dimCount, tileCount
    integer, pointer              :: Dof(:)
    real(r8), pointer             :: fldptr1(:)
    real(r8), pointer             :: fldptr2(:,:)
    real(r8), allocatable         :: ownedElemCoords(:), ownedElemCoords_x(:), ownedElemCoords_y(:)
    character(CS)                 :: cnumber
    character(CL)                 :: tmpstr
    type(ESMF_Field)              :: lfield
    integer                       :: rank
    integer                       :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fields
    integer                       :: gridToFieldMap(1)  ! currently the size must equal 1 for rank 2 fields
    logical                       :: tiles
    character(CL), allocatable    :: fieldNameList(:)
    character(*),parameter :: subName = '(med_io_write_FB) '
    !-------------------------------------------------------------------------------

    rc = ESMF_Success

    lfillvalue = fillvalue
    if (present(fillval)) lfillvalue = fillval
    lpre = ' '
    if (present(pre)) lpre = trim(pre)
    luse_float = .false.
    if (present(use_float)) luse_float = use_float

    tiles = .false.
    if (present(ntile)) then
      if (ntile > 0) tiles = .true.
    end if

    ! Error check
    if (.not. ESMF_FieldBundleIsCreated(FB, rc=rc)) then
       call ESMF_LogWrite(trim(subname)//" FB "//trim(lpre)//" not created", ESMF_LOGMSG_INFO)
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
       endif
       return
    endif

    ! Get number of fields
    if (present(flds)) then
       nf = size(flds)
    else
       call ESMF_FieldBundleGet(FB, fieldCount=nf, rc=rc)
       write(tmpstr,*) subname//' field count = '//trim(lpre), nf
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       if (nf < 1) then
          call ESMF_LogWrite(trim(subname)//" FB "//trim(lpre)//" empty", ESMF_LOGMSG_INFO)
          if (dbug_flag > 5) then
             call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
          endif
          rc = ESMF_Success
          return
       endif
       allocate(fieldNameList(nf))
       call ESMF_FieldBundleGet(FB, fieldNameList=fieldNameList, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Get field bundle mesh from first field
    call FB_getFieldN(FB, 1, field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, mesh=mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get mesh distgrid and number of elements
    call ESMF_MeshGet(mesh, elementDistgrid=distgrid, coordSys=coordsys, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh, spatialDim=ndims, numOwnedElements=nelements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(tmpstr,*) subname, 'ndims, nelements = ', ndims, nelements
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    ! Define coordinate attributes according to CoordSys
    if (coordsys == ESMF_COORDSYS_CART) then
       coordvarnames(1) = trim(lpre)//'_x'
       coordvarnames(2) = trim(lpre)//'_y'
       coordnames = (/'x-coordinate', 'y-coordinate'/)
       coordunits = (/'unitless','unitless'/)
    else
       coordvarnames(1) = trim(lpre)//'_lon'
       coordvarnames(2) = trim(lpre)//'_lat'
       coordnames = (/'longitude', 'latitude '/)
       if (coordsys == ESMF_COORDSYS_SPH_DEG) coordunits = (/'degrees_E', 'degrees_N'/)
       if (coordsys == ESMF_COORDSYS_SPH_RAD) coordunits = (/'radians  ', 'radians  '/)
    end if

    ! Set element coordinates
    if (.not. allocated(ownedElemCoords) .and. ndims > 0 .and. nelements > 0) then
       allocate(ownedElemCoords(ndims*nelements))
       allocate(ownedElemCoords_x(ndims*nelements/2))
       allocate(ownedElemCoords_y(ndims*nelements/2))
       call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ownedElemCoords_x = ownedElemCoords(1::2)
       ownedElemCoords_y = ownedElemCoords(2::2)
    end if

    ! Get tile info
    call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount))
    call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! write(tmpstr,*) subname,' counts = ',dimcount,tilecount,minindexptile,maxindexptile
    ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    ! TODO: this is not getting the global size correct for a FB coming in that does not have
    ! all the global grid values in the distgrid - e.g. CTSM

    ng = maxval(maxIndexPTile)
    if (tiles) then
      lnx = nx
      lny = ny
      lntile = ng/(lnx*lny)
      write(tmpstr,*) subname, 'ng,lnx,lny,lntile = ',ng,lnx,lny,lntile
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
      if (lntile /= ntile) then
         call ESMF_LogWrite(trim(subname)//' ERROR: grid2d size and ntile are not consistent ', ESMF_LOGMSG_INFO)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif
    else
      lnx = ng
      lny = 1
      if (nx > 0) lnx = nx
      if (ny > 0) lny = ny
      if (lnx*lny /= ng) then
         write(tmpstr,*) subname,' WARNING: grid2d size not consistent ',ng,lnx,lny
         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
      endif
    end if
    deallocate(minIndexPTile, maxIndexPTile)

    if (present(nt)) then
       frame = nt
    else
       frame = -1
    end if

    ! Write header
    if (whead) then
      if (tiles) then
       rcode = pio_def_dim(io_file, trim(lpre)//'_nx', lnx, dimid3(1))
       rcode = pio_def_dim(io_file, trim(lpre)//'_ny', lny, dimid3(2))
       rcode = pio_def_dim(io_file, trim(lpre)//'_ntile', ntile, dimid3(3))
       if (present(nt)) then
          dimid4(1:3) = dimid3
          rcode = pio_inq_dimid(io_file, 'time', dimid4(4))
          dimid => dimid4
       else
          dimid => dimid3
       endif
      else
       rcode = pio_def_dim(io_file, trim(lpre)//'_nx', lnx, dimid2(1))
       rcode = pio_def_dim(io_file, trim(lpre)//'_ny', lny, dimid2(2))
       if (present(nt)) then
          dimid3(1:2) = dimid2
          rcode = pio_inq_dimid(io_file, 'time', dimid3(3))
          dimid => dimid3
       else
          dimid => dimid2
       endif
      endif
      write(tmpstr,*) subname,' dimid = ',dimid
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       do k = 1,nf
          ! Determine field name
          if (present(flds)) then
             itemc = trim(flds(k))
          else
             itemc = trim(fieldNameList(k))
          end if

          ! Determine rank of field with name itemc
          call ESMF_FieldBundleGet(FB, itemc,  field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, rank=rank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! TODO (mvertens, 2019-03-13): this is a temporary mod to NOT write hgt
          if (trim(itemc) /= "hgt") then
             if (rank == 2) then
                call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                write(cnumber,'(i0)') ungriddedUbound(1)
                call ESMF_LogWrite(trim(subname)//':'//'field '//trim(itemc)// &
                     ' has an griddedUBound of  '//trim(cnumber), ESMF_LOGMSG_INFO)

                ! Create a new output variable for each element of the undistributed dimension
                do n = 1,ungriddedUBound(1)
                   if (trim(itemc) /= "hgt") then
                      write(cnumber,'(i0)') n
                      name1 = trim(lpre)//'_'//trim(itemc)//trim(cnumber)
                      call ESMF_LogWrite(trim(subname)//': defining '//trim(name1), ESMF_LOGMSG_INFO)
                      if (luse_float) then
                         rcode = pio_def_var(io_file, trim(name1), PIO_REAL, dimid, varid)
                         rcode = pio_put_att(io_file, varid,"_FillValue",real(lfillvalue,r4))
                      else
                         rcode = pio_def_var(io_file, trim(name1), PIO_DOUBLE, dimid, varid)
                         rcode = pio_put_att(io_file,varid,"_FillValue",lfillvalue)
                      end if
                      if (NUOPC_FieldDictionaryHasEntry(trim(itemc))) then
                         call NUOPC_FieldDictionaryGetEntry(itemc, canonicalUnits=cunit, rc=rc)
                         if (chkerr(rc,__LINE__,u_FILE_u)) return
                         rcode = pio_put_att(io_file, varid, "units"        , trim(cunit))
                      end if
                      rcode = pio_put_att(io_file, varid, "standard_name", trim(name1))
                      if (present(tavg)) then
                         if (tavg) then
                            rcode = pio_put_att(io_file, varid, "cell_methods", "time: mean")
                         endif
                      endif
                   end if
                end do
             else
                name1 = trim(lpre)//'_'//trim(itemc)
                call ESMF_LogWrite(trim(subname)//':'//trim(itemc)//':'//trim(name1),ESMF_LOGMSG_INFO)
                if (luse_float) then
                   rcode = pio_def_var(io_file, trim(name1), PIO_REAL, dimid, varid)
                   rcode = pio_put_att(io_file, varid, "_FillValue", real(lfillvalue, r4))
                else
                   rcode = pio_def_var(io_file, trim(name1), PIO_DOUBLE, dimid, varid)
                   rcode = pio_put_att(io_file, varid, "_FillValue", lfillvalue)
                end if
                if (NUOPC_FieldDictionaryHasEntry(trim(itemc))) then
                   call NUOPC_FieldDictionaryGetEntry(itemc, canonicalUnits=cunit, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return
                   rcode = pio_put_att(io_file, varid, "units", trim(cunit))
                end if
                rcode = pio_put_att(io_file, varid, "standard_name", trim(name1))
                if (present(tavg)) then
                   if (tavg) then
                      rcode = pio_put_att(io_file, varid, "cell_methods", "time: mean")
                   endif
                end if
             end if
          end if
       end do

       ! Add coordinate information to file
       do n = 1,ndims
          if (luse_float) then
             rcode = pio_def_var(io_file, trim(coordvarnames(n)), PIO_REAL, dimid, varid)
          else
             rcode = pio_def_var(io_file, trim(coordvarnames(n)), PIO_DOUBLE, dimid, varid)
          end if
          rcode = pio_put_att(io_file, varid, "long_name", trim(coordnames(n)))
          rcode = pio_put_att(io_file, varid, "units", trim(coordunits(n)))
          rcode = pio_put_att(io_file, varid, "standard_name", trim(coordnames(n)))
       end do
    end if

    if (wdata) then
       ! use distgrid extracted from field 1 above
       call ESMF_DistGridGet(distgrid, localDE=0, elementCount=ns, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(dof(ns))
       call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
       write(tmpstr,*) subname,' dof = ',ns,size(dof),dof(1),dof(ns)  !,minval(dof),maxval(dof)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       if (tiles) then
          call pio_initdecomp(io_subsystem, pio_double, (/lnx,lny,ntile/), dof, iodesc)
       else
          call pio_initdecomp(io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
         !call pio_writedof(lpre, (/lnx,lny/), int(dof,kind=PIO_OFFSET_KIND), mpicom)
       end if
       deallocate(dof)

       do k = 1,nf
          ! Determine field name
          if (present(flds)) then
             itemc = trim(flds(k))
          else
             itemc = trim(fieldNameList(k))
          end if

          call FB_getFldPtr(FB, itemc, &
               fldptr1=fldptr1, fldptr2=fldptr2, rank=rank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! TODO (mvertens, 2019-03-13): this is a temporary mod to NOT write hgt
          if (trim(itemc) /= "hgt") then
             if (rank == 2) then

                ! Determine the size of the ungridded dimension and the index where the undistributed dimension is located
                call ESMF_FieldBundleGet(FB, itemc,  field=lfield, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, gridToFieldMap=gridToFieldMap, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! Output for each ungriddedUbound index
                do n = 1,ungriddedUBound(1)
                   write(cnumber,'(i0)') n
                   name1 = trim(lpre)//'_'//trim(itemc)//trim(cnumber)
                   rcode = pio_inq_varid(io_file, trim(name1), varid)
                   call pio_setframe(io_file,varid,frame)

                   if (gridToFieldMap(1) == 1) then
                      call pio_write_darray(io_file, varid, iodesc, fldptr2(:,n), rcode, fillval=lfillvalue)
                   else if (gridToFieldMap(1) == 2) then
                      call pio_write_darray(io_file, varid, iodesc, fldptr2(n,:), rcode, fillval=lfillvalue)
                   end if
                end do
             else if (rank == 1 .or. rank == 0) then
                name1 = trim(lpre)//'_'//trim(itemc)
                rcode = pio_inq_varid(io_file, trim(name1), varid)
                call pio_setframe(io_file,varid,frame)
                ! fix for writing data on exchange grid, which has no data in some PETs
                if (rank == 0) nullify(fldptr1)
                call pio_write_darray(io_file, varid, iodesc, fldptr1, rcode, fillval=lfillvalue)
             end if  ! end if rank is 2 or 1 or 0

          end if ! end if not "hgt"
       end do  ! end loop over fields in FB

       ! Fill coordinate variables - why is this being done each time?
       rcode = pio_inq_varid(io_file, trim(coordvarnames(1)), varid)
       call pio_setframe(io_file,varid,frame)
       call pio_write_darray(io_file, varid, iodesc, ownedElemCoords_x, rcode, fillval=lfillvalue)

       rcode = pio_inq_varid(io_file, trim(coordvarnames(2)), varid)
       call pio_setframe(io_file,varid,frame)
       call pio_write_darray(io_file, varid, iodesc, ownedElemCoords_y, rcode, fillval=lfillvalue)

       call pio_syncfile(io_file)
       call pio_freedecomp(io_file, iodesc)
    endif
    if(allocated(ownedElemCoords)) then
       deallocate(ownedElemCoords)
    endif
    if(allocated(ownedElemCoords_x)) then
       deallocate(ownedElemCoords_x)
    endif
    if(allocated(ownedElemCoords_y)) then
       deallocate(ownedElemCoords_y)
    endif

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_io_write_FB

  !===============================================================================
  subroutine med_io_write_int(io_file, idata, dname, whead, wdata, rc)

    use pio, only : var_desc_t, pio_def_var, pio_put_att, pio_int, pio_inq_varid, pio_put_var

    !---------------
    ! Write scalar integer to netcdf file
    !---------------

    ! intput/output variables
    type(file_desc_t)            :: io_file
    integer          ,intent(in) :: idata    ! data to be written
    character(len=*) ,intent(in) :: dname    ! name of data
    logical          ,intent(in) :: whead    ! write header
    logical          ,intent(in) :: wdata    ! write data
    integer          ,intent(out):: rc

    ! local variables
    integer          :: rcode
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(*),parameter :: subName = '(med_io_write_int) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (whead) then
       if (NUOPC_FieldDictionaryHasEntry(trim(dname))) then
          call NUOPC_FieldDictionaryGetEntry(dname, canonicalUnits=cunit, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          rcode = pio_put_att(io_file,varid,"units",trim(cunit))
       end if
       rcode = pio_def_var(io_file,trim(dname),PIO_INT,varid)
       rcode = pio_put_att(io_file,varid,"standard_name",trim(dname))
    endif
    if (wdata) then
       rcode = pio_inq_varid(io_file,trim(dname),varid)
       rcode = pio_put_var(io_file,varid,idata)
    endif

  end subroutine med_io_write_int

  !===============================================================================
  subroutine med_io_write_int1d(io_file, idata, dname, whead, wdata, file_ind, rc)

    !---------------
    ! Write 1d integer array to netcdf file
    !---------------

    use pio     , only : var_desc_t, pio_def_dim, pio_def_var
    use pio     , only : pio_put_att, pio_inq_varid, pio_put_var
    use pio     , only : pio_int, pio_def_var

    ! input/output arguments
    type(file_desc_t)            :: io_file
    integer          ,intent(in) :: idata(:) ! data to be written
    character(len=*) ,intent(in) :: dname    ! name of data
    logical          ,intent(in) :: whead    ! write header
    logical          ,intent(in) :: wdata    ! write data
    integer,optional ,intent(in) :: file_ind
    integer         , intent(out):: rc

    ! local variables
    integer          :: rcode
    integer          :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    integer          :: lnx
    integer          :: lfile_ind
    character(*),parameter :: subName = '(med_io_write_int1d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    if (whead) then
       if (NUOPC_FieldDictionaryHasEntry(trim(dname))) then
          call NUOPC_FieldDictionaryGetEntry(dname, canonicalUnits=cunit, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          rcode = pio_put_att(io_file,varid,"units",trim(cunit))
       end if
       lnx = size(idata)
       rcode = pio_def_dim(io_file,trim(dname),lnx,dimid(1))
       rcode = pio_def_var(io_file,trim(dname),PIO_INT,dimid,varid)
       rcode = pio_put_att(io_file,varid,"standard_name",trim(dname))
    else if (wdata) then
       rcode = pio_inq_varid(io_file,trim(dname),varid)
       rcode = pio_put_var(io_file,varid,idata)
    endif

  end subroutine med_io_write_int1d

  !===============================================================================
  subroutine med_io_write_r8(io_file, rdata, dname, whead, wdata, rc)

    !---------------
    ! Write scalar double to netcdf file
    !---------------

    use pio , only : var_desc_t, pio_def_var, pio_put_att
    use pio , only : pio_double, pio_noerr, pio_inq_varid, pio_put_var

    ! input/output arguments
    type(file_desc_T)            :: io_file
    real(r8)         ,intent(in) :: rdata    ! data to be written
    character(len=*) ,intent(in) :: dname    ! name of data
    logical          ,intent(in) :: whead    ! write header
    logical          ,intent(in) :: wdata    ! write data
    integer          ,intent(out):: rc

    ! local variables
    integer          :: rcode
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(*),parameter :: subName = '(med_io_write_r8) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (whead) then
       rcode = pio_def_var(io_file,trim(dname),PIO_DOUBLE,varid)
       if (rcode==PIO_NOERR) then
          if (NUOPC_FieldDictionaryHasEntry(trim(dname))) then
             call NUOPC_FieldDictionaryGetEntry(dname, canonicalUnits=cunit, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             rcode = pio_put_att(io_file,varid,"units",trim(cunit))
          end if
          rcode = pio_put_att(io_file,varid,"standard_name",trim(dname))
       end if
    else if (wdata) then
       rcode = pio_inq_varid(io_file,trim(dname),varid)
       rcode = pio_put_var(io_file,varid,rdata)
    endif

  end subroutine med_io_write_r8

  !===============================================================================
  subroutine med_io_write_r81d(io_file, rdata, dname, whead, wdata, rc)

    !---------------
    ! Write 1d double array to netcdf file
    !---------------

    use pio , only : var_desc_t, pio_def_dim, pio_def_var
    use pio , only : pio_inq_varid, pio_put_var, pio_double, pio_put_att

    ! !INPUT/OUTPUT PARAMETERS:
    type(file_desc_t)            :: io_file
    real(r8)         ,intent(in) :: rdata(:) ! data to be written
    character(len=*) ,intent(in) :: dname    ! name of data
    logical          ,intent(in) :: whead    ! write header
    logical          ,intent(in) :: wdata    ! write data
    integer          ,intent(out):: rc

    ! local variables
    integer          :: rcode
    integer          :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    integer          :: lnx
    character(*),parameter :: subName = '(med_io_write_r81d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (whead) then
       lnx = size(rdata)
       rcode = pio_def_dim(io_file,trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(io_file,trim(dname),PIO_DOUBLE,dimid,varid)
       if (NUOPC_FieldDictionaryHasEntry(trim(dname))) then
          call NUOPC_FieldDictionaryGetEntry(dname, canonicalUnits=cunit, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          rcode = pio_put_att(io_file,varid,"units",trim(cunit))
       end if
       rcode = pio_put_att(io_file,varid,"standard_name",trim(dname))
    endif

    if (wdata) then
       rcode = pio_inq_varid(io_file,trim(dname),varid)
       rcode = pio_put_var(io_file,varid,rdata)
    endif

  end subroutine med_io_write_r81d

  !===============================================================================
  subroutine med_io_write_char(io_file, rdata, dname, whead, wdata, rc)

    !---------------
    ! Write char string to netcdf file
    !---------------

    use pio , only : var_desc_t, pio_def_dim, pio_put_att, pio_def_var, pio_inq_varid
    use pio , only : pio_char, pio_put_var

    ! input/output arguments
    type(file_desc_t)            :: io_file
    character(len=*) ,intent(in) :: rdata    ! data to be written
    character(len=*) ,intent(in) :: dname    ! name of data
    logical          ,intent(in) :: whead    ! write header
    logical          ,intent(in) :: wdata    ! write data
    integer          ,intent(out):: rc

    ! local variables
    integer          :: rcode
    integer          :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    integer          :: lnx
    character(CL)    :: charvar   ! buffer for string read/write
    character(*),parameter :: subName = '(med_io_write_char) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (whead) then
       lnx = len(charvar)
       rcode = pio_def_dim(io_file,trim(dname)//'_len',lnx,dimid(1))
       rcode = pio_def_var(io_file,trim(dname),PIO_CHAR,dimid,varid)
       if (NUOPC_FieldDictionaryHasEntry(trim(dname))) then
          call NUOPC_FieldDictionaryGetEntry(dname, canonicalUnits=cunit, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       rcode = pio_put_att(io_file,varid,"standard_name",trim(dname))
    else if (wdata) then
       charvar = ''
       charvar = trim(rdata)
       rcode = pio_inq_varid(io_file,trim(dname),varid)
       rcode = pio_put_var(io_file,varid,charvar)
    endif

  end subroutine med_io_write_char

  !===============================================================================
  subroutine med_io_define_time(io_file, time_units, calendar, rc)

    use ESMF, only : operator(==), operator(/=)
    use ESMF, only : ESMF_Calendar, ESMF_CalendarIsCreated
    use ESMF, only : ESMF_CALKIND_360DAY, ESMF_CALKIND_GREGORIAN
    use ESMF, only : ESMF_CALKIND_JULIAN, ESMF_CALKIND_JULIANDAY, ESMF_CALKIND_MODJULIANDAY
    use ESMF, only : ESMF_CALKIND_NOCALENDAR, ESMF_CALKIND_NOLEAP
    use ESMF, only : ESMF_LOGMSG_ERROR, ESMF_FAILURE
    use pio , only : var_desc_t, PIO_UNLIMITED
    use pio , only : pio_double, pio_def_dim, pio_def_var, pio_put_att
    use pio , only : pio_inq_varid, pio_put_var

    ! input/output variables
    type(file_desc_t)                :: io_file
    character(len=*)    , intent(in) :: time_units ! units of time
    type(ESMF_Calendar) , intent(in) :: calendar   ! calendar
    integer             , intent(out):: rc

    ! local variables
    integer          :: rcode
    integer          :: dimid(1)
    integer          :: dimid2(2)
    type(var_desc_t) :: varid
    character(CL)    :: calname        ! calendar name
    character(*),parameter :: subName = '(med_io_define_time) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. ESMF_CalendarIsCreated(calendar)) then
       call ESMF_LogWrite(trim(subname)//' ERROR: calendar is not created ', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    end if

    ! define time and add calendar attribute
    rcode = pio_def_dim(io_file, 'time', PIO_UNLIMITED, dimid(1))
    rcode = pio_def_var(io_file, 'time', PIO_DOUBLE, dimid, varid)
    rcode = pio_put_att(io_file, varid, 'units', trim(time_units))
    if (calendar == ESMF_CALKIND_360DAY) then
       calname = '360_day'
    else if (calendar == ESMF_CALKIND_GREGORIAN) then
       calname = 'gregorian'
    else if (calendar == ESMF_CALKIND_JULIAN) then
       calname = 'julian'
    else if (calendar == ESMF_CALKIND_JULIANDAY) then
       calname = 'ESMF_CALKIND_JULIANDAY'
    else if (calendar == ESMF_CALKIND_MODJULIANDAY) then
       calname = 'ESMF_CALKIND_MODJULIANDAY'
    else if (calendar == ESMF_CALKIND_NOCALENDAR) then
       calname = 'none'
    else if (calendar == ESMF_CALKIND_NOLEAP) then
       calname = 'noleap'
    end if
    rcode = pio_put_att(io_file, varid, 'calendar', trim(calname))

    ! define time bounds
    dimid2(2) = dimid(1)
    rcode = pio_def_dim(io_file, 'ntb', 2, dimid2(1))
    rcode = pio_def_var(io_file, 'time_bnds', PIO_DOUBLE, dimid2, varid)
    rcode = pio_put_att(io_file, varid, 'bounds', 'time_bnds')

  end subroutine med_io_define_time

  !===============================================================================
  subroutine med_io_write_time(io_file, time_val, tbnds, nt, rc)

    !---------------
    ! Write time variable to netcdf file
    !---------------

    use pio, only : pio_put_att, pio_inq_varid, pio_put_var

    ! input/output variables
    type(file_desc_t)               :: io_file
    real(r8) ,           intent(in) :: time_val   ! data to be written
    real(r8) ,           intent(in) :: tbnds(2)   ! time bounds
    integer  ,           intent(in) :: nt
    integer  ,           intent(out):: rc

    ! local variables
    integer :: rcode
    integer :: varid
    integer :: start(2),count(2)
    character(*),parameter :: subName = '(med_io_write_time) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! write time
    count = 1; start = nt
    rcode = pio_inq_varid(io_file, 'time', varid)
    rcode = pio_put_var(io_file, varid, start(1:1), count(1:1), (/time_val/))

    ! write time bounds
    rcode = pio_inq_varid(io_file, 'time_bnds', varid)
    start(1) = 1; start(2) = nt
    count(1) = 2; count(2) = 1
    rcode = pio_put_var(io_file, varid, start(1:2), count(1:2), tbnds)

  end subroutine med_io_write_time

  !===============================================================================
  subroutine med_io_read_FB(filename, vm, FB, pre, frame, rc)

    !---------------
    ! Read FB from netcdf file
    !---------------

    use ESMF , only : ESMF_FieldBundle, ESMF_Field, ESMF_Mesh, ESMF_DistGrid
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE
    use ESMF , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF , only : ESMF_FieldGet, ESMF_MeshGet, ESMF_DistGridGet
    use pio  , only : file_desc_T, var_desc_t, io_desc_t, pio_nowrite, pio_openfile
    use pio  , only : pio_noerr, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR
    use pio  , only : pio_inq_varid
    use pio  , only : pio_double, pio_get_att, pio_seterrorhandling, pio_freedecomp, pio_closefile
    use pio  , only : pio_read_darray, pio_offset_kind, pio_setframe

    ! input/output arguments
    character(len=*)                        ,intent(in)  :: filename
    type(ESMF_VM)                           ,intent(in)  :: vm
    type(ESMF_FieldBundle)                  ,intent(in)  :: FB       ! data to be read
    character(len=*)              ,optional ,intent(in)  :: pre      ! prefix to variable name
    integer(kind=PIO_OFFSET_KIND) ,optional ,intent(in)  :: frame
    integer                                 ,intent(out) :: rc

    ! local variables
    type(ESMF_Field)              :: lfield
    integer                       :: rcode
    integer                       :: nf
    integer                       :: k,n,l
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    type(io_desc_t)               :: iodesc
    character(CL)                 :: itemc       ! string converted to char
    character(CL)                 :: name1       ! var name
    character(CL)                 :: lpre        ! local prefix
    real(r8)                      :: lfillvalue
    integer                       :: rank, lsize
    real(r8), pointer             :: fldptr1(:), fldptr1_tmp(:)
    real(r8), pointer             :: fldptr2(:,:)
    character(CL)                 :: tmpstr
    character(len=16)             :: cnumber
    integer(kind=Pio_Offset_Kind) :: lframe
    integer                       :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fieldds
    integer                       :: gridToFieldMap(1)  ! currently the size must equal 1 for rank 2 fieldds
    character(*),parameter :: subName = '(med_io_read_FB) '
    !-------------------------------------------------------------------------------
    rc = ESMF_Success
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lpre = ' '
    if (present(pre)) then
       lpre = trim(pre)
    endif
    if (present(frame)) then
       lframe = frame
    else
       lframe = 1
    endif
    if (.not. ESMF_FieldBundleIsCreated(FB,rc=rc)) then
       call ESMF_LogWrite(trim(subname)//" FB "//trim(lpre)//" not created", ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       endif
       return
    endif

    call ESMF_FieldBundleGet(FB, fieldCount=nf, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(tmpstr,*) subname//' field count = '//trim(lpre),nf
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (nf < 1) then
       call ESMF_LogWrite(trim(subname)//" FB "//trim(lpre)//" empty", ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       endif
       return
    endif

    if (med_io_file_exists(vm, trim(filename))) then
       rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename),pio_nowrite)
       call ESMF_LogWrite(trim(subname)//' open file '//trim(filename), ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//' ERROR: file invalid '//trim(filename), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    endif

    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)

    do k = 1,nf
       ! Get name of field
       call FB_getNameN(FB, k, itemc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Get iodesc for all fields based on iodesc of first field (assumes that all fields have
       ! the same iodesc)
       if (k == 1) then
          call ESMF_FieldBundleGet(FB, itemc,  field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, rank=rank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (rank == 2) then
             name1 = trim(lpre)//'_'//trim(itemc)//'1'
          else if (rank == 1) then
             name1 = trim(lpre)//'_'//trim(itemc)
          end if
          call med_io_read_init_iodesc(FB, name1, pioid, iodesc, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       call ESMF_LogWrite(trim(subname)//' reading field '//trim(itemc), ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Get pointer to field bundle field
       ! Field bundle might be 2d or 1d - but field on mediator history or restart file will always be 1d
       call FB_getFldPtr(FB, itemc, &
            fldptr1=fldptr1, fldptr2=fldptr2, rank=rank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (rank == 2) then

          ! Determine the size of the ungridded dimension and the
          ! index where the undistributed dimension is located
          call ESMF_FieldBundleGet(FB, itemc,  field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, gridToFieldMap=gridToFieldMap, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (gridToFieldMap(1) == 1) then
             lsize = size(fldptr2, dim=1)
          else if (gridToFieldMap(1) == 2) then
             lsize = size(fldptr2, dim=2)
          end if
          allocate(fldptr1_tmp(lsize))

          do n = 1,ungriddedUBound(1)
             ! Create a name for the 1d field on the mediator history or restart file based on the
             ! ungridded dimension index of the field bundle 2d field
             write(cnumber,'(i0)') n
             name1 = trim(lpre)//'_'//trim(itemc)//trim(cnumber)

             rcode = pio_inq_varid(pioid, trim(name1), varid)
             if (rcode == pio_noerr) then
                call ESMF_LogWrite(trim(subname)//' read field '//trim(name1), ESMF_LOGMSG_INFO)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call pio_setframe(pioid, varid, lframe)
                call pio_read_darray(pioid, varid, iodesc, fldptr1_tmp, rcode)
                rcode = pio_get_att(pioid, varid, "_FillValue", lfillvalue)
                if (rcode /= pio_noerr) then
                   lfillvalue = fillvalue
                endif
                do l = 1,size(fldptr1_tmp)
                   if (fldptr1_tmp(l) == lfillvalue) fldptr1_tmp(l) = 0.0_r8
                enddo
             else
                fldptr1_tmp = 0.0_r8
             endif
             if (gridToFieldMap(1) == 1) then
                fldptr2(:,n) = fldptr1_tmp(:)
             else if (gridToFieldMap(1) == 2) then
                fldptr2(n,:) = fldptr1_tmp(:)
             end if
          end do

          deallocate(fldptr1_tmp)

       else if (rank == 1) then
          name1 = trim(lpre)//'_'//trim(itemc)

          rcode = pio_inq_varid(pioid, trim(name1), varid)
          if (rcode == pio_noerr) then
             call ESMF_LogWrite(trim(subname)//' read field '//trim(name1), ESMF_LOGMSG_INFO)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call pio_setframe(pioid,varid,lframe)
             call pio_read_darray(pioid, varid, iodesc, fldptr1, rcode)
             rcode = pio_get_att(pioid,varid,"_FillValue",lfillvalue)
             if (rcode /= pio_noerr) then
                lfillvalue = fillvalue
             endif
             do n = 1,size(fldptr1)
                if (fldptr1(n) == lfillvalue) fldptr1(n) = 0.0_r8
             enddo
          else
             fldptr1 = 0.0_r8
          endif
       end if

    enddo ! end of loop over fields
    call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)

    call pio_freedecomp(pioid, iodesc)
    call pio_closefile(pioid)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_io_read_FB

  !===============================================================================
  subroutine med_io_read_init_iodesc(FB, name1, pioid, iodesc, rc)

    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundle, ESMF_Mesh, ESMF_DistGrid
    use ESMF , only : ESMF_FieldBundleGet, ESMF_FieldGet, ESMF_MeshGet, ESMF_DistGridGet
    use ESMF , only : ESMF_Field, ESMF_FieldGet, ESMF_AttributeGet
    use pio  , only : file_desc_T, var_desc_t, io_desc_t, pio_nowrite, pio_openfile
    use pio  , only : pio_noerr, pio_inq_varndims
    use pio  , only : pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, pio_inq_vardimid
    use pio  , only : pio_double, pio_seterrorhandling, pio_initdecomp

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: name1
    type(file_desc_t)      , intent(in)    :: pioid
    type(io_desc_t)        , intent(inout) :: iodesc
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)    :: field
    type(ESMF_Mesh)     :: mesh
    type(ESMF_Distgrid) :: distgrid
    integer             :: rcode
    integer             :: ns,ng
    integer             :: ndims
    integer, pointer    :: dimid(:)
    type(var_desc_t)    :: varid
    integer             :: lnx,lny
    integer, pointer    :: minIndexPTile(:,:)
    integer, pointer    :: maxIndexPTile(:,:)
    integer             :: dimCount, tileCount
    integer, pointer    :: Dof(:)
    character(CL)       :: tmpstr
    character(*),parameter :: subName = '(med_io_read_init_iodesc) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    rcode = pio_inq_varid(pioid, trim(name1), varid)
    if (rcode == pio_noerr) then

       rcode = pio_inq_varndims(pioid, varid, ndims)
       write(tmpstr,*) trim(subname),' ndims = ',ndims
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       allocate(dimid(ndims))
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
       write(tmpstr,*) trim(subname),' lnx = ',lnx
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       if (ndims>=2) then
          rcode = pio_inq_dimlen(pioid, dimid(2), lny)
       else
          lny = 1
       end if
       deallocate(dimid)

       write(tmpstr,*) trim(subname),' lny = ',lny
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       ng = lnx * lny
       call FB_getFieldN(FB, 1, field, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, mesh=mesh, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       allocate(minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount))
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (ng > maxval(maxIndexPTile)) then
          write(tmpstr,*) subname,' WARNING: dimensions do not match', lnx, lny, maxval(maxIndexPTile)
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
          !TODO: this should not be an error for say CTSM which does not send a global grid
          !rc = ESMF_Failure
          !return
       endif

       call ESMF_DistGridGet(distgrid, localDE=0, elementCount=ns, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       allocate(dof(ns))
       call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
       write(tmpstr,*) subname,' dof = ',ns,size(dof),dof(1),dof(ns)  !,minval(dof),maxval(dof)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       call pio_initdecomp(io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
       deallocate(dof)

       deallocate(minIndexPTile, maxIndexPTile)
    else
       if(maintask) write(logunit,'(a)') trim(subname)//' ERROR: '//trim(name1)//' is not present, aborting '
       call ESMF_LogWrite(trim(subname)//' ERROR: '//trim(name1)//' is not present, aborting ', ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
    end if ! end if rcode check

  end subroutine med_io_read_init_iodesc

  !===============================================================================
  subroutine med_io_read_int(filename, vm, idata, dname, rc)

    !---------------
    ! Read scalar integer from netcdf file
    !---------------

    ! input/output arguments
    character(len=*) , intent(in)    :: filename ! file
    type(ESMF_VM)                    :: vm
    integer          , intent(inout) :: idata    ! integer data
    character(len=*) , intent(in)    :: dname    ! name of data
    integer          , intent(out)   :: rc

    ! local variables
    integer :: i1d(1)
    character(*),parameter :: subName = '(med_io_read_int) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call med_io_read_int1d(filename, vm, i1d, dname, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    idata = i1d(1)

  end subroutine med_io_read_int

  !===============================================================================
  subroutine med_io_read_int1d(filename, vm, idata, dname, rc)

    !---------------
    ! Read 1d integer array from netcdf file
    !---------------

    use pio , only : var_desc_t, file_desc_t, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, pio_seterrorhandling
    use pio , only : pio_get_var, pio_inq_varid, pio_get_att, pio_openfile
    use pio , only : pio_nowrite, pio_openfile, pio_global
    use pio , only : pio_closefile

    ! input/output arguments
    character(len=*), intent(in)    :: filename ! file
    type(ESMF_VM)                   :: vm
    integer         , intent(inout) :: idata(:) ! integer data
    character(len=*), intent(in)    :: dname    ! name of data
    integer         , intent(out)   :: rc

    ! local variables
    integer           :: rcode
    type(file_desc_t) :: pioid
    type(var_desc_t)  :: varid
    character(CL)     :: lversion
    character(CL)     :: name1
    integer           :: iam
    character(*),parameter :: subName = '(med_io_read_int1d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lversion=trim(version)

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (med_io_file_exists(vm, filename)) then
       rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename),pio_nowrite)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call ESMF_LogWrite(trim(subname)//'ERROR: file invalid '//trim(filename)//' '//trim(dname), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,idata)

    call pio_closefile(pioid)
  end subroutine med_io_read_int1d

  !===============================================================================
  subroutine med_io_read_r8(filename, vm, rdata, dname, rc)

    !---------------
    ! Read scalar double from netcdf file
    !---------------

    ! input/output arguments
    character(len=*) , intent(in)    :: filename ! file
    type(ESMF_VM)                    :: vm
    real(r8)         , intent(inout) :: rdata    ! real data
    character(len=*) , intent(in)    :: dname    ! name of data
    integer          , intent(out)   :: rc

    ! local variables
    real(r8) :: r1d(1)
    character(*),parameter :: subName = '(med_io_read_r8) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call med_io_read_r81d(filename, vm, r1d,dname, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    rdata = r1d(1)

  end subroutine med_io_read_r8

  !===============================================================================
  subroutine med_io_read_r81d(filename, vm, rdata, dname, rc)

    !---------------
    ! Read 1d double array from netcdf file
    !---------------

    use pio , only : file_desc_t, var_desc_t, pio_openfile, pio_closefile, pio_seterrorhandling
    use pio , only : PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, pio_inq_varid, pio_get_var
    use pio , only : pio_nowrite, pio_openfile, pio_global, pio_get_att

    ! input/output arguments
    character(len=*), intent(in)    :: filename ! file
    type(ESMF_VM)                   :: vm
    real(r8)        , intent(inout) :: rdata(:) ! real data
    character(len=*), intent(in)    :: dname    ! name of data
    integer         , intent(out)   :: rc

    ! local variables
    integer           :: rcode
    type(file_desc_T) :: pioid
    type(var_desc_t)  :: varid
    character(CL)     :: lversion
    character(CL)     :: name1
    integer           :: iam
    character(*),parameter :: subName = '(med_io_read_r81d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lversion=trim(version)

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (med_io_file_exists(vm, filename)) then
       rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename),pio_nowrite)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call ESMF_LogWrite(trim(subname)//'ERROR: file invalid '//trim(filename)//' '//trim(dname), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,rdata)

    call pio_closefile(pioid)
  end subroutine med_io_read_r81d

  !===============================================================================
  subroutine med_io_read_char(filename, vm, rdata, dname, rc)

    !---------------
    ! Read char string from netcdf file
    !---------------

    use pio , only : file_desc_t, var_desc_t, pio_seterrorhandling, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR
    use pio , only : pio_closefile, pio_inq_varid, pio_get_var
    use pio , only : pio_openfile, pio_global, pio_get_att, pio_nowrite

    ! input/output arguments
    character(len=*), intent(in)    :: filename ! file
    type(ESMF_VM)                   :: vm
    character(len=*), intent(inout) :: rdata    ! character data
    character(len=*), intent(in)    :: dname    ! name of data
    integer         , intent(out)   :: rc

    ! local variables
    integer           :: rcode
    type(file_desc_T) :: pioid
    type(var_desc_t)  :: varid
    character(CL)     :: lversion
    character(CL)     :: name1
    integer           :: iam
    character(CL)     :: charvar   ! buffer for string read/write
    character(*),parameter :: subName = '(med_io_read_char) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lversion=trim(version)

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (med_io_file_exists(vm, filename)) then
       rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename),pio_nowrite)
       ! write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call ESMF_LogWrite(trim(subname)//'ERROR: file invalid '//trim(filename)//' '//trim(dname), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,charvar)
    rdata = trim(charvar)

    call pio_closefile(pioid)
  end subroutine med_io_read_char

  !===============================================================================
  subroutine med_io_date2ymd_int (date,year,month,day)
    ! Converts coded-date (yyyymmdd) to year/month/day.
    ! input/output variables
    integer,intent(in)  :: date             ! coded-date (yyyymmdd)
    integer,intent(out) :: year,month,day   ! calendar year,month,day
    ! local variables
    integer :: tdate   ! temporary date
    !-------------------------------------------------------------------------------

    tdate = abs(date)
    year =int(tdate/10000)
    if (date < 0) year = -year
    month = int( mod(tdate,10000)/  100)
    day = mod(tdate,  100)
  end subroutine med_io_date2ymd_int

  subroutine med_io_date2ymd_long (date,year,month,day)
    ! Converts coded-date (yyyymmdd) to year/month/day.
    ! input/output variables
    integer(I8),intent(in)  :: date             ! coded-date ([yy]yyyymmdd)
    integer    ,intent(out) :: year,month,day   ! calendar year,month,day
    ! local variables
    integer(I8) :: tdate   ! temporary date
    character(*),parameter :: subName = "(med_io_date2ymd_long)"
    !-------------------------------------------------------------------------------

    tdate = abs(date)
    year =int(tdate/10000)
    if (date < 0) year = -year
    month = int( mod(tdate,10000_I8)/  100)
    day = int(mod(tdate,  100_I8))
  end subroutine med_io_date2ymd_long

  !===============================================================================
  subroutine med_io_datetod2string_int(date_str, ymd, tod)
    ! Converts coded date (yyyymmdd) and optional time of day to a string like
    ! 'yyyy-mm-dd-ttttt' (if tod is present) or 'yyyy-mm-dd' (if tod is absent).
    ! yyyy in the output string will have at least 4 but no more than 6 characters (with
    ! leading zeroes if necessary).

    ! input/output variables
    character(len=*) , intent(out) :: date_str
    integer          , intent(in)  :: ymd
    integer, optional, intent(in)  :: tod

    ! local variables
    integer          :: yy, mm, dd
    character(len=6) :: year_str
    character(len=3) :: month_str
    character(len=3) :: day_str
    character(len=6) :: time_str
    !---------------------------------------

    call med_io_date2ymd(ymd, yy, mm, dd)

    ! Convert year, month, day and time of day to a string like 'yyyy-mm-dd-ttttt'.
    ! yyyy in the output string will have at least 4 but no more than 6 characters (with
    ! leading zeroes if necessary).
    write(year_str,'(i6.4)') yy
    year_str = adjustl(year_str)
    write(month_str,'(a,i2.2)') '-',mm
    write(day_str  ,'(a,i2.2)') '-',dd
    if (present(tod)) then
       write(time_str,'(a,i5.5)') '-',tod
    else
       time_str = ' '
    end if
    date_str = trim(year_str) // trim(month_str) // trim(day_str) // trim(time_str)

  end subroutine med_io_datetod2string_int

  subroutine med_io_datetod2string_long(date_str, ymd, tod)
    ! Converts coded date (yyyymmdd) and optional time of day to a string like
    ! 'yyyy-mm-dd-ttttt' (if tod is present) or 'yyyy-mm-dd' (if tod is absent).
    ! yyyy in the output string will have at least 4 but no more than 6 characters (with
    ! leading zeroes if necessary).

    ! input/output variables
    character(len=*) , intent(out) :: date_str
    integer(i8)      , intent(in)  :: ymd
    integer, optional, intent(in)  :: tod

    ! local variables
    integer          :: yy, mm, dd
    character(len=6) :: year_str
    character(len=3) :: month_str
    character(len=3) :: day_str
    character(len=6) :: time_str
    !---------------------------------------

    call med_io_date2ymd(ymd, yy, mm, dd)

    ! Convert year, month, day and time of day to a string like 'yyyy-mm-dd-ttttt'.
    ! yyyy in the output string will have at least 4 but no more than 6 characters (with
    ! leading zeroes if necessary).
    write(year_str,'(i6.4)') yy
    year_str = adjustl(year_str)
    write(month_str,'(a,i2.2)') '-',mm
    write(day_str  ,'(a,i2.2)') '-',dd
    if (present(tod)) then
       write(time_str,'(a,i5.5)') '-',tod
    else
       time_str = ' '
    end if
    date_str = trim(year_str) // trim(month_str) // trim(day_str) // trim(time_str)

  end subroutine med_io_datetod2string_long

  !===============================================================================
  subroutine med_io_ymd2date_int(year,month,day,date)
    ! Converts  year, month, day to coded-date

    ! input/output variables
    integer,intent(in ) :: year,month,day  ! calendar year,month,day
    integer,intent(out) :: date            ! coded (yyyymmdd) calendar date
    !---------------------------------------

    ! NOTE: this calendar has a year zero (but no day or month zero)
    date = abs(year)*10000 + month*100 + day  ! coded calendar date
    if (year < 0) date = -date
  end subroutine med_io_ymd2date_int

  subroutine med_io_ymd2date_long(year,month,day,date)
    ! Converts  year, month, day to coded-date

    ! input/output variables
    integer    ,intent(in ) :: year,month,day  ! calendar year,month,day
    integer(I8),intent(out) :: date            ! coded ([yy]yyyymmdd) calendar date
    !---------------------------------------

    ! NOTE: this calendar has a year zero (but no day or month zero)
    date = abs(year)*10000_I8 + month*100 + day  ! coded calendar date
    if (year < 0) date = -date
  end subroutine med_io_ymd2date_long

end module med_io_mod
