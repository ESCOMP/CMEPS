module shr_lnd2rof_tracers_mod

  !========================================================================
  ! read lnd2rof_tracers_inparm namelist and sets up driver list of fields for
  ! lnd -> river communications
  !========================================================================

  use ESMF         , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMGet, ESMF_VMBroadcast
  use ESMF         , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_SUCCESS
  use shr_sys_mod  , only : shr_sys_abort
  use shr_log_mod  , only : shr_log_getLogUnit
  use shr_kind_mod , only : r8 => shr_kind_r8, cs => shr_kind_cs
  use shr_nl_mod   , only : shr_nl_find_group_name

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS
  public :: shr_lnd2rof_tracers_readnl       ! Read namelist

  character(len=*), parameter :: &
       u_FILE_u=__FILE__

!====================================================================================
CONTAINS
!====================================================================================

  subroutine shr_lnd2rof_tracers_readnl(NLFilename, lnd2rof_tracer_list)

    ! input/output variables
    character(len=*), intent(in)  :: NLFilename          ! Namelist filename
    character(len=*), intent(out) :: lnd2rof_tracer_list ! Colon delimited string of liquid lnd2rof tracers

    !----- local -----
    type(ESMF_VM)      :: vm
    integer            :: unitn                  ! namelist unit number
    integer            :: ierr                   ! error code
    logical            :: exists                 ! if file exists or not
    integer            :: rc
    integer            :: localpet
    integer            :: mpicom
    integer            :: logunit
    character(len=CS)  :: lnd2rof_tracers
    character(*),parameter :: subName = '(shr_lnd2rof_tracers_readnl) '
    ! ------------------------------------------------------------------

    namelist /lnd2rof_tracers_inparm/ lnd2rof_tracers

    !-----------------------------------------------------------------------------
    ! Read namelist and figure out the lnd2rof_tracers field list to pass
    ! First check if file exists and if not, n_lnd2rof_tracers will be zero
    !-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0 ) then
       call shr_sys_abort( subName//'ERROR: nlfilename not set' )
    end if
    call shr_log_getLogUnit(logunit)

    lnd2rof_tracers = ' '
    lnd2rof_tracer_list = ' '

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=mpicom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (localpet==0) then
       inquire(file=trim(NLFileName), exist=exists)
       if ( exists ) then
          open(newunit=unitn, file=trim(NLFilename), status='old' )
          call shr_nl_find_group_name(unitn, 'lnd2rof_tracers_inparm', ierr)
          if (ierr == 0) then
             ! Note that if ierr /= 0, no namelist is present.
             read(unitn, lnd2rof_tracers_inparm, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort(trim(subName) //'problem of read of lnd2rof_tracers_inparm ')
             endif
          endif
          close( unitn )
       end if
    end if
    call ESMF_VMBroadcast(vm, lnd2rof_tracers, CS, 0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    
    if (lnd2rof_tracers /= ' ') then
       lnd2rof_tracer_list = trim(lnd2rof_tracers)
    end if

  end subroutine shr_lnd2rof_tracers_readnl

end module shr_lnd2rof_tracers_mod
