module shr_lightning_coupling_mod

  !========================================================================
  ! Module for handling namelist variables related to lightning coupling
  !========================================================================

  use ESMF         , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMGet
  use ESMF         , only : ESMF_LOGERR_PASSTHRU, ESMF_SUCCESS
  use ESMF         , only : ESMF_VMBroadCast, ESMF_Logical, assignment(=)
  use shr_sys_mod  , only : shr_sys_abort
  use shr_log_mod  , only : shr_log_getLogUnit
  use shr_nl_mod   , only : shr_nl_find_group_name
  use nuopc_shr_methods, only : chkerr

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS
  public shr_lightning_coupling_readnl  ! Read namelist

  character(len=*), parameter :: &
       u_FILE_u=__FILE__

  !====================================================================================
CONTAINS
  !====================================================================================

  subroutine shr_lightning_coupling_readnl(NLFilename, atm_provides_lightning_out)

    !========================================================================
    ! reads lightning_coupling_nl namelist and returns a variable specifying
    ! if atmosphere model provides lightning flash frequency field to mediator
    !========================================================================

    ! input/output variables
    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    logical, intent(out) :: atm_provides_lightning_out ! if TRUE atm will provide lightning flash frequency

    !----- local -----
    logical :: atm_provides_lightning
    type(ESMF_VM)     :: vm
    integer           :: unitn                  ! namelist unit number
    integer           :: ierr                   ! error code
    logical           :: exists                 ! if file exists or not
    type(ESMF_Logical):: ltmp(1)
    integer           :: rc
    integer           :: localpet
    integer           :: mpicom
    integer           :: s_logunit
    character(len=*), parameter :: atm_ozone_frequency_not_present = 'NOT_PRESENT'
    character(len=*), parameter :: subname = '(shr_lightning_coupling_readnl) '
    ! ------------------------------------------------------------------

    namelist /lightning_coupling_nl/ atm_provides_lightning

    rc = ESMF_SUCCESS

    atm_provides_lightning_out = .false.
    ltmp(1) = .false.

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0 ) then
       call shr_sys_abort( subname//'ERROR: nlfilename not set' )
    end if
    call shr_log_getLogUnit(s_logunit)
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localpet, mpiCommunicator=mpicom, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (localpet==0) then
       ! ------------------------------------------------------------------------
       ! Set default values in case namelist file doesn't exist, lightning_coupling_nl group
       ! doesn't exist within the file, or a given variable isn't present in the namelist
       ! group in the file.
       ! ------------------------------------------------------------------------
       atm_provides_lightning = .false.

       ! ------------------------------------------------------------------------
       ! Read namelist file
       ! ------------------------------------------------------------------------
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          open(newunit=unitn, file=trim(NLFilename), status='old' )
          write(s_logunit,'(a)') subname,'Read in lightning_coupling_nl namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'lightning_coupling_nl', ierr)
          if (ierr == 0) then
             ! Note that ierr /= 0 means no namelist is present.
             read(unitn, lightning_coupling_nl, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort(subname//'problem reading lightning_coupling_nl')
             end if
          end if
          close( unitn )
       end if

       ltmp(1) = atm_provides_lightning

    end if

    ! ------------------------------------------------------------------------
    ! Broadcast values to all tasks
    ! ------------------------------------------------------------------------
    call ESMF_VMBroadcast(vm,  ltmp, count=1, rootPet=0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    atm_provides_lightning_out = ltmp(1)

  end subroutine shr_lightning_coupling_readnl

end module shr_lightning_coupling_mod
