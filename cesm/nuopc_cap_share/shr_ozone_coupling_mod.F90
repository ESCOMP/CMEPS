module shr_ozone_coupling_mod

  !========================================================================
  ! Module for handling namelist variables related to ozone coupling
  !========================================================================

  use ESMF         , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMGet
  use ESMF         , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_SUCCESS
  use shr_sys_mod  , only : shr_sys_abort
  use shr_log_mod  , only : s_logunit => shr_log_Unit
  use shr_nl_mod   , only : shr_nl_find_group_name
  use shr_mpi_mod  , only : shr_mpi_bcast

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS
  public shr_ozone_coupling_readnl  ! Read namelist

  ! !PUBLIC DATA MEMBERS
  ! atm_ozone_frequency can be one of the following values
  integer, parameter, public :: atm_ozone_frequency_unset = 0
  integer, parameter, public :: atm_ozone_frequency_subdaily = 1
  integer, parameter, public :: atm_ozone_frequency_multiday_average = 2

  character(len=*), parameter :: &
       u_FILE_u=__FILE__

  !====================================================================================
CONTAINS
  !====================================================================================

  subroutine shr_ozone_coupling_readnl(NLFilename, atm_ozone_frequency_val)

    !========================================================================
    ! reads ozone_coupling_nl namelist and returns a variable specifying the frequency at
    ! which the atmosphere model computes surface ozone
    !========================================================================

    ! input/output variables
    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    ! atm_ozone_frequency will be one of the above constants (atm_ozone_frequency_*),
    ! specifying the frequency at which the atmosphere model computes surface ozone
    integer         , intent(out) :: atm_ozone_frequency_val

    !----- local -----
    character(len=64) :: atm_ozone_frequency
    type(ESMF_VM)     :: vm
    integer           :: unitn                  ! namelist unit number
    integer           :: ierr                   ! error code
    logical           :: exists                 ! if file exists or not
    integer           :: rc
    integer           :: localpet
    integer           :: mpicom

    character(len=*), parameter :: atm_ozone_frequency_not_present = 'NOT_PRESENT'
    character(len=*), parameter :: subname = '('//__FILE__//':shr_ozone_coupling_readnl)'
    ! ------------------------------------------------------------------

    namelist /ozone_coupling_nl/ atm_ozone_frequency

    rc = ESMF_SUCCESS

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0 ) then
       call shr_sys_abort( subname//'ERROR: nlfilename not set' )
    end if

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm, localPet=localpet, mpiCommunicator=mpicom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (localpet==0) then
       ! ------------------------------------------------------------------------
       ! Set default values in case namelist file doesn't exist, ozone_coupling_nl group
       ! doesn't exist within the file, or a given variable isn't present in the namelist
       ! group in the file.
       ! ------------------------------------------------------------------------
       atm_ozone_frequency = atm_ozone_frequency_not_present

       ! ------------------------------------------------------------------------
       ! Read namelist file
       ! ------------------------------------------------------------------------
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          open(newunit=unitn, file=trim(NLFilename), status='old' )
          write(s_logunit,'(a)') '(shr_ozone_coupling_readnl) Read in ozone_coupling_nl namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'ozone_coupling_nl', ierr)
          if (ierr == 0) then
             ! Note that ierr /= 0 means no namelist is present.
             read(unitn, ozone_coupling_nl, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort(trim(subname)//'problem reading ozone_coupling_nl ')
             end if
          end if
          close( unitn )
       end if

       ! ------------------------------------------------------------------------
       ! Translate read-in values to appropriate return values
       ! ------------------------------------------------------------------------
       select case(atm_ozone_frequency)
       case(atm_ozone_frequency_not_present)
          atm_ozone_frequency_val = atm_ozone_frequency_unset
       case("subdaily")
          atm_ozone_frequency_val = atm_ozone_frequency_subdaily
       case("multiday_average")
          atm_ozone_frequency_val = atm_ozone_frequency_multiday_average
       case default
          call shr_sys_abort(trim(subname)//'unknown value for atm_ozone_frequency: '// &
               trim(atm_ozone_frequency))
       end select
    end if

    ! ------------------------------------------------------------------------
    ! Broadcast values to all processors
    ! ------------------------------------------------------------------------
    call shr_mpi_bcast(atm_ozone_frequency_val, mpicom)

  end subroutine shr_ozone_coupling_readnl

end module shr_ozone_coupling_mod
