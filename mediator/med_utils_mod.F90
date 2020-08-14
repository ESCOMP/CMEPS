module med_utils_mod

  use med_kind_mod, only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8

  implicit none
  private

  public :: med_memcheck
  public :: med_utils_ChkErr
  public :: med_log_clock_advance

  integer     , parameter :: memdebug_level=1
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_memcheck(string, level, mastertask)
    character(len=*), intent(in) :: string
    integer, intent(in) :: level
    logical, intent(in) :: mastertask
    integer :: ierr
#ifndef INTERNAL_PIO_INIT
    integer, external :: GPTLprint_memusage
    if((mastertask .and. memdebug_level > level) .or. memdebug_level > level+1) then
       ierr = GPTLprint_memusage(string)
    endif
#endif
  end subroutine med_memcheck

!===============================================================================

  logical function med_utils_ChkErr(rc, line, file, mpierr)
#ifdef USE_MPI2
    use mpi , only : MPI_ERROR_STRING, MPI_MAX_ERROR_STRING, MPI_SUCCESS
#else
    use mpi, only : MPI_SUCCESS
#endif
    use ESMF, only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO
    use ESMF, only : ESMF_FAILURE, ESMF_LogWrite

    integer, intent(in) :: rc
    integer, intent(in) :: line

    character(len=*), intent(in) :: file
    logical, optional, intent(in) :: mpierr
#ifndef USE_MPI2
    integer, parameter :: MPI_MAX_ERROR_STRING=80
#endif
    character(MPI_MAX_ERROR_STRING) :: lstring
    integer :: lrc, len, ierr

    med_utils_ChkErr = .false.
    lrc = rc
    if (present(mpierr) .and. mpierr) then
       if (rc == MPI_SUCCESS) return
#ifdef USE_MPI2
       call MPI_ERROR_STRING(rc, lstring, len, ierr)
#else
       write(lstring,*) "ERROR in mct mpi-serial library rc=",rc
#endif
       call ESMF_LogWrite("ERROR: "//trim(lstring), ESMF_LOGMSG_INFO, line=line, file=file)
       lrc = ESMF_FAILURE
    endif

    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
      med_utils_ChkErr = .true.
    endif

  end function med_utils_ChkErr

!===============================================================================

  subroutine med_log_clock_advance(clock, component, logunit)
    use ESMF, only : ESMF_Clock, ESMF_ClockPrint

    type(ESMF_Clock) :: clock
    character(len=*), intent(in) :: component
    integer, intent(in) :: logunit

    character(len=CL) :: cvalue, prestring
    integer :: rc

    write(prestring, *) "------>Advancing ",trim(component)," from: "
    call ESMF_ClockPrint(clock, options="currTime", unit=cvalue, &
         preString=trim(prestring), rc=rc)
    if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

    call ESMF_ClockPrint(clock, options="stopTime", unit=cvalue, &
         preString="--------------------------------> to: ", rc=rc)
    if (med_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

  end subroutine med_log_clock_advance

end module med_utils_mod
