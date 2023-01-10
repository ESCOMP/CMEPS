module esm_utils_mod

  implicit none
  public

  logical :: maintask
  integer :: logunit
  integer :: dbug_flag = 0

  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  logical function ChkErr(rc, line, file, mpierr)
#ifndef NO_MPI2
    use mpi, only : MPI_ERROR_STRING, MPI_MAX_ERROR_STRING, MPI_SUCCESS
#else
    use mpi, only : MPI_SUCCESS
#endif
    use ESMF, only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO
    use ESMF, only : ESMF_FAILURE, ESMF_LogWrite

    integer, intent(in) :: rc
    integer, intent(in) :: line

    character(len=*), intent(in) :: file
    logical, optional, intent(in) :: mpierr
#ifdef NO_MPI2
    integer, parameter :: MPI_MAX_ERROR_STRING=80
#endif
    character(MPI_MAX_ERROR_STRING) :: lstring
    integer :: dbrc, lrc, len, ierr

    ChkErr = .false.
    lrc = rc
    if (present(mpierr)) then
       if(mpierr) then
          if (rc == MPI_SUCCESS) return
#ifdef USE_MPI2
          call MPI_ERROR_STRING(rc, lstring, len, ierr)
#else
          write(lstring,*) "ERROR in mct mpi-serial library rc=",rc
#endif
          call ESMF_LogWrite("ERROR: "//trim(lstring), ESMF_LOGMSG_INFO, line=line, file=file, rc=dbrc)
          lrc = ESMF_FAILURE
       endif
    endif

    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       ChkErr = .true.
    endif
  end function ChkErr

end module esm_utils_mod
