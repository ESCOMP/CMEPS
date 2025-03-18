module shr_abort_mod
  ! This module defines procedures that can be used to abort the model cleanly in a
  ! system-specific manner
  !
  ! The public routines here are only meant to be used directly by shr_sys_mod. Other code
  ! that wishes to use these routines should use the republished names from shr_sys_mod
  ! (shr_sys_abort, shr_sys_backtrace). (This is for consistency with older code, from
  ! when these routines were defined in shr_sys_mod.)

  use shr_kind_mod, only : shr_kind_in, shr_kind_cx

#ifdef CPRNAG
  ! NAG does not provide this as an intrinsic, but it does provide modules
  ! that implement commonly used POSIX routines.
  use f90_unix_proc, only: abort
#endif

  implicit none

  ! PUBLIC: Public interfaces

  private

  ! The public routines here are only meant to be used directly by shr_sys_mod. Other code
  ! that wishes to use these routines should use the republished names from shr_sys_mod
  ! (shr_sys_abort, shr_sys_backtrace). (This is for consistency with older code, from
  ! when these routines were defined in shr_sys_mod.)
  public :: shr_abort_abort     ! abort a program
  public :: shr_abort_backtrace ! print a backtrace, if possible

contains

  !===============================================================================
  subroutine shr_abort_abort(string,rc, line, file)
    use esmf, only : ESMF_LOGWRITE, ESMF_LOGMSG_ERROR, ESMF_FINALIZE, ESMF_END_ABORT
    use shr_log_mod, only : shr_log_error
    ! Consistent stopping mechanism

    !----- arguments -----
    character(len=*)    , intent(in), optional :: string  ! error message string
    integer(shr_kind_in), intent(in), optional :: rc      ! error code
    integer(shr_kind_in), intent(in), optional :: line
    character(len=*), intent(in), optional :: file
   
    ! Local version of the string.
    ! (Gets a default value if string is not present.)
    character(len=shr_kind_cx) :: local_string
    integer :: lrc
    !-------------------------------------------------------------------------------

    if (present(string)) then
       local_string = trim(string)
    else
       local_string = "Unknown error submitted to shr_abort_abort."
    end if
    if(present(rc)) then
       write(local_string, *) trim(local_string), ' rc=',rc
       lrc = rc
    else
       lrc = 0
    endif

    call shr_log_error(local_string, rc=lrc, line=line, file=file)
    
    call shr_abort_backtrace()

    call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! A compiler's abort method may print a backtrace or do other nice
    ! things, but in fact we can rarely leverage this, because MPI_Abort
    ! usually sends SIGTERM to the process, and we don't catch that signal.
    call abort()

  end subroutine shr_abort_abort
  !===============================================================================

  !===============================================================================
  subroutine shr_abort_backtrace()
    ! This routine uses compiler-specific facilities to print a backtrace to
    ! error_unit (standard error, usually unit 0).

#if defined(CPRIBM)

    ! This theoretically should be in xlfutility, but using it from that
    ! module doesn't seem to always work.
    interface
       subroutine xl_trbk()
       end subroutine xl_trbk
    end interface

    call xl__trbk()

#elif defined(CPRGNU) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8 ))

    ! gfortran 4.8 and later implement this intrinsic. We explicitly call it
    ! out as such to make sure that it really is available, just in case the
    ! CPP logic above screws up.
    intrinsic :: backtrace

    call backtrace()

#elif defined(CPRINTEL)

    ! tracebackqq uses optional arguments, so *must* have an explicit
    ! interface.
    use ifcore, only: tracebackqq

    ! An exit code of -1 is a special value that prevents this subroutine
    ! from aborting the run.
    call tracebackqq(user_exit_code=-1)

#else

    ! Currently we have no means to request a backtrace from the NAG runtime,
    ! even though it is capable of emitting backtraces itself, if you use the
    ! "-gline" option.

    ! Similarly, PGI has a -traceback option, but no user interface for
    ! requesting a backtrace to be printed.

#endif

  end subroutine shr_abort_backtrace
  !===============================================================================

end module shr_abort_mod
