!BOP ===========================================================================
!
! !MODULE: shr_log_mod -- variables and methods for logging
!
! !DESCRIPTION:
!    Low-level shared variables for logging.
!
!    Also, routines for generating log file messages.
!
! !INTERFACE: ------------------------------------------------------------------

module shr_log_mod

! !USES:

  use shr_kind_mod, only: shr_kind_in, shr_kind_cx
  use shr_strconvert_mod, only: toString

  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

  implicit none
  private

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_log_errMsg
  public :: shr_log_OOBMsg
  public :: shr_log_setLogUnit
  public :: shr_log_getLogUnit
  public :: shr_log_error

! !PUBLIC DATA MEMBERS:

  public :: shr_log_Level
  public :: shr_log_Unit

!EOP

  ! low-level shared variables for logging, these may not be parameters
  integer(SHR_KIND_IN) :: shr_log_Level = 0
  integer(SHR_KIND_IN) :: shr_log_Unit  = output_unit

contains

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_log_errMsg -- Return an error message containing file & line info
!
! !DESCRIPTION:
!     Return an error message containing file & line info
!     \newline
!     errMsg = shr\_log\_errMsg(__FILE__, __LINE__)
!
! This is meant to be used when a routine expects a string argument for some message,
! but you want to provide file and line information.
!
! However: Note that the performance of this function can be very bad. It is currently
! maintained because it is used by old code, but you should probably avoid using this
! in new code if possible.
!
! !REVISION HISTORY:
!     2013-July-23 - Bill Sacks
!
! !INTERFACE: ------------------------------------------------------------------

  pure function shr_log_errMsg(file, line)

    ! !INPUT/OUTPUT PARAMETERS:

    character(len=SHR_KIND_CX)   :: shr_log_errMsg
    character(len=*), intent(in) :: file
    integer         , intent(in) :: line

    !EOP

    shr_log_errMsg = 'ERROR in '//trim(file)//' at line '//toString(line)

  end function shr_log_errMsg

  ! Create a message for an out of bounds error.
  pure function shr_log_OOBMsg(operation, bounds, idx) result(OOBMsg)

    ! A name for the operation being attempted when the bounds error
    ! occurred. A string containing the subroutine name is ideal, but more
    ! generic descriptions such as "read", "modify", or "insert" could be used.
    character(len=*), intent(in) :: operation

    ! Upper and lower bounds allowed for the operation.
    integer, intent(in) :: bounds(2)

    ! Index at which access was attempted.
    integer, intent(in) :: idx

    ! Output message
    character(len=:), allocatable :: OOBMsg

    allocate(OOBMsg, source=(operation//": "//toString(idx)//" not in range ["//&
         toString(bounds(1))//", "//toString(bounds(2))//"]."))

  end function shr_log_OOBMsg

  subroutine shr_log_setLogUnit(unit)
    integer, intent(in) :: unit

    shr_log_unit = unit

  end subroutine shr_log_setLogUnit

  subroutine shr_log_getLogUnit(unit)
    integer, intent(out) :: unit

     unit = shr_log_unit

  end subroutine shr_log_getLogUnit

  subroutine shr_log_error(string, rc, line, file)
    use esmf, only : ESMF_LOGWRITE, ESMF_LOGMSG_ERROR, ESMF_FINALIZE, ESMF_END_ABORT, ESMF_FAILURE, ESMF_SUCCESS
    ! This routine prints error messages to shr_log_unit (which is standard output
    ! for most tasks in CESM), to the ESMF PET files and to standard error if shr_log_unit is a
    ! file.  Sets rc to ESMF_FAILURE on return.

    !----- arguments -----
    character(len=*)    , intent(in) :: string  ! error message string
    integer(shr_kind_in), intent(inout), optional :: rc      ! error code
    integer(shr_kind_in), intent(in), optional :: line
    character(len=*), intent(in), optional :: file

    ! Local version of the string.
    ! (Gets a default value if string is not present.)
    character(len=shr_kind_cx) :: local_string
    integer, allocatable :: log_units(:)
    integer :: i
    !-------------------------------------------------------------------------------

    local_string = trim(string)
    if(present(rc)) then
       if (rc /= ESMF_SUCCESS) then
          write(local_string, *) trim(local_string), ' rc=',rc
       endif
       rc = ESMF_FAILURE
    endif

    call ESMF_LogWrite(local_string, ESMF_LOGMSG_ERROR, line=line, file=file)
    if (shr_log_unit == output_unit .or. shr_log_unit == error_unit) then
       ! If the log unit number is standard output or standard error, just
       ! print to that.
       allocate(log_units(1), source=[shr_log_unit])
    else
       ! Otherwise print the same message to both the log unit and standard
       ! error.
       allocate(log_units(2), source=[error_unit, shr_log_unit])
    end if

    do i = 1, size(log_units)
       write(log_units(i),*) trim(local_string)
       flush(log_units(i))
    end do

  end subroutine shr_log_error

end module shr_log_mod
