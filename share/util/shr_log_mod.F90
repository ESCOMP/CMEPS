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

  use shr_kind_mod
  use, intrinsic :: iso_fortran_env, only: output_unit

  implicit none
  private

  public :: shr_log_Level
  public :: shr_log_Unit

  ! low-level shared variables for logging, these may not be parameters
  integer(SHR_KIND_IN) :: shr_log_Level = 0
  integer(SHR_KIND_IN) :: shr_log_Unit  = output_unit

end module shr_log_mod
