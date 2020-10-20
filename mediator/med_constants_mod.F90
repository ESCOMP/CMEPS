module med_constants_mod

  use med_kind_mod, only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8

  implicit none
  public

  logical,  parameter :: med_constants_statewrite_flag = .false.
  real(R8), parameter :: med_constants_spval_init      = 0.0_R8  ! spval for initialization
  real(R8), parameter :: med_constants_spval           = 0.0_R8  ! spval
  real(R8), parameter :: med_constants_czero           = 0.0_R8  ! spval
  integer,  parameter :: med_constants_ispval_mask     = -987987 ! spval for RH mask values
  integer,  parameter :: med_constants_SecPerDay       = 86400   ! Seconds per day
  integer             :: med_constants_dbug_flag       = 2

end module med_constants_mod
