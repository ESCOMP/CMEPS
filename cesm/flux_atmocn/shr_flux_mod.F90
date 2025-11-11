module shr_flux_mod

  ! constants for atm/ocn/flux calculations

  use shr_kind_mod,  only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN ! shared kinds
  use shr_const_mod, only : shr_const_zvir, shr_const_cpdair, shr_const_cpvir, shr_const_karman, shr_const_g ! shared constants
  use shr_const_mod, only : shr_const_latvap, shr_const_latice, shr_const_stebol, shr_const_tkfrz
  use shr_const_mod, only : shr_const_ocn_ref_sal, shr_const_zsrflyr, shr_const_rgas
  use shr_sys_mod,   only : shr_sys_abort   ! shared system routines

  implicit none
  public

  public :: shr_flux_adjust_constants ! adjust constant values used in flux calculations. (used by CAM as well)

  integer, parameter :: debug = 0 ! internal debug level

  ! The follow variables are not declared as parameters so that they can be
  ! adjusted to support aquaplanet and potentially other simple model modes.
  ! The flux_adjust_constants subroutine is called to set the desired
  ! values.  The default values are from shr_const_mod.  Currently they are
  ! only used by the flux_atmocn routine.

  real(R8) :: loc_zvir   = shr_const_zvir
  real(R8) :: loc_cpdair = shr_const_cpdair
  real(R8) :: loc_cpvir  = shr_const_cpvir
  real(R8) :: loc_karman = shr_const_karman
  real(R8) :: loc_g      = shr_const_g
  real(R8) :: loc_latvap = shr_const_latvap
  real(R8) :: loc_latice = shr_const_latice
  real(R8) :: loc_stebol = shr_const_stebol
  real(R8) :: loc_tkfrz  = shr_const_tkfrz

  ! These control convergence of the iterative flux calculation
  ! (For Large and Pond scheme only; not UA or COARE).
  real(r8) :: flux_con_tol = 0.0_R8
  integer  :: flux_con_max_iter = 2

  !--- cold air outbreak parameters  (Mahrt & Sun 1995,MWR) -------------
  logical :: use_coldair_outbreak_mod = .false.

  real(R8),parameter    :: alpha = 1.4_R8
  real(R8),parameter    :: maxscl =2._R8  ! maximum wind scaling for flux
  real(R8),parameter    :: td0 = -10._R8  ! start t-ts for scaling

!===============================================================================
contains
!===============================================================================

  subroutine shr_flux_adjust_constants( &
       zvir, cpair, cpvir, karman, gravit, &
       latvap, latice, stebol, &
       flux_convergence_tolerance, &
       flux_convergence_max_iteration, &
       coldair_outbreak_mod)

    ! Adjust local constants.  Used to support simple models.

    real(R8)    , optional, intent(in) :: zvir
    real(R8)    , optional, intent(in) :: cpair
    real(R8)    , optional, intent(in) :: cpvir
    real(R8)    , optional, intent(in) :: karman
    real(R8)    , optional, intent(in) :: gravit
    real(R8)    , optional, intent(in) :: latvap
    real(R8)    , optional, intent(in) :: latice
    real(R8)    , optional, intent(in) :: stebol
    real(r8)    , optional, intent(in) :: flux_convergence_tolerance
    integer(in) , optional, intent(in) :: flux_convergence_max_iteration
    logical     , optional, intent(in) :: coldair_outbreak_mod
    !----------------------------------------------------------------------------

    if (present(zvir))   loc_zvir   = zvir
    if (present(cpair))  loc_cpdair = cpair
    if (present(cpvir))  loc_cpvir  = cpvir
    if (present(karman)) loc_karman = karman
    if (present(gravit)) loc_g      = gravit
    if (present(latvap)) loc_latvap = latvap
    if (present(latice)) loc_latice = latice
    if (present(stebol)) loc_stebol = stebol
    if (present(flux_convergence_tolerance)) flux_con_tol = flux_convergence_tolerance
    if (present(flux_convergence_max_iteration)) flux_con_max_iter = flux_convergence_max_iteration
    if (present(coldair_outbreak_mod)) use_coldair_outbreak_mod = coldair_outbreak_mod

  end subroutine shr_flux_adjust_constants

end module shr_flux_mod
