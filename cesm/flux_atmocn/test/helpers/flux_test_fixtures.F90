module flux_test_fixtures

  !-----------------------------------------------------------------------------
  ! Shared setup, teardown, canned input states and tolerances for the
  ! cesm/flux_atmocn unit tests.
  !
  ! Every test's setUp should call fixture_init and every tearDown should call
  ! fixture_reset, so that tests cannot leak module state into each other.
  !-----------------------------------------------------------------------------

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use shr_const_mod,  only : shr_const_tkfrz, shr_const_tktrip
  use shr_const_mod,  only : shr_const_mwwv, shr_const_mwdair
  use shr_const_mod,  only : shr_const_zvir, shr_const_cpdair, shr_const_cpvir
  use shr_const_mod,  only : shr_const_karman, shr_const_g, shr_const_latvap
  use shr_const_mod,  only : shr_const_latice, shr_const_stebol
  use shr_wv_sat_mod, only : shr_wv_sat_init, shr_wv_sat_final
  use shr_flux_mod,   only : shr_flux_adjust_constants
  use funit,          only : throw

  implicit none
  private

  public :: fixture_init
  public :: fixture_reset
  public :: flux_inputs_t
  public :: make_unstable_state

  ! Tolerances
  real(r8), parameter, public :: tol_exact    = 1.e-13_r8 ! same computation, different path
  real(r8), parameter, public :: tol_ref      = 1.e-12_r8 ! known-answer reference values
  real(r8), parameter, public :: tol_physical = 1.e-6_r8  ! loose physical checks

  ! Missing value passed to the schemes for masked points
  real(r8), parameter, public :: test_spval = 1.e30_r8

  ! Minimum wind speed, as passed by med_phases_aofluxes_mod
  real(r8), parameter, public :: test_minwind = 0.5_r8

  ! Width of the liquid-ice transition range; driver default (wv_sat_transition_start)
  real(r8), parameter :: wv_sat_transition_start = 20.0_r8

  ! Inputs common to the bulk flux schemes, for nMax points
  type :: flux_inputs_t
     integer               :: nMax
     integer , allocatable :: mask(:)
     real(r8), allocatable :: zbot(:), ubot(:), vbot(:), thbot(:), qbot(:)
     real(r8), allocatable :: rainc(:), rbot(:), tbot(:), pslv(:)
     real(r8), allocatable :: us(:), vs(:), ts(:)
  end type flux_inputs_t

contains

  !-----------------------------------------------------------------------------
  subroutine fixture_init()
    ! Initialize shared module state the same way the CESM driver does (esm.F90).
    ! Without this, shr_wv_sat_mod's epsilo is uninitialized and qsat is bad.
    character(len=256) :: errstring

    call shr_wv_sat_init(shr_const_tkfrz, shr_const_tktrip, wv_sat_transition_start, &
         shr_const_mwwv/shr_const_mwdair, errstring)
    if (len_trim(errstring) > 0) call throw('fixture_init: '//trim(errstring))

  end subroutine fixture_init

  !-----------------------------------------------------------------------------
  subroutine fixture_reset()
    ! Restore water vapor saturation variables back to their default values,
    ! in case the test modifies them.

    call shr_wv_sat_final()
    call shr_flux_adjust_constants( &
         zvir=shr_const_zvir, cpair=shr_const_cpdair, cpvir=shr_const_cpvir, &
         karman=shr_const_karman, gravit=shr_const_g, latvap=shr_const_latvap, &
         latice=shr_const_latice, stebol=shr_const_stebol, &
         flux_convergence_tolerance=0.0_r8, flux_convergence_max_iteration=2, &
         coldair_outbreak_mod=.false.)

  end subroutine fixture_reset

  !-----------------------------------------------------------------------------
  pure function make_unstable_state(nMax) result(state)
    ! Warm ocean under cooler, dry air with moderate wind, repeated at nMax points.
    ! Expected: upward (negative, since fluxes are positive downward) sensible and
    ! latent heat flux, and evaporation.
    integer, intent(in) :: nMax
    type(flux_inputs_t) :: state

    state%nMax = nMax
    allocate(state%mask(nMax), state%zbot(nMax), state%ubot(nMax), state%vbot(nMax), &
         state%thbot(nMax), state%qbot(nMax), state%rainc(nMax), state%rbot(nMax), &
         state%tbot(nMax), state%pslv(nMax), state%us(nMax), state%vs(nMax), state%ts(nMax))

    state%mask  = 1
    state%zbot  = 10.0_r8      ! m
    state%ubot  = 8.0_r8       ! m/s
    state%vbot  = 3.0_r8       ! m/s
    state%tbot  = 290.0_r8     ! K
    state%thbot = 290.1_r8     ! K
    state%qbot  = 0.008_r8     ! kg/kg
    state%rainc = 0.0_r8       ! kg/m^2/s
    state%rbot  = 1.2_r8       ! kg/m^3
    state%pslv  = 101325.0_r8  ! Pa
    state%us    = 0.0_r8       ! m/s
    state%vs    = 0.0_r8       ! m/s
    state%ts    = 295.0_r8     ! K

  end function make_unstable_state

end module flux_test_fixtures
