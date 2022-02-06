module flux_atmocn_ccpp_mod

  use med_kind_mod,    only : R8=>SHR_KIND_R8
  use physcons,        only : p0 => con_p0
  use physcons,        only : cappa => con_rocp
  use MED_data,        only : physics 
  use med_ccpp_driver, only : med_ccpp_driver_init
  use med_ccpp_driver, only : med_ccpp_driver_run
  use med_ccpp_driver, only : med_ccpp_driver_finalize

  implicit none

  private ! default private

  public :: flux_atmOcn_ccpp ! computes atm/ocn fluxes

!===============================================================================
contains
!===============================================================================

  subroutine flux_atmOcn_ccpp(nMax, mask, psfc, pbot, tbot, qbot, zbot, &
             garea, ubot, usfc, vbot, vsfc, rbot, ts, lwdn, sen, lat, &
             lwup, evp, taux, tauy, qref, missval)

    implicit none

    !--- input arguments --------------------------------
    integer , intent(in)  :: nMax        ! data vector length
    integer , intent(in)  :: mask (nMax) ! ocn domain mask
    real(r8), intent(in)  :: psfc(nMax)  ! atm P (surface)                (Pa)
    real(r8), intent(in)  :: pbot(nMax)  ! atm P (bottom)                 (Pa)
    real(r8), intent(in)  :: tbot(nMax)  ! atm T (bottom)                 (K)
    real(r8), intent(in)  :: qbot(nMax)  ! atm specific humidity (bottom) (kg/kg)
    real(r8), intent(in)  :: zbot(nMax)  ! atm level height               (m)
    real(r8), intent(in)  :: garea(nMax) ! grid area                      (m^2)
    real(r8), intent(in)  :: ubot(nMax)  ! atm u wind (bottom)            (m/s)
    real(r8), intent(in)  :: usfc(nMax)  ! atm u wind (surface)           (m/s)
    real(r8), intent(in)  :: vbot(nMax)  ! atm v wind (bottom)            (m/s)    
    real(r8), intent(in)  :: vsfc(nMax)  ! atm v wind (surface)           (m/s)    
    real(r8), intent(in)  :: rbot(nMax)  ! atm density                    (kg/m^3)    
    real(r8), intent(in)  :: lwdn(nMax)  ! atm lw downward                (W/m^2)
    real(r8), intent(in)  :: ts(nMax)    ! ocn surface temperature        (K)
    real(r8), intent(in), optional :: missval ! masked value

    !--- output arguments -------------------------------
    real(r8), intent(out) :: sen(nMax)   ! heat flux: sensible            (W/m^2)
    real(r8), intent(out) :: lat(nMax)   ! heat flux: latent              (W/m^2)
    real(r8), intent(out) :: lwup(nMax)  ! heat flux: lw upward           (W/m^2)
    real(r8), intent(out) :: evp(nMax)   ! heat flux: evap                ((kg/s)/m^2)
    real(r8), intent(out) :: taux(nMax)  ! surface stress, zonal          (N)
    real(r8), intent(out) :: tauy(nMax)  ! surface stress, maridional     (N)
    real(r8), intent(out) :: qref(nMax)  ! diag: 2m ref humidity          (kg/kg)

    !--- local variables --------------------------------
    logical, save :: first_call = .true.
    character(len=*),parameter :: subname=' (flux_atmOcn_ccpp) '
    !---------------------------------------

    if (first_call) then
       ! allocate and initalize data structures
       call physics%statein%create(nMax)
       call physics%interstitial%create(nMax)
       call physics%coupling%create(nMax)
       call physics%grid%create(nMax)
       call physics%sfcprop%create(nMax)

       ! initalize dimension 
       physics%init%im = nMax

       ! initalize model related parameters
       ! TODO: part of these need to be ingested from FV3 input.nml or configured through ESMF config file
       call physics%model%init()

       ! call CCPP init
       ! TODO: suite name need to be provided by ESMF config file
       call med_ccpp_driver_init('FV3_sfc_ocean')
       first_call = .false.
    end if

    ! fill in atmospheric forcing
    physics%statein%pgr(:)   = psfc(:)
    physics%statein%ugrs(:)  = ubot(:)
    physics%statein%vgrs(:)  = vbot(:)
    physics%statein%qgrs(:)  = qbot(:)
    physics%statein%prsl(:)  = pbot(:)
    physics%statein%zlvl(:)  = zbot(:)
    physics%statein%prsik(:) = (psfc(:)/p0)**cappa
    physics%statein%prslk(:) = (pbot(:)/p0)**cappa
    physics%statein%u10m(:)  = usfc(:)
    physics%statein%v10m(:)  = vsfc(:)

    ! fill in grid related variables
    physics%grid%area(:) = garea(:)

    ! customization of host model options to calculate the fluxes
    physics%model%lseaspray = .true.
    physics%model%ivegsrc = 1
    physics%model%redrag = .true.

  end subroutine flux_atmOcn_ccpp

end module flux_atmocn_ccpp_mod
