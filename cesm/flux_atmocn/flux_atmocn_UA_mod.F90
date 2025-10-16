module flux_atmocn_UA_mod

  !===============================================================================
  ! !DESCRIPTION:
  !
  !     Internal atm/ocn flux calculation
  !     using University of Arizona method.
  !
  !     Reference:
  !         Zeng, X., M. Zhao, and R.E. Dickinson, 1998: Intercomparison of Bulk
  !             Aerodynamic Algorithms for the Computation of Sea Surface Fluxes
  !             Using TOGA COARE and TAO Data. J. Climate, 11, 2628–2644,
  !             https://doi.org/10.1175/1520-0442(1998)011<2628%3AIOBAAF>2.0.CO%3B2
  !
  !     Equation numbers are from this paper.
  !
  ! !REVISION HISTORY:
  !     2017-Aug-28 - J. Reeves Eyre - code re-written for E3SM
  !     2018-Oct-30 - J. Reeves Eyre - bug fix and add convective gustiness.
  !     2019-May-08 - J. Reeves Eyre - remove convective gustiness
  !                   and add cold air outbreak modification.
  !===============================================================================

  use shr_kind_mod,   only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN ! shared kinds
  use shr_flux_mod,   only : td0, maxscl, alpha
  use shr_flux_mod,   only : loc_zvir, loc_tkfrz, loc_cpdair, loc_cpvir, loc_g
  use shr_flux_mod,   only : use_coldair_outbreak_mod, loc_karman, loc_stebol

  implicit none
  private

  public :: flux_atmOcn_UA

  ! private member functions:
  private :: psi_ua
  private :: qsat_ua
  private :: rough_ua

  integer, private :: debug = 0

contains

  subroutine flux_atmOcn_UA(                &
       logunit,  spval, nMax,               &
       zbot, ubot, vbot, thbot,             &
       qbot, rbot, tbot, us, vs, pslv,      &
       ts, mask, sen, lat, lwup, evap,      &
       taux, tauy, tref, qref,              &
       duu10n,  ustar_sv, re_sv, ssq_sv)

    !--- input arguments --------------------------------
    integer    ,intent(in) :: logunit
    real(R8)   ,intent(in) :: spval
    integer    ,intent(in) :: nMax        ! data vector length
    integer    ,intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
    real(R8)   ,intent(in) :: zbot (nMax) ! atm level height (m)
    real(R8)   ,intent(in) :: ubot (nMax) ! atm u wind (m/s)
    real(R8)   ,intent(in) :: vbot (nMax) ! atm v wind (m/s)
    real(R8)   ,intent(in) :: thbot(nMax) ! atm potential T (K)
    real(R8)   ,intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
    real(R8)   ,intent(in) :: rbot (nMax) ! atm air density (kg/m^3)
    real(R8)   ,intent(in) :: tbot (nMax) ! atm T (K)
    real(R8)   ,intent(in) :: pslv (nMax) ! sea level pressure (Pa)
    real(R8)   ,intent(in) :: us (nMax) ! ocn u-velocity (m/s)
    real(R8)   ,intent(in) :: vs (nMax) ! ocn v-velocity (m/s)
    real(R8)   ,intent(in) :: ts (nMax) ! ocn temperature (K)

    !--- output arguments -------------------------------
    real(R8),intent(out)  ::  sen (nMax)     ! heat flux: sensible (W/m^2)
    real(R8),intent(out)  ::  lat (nMax)     ! heat flux: latent (W/m^2)
    real(R8),intent(out)  ::  lwup (nMax)     ! heat flux: lw upward (W/m^2)
    real(R8),intent(out)  ::  evap (nMax)     ! water flux: evap ((kg/s)/m^2)
    real(R8),intent(out)  ::  taux (nMax)     ! surface stress, zonal (N)
    real(R8),intent(out)  ::  tauy (nMax)     ! surface stress, maridional (N)
    real(R8),intent(out)  ::  tref (nMax)     ! diag:  2m ref height T (K)
    real(R8),intent(out)  ::  qref (nMax)     ! diag:  2m ref humidity (kg/kg)
    real(R8),intent(out)  :: duu10n(nMax)     ! diag: 10m wind speed squared (m/s)^2

    real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
    real(R8),intent(out),optional :: re_sv (nMax) ! diag: sqrt of exchange coefficient (water)
    real(R8),intent(out),optional :: ssq_sv (nMax) ! diag: sea surface humidity (kg/kg)

    !--- local constants --------------------------------
    real(R8),parameter :: zetam = -1.574_R8 ! Very unstable zeta cutoff for momentum (-)
    real(R8),parameter :: zetat = -0.465_R8 ! Very unstable zeta cutoff for T/q (-)
    real(R8),parameter :: umin  = 0.1_R8    ! minimum wind speed (m/s)
    real(R8),parameter :: zref  = 10.0_R8   ! reference height (m)
    real(R8),parameter :: ztref = 2.0_R8    ! reference height for air T (m)
    real(R8),parameter :: beta = 1.0_R8     ! constant used in W* calculation (-)
    real(R8),parameter :: zpbl = 1000.0_R8  ! PBL height used in W* calculation (m)
    real(R8),parameter :: gamma = 0.0098_R8         ! Dry adiabatic lapse rate (K/m)
    real(R8),parameter :: onethird = 1.0_R8/3.0_R8  ! Used repeatedly.

    !--- local variables --------------------------------
    integer  :: n           ! vector loop index
    integer  :: i           ! iteration loop index
    real(R8) :: vmag_abs    ! surface wind magnitude (m s-1)
    real(R8) :: vmag_rel    ! surface wind magnitude relative to surface current (m s-1)
    real(R8) :: vmag        ! surface wind magnitude with large eddy correction and minimum value (m s-1)
                            ! (This can change on each iteration.)
    real(R8) :: thv         ! virtual temperature (K)
    real(R8) :: ssq         ! sea surface humidity (kg/kg)
    real(R8) :: delth       ! potential T difference (K)
    real(R8) :: delthv      ! virtual potential T difference (K)
    real(R8) :: delq        ! humidity difference (kg/kg)
    real(R8) :: ustar       ! friction velocity (m s-1)
    real(R8) :: qstar       ! humidity scaling parameter (kg/kg)
    real(R8) :: tstar       ! temperature scaling parameter (K)
    real(R8) :: thvstar     ! virtual temperature scaling parameter (K)
    real(R8) :: wstar       ! convective velocity scale (m s-1)
    real(R8) :: zeta        ! dimensionless height (z / Obukhov length)
    real(R8) :: obu         ! Obukhov length (m)
    real(R8) :: tau         ! magnitude of wind stress (N m-2)
    real(R8) :: cp          ! specific heat of moist air (J kg-1 K-1)
    real(R8) :: xlv         ! Latent heat of vaporization (J kg-1)
    real(R8) :: visa        ! Kinematic viscosity of dry air (m2 s-1)
    real(R8) :: tbot_oC     ! Temperature used in visa (deg C)
    real(R8) :: rb          ! Bulk Richardson number (-)
    real(R8) :: zo          ! Roughness length for momentum (m)
    real(R8) :: zoq         ! Roughness length for moisture (m)
    real(R8) :: zot         ! Roughness length for heat (m)
    real(R8) :: u10         ! 10-metre wind speed (m s-1)
    real(R8) :: re          ! Moisture exchange coefficient for compatibility with default algorithm.
    real(R8) :: loc_epsilon ! Ratio of gas constants (-)

    !--- for cold air outbreak calc --------------------------------
    real(R8)    :: tdiff(nMax)  ! tbot - ts
    real(R8)    :: vscl

    !--- formats ----------------------------------------
    character(*),parameter :: subName = '(flux_atmOcn) '
    character(*),parameter ::   F00 = "('(flux_atmOcn) ',4a)"
    !---------------------------------------------------------------------------

    !-----
    ! Straight from original subroutine.
    if (debug > 0) write(logunit,F00) "enter"

    ! Evaluate loc_epsilon.
    loc_epsilon = 1.0_R8 / (1.0_R8 + loc_zvir)

    !--- for cold air outbreak calc --------------------------------
    tdiff = tbot - ts

    ! Loop over grid points.
    do n=1,nMax

       if (mask(n) /= 0) then

          !-----Calculate some required near surface variables.---------
          vmag_abs = sqrt( ubot(n)**2 + vbot(n)**2 )
          vmag_rel = sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2 )

          ! For Cold Air Outbreak Modification (based on Mahrt & Sun 1995,MWR):
          if (use_coldair_outbreak_mod) then
             ! Increase windspeed for negative tbot-ts
             if (tdiff(n).lt.td0) then
                vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag_rel))),maxscl)
                vmag_rel=vmag_rel*vscl
             endif
          endif

          delth = thbot(n) - ts(n)                     ! Pot. temp. difference with surface (K)
          ! Note this is equivalent to Zeng et al
          ! (1998) version = delt + 0.0098*zbot
          thv = thbot(n)*(1.0_R8+0.61_R8*qbot(n))      ! Virtual potential temperature (K)
          ! EQN (17):
          !ssq = 0.98_R8 * qsat_ua(ts(n),ps, &          ! Surface specific humidity (kg kg-1)
          !                        loc_epsilon)
          ssq = 0.98_R8 * qsat_ua(ts(n),pslv(n), &     ! Surface specific humidity (kg kg-1)
               loc_epsilon)
          delq = qbot(n) - ssq                         ! Difference to surface (kg kg-1)
          delthv = delth*(1.0_R8+0.61_R8*qbot(n)) + &  ! Difference of virtual potential
               & 0.61_R8*thbot(n)*delq               ! temperature with surface (K)

          xlv = 1.0e+6_R8 * &                           ! Latent heat of vaporization (J kg-1)
               & (2.501_R8 - 0.00237_R8 * (ts(n) - loc_tkfrz))
          tbot_oC = tbot(n) - loc_tkfrz
          visa = 1.326e-5_R8 * (1.0_R8 + &             ! Kinematic viscosity of dry
               & 6.542e-3_R8*tbot_oC + &               ! air (m2 s-1) from Andreas (1989)
               & 8.301e-6_R8*tbot_oC*tbot_oC - &       ! CRREL Rep. 89-11
               & 4.84e-9_R8*tbot_oC*tbot_oC*tbot_oC)
          cp = loc_cpdair*(1.0_R8 + loc_cpvir*ssq)     ! specific heat of moist air (J kg-1 K-1)

          !-----Initial values of u* and convective velocity.-----------
          ustar = 0.06_R8
          wstar = 0.5_R8
          ! Update wind speed if unstable regime.
          if (delthv.lt.0.0_R8) then
             ! EQN (19)
             vmag = sqrt( vmag_rel**2 + beta*beta*wstar*wstar )
          else
             ! EQN (18)
             vmag = max(umin,vmag_rel)
          endif

          !-----Iterate to compute new u* and z0.-----------------------
          do i = 1,5
             ! EQN (24)
             zo = 0.013_R8*ustar*ustar/loc_g + 0.11_R8*visa/ustar
             ! EQN (9) assuming neutral
             ustar = loc_karman*vmag/log(zbot(n)/zo)
          enddo

          !-----Assess stability.---------------------------------------
          rb = loc_g*zbot(n)*delthv / (thv*vmag*vmag)    ! bulk Richardson number

          if(rb.ge.0.0_R8) then
             ! Neutral or stable: EQNs (4), (9), (13) and definition of rb.
             zeta = rb*log(zbot(n)/zo) / &
                  & (1.0_R8 - 5.0_R8*min(rb,0.19_R8))
          else
             ! Unstable: EQNs (4), (8), (12) and definition of rb.
             zeta = rb*log(zbot(n)/zo)
          endif

          obu = zbot(n)/zeta                             ! Obukhov length
          obu = sign(max(zbot(n)/10.0_R8, abs(obu)), obu)

          !-----Main iterations (2-10 iterations would be fine).-------
          do i=1,10

             ! Update roughness lengths.
             call rough_ua(zo,zot,zoq,ustar,visa)

             ! Wind variables.
             zeta = zbot(n) / obu
             if (zeta.lt.zetam) then
                ! Very unstable regime
                ! EQN (7) with extra z0 term.
                ustar = loc_karman * vmag / (log(zetam*obu/zo) - &
                     & psi_ua(1_IN, zetam) + &
                     & psi_ua(1_IN, zo/obu) + &
                     & 1.14_R8 * ((-zeta)**onethird - (-zetam)**onethird) )
             else if (zeta.lt.0.0_R8) then
                ! Unstable regime
                ! EQN (8) with extra z0 term.
                ustar = loc_karman * vmag / (log(zbot(n)/zo) - &
                     & psi_ua(1_IN,zeta) + psi_ua(1_IN,zo/obu) )
             else if (zeta.le.1.0_R8) then
                ! Stable regime
                ! EQN (9) with extra z0 term.
                ustar = loc_karman * vmag / (log(zbot(n)/zo) + &
                     & 5.0_R8*zeta - 5.0_R8*zo/obu)
             else
                ! Very stable regime
                ! EQN (10) with extra z0 term.
                ustar = loc_karman * vmag / (log(obu/zo) + 5.0_R8 - &
                     &  5.0_R8*zo/obu + &
                     &  (5.0_R8*log(zeta) + zeta - 1.0_R8) )
             endif

             ! Temperature variables.
             if(zeta.lt.zetat) then
                ! Very unstable regime
                ! EQN (11) with extra z0 term.
                tstar = loc_karman * delth / (log(zetat*obu/zot) - &
                     & psi_ua(2_IN, zetat) + &
                     & psi_ua(2_IN, zot/obu) + &
                     & 0.8_R8*((-zetat)**(-onethird) - (-zeta)**(-onethird)) )
             else if (zeta.lt.0.0_R8) then
                ! Unstable regime
                ! EQN (12) with extra z0 term.
                tstar = loc_karman * delth / &
                     & (log(zbot(n)/zot) - psi_ua(2_IN,zeta) + psi_ua(2_IN,zot/obu))
             else if (zeta.le.1.0_R8) then
                ! Stable regime
                ! EQN (13) with extra z0 term.
                tstar = loc_karman * delth / (log(zbot(n)/zot) + &
                     &   5.0_R8*zeta - 5.0_R8*zot/obu)
             else
                ! Very stable regime
                ! EQN (14) with extra z0 term.
                tstar = loc_karman * delth / (log(obu/zot) + &
                     &   5.0_R8 - 5.0_R8*zot/obu  + &
                     &   (5.0_R8*log(zeta) + zeta - 1.0_R8) )
             endif

             ! Humidity variables.
             ! This is done with re to give variable to save out like
             ! in old algorithm.
             if (zeta.lt.zetat) then
                ! Very unstable regime
                ! EQN (11) with extra z0 term.
                re = loc_karman / (log(zetat*obu/zoq) - psi_ua(2_IN,zetat) + &
                     & psi_ua(2_IN,zoq/obu) + &
                     & 0.8_R8*((-zetat)**(-onethird) - (-zeta)**(-onethird)) )
             else if (zeta.lt.0.0_R8) then
                ! Unstable regime
                ! EQN (12) with extra z0 term.
                re = loc_karman / &
                     & (log(zbot(n)/zoq) - psi_ua(2_IN,zeta) + psi_ua(2_IN,zoq/obu))
             else if (zeta.le.1.0_R8) then
                ! Stable regime
                ! EQN (13) with extra z0 term.
                re = loc_karman / &
                     & (log(zbot(n)/zoq) + 5.0_R8*zeta - 5.0_R8*zoq/obu)
             else
                ! Very stable regime
                ! EQN (14) with extra z0 term.
                re = loc_karman / &
                     & (log(obu/zoq) + 5.0_R8 - 5.0_R8*zoq/obu + &
                     & (5.0_R8*log(zeta) + zeta - 1.0_R8) )
             endif
             qstar = re * delq

             ! Update Obukhov length.
             thvstar = tstar*(1.0_R8 + 0.61_R8*qbot(n)) + 0.61_R8*thbot(n)*qstar
             ! EQN (4)
             obu = ustar*ustar * thv / (loc_karman*loc_g*thvstar)
             obu = sign( max(zbot(n)/10.0_R8, abs(obu)) ,obu)

             ! Update wind speed if in unstable regime.
             if (delthv.lt.0.0_R8) then
                ! EQN (20)
                wstar = beta * (-loc_g*ustar*thvstar*zpbl/thv)**onethird
                ! EQN (19)
                vmag = sqrt(vmag_rel**2 + wstar*wstar)
             else
                ! EQN (18)
                vmag = max(umin,vmag_rel)
             endif

          enddo ! End of iterations for ustar, tstar, qstar etc.


          !-----Calculate fluxes and wind stress.---------------------

          !--- momentum flux ---
          ! This should ensure zero wind stress when (relative) wind speed is zero,
          ! components are consistent with total, and we don't ever divide by zero.
          ! EQN (21)
          tau = rbot(n) * ustar * ustar
          taux(n) = tau * (ubot(n)-us(n)) / max(umin, vmag_rel)
          tauy(n) = tau * (vbot(n)-vs(n)) / max(umin, vmag_rel)

          !--- heat flux ---
          ! EQNs (22) and (23)
          sen (n) =  cp * rbot(n) * tstar * ustar
          lat (n) = xlv * rbot(n) * qstar * ustar
          lwup(n) = -loc_stebol * ts(n)**4

          !--- water flux ---
          evap(n) = lat(n)/xlv

          !------------------------------------------------------------
          ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------

          zeta = zbot(n) / obu
          if (zeta.lt.zetat) then
             if (zeta.lt.zetam) then
                ! Very unstable regime for U.
                ! EQN (7)
                u10 = vmag_abs + (ustar/loc_karman) * &
                     & 1.14_R8 * ((-zref/obu)**onethird - (-zeta)**onethird)
             else
                ! Unstable regime for U.
                ! EQN (8)
                u10 = vmag_abs + (ustar/loc_karman) * &
                     & (log(zref/zbot(n)) - (psi_ua(1_IN,zref/obu) - psi_ua(1_IN,zeta)) )
             endif
             ! Very unstable regime for T and q.
             ! EQN (11)
             tref(n) = thbot(n) + (tstar/loc_karman) * &
                  & 0.8_R8 * ((-zeta)**(-onethird) - (-ztref/obu)**(-onethird))
             qref(n) = qbot(n) + (qstar/loc_karman) * &
                  & 0.8_R8 * ((-zeta)**(-onethird) - (-ztref/obu)**(-onethird))

          else if (zeta.lt.0.0_R8) then
             ! Unstable regime.
             ! EQN (8)
             u10 = vmag_abs + (ustar/loc_karman) * &
                  & (log(zref/zbot(n)) - (psi_ua(1_IN,zref/obu) - psi_ua(1_IN,zeta)) )
             ! EQN (12)
             tref(n) = thbot(n) + (tstar/loc_karman) * &
                  & (log(ztref/zbot(n)) - (psi_ua(2_IN,ztref/obu) - psi_ua(2_IN,zeta)) )
             qref(n) = qbot(n) + (qstar/loc_karman) * &
                  & (log(ztref/zbot(n)) - (psi_ua(2_IN,ztref/obu) - psi_ua(2_IN,zeta)) )
          else if (zeta.le.1.0_R8) then
             ! Stable regime.
             ! EQN (9)
             u10 = vmag_abs + (ustar/loc_karman) * &
                  & (log(zref/zbot(n)) + 5.0_R8*zref/obu - 5.0_R8*zeta)
             ! EQN (13)
             tref(n) = thbot(n) + (tstar/loc_karman) * &
                  & (log(ztref/zbot(n)) + 5.0_R8*ztref/obu - 5.0_R8*zeta)
             qref(n) = qbot(n) + (qstar/loc_karman) * &
                  & (log(ztref/zbot(n)) + 5.0_R8*ztref/obu - 5.0_R8*zeta)
          else
             ! Very stable regime.
             ! EQN (10)
             u10 = vmag_abs + (ustar/loc_karman) * &
                  & (5.0_R8*log(zref/zbot(n)) + zref/obu - zeta)
             ! EQN (14)
             tref(n) = thbot(n) + (tstar/loc_karman) * &
                  & (5.0_R8*log(ztref/zbot(n)) + ztref/obu - zeta)
             qref(n) = qbot(n) + (qstar/loc_karman) * &
                  & (5.0_R8*log(ztref/zbot(n)) + ztref/obu - zeta)

          endif

          tref(n) = tref(n) - gamma*ztref   ! pot. temp to temp correction
          duu10n(n) = u10*u10 ! 10m wind speed squared

          !------------------------------------------------------------
          ! optional diagnostics, needed for water tracer fluxes (dcn)
          !------------------------------------------------------------
          if (present(ustar_sv)) ustar_sv(n) = ustar
          if (present(ssq_sv  )) ssq_sv(n)   = ssq
          if (present(re_sv   )) re_sv(n)    = re

       else

          !------------------------------------------------------------
          ! no valid data here -- out of ocean domain
          !------------------------------------------------------------
          sen   (n) = spval  ! sensible         heat flux  (W/m^2)
          lat   (n) = spval  ! latent           heat flux  (W/m^2)
          lwup  (n) = spval  ! long-wave upward heat flux  (W/m^2)
          evap  (n) = spval  ! evaporative water flux ((kg/s)/m^2)
          evap_16O (n) = spval !water tracer flux (kg/s)/m^2)
          evap_HDO (n) = spval !HDO tracer flux  (kg/s)/m^2)
          evap_18O (n) = spval !H218O tracer flux (kg/s)/m^2)
          taux  (n) = spval  ! x surface stress (N)
          tauy  (n) = spval  ! y surface stress (N)
          tref  (n) = spval  !  2m reference height temperature (K)
          qref  (n) = spval  !  2m reference height humidity (kg/kg)
          duu10n(n) = spval  ! 10m wind speed squared (m/s)^2

          ! Optional diagnostics too:
          if (present(ustar_sv)) ustar_sv(n) = spval
          if (present(re_sv   )) re_sv   (n) = spval
          if (present(ssq_sv  )) ssq_sv  (n) = spval

       endif

    enddo ! loop over grid points

  end subroutine flux_atmOcn_UA


  !===============================================================================

  real(R8) function psi_ua(k,zeta)

    ! Stability function for rb < 0

    !-----Input variables.----------
    integer(IN), intent(in) :: k       ! Indicates whether this is for momentum (k=1)
    ! or for heat/moisture (k=2)
    real(R8), intent(in) :: zeta       ! Dimensionless height (=z/L)

    !-----Local variables.----------
    real(R8) :: chik                   ! Function of zeta.

    ! EQN (16)
    chik = (1.0_R8 - 16.0_R8*zeta)**0.25_R8

    if(k.eq.1) then
       ! EQN (15) for momentum
       psi_ua = 2.0_R8 * log((1.0_R8 + chik)*0.5_R8) + &
            &      log((1.0_R8 + chik*chik)*0.5_R8) - &
            & 2.0_R8 * atan(chik) + 2.0_R8 * atan(1.0_R8)
    else
       ! EQN (15) for heat/moisture
       psi_ua = 2.0_R8 * log((1.0_R8 + chik*chik)*0.5_R8)
    endif

  end function psi_ua

  !===============================================================================

  real(R8) function qsat_ua(t,p,loc_epsilon)

    ! Uses Tetens' formula for saturation vapor pressure from
    ! Buck(1981) JAM 20, 1527-1532

    !-----Input variables.----------
    real(R8), intent(in) :: t           ! temperature (K)
    real(R8), intent(in) :: p           ! pressure (Pa)
    real(R8), intent(in) :: loc_epsilon ! Ratio of gas constants (-)

    !-----Local variables.----------
    real(R8) :: esat                    ! saturated vapor pressure (hPa)

    ! Calculate saturated vapor pressure in hPa.
    esat = (1.0007_R8 + 0.00000346_R8 * (p/100.0_R8)) * 6.1121_R8 * &
         & exp(17.502_R8 * (t - loc_tkfrz) / (240.97_R8 + (t - loc_tkfrz)))

    ! Convert to specific humidity (kg kg-1).
    qsat_ua = loc_epsilon * esat / ((p/100.0_R8) - (1.0_R8 - loc_epsilon)*esat)

  end function qsat_ua

  !===============================================================================

  subroutine rough_ua(zo,zot,zoq,ustar,visa)

    ! Calculate roughness lengths: zo, zot, zoq.

    !-----Input variables.----------
    real(R8), intent(in) :: ustar      ! friction velocity (m s-1)
    real(R8), intent(in) :: visa       ! kinematic viscosity of dry air (m2 s-1)

    !-----Output variables.---------
    real(R8), intent(out) :: zo        ! roughness length for momentum (m)
    real(R8), intent(out) :: zot       ! roughness length for heat (m)
    real(R8), intent(out) :: zoq       ! roughness length for water vapor (m)

    !-----Local variables.----------
    real(R8) :: re_rough               ! Rougness Reynold's number (-)
    real(R8) :: xq                     ! Logarithm of roughness length ratios (moisture)
    real(R8) :: xt                     ! Logarithm of roughness length ratios (heat)

    zo = 0.013_R8*ustar*ustar/loc_g + 0.11_R8*visa/ustar      ! EQN (24)
    re_rough = ustar*zo/visa                                  ! By definition.
    xq = 2.67_R8*re_rough**0.25_R8 - 2.57_R8                  ! EQN (25)
    xt = xq                                                   ! EQN (26)
    zoq = zo/exp(xq)                                          ! By definition of xq
    zot = zo/exp(xt)                                          ! By definition of xt

  end subroutine rough_ua

end module flux_atmocn_UA_mod
