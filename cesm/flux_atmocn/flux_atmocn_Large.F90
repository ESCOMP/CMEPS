module flux_atmOcn_large_mod

  !-------------------------------------------------------------------------------
  ! PURPOSE:
  !   computes atm/ocn surface fluxes using Large and Pond
  !
  ! NOTES:
  !   o all fluxes are positive downward
  !   o net heat flux = net sw + lw up + lw down + sen + lat
  !   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
  !   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
  !
  ! ASSUMPTIONS:
  !  Large:
  !   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
  !   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
  !                                 ctn = .0180 sqrt(cdn), stable
  !   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
  !   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
  !-------------------------------------------------------------------------------

  use shr_kind_mod,   only: R8=>SHR_KIND_R8, IN=>SHR_KIND_IN ! shared kinds
  use shr_flux_mod,   only: loc_cpdair, loc_cpvir, loc_karman, loc_g, loc_zvir
  use shr_flux_mod,   only: loc_latvap, loc_stebol, use_coldair_outbreak_mod
  use shr_flux_mod,   only: flux_con_tol, flux_con_max_iter
  use shr_flux_mod,   only: alpha, maxscl, td0
  use shr_sys_mod,    only: shr_sys_abort

  implicit none
  public

  integer, private :: debug = 0

contains

  subroutine flux_atmOcn_large(              &
       logunit, spval, nMax,                 &
       zbot, ubot, vbot, thbot,              &
       qbot, rainc, rbot,                    &
       tbot, us, vs,  pslv,                  &
       ts, mask,  seq_flux_atmocn_minwind,   &
       sen, lat, lwup, evap,                 &
       taux, tauy, tref, qref,               &
       add_gusts, duu10n, ugust_out, u10res, &
       ustar_sv, re_sv, ssq_sv)

    !--- input arguments --------------------------------
    integer  ,intent(in) :: logunit
    real(R8) ,intent(in) :: spval       ! local missing value
    integer  ,intent(in) :: nMax        ! data vector length
    integer  ,intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
    logical  ,intent(in) :: add_gusts
    real(R8) ,intent(in) :: zbot (nMax) ! atm level height (m)
    real(R8) ,intent(in) :: ubot (nMax) ! atm u wind (m/s)
    real(R8) ,intent(in) :: vbot (nMax) ! atm v wind (m/s)
    real(R8) ,intent(in) :: thbot(nMax) ! atm potential T (K)
    real(R8) ,intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
    real(R8) ,intent(in) :: rainc(nMax) ! atm precip for convective gustiness (kg/m^3) - RBN 24Nov2008/MDF 31Jan2022
    real(R8) ,intent(in) :: rbot (nMax) ! atm air density (kg/m^3)
    real(R8) ,intent(in) :: tbot (nMax) ! atm T (K)
    real(R8) ,intent(in) :: pslv (nMax) ! atm sea level pressure(Pa)
    real(R8) ,intent(in) :: us (nMax)   ! ocn u-velocity (m/s)
    real(R8) ,intent(in) :: vs (nMax)   ! ocn v-velocity (m/s)
    real(R8) ,intent(in) :: ts (nMax)   ! ocn temperature (K)
    real(R8) ,intent(in) :: seq_flux_atmocn_minwind ! minimum wind speed for atmocn (m/s)

    !--- output arguments -------------------------------
    real(R8),intent(out)  :: sen (nMax)             ! heat flux: sensible (W/m^2)
    real(R8),intent(out)  :: lat (nMax)             ! heat flux: latent (W/m^2)
    real(R8),intent(out)  :: lwup (nMax)            ! heat flux: lw upward (W/m^2)
    real(R8),intent(out)  :: evap (nMax)            ! water flux: evap ((kg/s)/m^2)
    real(R8),intent(out)  :: taux (nMax)            ! surface stress, zonal (N)
    real(R8),intent(out)  :: tauy (nMax)            ! surface stress, maridional (N)
    real(R8),intent(out)  :: tref (nMax)            ! diag:  2m ref height T (K)
    real(R8),intent(out)  :: qref (nMax)            ! diag:  2m ref humidity (kg/kg)
    real(R8),intent(out)  :: duu10n(nMax)           ! diag: 10m wind speed squared (m/s)^2
    real(R8),intent(out)  :: ugust_out(nMax)        ! diag: gustiness addition to U10 (m/s)
    real(R8),intent(out)  :: u10res(nMax)           ! diag: gustiness addition to U10 (m/s)

    real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
    real(R8),intent(out),optional :: re_sv (nMax)   ! diag: sqrt of exchange coefficient (water)
    real(R8),intent(out),optional :: ssq_sv (nMax)  ! diag: sea surface humidity (kg/kg)

    !--- local constants --------------------------------
    real(R8),parameter :: zref  = 10.0_R8 ! reference height (m)
    real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)

    !real(R8),parameter :: cexcd  = 0.0346_R8 ! ratio Ch(water)/CD
    !real(R8),parameter :: chxcds = 0.018_R8  ! ratio Ch(heat)/CD for stable case
    !real(R8),parameter :: chxcdu = 0.0327_R8 ! ratio Ch(heat)/CD for unstable case

    !--- local variables --------------------------------
    integer  :: n      ! vector loop index
    integer  :: iter
    real(R8) :: vmag   ! surface wind magnitude (m/s)
    real(R8) :: ssq    ! sea surface humidity (kg/kg)
    real(R8) :: delt   ! potential T difference (K)
    real(R8) :: delq   ! humidity difference (kg/kg)
    real(R8) :: stable ! stability factor
    real(R8) :: rdn    ! sqrt of neutral exchange coeff (momentum)
    real(R8) :: rhn    ! sqrt of neutral exchange coeff (heat)
    real(R8) :: ren    ! sqrt of neutral exchange coeff (water)
    real(R8) :: rd     ! sqrt of exchange coefficient (momentum)
    real(R8) :: rh     ! sqrt of exchange coefficient (heat)
    real(R8) :: re     ! sqrt of exchange coefficient (water)
    real(R8) :: ustar  ! ustar
    real(r8) :: ustar_prev
    real(R8) :: qstar  ! qstar
    real(R8) :: tstar  ! tstar
    real(R8) :: hol    ! H (at zbot) over L
    real(R8) :: xsq    ! ?
    real(R8) :: xqq    ! ?
    real(R8) :: psimh  ! stability function at zbot (momentum)
    real(R8) :: psixh  ! stability function at zbot (heat and water)
    real(R8) :: psix2  ! stability function at ztref reference height
    real(R8) :: alz    ! ln(zbot/zref)
    real(R8) :: al2    ! ln(zref/ztref)
    real(R8) :: u10n   ! 10m neutral wind
    real(R8) :: tau    ! stress at zbot
    real(R8) :: cp     ! specific heat of moist air
    real(R8) :: fac    ! vertical interpolation factor
    real(R8) :: wind0  ! resolved large-scale 10m wind (no gust added)

    !--- local functions --------------------------------
    real(R8) :: qsat   ! function: the saturation humididty of air (kg/m^3)

    ! (formula v*=[c4/U10+c5+c6*U10]*U10 in Large et al. 1994)
    real(R8) :: cdn    ! function: neutral drag coeff at 10m

    ! Large only (stability functions)
    real(R8) :: psimhu ! function: unstable part of psimh
    real(R8) :: psixhu ! function: unstable part of psimx
    real(R8) :: Umps   ! dummy arg ~ wind velocity (m/s)
    real(R8) :: Tk     ! dummy arg ~ temperature (K)
    real(R8) :: xd     ! dummy arg ~ ?

    !--- for cold air outbreak calc --------------------------------
    real(R8) :: tdiff(nMax) ! tbot - ts
    real(R8) :: vscl

    real(R8) :: ugust ! function: gustiness as a function of convective rainfall.
    real(R8) :: gprec ! convective rainfall argument for ugust
    ! -------------------------------------------------------------------------

    qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)

    ! Large and Yeager 2009
    cdn(Umps)  =  0.0027_R8 / min(33.0000_R8,Umps) + 0.000142_R8 + &
         0.0000764_R8 * min(33.0000_R8,Umps) - 3.14807e-13_r8 * min(33.0000_R8,Umps)**6

    ! Capped Large and Pond by wind
    !   cdn(Umps)  =   0.0027_R8 / min(30.0_R8,Umps) + 0.000142_R8 + 0.0000764_R8 * min(30.0_R8,Umps)
    ! Capped Large and Pond by Cd
    !   cdn(Umps) = min(0.0025_R8, (0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps ))
    ! Large and Pond
    !   cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps

    psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
    psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

    ! Convective gustiness appropriate for input precipitation.
    ! Following Regelsperger et al. (2000, J. Clim)
    ! Ug = log(1.0+6.69R-0.476R^2)
    ! Coefficients X by 8640 for mm/s (from cam) -> cm/day (for above forumla)
    ugust(gprec) = log(1._R8+57801.6_r8*gprec-3.55332096e7_r8*(gprec**2))

    !--- formats ----------------------------------------
    character(*),parameter :: subName = '(flux_atmOcn) '
    character(*),parameter ::   F00 = "('(flux_atmOcn) ',4a)"
    ! --------------------------------------------------------------------------

    if (debug > 0) write(logunit,F00) "enter"

    u10n  = spval
    rh    = spval
    psixh = spval
    hol   = spval

    !--- for cold air outbreak calc --------------------------------
    tdiff= tbot - ts

    al2 = log(zref/ztref)

    DO n=1,nMax
       if (mask(n) /= 0) then

          !--- compute some needed quantities ---
          if (add_gusts) then
             vmag = max(seq_flux_atmocn_minwind, &
                        sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2 + (1.0_R8*ugust(min(rainc(n),6.94444e-4_r8))**2)) )
             ugust_out(n) = ugust(min(rainc(n),6.94444e-4_r8))
          else
             vmag = max(seq_flux_atmocn_minwind, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
             ugust_out(n) = 0.0_r8
          end if
          wind0 = max(seq_flux_atmocn_minwind, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )

          if (use_coldair_outbreak_mod) then
             ! Cold Air Outbreak Modification:
             ! Increase windspeed for negative tbot-ts
             ! based on Mahrt & Sun 1995,MWR

             if (tdiff(n).lt.td0) then
                ! if add_gusts wind0 and vmag are different, both need this factor.
                vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
                vmag=vmag*vscl
                vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(wind0))),maxscl)
                wind0=wind0*vscl
             endif
          endif

          ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)
          delt   = thbot(n) - ts(n)                  ! pot temp diff (K)
          delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
          alz    = log(zbot(n)/zref)
          cp     = loc_cpdair*(1.0_R8 + loc_cpvir*ssq)

          !------------------------------------------------------------
          ! first estimate of Z/L and ustar, tstar and qstar
          !------------------------------------------------------------

          !--- neutral coefficients, z/L = 0.0 ---
          stable = 0.5_R8 + sign(0.5_R8 , delt)
          rdn    = sqrt(cdn(vmag))
          rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8
          !(1.0_R8-stable) * chxcdu + stable * chxcds
          ren    = 0.0346_R8 !cexcd

          !--- ustar, tstar, qstar ---
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq
          ustar_prev = ustar*2.0_R8
          iter = 0
          do while( abs((ustar - ustar_prev)/ustar) > flux_con_tol .and. iter < flux_con_max_iter)
             iter = iter + 1
             ustar_prev = ustar
             !--- compute stability & evaluate all stability functions ---
             hol  = loc_karman*loc_g*zbot(n)*  &
                  (tstar/thbot(n)+qstar/(1.0_R8/loc_zvir+qbot(n)))/ustar**2
             hol  = sign( min(abs(hol),10.0_R8), hol )
             stable = 0.5_R8 + sign(0.5_R8 , hol)
             xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
             xqq    = sqrt(xsq)
             psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
             psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

             !--- shift wind speed using old coefficient ---
             rd   = rdn / (1.0_R8 + max(rdn/loc_karman*(alz-psimh), -0.5_r8))
             u10n = vmag * rd / rdn

             !--- update transfer coeffs at 10m and neutral stability ---
             rdn = sqrt(cdn(u10n))
             ren = 0.0346_R8 !cexcd
             rhn = (1.0_R8-stable)*0.0327_R8 + stable * 0.018_R8
             !(1.0_R8-stable) * chxcdu + stable * chxcds

             !--- shift all coeffs to measurement height and stability ---
             rd = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             rh = rhn / (1.0_R8 + rhn/loc_karman*(alz-psixh))
             re = ren / (1.0_R8 + ren/loc_karman*(alz-psixh))

             !--- update ustar, tstar, qstar using updated, shifted coeffs --
             ustar = rd * vmag
             tstar = rh * delt
             qstar = re * delq
          enddo
          if (iter < 1) then
             write(logunit,*) ustar,ustar_prev,flux_con_tol,flux_con_max_iter
             call shr_sys_abort('No iterations performed in flux_atmocn_mod')
          end if
          !------------------------------------------------------------
          ! compute the fluxes
          !------------------------------------------------------------

          tau = rbot(n) * ustar * ustar

          !--- momentum flux ---
          taux(n) = tau * (ubot(n)-us(n)) / vmag
          tauy(n) = tau * (vbot(n)-vs(n)) / vmag

          !--- heat flux ---
          sen (n) =          cp * tau * tstar / ustar
          lat (n) =  loc_latvap * tau * qstar / ustar
          lwup(n) = -loc_stebol * ts(n)**4

          !--- water flux ---
          evap(n) = lat(n)/loc_latvap

          !------------------------------------------------------------
          ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------
          hol = hol*ztref/zbot(n)
          xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
          xqq = sqrt(xsq)
          psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
          fac     = (rh/loc_karman) * (alz + al2 - psixh + psix2 )
          tref(n) = thbot(n) - delt*fac
          tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
          fac     = (re/loc_karman) * (alz + al2 - psixh + psix2 )
          qref(n) =  qbot(n) - delq*fac

          duu10n(n) = u10n*u10n ! 10m wind speed squared
          u10res(n) = u10n * (wind0/vmag)  ! resolved 10m wind

          !------------------------------------------------------------
          ! optional diagnostics, needed for water tracer fluxes (dcn)
          !------------------------------------------------------------
          if (present(ustar_sv)) ustar_sv(n) = ustar
          if (present(re_sv   )) re_sv(n)    = re
          if (present(ssq_sv  )) ssq_sv(n)   = ssq

       else
          !------------------------------------------------------------
          ! no valid data here -- out of domain
          !------------------------------------------------------------

          sen (n)      = spval ! sensible         heat flux (W/m^2)
          lat (n)      = spval ! latent           heat flux (W/m^2)
          lwup (n)     = spval ! long-wave upward heat flux (W/m^2)
          evap (n)     = spval ! evaporative water flux ((kg/s)/m^2)
          taux (n)     = spval ! x surface stress (N)
          tauy (n)     = spval ! y surface stress (N)
          tref (n)     = spval ! 2m reference height temperature (K)
          qref (n)     = spval ! 2m reference height humidity (kg/kg)
          duu10n(n)    = spval ! 10m wind speed squared (m/s)^2
          ugust_out(n) = spval ! gustiness addition (m/s)
          u10res(n)    = spval ! 10m resolved wind (no gusts) (m/s)

          if (present(ustar_sv)) ustar_sv(n) = spval
          if (present(re_sv   )) re_sv (n) = spval
          if (present(ssq_sv  )) ssq_sv (n) = spval
       endif
    enddo

  end subroutine flux_atmOcn_large

end module flux_atmOcn_large_mod
