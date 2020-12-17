module shr_flux_mod

  use shr_kind_mod    ! shared kinds
  use shr_const_mod   ! shared constants
  use shr_sys_mod     ! shared system routines
  use shr_log_mod, only: s_logunit => shr_log_Unit

  implicit none

  private ! default private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_flux_atmOcn           ! computes atm/ocn fluxes
  public :: shr_flux_adjust_constants ! adjust constant values used in flux calculations.

  !--- rename kinds for local readability only ---
  integer,parameter :: R8 = SHR_KIND_R8  ! 8 byte real
  integer,parameter :: IN = SHR_KIND_IN  ! native/default integer

  ! The follow variables are not declared as parameters so that they can be
  ! adjusted to support aquaplanet and potentially other simple model modes.
  ! The shr_flux_adjust_constants subroutine is called to set the desired
  ! values.  The default values are from shr_const_mod.  Currently they are
  ! only used by the shr_flux_atmocn and shr_flux_atmice routines.
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
  real(r8)    :: flux_con_tol = 0.0_R8
  integer(IN) :: flux_con_max_iter = 2

  character(len=*), parameter :: sourcefile = &
       __FILE__

  !--- cold air outbreak parameters  (Mahrt & Sun 1995,MWR) -------------
  logical :: use_coldair_outbreak_mod = .false.
  real(R8),parameter    :: alpha = 1.4_R8
  real(R8),parameter    :: maxscl =2._R8  ! maximum wind scaling for flux
  real(R8),parameter    :: td0 = -10._R8   ! start t-ts for scaling

!===============================================================================
contains
!===============================================================================

  subroutine shr_flux_adjust_constants( &
       zvir, cpair, cpvir, karman, gravit, &
       latvap, latice, stebol, flux_convergence_tolerance, &
       flux_convergence_max_iteration, &
       coldair_outbreak_mod)

    ! Adjust local constants.  Used to support simple models.

    real(R8), optional, intent(in) :: zvir
    real(R8), optional, intent(in) :: cpair
    real(R8), optional, intent(in) :: cpvir
    real(R8), optional, intent(in) :: karman
    real(R8), optional, intent(in) :: gravit
    real(R8), optional, intent(in) :: latvap
    real(R8), optional, intent(in) :: latice
    real(R8), optional, intent(in) :: stebol
    real(r8), optional, intent(in)  :: flux_convergence_tolerance
    integer(in), optional, intent(in) :: flux_convergence_max_iteration
    logical, optional, intent(in) :: coldair_outbreak_mod
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
    if(present(coldair_outbreak_mod)) use_coldair_outbreak_mod = coldair_outbreak_mod
  end subroutine shr_flux_adjust_constants

  !===============================================================================
  subroutine shr_flux_atmOcn(nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   &
       &               qbot  ,s16O  ,sHDO  ,s18O  ,rbot  ,   &
       &               tbot  ,us    ,vs    ,   &
       &               ts    ,mask  ,seq_flux_atmocn_minwind, &
       &               sen   ,lat   ,lwup  ,   &
       &               r16O, rhdo, r18O, &
       &               evap  ,evap_16O, evap_HDO, evap_18O, &
       &               taux  ,tauy  ,tref  ,qref  ,   &
       &               ocn_surface_flux_scheme, &
       &               duu10n,  ustar_sv   ,re_sv ,ssq_sv,   &
       &               missval    )

    implicit none

    !--- input arguments --------------------------------
    integer(IN),intent(in) ::       nMax  ! data vector length
    integer(IN),intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
    real(R8)   ,intent(in) :: zbot (nMax) ! atm level height                     (m)
    real(R8)   ,intent(in) :: ubot (nMax) ! atm u wind (bottom or 10m)           (m/s)
    real(R8)   ,intent(in) :: vbot (nMax) ! atm v wind (bottom or 10m)           (m/s)
    real(R8)   ,intent(in) :: thbot(nMax) ! atm potential T                      (K)
    real(R8)   ,intent(in) :: qbot (nMax) ! atm specific humidity (bottom or 2m) (kg/kg)
    real(R8)   ,intent(in) :: s16O (nMax) ! atm H216O tracer conc.               (kg/kg)
    real(R8)   ,intent(in) :: sHDO (nMax) ! atm HDO tracer conc.                 (kg/kg)
    real(R8)   ,intent(in) :: s18O (nMax) ! atm H218O tracer conc.               (kg/kg)
    real(R8)   ,intent(in) :: r16O (nMax) ! ocn H216O tracer ratio/Rstd
    real(R8)   ,intent(in) :: rHDO (nMax) ! ocn HDO tracer ratio/Rstd
    real(R8)   ,intent(in) :: r18O (nMax) ! ocn H218O tracer ratio/Rstd
    real(R8)   ,intent(in) :: rbot (nMax) ! atm air density                      (kg/m^3)
    real(R8)   ,intent(in) :: tbot (nMax) ! atm T (bottom or 2m)                 (K)
    real(R8)   ,intent(in) :: us   (nMax) ! ocn u-velocity                       (m/s)
    real(R8)   ,intent(in) :: vs   (nMax) ! ocn v-velocity                       (m/s)
    real(R8)   ,intent(in) :: ts   (nMax) ! ocn temperature                      (K)
    integer(IN),intent(in), optional :: ocn_surface_flux_scheme
    real(R8)   ,intent(in), optional :: seq_flux_atmocn_minwind ! minimum wind speed for atmocn (m/s)

    !--- output arguments -------------------------------
    real(R8),intent(out)  ::  sen  (nMax) ! heat flux: sensible    (W/m^2)
    real(R8),intent(out)  ::  lat  (nMax) ! heat flux: latent      (W/m^2)
    real(R8),intent(out)  ::  lwup (nMax) ! heat flux: lw upward   (W/m^2)
    real(R8),intent(out)  ::  evap (nMax) ! water flux: evap  ((kg/s)/m^2)
    real(R8),intent(out)  ::  evap_16O (nMax) ! water flux: evap ((kg/s/m^2)
    real(R8),intent(out)  ::  evap_HDO (nMax) ! water flux: evap ((kg/s)/m^2)
    real(R8),intent(out)  ::  evap_18O (nMax) ! water flux: evap ((kg/s/m^2)
    real(R8),intent(out)  ::  taux (nMax) ! surface stress, zonal      (N)
    real(R8),intent(out)  ::  tauy (nMax) ! surface stress, maridional (N)
    real(R8),intent(out)  ::  tref (nMax) ! diag:  2m ref height T     (K)
    real(R8),intent(out)  ::  qref (nMax) ! diag:  2m ref humidity (kg/kg)
    real(R8),intent(out)  :: duu10n(nMax) ! diag: 10m wind speed squared (m/s)^2

    real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
    real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
    real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)

    real(R8),intent(in) ,optional :: missval        ! masked value

    ! !EOP

    !--- local constants --------------------------------
    real(R8),parameter :: umin  =  0.5_R8 ! minimum wind speed       (m/s)
    real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
    real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)
    !!++ Large only
    !real(R8),parameter :: cexcd  = 0.0346_R8 ! ratio Ch(water)/CD
    !real(R8),parameter :: chxcds = 0.018_R8  ! ratio Ch(heat)/CD for stable case
    !real(R8),parameter :: chxcdu = 0.0327_R8 ! ratio Ch(heat)/CD for unstable case
    !!++ COARE only
    real(R8),parameter :: zpbl =700.0_R8 ! PBL depth [m] for gustiness parametriz.

    !--- local variables --------------------------------
    integer(IN) :: n      ! vector loop index
    integer(IN) :: iter
    real(R8)    :: vmag   ! surface wind magnitude   (m/s)
    real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
    real(R8)    :: delt   ! potential T difference   (K)
    real(R8)    :: delq   ! humidity difference      (kg/kg)
    real(R8)    :: stable ! stability factor
    real(R8)    :: rdn    ! sqrt of neutral exchange coeff (momentum)
    real(R8)    :: rhn    ! sqrt of neutral exchange coeff (heat)
    real(R8)    :: ren    ! sqrt of neutral exchange coeff (water)
    real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)
    real(R8)    :: rh     ! sqrt of exchange coefficient (heat)
    real(R8)    :: re     ! sqrt of exchange coefficient (water)
    real(R8)    :: ustar  ! ustar
    real(r8)     :: ustar_prev
    real(R8)    :: qstar  ! qstar
    real(R8)    :: tstar  ! tstar
    real(R8)    :: hol    ! H (at zbot) over L
    real(R8)    :: xsq    ! ?
    real(R8)    :: xqq    ! ?
    !!++ Large only
    real(R8)    :: psimh  ! stability function at zbot (momentum)
    real(R8)    :: psixh  ! stability function at zbot (heat and water)
    real(R8)    :: psix2  ! stability function at ztref reference height
    real(R8)    :: alz    ! ln(zbot/zref)
    real(R8)    :: al2    ! ln(zref/ztref)
    real(R8)    :: u10n   ! 10m neutral wind
    real(R8)    :: tau    ! stress at zbot
    real(R8)    :: cp     ! specific heat of moist air
    real(R8)    :: fac    ! vertical interpolation factor
    real(R8)    :: spval  ! local missing value
    !!++ COARE only
    real(R8)    :: zo,zot,zoq      ! roughness lengths
    real(R8)    :: hsb,hlb         ! sens & lat heat flxs at zbot
    real(R8) :: trf,qrf,urf,vrf ! reference-height quantities

    !--- local functions --------------------------------
    real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
    !!++ Large only (formula v*=[c4/U10+c5+c6*U10]*U10 in Large et al. 1994)
    real(R8)    :: cdn    ! function: neutral drag coeff at 10m
    !!++ Large only (stability functions)
    real(R8)    :: psimhu ! function: unstable part of psimh
    real(R8)    :: psixhu ! function: unstable part of psimx
    real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
    real(R8)    :: Tk     ! dummy arg ~ temperature (K)
    real(R8)    :: xd     ! dummy arg ~ ?
    !--- for cold air outbreak calc --------------------------------
    real(R8)    :: tdiff(nMax)               ! tbot - ts
    real(R8)    :: vscl

    qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)
    cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps
    psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
    psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

    !--- formats ----------------------------------------
    character(*),parameter :: subName = '(shr_flux_atmOcn) '
    character(*),parameter ::   F00 = "('(shr_flux_atmOcn) ',4a)"

    !-------------------------------------------------------------------------------
    ! PURPOSE:
    !   computes atm/ocn surface fluxes
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
    !  COARE:
    !   o use COAREv3.0 function (tht 22/11/2013)
    !-------------------------------------------------------------------------------

    if (present(missval)) then
       spval = missval
    else
       spval = shr_const_spval
    endif
    u10n = spval
    rh = spval
    psixh = spval
    hol=spval

    !--- for cold air outbreak calc --------------------------------
    tdiff= tbot - ts

    al2 = log(zref/ztref)
    DO n=1,nMax
       if (mask(n) /= 0) then

          !--- compute some needed quantities ---
          vmag   = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
          if (use_coldair_outbreak_mod) then
             ! Cold Air Outbreak Modification:
             ! Increase windspeed for negative tbot-ts
             ! based on Mahrt & Sun 1995,MWR

             if (tdiff(n).lt.td0) then
                vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
                vmag=vmag*vscl
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
             rd   = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             if (ocn_surface_flux_scheme == -1)then
                u10n = vmag
             else
                u10n = vmag * rd / rdn
             end if

             !--- update transfer coeffs at 10m and neutral stability ---
             rdn = sqrt(cdn(u10n))
             ren = 0.0346_R8 !cexcd
             rhn = (1.0_R8-stable)*0.0327_R8 + stable * 0.018_R8
             !(1.0_R8-stable) * chxcdu + stable * chxcds

             !--- shift all coeffs to measurement height and stability ---
             if (ocn_surface_flux_scheme == -1)then
               rd = rdn
             else
               rd = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             end if
             rh = rhn / (1.0_R8 + rhn/loc_karman*(alz-psixh))
             re = ren / (1.0_R8 + ren/loc_karman*(alz-psixh))

             !--- update ustar, tstar, qstar using updated, shifted coeffs --
             ustar = rd * vmag
             tstar = rh * delt
             qstar = re * delq
          enddo
          if (iter < 1) then
             write(s_logunit,*) ustar,ustar_prev,flux_con_tol,flux_con_max_iter
             call shr_sys_abort('shr_flux_mod: No iterations performed ')
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

          if (present(ustar_sv)) ustar_sv(n) = spval
          if (present(re_sv   )) re_sv   (n) = spval
          if (present(ssq_sv  )) ssq_sv  (n) = spval
       endif
    end DO

  end subroutine shr_flux_atmOcn

end module shr_flux_mod
