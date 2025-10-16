module flux_atmocn_diurnal_mod

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
   !   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
   !   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
   !                                 ctn = .0180 sqrt(cdn), stable
   !   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
   !   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
   !-------------------------------------------------------------------------------

   use shr_kind_mod,          only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN ! shared kinds
   use shr_flux_mod,          only : td0, maxscl, alpha, use_coldair_outbreak_mod
   use shr_const_mod,         only : shr_const_zvir, shr_const_cpdair, shr_const_karman, shr_const_g
   use shr_const_mod,         only : shr_const_latvap, shr_const_latice, shr_const_stebol, shr_const_tkfrz
   use shr_const_mod,         only : shr_const_pi, shr_const_spval, shr_const_cpvir
   use shr_const_mod,         only : shr_const_ocn_ref_sal, shr_const_zsrflyr, shr_const_rgas
   use shr_sys_mod,           only : shr_sys_abort
   use flux_atmocn_COARE_mod, only : cor30a

   implicit none
   private

   public :: flux_atmOcn_Diurnal

   private :: cuberoot

   integer  :: flux_con_max_iter = 2
   real(r8) :: flux_con_tol = 0.0_R8
   integer  :: debug = 0

contains

  subroutine flux_atmOcn_diurnal(                      &
       logunit, spval, ocn_surface_flux_scheme,        &
       nMax, zbot, ubot, vbot, thbot,                  &
       qbot,  rbot, tbot, us, vs,                      &
       ts, mask,  seq_flux_atmocn_minwind,             &
       sen, lat, lwup,                                 &
       evap,  taux, tauy, tref, qref,                  &
       uGust, lwdn,  swdn,  swup, prec,                &
       swpen, ocnsal, ocn_prognostic,                  &
       latt, long,  warm,  salt,  speed, regime,       &
       warmMax, windMax, qSolAvg, windAvg,             &
       warmMaxInc, windMaxInc, qSolInc, windInc, nInc, &
       tBulk, tSkin, tSkin_day, tSkin_night,           &
       cSkin, cSkin_night, secs, dt,                   &
       duu10n,  ustar_sv, re_sv, ssq_sv, cold_start)

    ! Arguments
    !
    integer  ,intent(in)    :: logunit
    real(r8) ,intent(in)    :: spval
    integer  ,intent(in)    :: ocn_surface_flux_scheme
    integer  ,intent(in)    :: nMax                    ! data vector length
    integer  ,intent(in)    :: mask (nMax)             ! ocn domain mask 0 <=> out of domain
    real(R8) ,intent(in)    :: zbot (nMax)             ! atm level height(m)
    real(R8) ,intent(in)    :: ubot (nMax)             ! atm u wind(m/s)
    real(R8) ,intent(in)    :: vbot (nMax)             ! atm v wind(m/s)
    real(R8) ,intent(in)    :: thbot(nMax)             ! atm potential T (K)
    real(R8) ,intent(in)    :: qbot (nMax)             ! atm specific humidity(kg/kg)
    real(R8) ,intent(in)    :: rbot (nMax)             ! atm air density(kg/m^3)
    real(R8) ,intent(in)    :: tbot (nMax)             ! atm T(K)
    real(R8) ,intent(in)    :: us   (nMax)             ! ocn u-velocity (m/s)
    real(R8) ,intent(in)    :: vs   (nMax)             ! ocn v-velocity (m/s)
    real(R8) ,intent(in)    :: ts   (nMax)             ! ocn temperature(K)
    real(R8) ,intent(in)    :: uGust (nMax)            ! NEW not used
    real(R8) ,intent(in)    :: lwdn  (nMax)            ! NEW
    real(R8) ,intent(in)    :: swdn  (nMax)            ! NEW
    real(R8) ,intent(in)    :: swup  (nMax)            ! NEW
    real(R8) ,intent(in)    :: prec  (nMax)            ! NEW
    real(R8) ,intent(in)    :: latt  (nMax)            ! NEW
    real(R8) ,intent(in)    :: long  (nMax)            ! NEW
    logical  ,intent(in)    :: ocn_prognostic          ! NEW
    integer  ,intent(in)    :: secs                    ! NEW  elsapsed seconds in day (GMT)
    integer  ,intent(in)    :: dt                      ! NEW
    real(R8) ,intent(inout) :: swpen (nMax)            ! NEW
    real(R8) ,intent(inout) :: ocnsal(nMax)            ! NEW (kg/kg)
    real(R8) ,intent(inout) :: warm  (nMax)            ! NEW
    real(R8) ,intent(inout) :: salt  (nMax)            ! NEW
    real(R8) ,intent(inout) :: speed (nMax)            ! NEW
    real(R8) ,intent(inout) :: regime(nMax)            ! NEW
    real(R8) ,intent(out)   :: warmMax(nMax)           ! NEW
    real(R8) ,intent(out)   :: windMax(nMax)           ! NEW
    real(R8) ,intent(inout) :: qSolAvg(nMax)           ! NEW
    real(R8) ,intent(inout) :: windAvg(nMax)           ! NEW
    real(R8) ,intent(inout) :: warmMaxInc(nMax)        ! NEW
    real(R8) ,intent(inout) :: windMaxInc(nMax)        ! NEW
    real(R8) ,intent(inout) :: qSolInc(nMax)           ! NEW
    real(R8) ,intent(inout) :: windInc(nMax)           ! NEW
    real(R8) ,intent(inout) :: nInc(nMax)              ! NEW
    real(R8) ,intent(out)   :: tBulk (nMax)            ! NEW
    real(R8) ,intent(out)   :: tSkin (nMax)            ! NEW
    real(R8) ,intent(out)   :: tSkin_day (nMax)        ! NEW
    real(R8) ,intent(out)   :: tSkin_night (nMax)      ! NEW
    real(R8) ,intent(out)   :: cSkin (nMax)            ! NEW
    real(R8) ,intent(out)   :: cSkin_night (nMax)      ! NEW
    logical  ,intent(in)    :: cold_start              ! cold start flag
    real(R8) ,intent(in)    :: seq_flux_atmocn_minwind ! minimum wind speed for atmocn      (m/s)
    real(R8) ,intent(out)   :: sen  (nMax)             ! heat flux: sensible    (W/m^2)
    real(R8) ,intent(out)   :: lat  (nMax)             ! heat flux: latent      (W/m^2)
    real(R8) ,intent(out)   :: lwup (nMax)             ! heat flux: lw upward   (W/m^2)
    real(R8) ,intent(out)   :: evap (nMax)             ! water flux: evap  ((kg/s)/m^2)
    real(R8) ,intent(out)   :: taux (nMax)             ! surface stress, zonal      (N)
    real(R8) ,intent(out)   :: tauy (nMax)             ! surface stress, maridional (N)
    real(R8) ,intent(out)   :: tref (nMax)             ! diag:  2m ref height T     (K)
    real(R8) ,intent(out)   :: qref (nMax)             ! diag:  2m ref humidity (kg/kg)
    real(R8) ,intent(out)   :: duu10n(nMax)            ! diag: 10m wind speed squared (m/s)^2
    real(R8) ,intent(out),optional :: ustar_sv(nMax)   ! diag: ustar
    real(R8) ,intent(out),optional :: re_sv   (nMax)   ! diag: sqrt of exchange coefficient (water)
    real(R8) ,intent(out),optional :: ssq_sv  (nMax)   ! diag: sea surface humidity  (kg/kg)

    !--- local constants --------------------------------
    real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
    real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)

    real(R8),parameter :: lambdaC  = 6.0_R8
    real(R8),parameter :: lambdaL  = 0.0_R8
    real(R8),parameter :: doLMax   = 1.0_R8
    real(R8),parameter :: pwr      = 0.2_R8
    real(R8),parameter :: Rizero   = 1.0_R8
    real(R8),parameter :: NUzero   = 40.0e-4_R8
    real(R8),parameter :: Prandtl  = 1.0_R8
    real(R8),parameter :: kappa0   = 0.2e-4_R8

    real(R8),parameter :: F0       = 0.5_R8
    real(R8),parameter :: F1       = 0.15_R8
    real(R8),parameter :: R1       = 10.0_R8

    real(R8),parameter :: Ricr     = 0.30_R8
    real(R8),parameter :: tiny     = 1.0e-12_R8
    real(R8),parameter :: tiny2    = 1.0e-6_R8
    real(R8),parameter :: pi       = SHR_CONST_PI

    !!++ COARE only
    real(R8),parameter :: zpbl =700.0_R8 ! PBL depth [m] for gustiness parametriz.

    !--- local variables --------------------------------
    integer(IN) :: n       ! vector loop index
    integer(IN) :: iter       ! iteration loop index
    integer(IN) :: lsecs   ! local seconds elapsed
    integer(IN) :: lonsecs ! incrememnt due to lon offset
    real(R8)    :: vmag    ! surface wind magnitude   (m/s)
    real(R8)    :: ssq     ! sea surface humidity     (kg/kg)
    real(R8)    :: delt    ! potential T difference   (K)
    real(R8)    :: delq    ! humidity difference      (kg/kg)
    real(R8)    :: stable  ! stability factor
    real(R8)    :: rdn     ! sqrt of neutral exchange coeff (momentum)
    real(R8)    :: rhn     ! sqrt of neutral exchange coeff (heat)
    real(R8)    :: ren     ! sqrt of neutral exchange coeff (water)
    real(R8)    :: rd      ! sqrt of exchange coefficient (momentum)
    real(R8)    :: rh      ! sqrt of exchange coefficient (heat)
    real(R8)    :: re      ! sqrt of exchange coefficient (water)
    real(R8)    :: ustar   ! ustar
    real(R8)    :: ustar_prev   ! ustar
    real(R8)    :: qstar   ! qstar
    real(R8)    :: tstar   ! tstar
    real(R8)    :: hol     ! H (at zbot) over L
    real(R8)    :: xsq     ! ?
    real(R8)    :: xqq     ! ?
    real(R8)    :: psimh   ! stability function at zbot (momentum)
    real(R8)    :: psixh   ! stability function at zbot (heat and water)
    real(R8)    :: psix2   ! stability function at ztref reference height
    real(R8)    :: alz     ! ln(zbot/zref)
    real(R8)    :: al2     ! ln(zref/ztref)
    real(R8)    :: u10n    ! 10m neutral wind
    real(R8)    :: tau     ! stress at zbot
    real(R8)    :: cp      ! specific heat of moist air
    real(R8)    :: fac     ! vertical interpolation factor
    real(R8)    :: DTiter  !
    real(R8)    :: DSiter  !
    real(R8)    :: DViter  !

    real(R8)    :: Dcool   !
    real(R8)    :: Qdel    ! net cool skin heating
    real(R8)    :: Hd      ! net heating above -z=d
    real(R8)    :: Hb      ! net kinematic heating above -z = delta
    real(R8)    :: lambdaV !
    real(R8)    :: Fd      ! net fresh water forcing above -z=d
    real(R8)    :: ustarw  ! surface wind forcing of layer above -z=d

    real(R8)    :: Qsol   ! solar heat flux (W/m2)
    real(R8)    :: Qnsol  ! non-solar heat flux (W/m2)

    real(R8)    :: SSS  ! sea surface salinity
    real(R8)    :: alphaT  !
    real(R8)    :: betaS  !

    real(R8)    :: doL     ! ocean forcing stablity parameter
    real(R8)    :: Rid     ! Richardson number at depth d
    real(R8)    :: Ribulk  ! Bulk  Richardson number at depth d
    real(R8)    :: FofRi   ! Richardon number dependent diffusivity
    real(R8)    :: Smult   ! multiplicative term based on regime
    real(R8)    :: Sfact   ! multiplicative term based on regime
    real(R8)    :: Kdiff   ! diffusive term based on regime
    real(R8)    :: Kvisc   ! viscosity term based on regime
    real(R8)    :: rhocn   !
    real(R8)    :: rcpocn  !
    real(R8)    :: Nreset  ! value for multiplicative reset factor
    logical     :: lmidnight
    logical     :: ltwopm
    logical     :: ltwoam
    logical     :: lfullday
    integer     :: nsum
    real(R8)    :: pexp   ! eqn 19
    real(R8)    :: AMP    ! eqn 18
    real(R8)    :: dif3
    real(R8)    :: phid

    !!++ COARE only
    real(R8)    :: zo,zot,zoq      ! roughness lengths
    real(R8)    :: hsb,hlb         ! sens & lat heat flxs at zbot
    real(R8)    :: trf,qrf,urf,vrf ! reference-height quantities

    !--- local functions --------------------------------
    real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
    real(R8)    :: cdn    ! function: neutral drag coeff at 10m
    real(R8)    :: psimhu ! function: unstable part of psimh
    real(R8)    :: psixhu ! function: unstable part of psimx
    real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
    real(R8)    :: Tk     ! dummy arg ~ temperature (K)
    real(R8)    :: xd     ! dummy arg ~ ?
    real(R8)    :: molvisc ! molecular viscosity
    real(R8)    :: molPr   ! molecular Prandtl number

    !--- for cold air outbreak calc --------------------------------
    real(R8)    :: tdiff(nMax)               ! tbot - ts
    real(R8)    :: vscl

    qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)
    cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps
    psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
    psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)
    molvisc(Tk)  = 1.623e-6_R8 * exp((-1.0_R8*(Tk-273.15_R8))/45.2_R8)
    molPr(Tk)    = 11.64_R8 * exp((-1.0_R8*(Tk-273.15_R8))/40.7_R8)

    !--- formats ----------------------------------------
    character(*),parameter :: subName = '(flux_atmOcn_diurnal) '
    character(*),parameter ::   F00 = "('(flux_atmOcn_diurnal) ',4a)"

    if (debug > 0) write(logunit,F00) "enter"

    rh = spval
    dviter = spval
    dtiter = spval
    dsiter = spval
    al2 = log(zref/ztref)

    !--- for cold air outbreak calc --------------------------------
    tdiff= tbot - ts

    ! equations 18 and 19
    AMP = 1.0_R8/F0-1.0_R8
    pexp = log( (1.0_R8/F1-F0) / (1.0_R8-F0) ) / log(R1)

    if (.not. ocn_prognostic) then
       ! Set swpen and ocean salinity from following analytic expressions
       swpen(:) = 0.67_R8*(exp((-1._R8*shr_const_zsrflyr)/1.0_R8)) + &
            0.33_R8*exp((-1._R8*shr_const_zsrflyr)/17.0_R8)
       ocnsal(:) = shr_const_ocn_ref_sal/1000.0_R8
    else
       ! use swpen and ocnsal from input argument
    endif

    if (cold_start) then
       write(logunit,F00) "Initialize diurnal cycle fields"
       warm       (:) = 0.0_R8
       salt       (:) = 0.0_R8
       speed      (:) = 0.0_R8
       regime     (:) = 0.0_R8
       qSolAvg    (:) = 0.0_R8
       windAvg    (:) = 0.0_R8
       warmMax    (:) = 0.0_R8
       windMax    (:) = 0.0_R8
       warmMaxInc (:) = 0.0_R8
       windMaxInc (:) = 0.0_R8
       qSolInc    (:) = 0.0_R8
       windInc    (:) = 0.0_R8
       nInc       (:) = 0.0_R8
       tSkin_day  (:) = ts(:)
       tSkin_night(:) = ts(:)
       cSkin_night(:) = 0.0_R8
    endif
    u10n = 0.0_r8
    stable = 0.0_r8
    DO n=1,nMax

       if (mask(n) /= 0) then

          !--- compute some initial and useful flux quantities ---

          vmag = max(seq_flux_atmocn_minwind, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
          if (use_coldair_outbreak_mod) then
             ! Cold Air Outbreak Modification:
             ! Increase windspeed for negative tbot-ts
             ! based on Mahrt & Sun 1995,MWR

             if (tdiff(n).lt.td0) then
                vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
                vmag=vmag*vscl
             endif
          endif
          alz      = log(zbot(n)/zref)
          hol      = 0.0
          psimh    = 0.0
          psixh    = 0.0
          rdn      = sqrt(cdn(vmag))

          tBulk(n) = ts(n)+warm(n)    ! first guess for tBulk from read in ts,warm
          tSkin(n) = tBulk(n)
          Qsol     = swdn(n) + swup(n)
          SSS      = 1000.0_R8*ocnsal(n)+salt(n)
          lambdaV = lambdaC

          alphaT   = 0.000297_R8*(1.0_R8+0.0256_R8*(ts(n)-298.15_R8)+0.003_R8*(SSS - 35.0_R8))
          betaS    = 0.000756_R8*(1.0_R8-0.0016_R8*(ts(n)-298.15_R8))
          rhocn    = 1023.342_R8*(1.0_R8-0.000297_R8*(ts(n)-298.15_R8)+0.000756_R8 * (SSS - 35.0_R8))
          rcpocn   = rhocn * 3990.0_R8*(1.0_R8-0.0012_R8*(SSS - 35.0_R8))

          Rid =  shr_const_g * (alphaT*warm(n) - betaS*salt(n)) *pwr*shr_const_zsrflyr  / &
               ( pwr*MAX(tiny,speed(n)) )**2

          Ribulk = 0.0

          !----------------------------------------------------------
          ! convert elapsed time from GMT to local &
          ! check elapsed time. reset warm if near lsecs = reset_sec
          !----------------------------------------------------------
          Nreset = 1.0_R8

          lonsecs   = ceiling(long(n)/360.0_R8*86400.0)
          lsecs     = mod(secs + lonsecs,86400)

          lmidnight = (lsecs >= 0     .and. lsecs < dt)        ! 0 = midnight
          ltwopm    = (lsecs >= 48600 .and. lsecs < 48600+dt)  ! 48600 = 1:30pm
          ltwoam    = (lsecs >= 5400  .and. lsecs < 5400 +dt)  ! 5400 = 1:30am
          lfullday  = (lsecs > 86400-dt .and. lsecs <= 86400)
          nsum = nint(nInc(n))

          if ( lmidnight ) then
             Regime(n)  = 1.0_R8               !  RESET DIURNAL
             warm(n)    = 0.0_R8
             salt(n)    = 0.0_R8
             speed(n)   = 0.0_R8
          endif

          ssq    = 0.98_R8 * qsat(tBulk(n)) / rbot(n)   ! sea surf hum (kg/kg)
          delt   = thbot(n) - tBulk(n)                  ! pot temp diff (K)
          delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
          cp     = shr_const_cpdair*(1.0_R8 + shr_const_cpvir*ssq)

          !!.................................................................
          !! ocn_surface_flux_scheme = 0 : Default E3SMv1
          !!                         = 1 : COARE algorithm
          !!.................................................................
          if (ocn_surface_flux_scheme .eq. 0) then! use Large algorithm
             stable = 0.5_R8 + sign(0.5_R8 , delt)


             !--- shift wind speed using old coefficient  and stability function

             rd   = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
             u10n = vmag * rd / rdn

             !--- initial neutral  transfer coeffs at 10m
             rdn    = sqrt(cdn(u10n))
             rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8
             ren    = 0.0346_R8

             !--- initial ustar, tstar, qstar ---
             ustar = rdn * vmag
             tstar = rhn * delt
             qstar = ren * delq

          else if (ocn_surface_flux_scheme .eq. 1) then! use COARE algorithm

             call cor30a(ubot(n),vbot(n),tbot(n),qbot(n),rbot(n) &  ! in atm params
                  & ,us(n),vs(n),tBulk(n),ssq                &  ! in surf params (NB ts -> tBulk)
                  & ,zpbl,zbot(n),zbot(n),zref,ztref,ztref   &  ! in heights
                  & ,tau,hsb,hlb                             &  ! out: fluxes
                  & ,zo,zot,zoq,hol,ustar,tstar,qstar        &  ! out: ss scales
                  & ,rd,rh,re                                &  ! out: exch. coeffs
                  & ,trf,qrf,urf,vrf)                             ! out: reference-height params

             ! for the sake of maintaining same defs
             hol=zbot(n)/hol
             rd=sqrt(rd)
             rh=sqrt(rh)
             re=sqrt(re)

          ELSE  ! N.B.: *no* valid ocn_surface_flux_scheme=2 option if diurnal=.true.

             call shr_sys_abort(subName//" flux_atmOcn_diurnal requires ocn_surface_flux_scheme = 0 or 1")

          ENDIF

          ustar_prev = ustar * 2.0_R8
          iter = 0
          ! --- iterate ---
          ! Originally this code did three iterations while the non-diurnal version did two
          ! So in the new loop this is <= flux_con_max_iter instead of < so that the same defaults
          ! will give the same answers in both cases.
          do while( abs((ustar - ustar_prev)/ustar) > flux_con_tol .and. iter <= flux_con_max_iter)
             iter = iter + 1
             ustar_prev = ustar
             !------------------------------------------------------------
             ! iterate to converge on FLUXES  Z/L, ustar, tstar and qstar
             ! and on Rid  in the DIURNAL CYCLE
             !------------------------------------------------------------
             Smult = 0.0_R8
             Sfact = 0.0_R8
             Kdiff = 0.0_R8
             Kvisc = 0.0_R8
             dif3 = 0.0_R8

             ustarw  = ustar*sqrt(max(tiny,rbot(n)/rhocn))
             Qnsol   = lwdn(n) - shr_const_stebol*(tSkin(n))**4 + &
                  rbot(n)*ustar*(cp*tstar + shr_const_latvap*qstar)
             Hd      = (Qnsol   + Qsol*(1.0_R8-swpen(n)) ) / rcpocn
             Fd      = (prec(n) + rbot(n)*ustar*qstar ) * SSS / rhocn

             !--- COOL SKIN EFFECT ---
             Dcool  = lambdaV*molvisc(tBulk(n)) / ustarw
             Qdel   = Qnsol + Qsol * &
                  (0.137_R8 + 11.0_R8*Dcool - 6.6e-5/Dcool *(1.0_R8 - exp((-1.0_R8*Dcool)/8.0e-4)))
             Hb = (Qdel/rcpocn)+(Fd*betaS/alphaT)
             Hb = min(Hb , 0.0_R8)

             !            lambdaV = lambdaC*(1.0_R8 + ( (0.0_R8-Hb)*16.0_R8*molvisc(tBulk(n))* &
             !                 shr_const_g*alphaT*molPr(tBulk(n))**2/ustarw**4)**0.75)**(-1._R8/3._R8)
             lambdaV = 6.5_R8
             cSkin(n) =  MIN(0.0_R8, lambdaV * molPr(tBulk(n)) * Qdel / ustarw / rcpocn )

             !--- REGIME ---
             doL = shr_const_zsrflyr*shr_const_karman*shr_const_g* &
                  (alphaT*Hd + betaS*Fd ) / ustarw**3
             Rid = MAX(0.0_R8,Rid)
             Smult = dt * (pwr+1.0_R8) / (shr_const_zsrflyr*pwr)
             Sfact = dt * (pwr+1.0_R8) / (shr_const_zsrflyr)**2
             FofRi = 1.0_R8/(1.0_R8 + AMP*(Rid/Rizero)**pexp)

             if ( (doL.gt.0.0_R8) .and. (Qsol.gt.0.0)  ) then
                phid  = MIN(1.0_R8 + 5.0_R8 * doL, 5.0_R8 + doL)
                FofRi = 1.0_R8/(1.0_R8 + AMP*(Rid/Rizero)**pexp)
                dif3 = (kappa0 + NUzero *FofRi)

                if ((doL.le.lambdaL).and.(NINT(regime(n)).le.2)) then
                   regime(n) = 2.0_R8
                   Kdiff =  shr_const_karman * ustarw * shr_const_zsrflyr / phid
                   Kvisc = Kdiff * (1.0_R8 - doL/lambdaL)**2 + &
                        dif3 * (doL/lambdaL)**2 * (3.0_R8 - 2.0_R8 * doL/lambdaL)
                   Kdiff = Kvisc
                else
                   regime(n) = 3.0_R8
                   Kdiff =          kappa0 + NUzero * FofRi
                   Kvisc = Prandtl* kappa0 + NUzero * FofRi
                endif
             else
                if (regime(n).eq.1.0_R8) then
                   Smult      = 0.0_R8
                else
                   if (Ribulk .gt. Ricr) then
                      regime(n) = 3.0_R8
                      Kdiff =          kappa0 + NUzero * FofRi
                      Kvisc = Prandtl* kappa0 + NUzero * FofRi
                   else
                      regime(n) = 4.0_R8
                      Kdiff = shr_const_karman*ustarw*shr_const_zsrflyr *cuberoot(1.0_R8-7.0_R8*doL)
                      Kvisc = Kdiff
                   endif
                endif

             endif

             !--- IMPLICIT INTEGRATION ---

             DTiter = (warm(n)  +(Smult*Hd))               /(1.+ Sfact*Kdiff)
             DSiter = (salt(n)  -(Smult*Fd))               /(1.+ Sfact*Kdiff)
             DViter = (speed(n) +(Smult*ustarw*ustarw))    /(1.+ Sfact*Kvisc)
             DTiter = MAX( 0.0_R8, DTiter)
             DViter = MAX( 0.0_R8, DViter)

             Rid =(shr_const_g*(alphaT*DTiter-betaS*DSiter)*pwr*shr_const_zsrflyr)  / &
                  (pwr*MAX(tiny,DViter))**2
             Ribulk = Rid * pwr
             Ribulk = 0.0_R8
             tBulk(n) = ts(n) + DTiter
             tSkin(n) = tBulk(n) + cskin(n)

             !--need to update ssq,delt,delq as function of tBulk ----

             ssq    = 0.98_R8 * qsat(tBulk(n)) / rbot(n)   ! sea surf hum (kg/kg)
             delt   = thbot(n) - tBulk(n)                  ! pot temp diff (K)
             delq   = qbot(n) - ssq                        ! spec hum dif (kg/kg)

             !--- UPDATE FLUX ITERATION ---

             !!.................................................................
             !! ocn_surface_flux_scheme = 0 : Default CESM1.2
             !!                         = 1 : COARE algorithm
             !!.................................................................
             if (ocn_surface_flux_scheme .eq. 0) then! use Large algorithm

                !--- compute stability & evaluate all stability functions ---
                hol  = shr_const_karman*shr_const_g*zbot(n)*  &
                     (tstar/thbot(n)+qstar/(1.0_R8/shr_const_zvir+qbot(n)))/ustar**2
                hol  = sign( min(abs(hol),10.0_R8), hol )
                stable = 0.5_R8 + sign(0.5_R8 , hol)
                xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
                xqq    = sqrt(xsq)
                psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
                psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

                !--- shift wind speed using old coefficient  and stability function  ---
                rd   = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
                u10n = vmag * rd / rdn

                !--- update neutral  transfer coeffs at 10m
                rdn    = sqrt(cdn(u10n))
                rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8
                ren    = 0.0346_R8

                !--- shift all coeffs to measurement height and stability ---
                rd = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
                rh = rhn / (1.0_R8 + rhn/shr_const_karman*(alz-psixh))
                re = ren / (1.0_R8 + ren/shr_const_karman*(alz-psixh))

                ustar = rd * vmag
                tstar = rh * delt
                qstar = re * delq

                !--- heat flux ---

                tau     = rbot(n) * ustar * ustar
                sen (n) =                cp * tau * tstar / ustar
                lat (n) = shr_const_latvap * tau * qstar / ustar

             else if (ocn_surface_flux_scheme .eq. 1) then! use COARE algorithm

                call cor30a(ubot(n),vbot(n),tbot(n),qbot(n),rbot(n) &  ! in atm params
                     & ,us(n),vs(n),tBulk(n),ssq                &  ! in surf params (NB ts -> tBulk)
                     & ,zpbl,zbot(n),zbot(n),zref,ztref,ztref   &  ! in heights
                     & ,tau,hsb,hlb                             &  ! out: fluxes
                     & ,zo,zot,zoq,hol,ustar,tstar,qstar        &  ! out: ss scales
                     & ,rd,rh,re                                &  ! out: exch. coeffs
                     & ,trf,qrf,urf,vrf)                               ! out: reference-height params

                ! for the sake of maintaining same defs
                hol=zbot(n)/hol
                rd=sqrt(rd)
                rh=sqrt(rh)
                re=sqrt(re)

                !--- heat flux ---
                sen (n) =  hsb
                lat (n) =  hlb

             else ! N.B.: NO ocn_surface_flux_scheme=2 option
                call shr_sys_abort(subName//", flux_diurnal requires ocn_surface_flux_scheme = 0 or 1")
             endif

          enddo   ! end iteration loop
          if (iter < 1) then
             call shr_sys_abort('No iterations performed ')
          end if

          !--- COMPUTE FLUXES TO ATMOSPHERE AND OCEAN ---

          !--- momentum flux ---
          taux(n) = tau * (ubot(n)-us(n)) / vmag
          tauy(n) = tau * (vbot(n)-vs(n)) / vmag

          !--- LW radiation ---
          lwup(n) = -shr_const_stebol * Tskin(n)**4

          !--- water flux ---
          evap(n) = lat(n)/shr_const_latvap

          !------------------------------------------------------------
          ! compute diagnostics: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------

          if (ocn_surface_flux_scheme .eq. 0) then ! use Large algorithm

             hol = hol*ztref/zbot(n)
             xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
             xqq = sqrt(xsq)
             psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
             fac     = (rh/shr_const_karman) * (alz + al2 - psixh + psix2 )
             tref(n) = thbot(n) - delt*fac
             tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
             fac     = (re/shr_const_karman) * (alz + al2 - psixh + psix2 )
             qref(n) =  qbot(n) - delq*fac

             duu10n(n) = u10n*u10n ! 10m wind speed squared

          else if (ocn_surface_flux_scheme .eq. 1) then! use COARE algorithm

             tref(n) = trf
             qref(n) = qrf
             duu10n(n) = urf**2+vrf**2
             u10n = sqrt(duu10n(n))
          endif

          !------------------------------------------------------------
          ! update new prognostic variables
          !------------------------------------------------------------

          warm  (n) = DTiter
          salt  (n) = DSiter
          speed (n) = DViter

          if (ltwopm) then
             tSkin_day(n) = tSkin(n)
             warmmax(n) = max(DTiter,0.0_R8)
          endif

          if (ltwoam) then
             tSkin_night(n) = tSkin(n)
             cSkin_night(n) = cSkin(n)
          endif

          if ((lmidnight).and.(lfullday)) then
             qSolAvg(n) = qSolInc(n)/real(nsum+1,R8)
             windAvg(n) = windInc(n)/real(nsum+1,R8)
             ! warmMax(n) = max(DTiter,warmMaxInc(n))
             windMax(n) = max(u10n,windMaxInc(n))

             nsum = 0

             qSolInc(n) = Qsol
             windInc(n) = u10n

             ! warmMaxInc(n) = 0.0_R8
             windMaxInc(n) = 0.0_R8
          endif

          nInc(n) = real(nsum,R8) ! set nInc to incremented or reset nsum

          if (present(ustar_sv)) ustar_sv(n) = ustar
          if (present(re_sv   )) re_sv   (n) = re
          if (present(ssq_sv  )) ssq_sv  (n) = ssq

       else              ! mask = 0

          !------------------------------------------------------------
          ! no valid data here -- out of domain
          !------------------------------------------------------------
          warm       (n) = spval
          salt       (n) = spval
          speed      (n) = spval
          regime     (n) = spval
          tBulk      (n) = spval
          tSkin      (n) = spval
          tSkin_night(n) = spval
          tSkin_day  (n) = spval
          cSkin      (n) = spval
          cSkin_night(n) = spval
          warmMax    (n) = spval
          windMax    (n) = spval
          qSolAvg    (n) = spval
          windAvg    (n) = spval
          warmMaxInc (n) = spval
          windMaxInc (n) = spval
          qSolInc    (n) = spval
          windInc    (n) = spval
          nInc       (n) = 0.0_R8
          sen   (n)    = spval  ! sensible         heat flux  (W/m^2)
          lat   (n)    = spval  ! latent           heat flux  (W/m^2)
          lwup  (n)    = spval  ! long-wave upward heat flux  (W/m^2)
          evap  (n)    = spval  ! evaporative water flux ((kg/s)/m^2)
          taux  (n)    = spval  ! x surface stress (N)
          tauy  (n)    = spval  ! y surface stress (N)
          tref  (n)    = spval  ! 2m reference height temperature (K)
          qref  (n)    = spval  ! 2m reference height humidity (kg/kg)
          duu10n(n)    = spval  ! 10m wind speed squared (m/s)^2

          if (present(ustar_sv)) ustar_sv(n) = spval
          if (present(re_sv   )) re_sv   (n) = spval
          if (present(ssq_sv  )) ssq_sv  (n) = spval

       endif   ! mask
    end DO     ! loop over n

  end subroutine flux_atmOcn_diurnal

   ! ===================================================================

   real(R8) elemental function cuberoot(a)
      real(R8), intent(in) :: a
      real(R8), parameter :: one_third = 1._R8/3._R8
      cuberoot = sign(abs(a)**one_third, a)
   end function cuberoot


end module flux_atmocn_diurnal_mod
