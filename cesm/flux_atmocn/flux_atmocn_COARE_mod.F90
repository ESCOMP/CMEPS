module flux_atmocn_COARE_mod

  !-------------------------------------------------------------------------------
  ! PURPOSE:
  !   computes atm/ocn surface fluxes using  COARE v3.0 parametrisation
  !
  ! NOTES:
  !   o all fluxes are positive downward
  !   o net heat flux = net sw + lw up + lw down + sen + lat
  !   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
  !   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
  !
  ! !REVISION HISTORY:
  !   2013-Nov-22: Thomas Toniazzo's adaptation of Chris Fairall's code,
  !    downloaded from
  !    ftp://ftp1.esrl.noaa.gov/users/cfairall/wcrp_wgsf/computer_programs/cor3_0/
  !     * no wave, standard coare 2.6 charnock
  !     * skin parametrisation also off (would require radiative fluxes and
  !      rainrate in input)
  !     * added diagnostics, comments and references
  !-------------------------------------------------------------------------------

  use shr_kind_mod,   only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN ! shared kinds
  use shr_const_mod,  only : shr_const_stebol, shr_const_latvap, shr_const_g
  use shr_const_mod,  only : shr_const_rgas, shr_const_cpdair
  use shr_flux_mod,   only : td0, maxscl, alpha
  use shr_flux_mod,   only : use_coldair_outbreak_mod
  use shr_wv_sat_mod, only : shr_wv_sat_qsat_liquid ! use saturation calculation consistent with CAM

  implicit none
  private

  public :: flux_atmOcn_COARE
  public :: cor30a

  private :: psiuo
  private :: psit_30

  integer :: debug = 0 ! internal debug level

contains

  subroutine flux_atmOcn_COARE(                       &
       logunit, spval, nMax, zbot, ubot, vbot, thbot, &
       qbot,  rainc, rbot, tbot ,us ,vs,  pslv,       &
       ts, mask, seq_flux_atmocn_minwind,             &
       sen, lat, lwup, evap,                          &
       taux ,tauy, tref, qref,                        &
       duu10n, ugust_out, u10res,                     &
       ustar_sv, re_sv, ssq_sv)

    !--- input arguments --------------------------------
    integer  , intent(in) :: logunit
    real(R8) , intent(in) :: spval
    integer  , intent(in) :: nMax        ! data vector length
    integer  , intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
    real(R8) , intent(in) :: zbot (nMax) ! atm level height      (m)
    real(R8) , intent(in) :: ubot (nMax) ! atm u wind            (m/s)
    real(R8) , intent(in) :: vbot (nMax) ! atm v wind            (m/s)
    real(R8) , intent(in) :: thbot(nMax) ! atm potential T       (K)
    real(R8) , intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
    real(R8) , intent(in) :: rainc(nMax) ! atm precip for convective gustiness (kg/m^3) - RBN 24Nov2008/MDF 31Jan2022
    real(R8) , intent(in) :: rbot (nMax) ! atm air density       (kg/m^3)
    real(R8) , intent(in) :: tbot (nMax) ! atm T                 (K)
    real(R8) , intent(in) :: pslv (nMax) ! atm sea level pressure(Pa)
    real(R8) , intent(in) :: us   (nMax) ! ocn u-velocity        (m/s)
    real(R8) , intent(in) :: vs   (nMax) ! ocn v-velocity        (m/s)
    real(R8) , intent(in) :: ts   (nMax) ! ocn temperature       (K)
    real(R8) , intent(in) :: seq_flux_atmocn_minwind        ! minimum wind speed for atmocn      (m/s)

    !--- output arguments -------------------------------
    real(R8),intent(out)  ::  sen  (nMax) ! heat flux: sensible    (W/m^2)
    real(R8),intent(out)  ::  lat  (nMax) ! heat flux: latent      (W/m^2)
    real(R8),intent(out)  ::  lwup (nMax) ! heat flux: lw upward   (W/m^2)
    real(R8),intent(out)  ::  evap (nMax) ! water flux: evap  ((kg/s)/m^2)
    real(R8),intent(out)  ::  taux (nMax) ! surface stress, zonal      (N)
    real(R8),intent(out)  ::  tauy (nMax) ! surface stress, maridional (N)
    real(R8),intent(out)  ::  tref (nMax) ! diag:  2m ref height T     (K)
    real(R8),intent(out)  ::  qref (nMax) ! diag:  2m ref humidity (kg/kg)
    real(R8),intent(out)  :: duu10n(nMax) ! diag: 10m wind speed squared (m/s)^2
    real(R8),intent(out)  :: ugust_out(nMax) ! diag: gustiness addition to U10 (m/s)
    real(R8),intent(out)  :: u10res(nMax) ! diag: gustiness addition to U10 (m/s)

    real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
    real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
    real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)

    !--- local constants --------------------------------
    real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
    real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)
    real(R8),parameter :: zpbl =700.0_R8 ! PBL depth [m] for gustiness parametriz.

    !--- local variables --------------------------------
    integer(IN) :: n               ! vector loop index
    real(R8)    :: vmag            ! surface wind magnitude   (m/s)
    real(R8)    :: ssq             ! sea surface humidity     (kg/kg)
    real(R8)    :: delt            ! potential T difference   (K)
    real(R8)    :: delq            ! humidity difference      (kg/kg)
    real(R8)    :: stable          ! stability factor
    real(R8)    :: rd              ! sqrt of exchange coefficient (momentum)
    real(R8)    :: rh              ! sqrt of exchange coefficient (heat)
    real(R8)    :: re              ! sqrt of exchange coefficient (water)
    real(R8)    :: ustar           ! ustar
    real(R8)    :: qstar           ! qstar
    real(R8)    :: tstar           ! tstar
    real(R8)    :: hol             ! H (at zbot) over L
    real(R8)    :: zo,zot,zoq      ! roughness lengths
    real(R8)    :: hsb,hlb         ! sens & lat heat flxs at zbot
    real(R8)    :: tau             ! stress at zbot
    real(R8)    :: trf,qrf,urf,vrf ! reference-height quantities

    !--- local functions --------------------------------
    real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
    real(R8)    :: Tk     ! dummy arg ~ temperature (K)

    !--- for cold air outbreak calc --------------------------------
    real(R8)    :: tdiff(nMax) ! tbot - ts
    real(R8)    :: vscl

    !--- functions ---
    qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)

    !--- formats ----------------------------------------
    character(*),parameter :: subName = '(flux_atmOcn_COARE) '
    character(*),parameter ::   F00 = "('(flux_atmOcn_COARE) ',4a)"

    if (debug > 0) write(logunit,F00) "enter"

    rh = spval
    hol= spval

    !--- for cold air outbreak calc --------------------------------
    tdiff= tbot - ts

    DO n=1,nMax
       if (mask(n) /= 0) then

          !--- compute some needed quantities ---
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

          call shr_wv_sat_qsat_liquid(ts(n), pslv(n), qsat, ssq)
          ssq = 0.98_R8 * ssq   ! sea surf hum (kg/kg)

          call cor30a(ubot(n),vbot(n),tbot(n),qbot(n),rbot(n), & ! in atm params
               us(n),vs(n),ts(n),ssq,                          & ! in surf params
               zpbl,zbot(n),zbot(n),zref,ztref,ztref,          & ! in heights
               tau,hsb,hlb,                                    & ! out: fluxes
               zo,zot,zoq,hol,ustar,tstar,qstar,               & ! out: ss scales
               rd,rh,re,                                       & ! out: exch. coeffs
               trf,qrf,urf,vrf)                                  ! out: reference-height params

          ! for the sake of maintaining same defs
          hol = zbot(n)/hol
          rd = sqrt(rd)
          rh = sqrt(rh)
          re = sqrt(re)

          !--- momentum flux ---
          taux(n) = tau * (ubot(n)-us(n)) / vmag
          tauy(n) = tau * (vbot(n)-vs(n)) / vmag

          !--- heat flux ---
          sen (n) = hsb
          lat (n) = hlb
          lwup(n) = -shr_const_stebol * ts(n)**4

          !--- water flux ---
          evap(n) = lat(n)/shr_const_latvap

          !------------------------------------------------------------
          ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------
          tref(n) = trf
          qref(n) = qrf
          duu10n(n) = urf**2+vrf**2

          !------------------------------------------------------------
          ! optional diagnostics, needed for water tracer fluxes (dcn)
          !------------------------------------------------------------
          if (present(ustar_sv)) ustar_sv(n) = ustar
          if (present(re_sv   )) re_sv(n) = re
          if (present(ssq_sv )) ssq_sv(n) = ssq

          u10res(n) = sqrt(duu10n(n))
          ugust_out(n) = 0._r8

       else

          !------------------------------------------------------------
          ! no valid data here -- out of domain
          !------------------------------------------------------------
          sen      (n) = spval  ! sensible         heat flux  (W/m^2)
          lat      (n) = spval  ! latent           heat flux  (W/m^2)
          lwup     (n) = spval  ! long-wave upward heat flux  (W/m^2)
          evap     (n) = spval  ! evaporative water flux ((kg/s)/m^2)
          taux     (n) = spval  ! x surface stress (N)
          tauy     (n) = spval  ! y surface stress (N)
          tref     (n) = spval  !  2m reference height temperature (K)
          qref     (n) = spval  !  2m reference height humidity (kg/kg)
          duu10n   (n) = spval  ! 10m wind speed squared (m/s)^2

          u10res   (n) = spval
          ugust_out(n) = spval

          if (present(ustar_sv)) ustar_sv(n) = spval
          if (present(re_sv   )) re_sv   (n) = spval
          if (present(ssq_sv  )) ssq_sv  (n) = spval
       endif
    enddo

  end subroutine flux_atmOcn_COARE

  !=================================================================

  subroutine cor30a(ubt,vbt,tbt,qbt,rbt, &    ! in atm params
        uss,vss,tss,qss,                 &    ! in surf params
        zbl,zbu,zbt,zrfu,zrfq,zrft,      &    ! in heights
        tau,hsb,hlb,                     &    ! out: fluxes
        zo,zot,zoq,L,usr,tsr,qsr,        &    ! out: ss scales
        Cd,Ch,Ce,                        &    ! out: exch. coeffs
        trf,qrf,urf,vrf)                      ! out: reference-height params

    ! Arguments
    real(R8), intent(in)  :: ubt,vbt,tbt,qbt,rbt
    real(R8), intent(in)  :: uss,vss,tss,qss
    real(R8), intent(in)  :: zbl,zbu,zbt,zrfu,zrfq,zrft
    real(R8), intent(out) :: tau,hsb,hlb
    real(R8), intent(out) :: zo,zot,zoq,L,usr,tsr,qsr
    real(R8), intent(out) :: Cd,Ch,Ce
    real(R8), intent(out) :: trf,qrf,urf,vrf

    ! Local variables
    real(R8):: ua,va,ta,q,rb,us,vs,ts,qs,zi,zu,zt,zq,zru,zrq,zrt ! internal vars

    real(R8):: cpa,rgas,grav,pi,von,beta ! phys. params
    real(R8):: le,rhoa,cpv               ! derived phys. params
    real(R8):: t,visa,du,dq,dt           ! params of problem

    real(R8):: u10,zo10,zot10,cd10,ch10,ct10,ct,cc,ribu,zetu,l10,charn ! init vars
    real(R8):: zet,rr,bf,ug,ut     ! loop iter vars
    real(R8):: cdn_10,chn_10,cen_10  ! aux. output vars

    integer(IN):: i,nits ! iter loop counters

    integer(IN):: jcool                  ! aux. cool-skin vars
    real(R8)   :: dter,wetc,dqer
    !----------------------------------------------------------------

    ua  = ubt  !wind components (m/s) at height zu (m)
    va  = vbt
    ta  = tbt  !bulk air temperature (K), height zt
    Q   = qbt  !bulk air spec hum (kg/kg), height zq
    rb  = rbt  !air density
    us  = uss  !surface current components (m/s)
    vs  = vss
    ts  = tss  !bulk water temperature (K) if jcool= 1, interface water T if jcool= 0
    qs  = qss  !bulk water spec hum (kg/kg) if jcool= 1 etc
    zi  = zbl  !PBL depth (m)
    zu  = zbu  !wind speed measurement height (m)
    zt  = zbt  !air T measurement height (m)
    zq  = zbt  !air q measurement height (m)
    zru = zrfu !reference height for st.diagn.U
    zrq = zrfq !reference height for st.diagn.T,q
    zrt = zrft !reference height for st.diagn.T,q

    !**** constants
    Beta= 1.2_R8
    von = 0.4_R8
    pi  = 3.141593_R8
    grav= SHR_CONST_G
    Rgas= SHR_CONST_RGAS
    cpa = SHR_CONST_CPDAIR

    !*** physical parameters
    Le  = SHR_CONST_LATVAP -.00237e6_R8*(ts-273.16_R8)

    ! cpv = shr_const_cpdair*(1.0_R8 + shr_const_cpvir*Qs) ! form in NCAR code
    cpv = cpa*(1.0_R8+0.84_R8*Q)

    ! rhoa= P/(Rgas*ta*(1+0.61*Q)) ! if input were pressure
    rhoa= rb

    ! parametrisation for air kinematic viscosity (Andreas 1989,p.31)
    t = ta-273.16_R8
    visa= 1.326e-5_R8*(1.0_R8+6.542e-3_R8*t+8.301e-6_R8*t*t-4.84e-9_R8*t*t*t)

    du = sqrt((ua-us)**2+(va-vs)**2)
    dt = ts-ta -.0098_R8*zt
    dq = Qs-Q

    !*** don't use cool-skin params for now, but assign values to Ter and Qer
    jcool=0_IN
    dter=0.3_R8
    wetc=0.622_R8*Le*Qs/(Rgas*ts**2)
    dqer=wetc*dter

    !***************** Begin bulk-model calculations ***************

    !*************** first guess
    ug=0.5_R8

    ut   = sqrt(du*du+ug*ug)
    u10  = ut*log(10.0_R8/1.0e-4_R8)/log(zu/1.0e-4_R8)
    usr  = .035_R8*u10
    zo10 = 0.011_R8*usr*usr/grav+0.11_R8*visa/usr
    Cd10 = (von/log(10.0_R8/zo10))**2
    Ch10 = 0.00115_R8
    Ct10 = Ch10/sqrt(Cd10)
    zot10= 10.0_R8/exp(von/Ct10)
    Cd   =(von/log(zu/zo10))**2
    Ct   = von/log(zt/zot10)
    CC   = von*Ct/Cd

    ! Bulk Richardson number
    Ribu=-grav*zu/ta*((dt-dter*jcool)+.61_R8*ta*dq)/ut**2

    ! initial guess for stability parameter...
    if (Ribu .LT. 0.0_R8) then
       ! pbl-height dependent
       zetu=CC*Ribu/( 1.0_R8 - (.004_R8*Beta**3*zi/zu) * Ribu )
    else
       zetu=CC*Ribu*(1.0_R8 + 27.0_R8/9.0_R8*Ribu/CC)
    endif

    ! ...and MO length
    L10=zu/zetu

    if (zetu .GT. 50.0_R8) then
       nits=1_IN
    else
       nits=3_IN
    endif

    usr =  ut*von/(log(zu/zo10)-psiuo(zu/L10))
    tsr = (dt-dter*jcool)*von/(log(zt/zot10)-psit_30(zt/L10))
    qsr = (dq-dqer*jcool)*von/(log(zq/zot10)-psit_30(zq/L10))

    ! parametrisation for Charney parameter (section 3c of Fairall et al. 2003)
    charn=0.011_R8
    if (ut .GT. 10.0_R8) then
       charn=0.011_R8+(ut-10.0_R8)/(18.0_R8-10.0_R8)*(0.018_R8-0.011_R8)
    endif
    if (ut .GT. 18.0_R8) then
       charn=0.018_R8
    endif
    !*************** end first guess ************

    !***************  iteration loop ************
    do i=1, nits

       ! stability parameter
       zet=-von*grav*zu/ta*(tsr*(1.0_R8+0.61_R8*Q)+.61_R8*ta*qsr)/(usr*usr)/(1.0_R8+0.61_R8*Q)

       ! momentum roughness length...
       zo = charn*usr*usr/grav+0.11_R8*visa/usr

       ! ...& MO length
       L = zu/zet

       ! tracer roughness length
       rr = zo*usr/visa
       zoq= min(1.15e-4_R8,5.5e-5_R8/rr**.6_R8)
       zot= zoq ! N.B. same for vapour and heat

       ! new surface-layer scales
       usr =  ut            *von/(log(zu/zo )-psiuo(zu/L))
       tsr = (dt-dter*jcool)*von/(log(zt/zot)-psit_30(zt/L))
       qsr = (dq-dqer*jcool)*von/(log(zq/zoq)-psit_30(zq/L))

       ! gustiness parametrisation
       Bf=-grav/ta*usr*(tsr+.61_R8*ta*qsr)
       if (Bf .GT. 0.0_R8) then
          ug=Beta*(Bf*zi)**.333_R8
       else
          ug=.2_R8
       endif
       ut=sqrt(du*du+ug*ug)

    enddo
    !***************     end loop    ************

    !******** fluxes @ measurement heights zu,zt,zq ********
    tau= rhoa*usr*usr*du/ut                !stress magnitude
    hsb=-rhoa*cpa*usr*tsr                  !heat downwards
    hlb=-rhoa*Le*usr*qsr                   !wv downwards

    !****** transfer coeffs relative to ut @meas. hts ******
    Cd= tau/rhoa/ut/max(.1_R8,du)
    if (tsr.ne.0._r8) then
       Ch= usr/ut*tsr/(dt-dter*jcool)
    else
       Ch= usr/ut* von/(log(zt/zot)-psit_30(zt/L))
    endif
    if (qsr.ne.0.0_R8) then
       Ce= usr/ut*qsr/(dq-dqer*jcool)
    else
       Ce= usr/ut* von/(log(zq/zoq)-psit_30(zq/L))
    endif

    !**********  10-m neutral coeff relative to ut *********
    Cdn_10=von*von/log(10.0_R8/zo)/log(10.0_R8/zo)
    Chn_10=von*von/log(10.0_R8/zo)/log(10.0_R8/zot)
    Cen_10=von*von/log(10.0_R8/zo)/log(10.0_R8/zoq)

    !**********  reference-height values for u,q,T *********
    urf=us+(ua-us)*(log(zru/zo)-psiuo(zru/L))/(log(zu/zo)-psiuo(zu/L))
    vrf=vs+(va-vs)*(log(zru/zo)-psiuo(zru/L))/(log(zu/zo)-psiuo(zu/L))
    qrf=qs-dq*(log(zrq/zoq)-psit_30(zrq/L))/(log(zq/zoq)-psit_30(zq/L))
    trf=ts-dt*(log(zrt/zot)-psit_30(zrt/L))/(log(zt/zot)-psit_30(zt/L))
    trf=trf+.0098_R8*zrt

  end subroutine cor30a

  !===============================================================================

  real (R8) function psiuo(zet)
     !======================================================================
     !   momentum stability functions adopted in COARE v3.0 parametrisation.
     !   Chris Fairall's code (see cor30a)
     !
     ! !REVISION HISTORY:
     !   22/11/2013: Thomas Toniazzo: comments added
     !======================================================================

    ! !INPUT/OUTPUT PARAMETERS:
    real(R8),intent(in)  :: zet
    real(R8) ::c,x,psik,psic,f
    !-----------------------------------------------------------------
    ! N.B.: z0/L always neglected compared to z/L and to 1
    !-----------------------------------------------------------------
    if(zet>0.0_R8)then
       ! Beljaars & Holtslag (1991)
       c=min(50._R8,.35_R8*zet)
       psiuo=-((1.0_R8+1.0_R8*zet)**1.0_R8+.667_R8*(zet-14.28_R8)/exp(c)+8.525_R8)
    else
       ! Dyer & Hicks (1974) for weak instability
       x=(1.0_R8-15.0_R8*zet)**.25_R8                   ! 15 instead of 16
       psik=2.0_R8*log((1.0_R8+x)/2.0_R8)+log((1.0_R8+x*x)/2.0_R8)-2.0_R8*atan(x)+2.0_R8*atan(1.0_R8)
       ! Fairall et al. (1996) for strong instability (Eq.(13))
       x=(1.0_R8-10.15_R8*zet)**.3333_R8
       psic= 1.5_R8*log((1.0_R8+x+x*x)/3.0_R8)-sqrt(3.0_R8)*atan((1.0_R8+2.0_R8*x)/sqrt(3.0_R8)) &
            & +4.0_R8*atan(1.0_R8)/sqrt(3.0_R8)
       f=zet*zet/(1.0_R8+zet*zet)
       psiuo=(1.0_R8-f)*psik+f*psic
    endif
  END FUNCTION psiuo

  real (R8) function psit_30(zet)
     !===============================================================================
     !   momentum stability functions adopted in COARE v3.0 parametrisation.
     !   Chris Fairall's code (see cor30a)
     !
     ! !REVISION HISTORY:
     !   22/11/2013: Thomas Toniazzo: comments added
     !===============================================================================

    ! !INPUT/OUTPUT PARAMETERS:
    real(R8),intent(in)  :: zet
    ! !EOP
    real(R8) ::c,x,psik,psic,f
    !-----------------------------------------------------------------
    ! N.B.: z0/L always neglected compared to z/L and to 1
    !-----------------------------------------------------------------
    if(zet>0.0_R8)then
       ! Beljaars & Holtslag (1991)
       c=min(50._R8,.35_R8*zet)
       psit_30=-((1.0_R8+2.0_R8/3.0_R8*zet)**1.5_R8+.667_R8*(zet-14.28_R8)/exp(c)+8.525_R8)
    else
       ! Dyer & Hicks (1974) for weak instability
       x=(1.0_R8-15.0_R8*zet)**.5_R8                    ! 15 instead of 16
       psik=2.0_R8*log((1.0_R8+x)/2.0_R8)
       ! Fairall et al. (1996) for strong instability
       x=(1.0_R8-(34.15_R8*zet))**.3333_R8
       psic= 1.5_R8*log((1.0_R8+x+x*x)/3.0_R8)-sqrt(3.0_R8)*atan((1.0_R8+2.0_R8*x)/sqrt(3.0_R8)) &
            & +4.0_R8*atan(1.0_R8)/sqrt(3.0_R8)
       f=zet*zet/(1.0_R8+zet*zet)
       psit_30=(1.0_R8-f)*psik+f*psic
    endif
  end FUNCTION psit_30

end module flux_atmocn_COARE_mod
