module flux_atmocn_driver_mod

  use shr_kind_mod,          only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN ! shared kinds
  use shr_const_mod,         only : shr_const_spval
  use flux_atmocn_Large_mod, only : flux_atmocn_Large
  use flux_atmocn_COARE_mod, only : flux_atmocn_COARE
  use flux_atmocn_UA_mod,    only : flux_atmocn_UA

  implicit none
  public

contains

  subroutine flux_atmOcn_driver(logunit, nMax, &
       zbot, ubot, vbot, thbot,                &
       qbot,  rainc, rbot,                     &
       tbot, us, vs, pslv,                     &
       ts, mask,  seq_flux_atmocn_minwind,     &
       sen, lat, lwup, evap,                   &
       taux, tauy, tref, qref,                 &
       ocn_surface_flux_scheme,                &
       add_gusts, duu10n, ugust_out, u10res,   &
       ustar_sv, re_sv, ssq_sv, missval)

    !--- input arguments --------------------------------
    integer  , intent(in) :: logunit
    integer  , intent(in) :: nMax        ! data vector length
    integer  , intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
    logical  , intent(in) :: add_gusts
    real(R8) , intent(in) :: zbot (nMax) ! atm level height (m)
    real(R8) , intent(in) :: ubot (nMax) ! atm u wind (m/s)
    real(R8) , intent(in) :: vbot (nMax) ! atm v wind (m/s)
    real(R8) , intent(in) :: thbot(nMax) ! atm potential T (K)
    real(R8) , intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
    real(R8) , intent(in) :: rainc(nMax) ! atm precip for convective gustiness (kg/m^3) - RBN 24Nov2008/MDF 31Jan2022
    real(R8) , intent(in) :: rbot (nMax) ! atm air density (kg/m^3)
    real(R8) , intent(in) :: tbot (nMax) ! atm T (K)
    real(R8) , intent(in) :: pslv (nMax) ! atm sea level pressure(Pa)
    real(R8) , intent(in) :: us   (nMax) ! ocn u-velocity (m/s)
    real(R8) , intent(in) :: vs   (nMax) ! ocn v-velocity (m/s)
    real(R8) , intent(in) :: ts   (nMax) ! ocn temperature (K)
    real(R8) , intent(in) :: seq_flux_atmocn_minwind        ! minimum wind speed for atmocn (m/s)
    integer  , intent(in) :: ocn_surface_flux_scheme

    !--- output arguments -------------------------------
    real(R8),intent(out)  ::  sen  (nMax)     ! heat flux: sensible (W/m^2)
    real(R8),intent(out)  ::  lat  (nMax)     ! heat flux: latent (W/m^2)
    real(R8),intent(out)  ::  lwup (nMax)     ! heat flux: lw upward (W/m^2)
    real(R8),intent(out)  ::  evap (nMax)     ! water flux: evap ((kg/s)/m^2)
    real(R8),intent(out)  ::  taux (nMax)     ! surface stress, zonal (N)
    real(R8),intent(out)  ::  tauy (nMax)     ! surface stress, maridional (N)
    real(R8),intent(out)  ::  tref (nMax)     ! diag:  2m ref height T (K)
    real(R8),intent(out)  ::  qref (nMax)     ! diag:  2m ref humidity (kg/kg)
    real(R8),intent(out)  :: duu10n(nMax)     ! diag: 10m wind speed squared (m/s)^2
    real(R8),intent(out)  :: ugust_out(nMax)  ! diag: gustiness addition to U10 (m/s)
    real(R8),intent(out)  :: u10res(nMax)     ! diag: gustiness addition to U10 (m/s)

    real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
    real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
    real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity (kg/kg)
    real(R8),intent(in) ,optional :: missval        ! masked value

    !--- local variables --------------------------------
    integer  :: n
    real(R8) :: spval  ! local missing value
    !--------------------------------------------------------------------------------

    !!.................................................................
    !! ocn_surface_flux_scheme = 0 : Large and Pond
    !!                         = 1 : COARE algorithm
    !!                         = 2 : UA algorithm
    !!.................................................................

    if (present(missval)) then
       spval = missval
    else
       spval = shr_const_spval
    endif

    ! Default flux scheme.
    if (ocn_surface_flux_scheme == 0) then

       call flux_atmOcn_Large(                    &
            logunit, spval, nMax,                 &
            zbot, ubot, vbot, thbot,              &
            qbot,  rainc,  rbot,                  &
            tbot, us, vs,  pslv,                  &
            ts, mask,  seq_flux_atmocn_minwind,   &
            sen, lat, lwup,  evap,                &
            taux, tauy, tref, qref,               &
            add_gusts, duu10n, ugust_out, u10res, &
            ustar_sv=ustar_sv, re_sv=re_sv, ssq_sv=ssq_sv)

    else if (ocn_surface_flux_scheme == 1) then

       call flux_atmOcn_COARE(                  &
            logunit, spval, nMax,               &
            zbot, ubot, vbot, thbot,            &
            qbot,  rainc, rbot,                 &
            tbot, us, vs,  pslv,                &
            ts, mask,  seq_flux_atmocn_minwind, &
            sen, lat, lwup, evap,               &
            taux, tauy, tref, qref,             &
            duu10n, ugust_out, u10res,          &
            ustar_sv=ustar_sv, re_sv=re_sv, ssq_sv=ssq_sv)

    else if (ocn_surface_flux_scheme == 2) then

       call flux_atmOcn_UA(                 &
            logunit, spval, nMax,           &
            zbot, ubot, vbot, thbot,        &
            qbot, rbot, tbot, us, vs, pslv, &
            ts, mask, sen, lat, lwup, evap, &
            taux, tauy, tref, qref,         &
            duu10n, ustar_sv=ustar_sv, re_sv=re_sv, ssq_sv=ssq_sv)

       do n = 1,nMax
          if (mask(n) /= 0) then
             u10res(n) = sqrt(duu10n(n))
             ugust_out(n) = 0._r8
          else
             u10res   (n) = spval
             ugust_out(n) = spval
          end if
       end do

    end if

  end subroutine flux_atmOcn_driver

end module flux_atmocn_driver_mod
