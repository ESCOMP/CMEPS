module flux_atmocn_driver_mod.F90

contains

  subroutine flux_atmOcn_driver(logunit, nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   &
       &               qbot,  rainc ,s16O  ,sHDO  ,s18O  ,rbot,  &
       &               tbot  ,us    ,vs,  pslv,  &
       &               ts    ,mask  , seq_flux_atmocn_minwind, &
       &               sen   ,lat   ,lwup  ,   &
       &               r16O, rhdo, r18O, &
       &               evap  ,evap_16O, evap_HDO, evap_18O, &
       &               taux  ,tauy  ,tref  ,qref  ,   &
       &               ocn_surface_flux_scheme, &
       &               add_gusts, &
       &               duu10n, &
       &               ugust_out, &
       &               u10res, &
       &               ustar_sv ,re_sv ,ssq_sv,   &
       &               missval)

    !!.................................................................
    !! ocn_surface_flux_scheme = 0 : Default CESM1.2
    !!                         = 1 : COARE algorithm
    !!                         = 2 : UA algorithm (separate subroutine)
    !!.................................................................

    ! Default flux scheme.
    if (ocn_surface_flux_scheme == 0) then

       call flux_atmOcn_large(&
            logunit, nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   &
            qbot,  rainc ,s16O  ,sHDO  ,s18O  ,rbot,  &
            tbot  ,us    ,vs,  pslv,  &
            ts    ,mask  , seq_flux_atmocn_minwind, &
            sen   ,lat   ,lwup  ,   &
            r16O, rhdo, r18O, &
            evap  ,evap_16O, evap_HDO, evap_18O, &
            taux  ,tauy  ,tref  ,qref  ,   &
            add_gusts, &
            duu10n, &
            ugust_out, &
            u10res, &
            ustar_sv ,re_sv ,ssq_sv,   &
            missval)

    else if (ocn_surface_flux_scheme == 1) then

       call flux_atmOcn_COARE(&
            logunit, nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   &
            qbot,  rainc ,s16O  ,sHDO  ,s18O  ,rbot,  &
            tbot  ,us    ,vs,  pslv,  &
            ts    ,mask  , seq_flux_atmocn_minwind, &
            sen   ,lat   ,lwup  ,   &
            r16O, rhdo, r18O, &
            evap  ,evap_16O, evap_HDO, evap_18O, &
            taux  ,tauy  ,tref  ,qref  ,   &
            add_gusts, &
            duu10n, &
            ugust_out, &
            u10res, &
            ustar_sv ,re_sv ,ssq_sv,   &
            missval)

    else if (ocn_surface_flux_scheme == 2) then

       call flux_atmOcn_UA(&
            logunit, nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   &
            qbot,  rainc ,s16O  ,sHDO  ,s18O  ,rbot,  &
            tbot  ,us    ,vs,  pslv,  &
            ts    ,mask  , seq_flux_atmocn_minwind, &
            sen   ,lat   ,lwup  ,   &
            r16O, rhdo, r18O, &
            evap  ,evap_16O, evap_HDO, evap_18O, &
            taux  ,tauy  ,tref  ,qref  ,   &
            add_gusts, &
            duu10n, &
            ugust_out, &
            u10res, &
            ustar_sv ,re_sv ,ssq_sv,   &
            missval)

       do n=1,nMax
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
