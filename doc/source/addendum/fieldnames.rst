.. _field_naming_convention:

CMEPS field names
=================

The following state names are currently supported. Note that each application might only use a subset of these fields.

.. csv-table:: "Atmospheric State Names (import to mediator)"
   :header: "stat name", "description"
   :widths: 20, 60

   "Sa_co2diag", "diagnostic CO2 at the lowest model level"
   "Sa_co2prog", "prognostic CO2 at the lowest model level"
   "Sa_dens", "air density at lowest model layer"
   "Sa_pbot", "air pressure at lowest model layer"
   "Sa_pslv", "air pressure at land and sea surface"
   "Sa_ptem", "potential temperature at lowest model layer"
   "Sa_shum", "air specific humidity at lowest model layer"
   "Sa_tbot", "air temperature at lowest model layer"
   "Sa_topo", "surface topographic height"
   "Sa_u", "air zonal wind at lowest model layer"
   "Sa_v", "air meridional wind at lowest model layer"
   "Sa_z", "air height wind at lowest model layer"

.. csv-table:: "Sea Ice State Names (import to mediator)"
   :header: "name", "description"
   :widths: 20, 60

   "Si_anidf", "sea ice near infrared diffuse albedo"
   "Si_anidr", "sea ice near infrared direct albedo"
   "Si_avsdf", "sea ice visible diffuse albedo"
   "Si_avsdr", "sea ice visible direct albedo"
   "Si_ifrac", "sea ice fraction"
   "Si_imask", "sea ice land mask"
   "Si_ifrac_n", "ice fraction by thickness category"
   "Si_qref", "reference height specific humidity"
   "Si_qref_wiso", "reference specific water isotope humidity at 2 meters"
   "Si_t", "sea ice surface temperature"
   "Si_tref", "reference height temperature"
   "Si_u10", "10m wind speed"
   "Si_vice", "volume of sea ice per unit area"
   "Si_snowh", "surface snow water equivalent"
   "Si_vsno", "volume of snow per unit area"

.. csv-table:: "Land State Names (import to mediator)"
   :header: "name", "description"
   :widths: 20, 60

   "Sl_anidf", ""
   "Sl_anidr", ""
   "Sl_avsdf", ""
   "Sl_avsdr", ""
   "Sl_ddvel", ""
   "Sl_fv", ""
   "Sl_fztop", ""
   "Sl_lfrac", ""
   "Sl_lfrin", ""
   "Sl_qref", ""
   "Sl_qref_wiso", ""
   "Sl_ram1", ""
   "Sl_snowh", ""
   "Sl_snowh_wiso", ""
   "Sl_t", ""
   "Sl_topo_elev", ""
   "Sl_topo", ""
   "Sl_tsrf_elev", ""
   "Sl_tsrf", ""
   "Sl_tref", ""
   "Sl_u10", ""

.. csv-table:: "Ocean State Names (import to mediator)"
   :header: "name", "description"
   :widths: 20, 60

   "So_blddepth", "ocean boundary layer depth"
   "So_anidf", "ocean near infrared diffuse albedo"
   "So_anidr", "ocean near infrared direct albedo"
   "So_avsdf", "ocean visible diffuse albedo"
   "So_avsdr", "ocean visible direct albedo"
   "So_bldepth", "ocean mixed layer depth"
   "So_dhdx", "sea surface slope in meridional direction"
   "So_dhdy", "sea surface slope in zonal direction"
   "So_duu10n", "10m wind speed"
   "So_fswpen", "shortwave penetration through sea ice (all bands)"
   "So_ofrac", "ocean fraction"
   "So_omask", "ocean land mask"
   "So_qref", "reference specific humidity at 2 meters"
   "So_re", "square of exchange coefficient for tracers (mediator aoflux)"
   "So_s", "sea surface salinity"
   "So_ssq", "surface saturation specific humidity in ocean (mediator aoflux)"
   "So_t", "sea surface temperature"
   "So_tref", "reference temperature at 2 meters"
   "So_u", "ocean current in zonal direction"
   "So_u10", "10m wind speed"
   "So_ustar", "friction velocity (mediator aoflux)"
   "So_v", "ocean current in meridional direction"

.. csv-table:: "Land Ice State Names (import to mediator)"
   :header: "name", "description"
   :widths: 20, 60

   "Sg_ice_covered", ""
   "Sg_ice_covered_elev", ""
   "Sg_icemask", ""
   "Sg_icemask_coupled_fluxes", ""
   "Sg_topo", ""
   "Sg_topo_elev", ""

.. csv-table:: "Wave State Names (import to mediator) "
   :header: "name", "description"
   :widths: 20, 60

   "Sw_hstokes", "Stokes drift depth"
   "Sw_lamult", "Langmuir multiplier"
   "Sw_ustokes", "Stokes drift u-component"
   "Sw_vstokes", "Stokes drift v-component"

.. csv-table:: "Mediator State Names (export from mediator)"
   :header: "name", "description"
   :widths: 20, 60

   "Sx_anidf", ""
   "Sx_anidr", ""
   "Sx_avsdf", ""
   "Sx_avsdr", ""
   "Sx_qref", "merged reference specific humidity at 2 meters"
   "Sx_t", "merged ice and ocean surface temperature"
   "Sx_tref", "merged reference temperature at 2 meters"
   "Sx_u10", "merged 10m wind speed"

State Variables
~~~~~~~~~~~~~~~

The following flux prefixes are used:

.. csv-table::
   :header: "flux prefix", "description"
   :widths: 20, 60

   "Faxa\_", "atm flux computed by atm"
   "Fall\_", "lnd-atm flux computed by lnd"
   "Fioi\_", "ice-ocn flux computed by ice"
   "Faii\_", "ice_atm flux computed by ice"
   "Flrr\_", "lnd-rof flux computed by rof"
   "Firr\_", "rof-ice flux computed by rof"
   "Faxx\_", "mediator merged fluxes sent to the atm"
   "Foxx\_", "mediator merged fluxes sent to the ocn"
   "Fixx\_", "mediator merged fluxes sent to the ice"

The following flux-names are used:

.. csv-table::
   :header: "flux name", "description"
   :widths: 20, 60

   "_evap", "air-ice evaporative water flux, positive downwards"
   "_lat", "air-ice latent heat, positive downwards"
   "_lwup", "air-ice surface longwave flux, positive downwards"
   "_sen", "air-ice sensible heat, positive downwards"
   "_swnet", "net short wave, positive downwards"
   "_melth", "net heat flux to ocean from ice"
   "_meltw", "fresh water flux to ocean from ice"
   "_salt", "salt to ocean from ice"
   "_swpen", "flux of shortwave through ice to ocean"
   "_swpen_vdr", "flux of visible direct shortwave through ice to ocean"
   "_swpen_vdf", "flux of visible diffuse shortwave through ice to ocean"
   "_swpen_idr", "flux of near infrared direct through ice to ocean"
   "_swpen_idf", "flux of near infrared diffuse through ice to ocean"
   "_taux", "zonal stress, positive downwards"
   "_tauy", "air-ice meridional stress, positive downwards"
   "_q", "ice-ocn freezing melting potential"
