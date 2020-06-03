.. _field_naming_convention:

Application Specific Field Exchange Specification
=================================================

CMEPS is a community mediator in that it supports multiple coupling systems.
CMEPS contains two coupled model specific files that determine:
- the fields that are exchanged between components
- how source fields are mapped to destination fields
- how source fields are merged after mapping to destination fields

This occurs via the following files:

- ``esmFldsExchange_cesm_mod.F90``
- ``esmFldsExchange_nems_mod.F90``
- ``esmFldsExchange_hafs_mod.F90``

CMEPS advertises **all possible fields** that can be imported to and
exported by the mediator for the target coupled system. Not all of
these fields will be connected to the various components. The
connections will be determined by what the components advertise in
their respective advtise phase.

The mediator variable names can be seen in the application specific
YAML field dictionary. Currently, three field dictionaries are
supported for the target coupled model applications:

- ``fd_cesm.yaml`` for CESM
- ``fd_nems.yaml`` for UFS-S2S
- ``fd_hafs.yaml`` for UFS-HAFS

Field Naming Convention
-----------------------

The CMEPS field name convention in these YAML files is independent of the model components.
The convention differentiates between variables that are state fields versus flux fields.

The CMEPS naming convention assumes the following one letter designation for the various components as
well as the mediator ::

  import to mediator
  ==================
  a => atmosphere
  i => sea-ice
  l => land
  g => land-ice
  o => ocean
  r => river
  w => wave

  export from mediator (after  mapping and merging)
  ==================
  x => mediator

State Variables
~~~~~~~~~~~~~~~

State variables have a prefix that always start with an ``S`` followned by a two character string::

  state-prefix
    first 3 characters: Sx_, Sa_, Si_, Sl_, So_
    one letter indices: x,a,i,g,l,o,g,r,w
  state-name
    what follows state prefix

As an example, ``Sx_t`` is the merged surface temperature from land, ice and ocean sent to the atmopshere for CESM.

The following state names are currently supported. Note that each application might only use a subset of these fields.

.. csv-table:: "Atmospheric State Names (import to mediator)"
   :header: "stat name", "description"
   :widths: 20, 60

   "Sa_co2diag", ""
   "Sa_co2prog", ""
   "Sa_dens", "air density at lowest model layer"
   "Sa_pbot", "air pressure at lowest model layer"
   "Sa_pslv", "air pressure at land and sea surface"
   "Sa_ptem", ""
   "Sa_shum", "air specific humidity at lowest model layer"
   "Sa_tbot", "air temperature at lowest model layer"
   "Sa_topo", ""
   "Sa_u", "air zonal wind at lowest model layer"
   "Sa_v", "air meridioinal wind at lowest model layer"
   "Sa_z", "air height wind at lowest model layer"

.. csv-table:: "Sea Ice State Names (import to mediator)"
   :header: "name", "description"
   :widths: 20, 60

   "Si_anidf", "ice near infrared diffuse albedo"
   "Si_anidr", "ice near infrared direct albedo"
   "Si_avsdf", "ice vissible diffuse albedo"
   "Si_avsdr", "ice vissible direct albedo"
   "Si_ifrac", "ice fraction"
   "Si_imask", "ice mask"
   "Si_ifrac_n", ""
   "Si_imask", ""
   "Si_qref", ""
   "Si_qref_wiso", ""
   "Si_t", ""
   "Si_tref", ""
   "Si_u10", ""
   "Si_vice", ""
   "Si_snowh", ""
   "Si_vsno", ""

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
   "So_anidr", ""
   "So_avsdf", ""
   "So_avsdr", ""
   "So_bldepth", ""
   "So_dhdx", ""
   "So_dhdy", ""
   "So_duu10n", ""
   "So_fswpen", ""
   "So_ofrac", ""
   "So_omask", ""
   "So_qref", ""
   "So_re", ""
   "So_s", ""
   "So_ssq", ""
   "So_t", ""
   "So_tref", ""
   "So_u", ""
   "So_u10", ""
   "So_ustar", ""
   "So_v", ""

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

   "Sw_hstokes", ""
   "Sw_lamult", ""
   "Sw_ustokes", ""
   "Sw_vstokes", ""

.. csv-table:: "Mediator State Names (export from mediator)"
   :header: "name", "description"
   :widths: 20, 60

   "Sx_anidf", ""
   "Sx_anidr", ""
   "Sx_avsdf", ""
   "Sx_avsdr", ""
   "Sx_qref", ""
   "Sx_t", ""
   "Sx_tref", ""
   "Sx_u10", ""

State Variables
~~~~~~~~~~~~~~~

Flux variables specify both source and destination components and have a 5 character prefix::

  flux-prefix
    first 5 characters: Flmn_
    lm => between components l and m
    n  => computed by component n
  flux-name
     what follows the flux-prefix

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
   "_swpen_vdf", "flux of visible diffuse shortwave thrrouh ice to ocean"
   "_swpen_idr", "flux of near infrared direct through ice to ocean"
   "_swpen_idf", "flux of near infrared diffuse through ice to ocean"
   "_taux", "zonal stress, positive downwards"
   "_tauy", "air-ice meridional stress, positive downwards"
   "_q", "ice-ocn freezing melting potential"
