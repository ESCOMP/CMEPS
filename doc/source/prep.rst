.. _prep_modules:

======================
 CMEPS `prep` modules
======================

The following modules comprise the "prep phase" CMEPS code:

`med_phases_prep_atm_mod.F90`
  prepares the mediator atmosphere export state

`med_phases_prep_ice_mod.F90`
  prepares the mediator ice export state

`med_phases_prep_glc_mod.F90`
  prepares the mediator land-ice export state

`med_phases_prep_lnd_mod.F90`
  prepares the mediator land export state

`med_phases_prep_ocn_mod.F90`
  prepares the mediator ocean export state

`med_phases_prep_rof_mod.F90`
  prepares the mediator river export state

`med_phases_prep_wav_mod.F90`
  prepares the mediator wave export state

Each prep phase module has several sections:

1. Mapping each source field that needs to be mapped to the destination mesh.
   This is obtained from the ``addmap`` calls in ``esmFldsExchange_cesm_mod.F90``.
   Each `prep` module will call the generic routine  ``med_map_FB_Regrid_Norm`` to do this mapping.

2. Merging the set of source fields that have been mapped to the destination mesh.
   This is obtained from the ``addmrg`` calls in ``esmFldsExchange_cesm_mod.F90``.

3. Carrying out optional custom calculations that cannot be specified
   via ``addmap`` or ``addmrg`` calls. Custom calculations are the
   only part of the CMEPS prep phases that can be can be application
   specific. The attribute ``coupling_mode`` is utilized to by the
   prep phases to determine if a particular customization is targeted
   for only one application. Currently prep phase customization
   encompasses the following:

   * ``med_phases_prep_atm``:

     * Calculation of ocean albedos and atmosphere/ocean fluxes (for CESM).
     * Calculation of land, ice and ocean fractions to send to the atmosphere if those components are present.
   * ``med_phases_pre_ice``:

     * Update the scalar data for the time of the next short wave calculation caried out by the atmosphere (this is needed to the
       ice component to determine the zenith angle) (for CESM)
     * applicate of precipitation factor received from the ocean component (for CESM)

   * ``med_phases_prep_glc``:

     * the land-ice component prep phase `ONLY` uses custom code. Land
       import fields that are destined for the land-ice component are
       in elevation classes, whereas the land-ice components requires
       import data that is not in elevation classes. In addition, the
       land-ice component couples at a much longer time scale than the
       land component. The custom code in this module carries out the
       mapping and merged to take data from the land component,
       accumulate it and map the data both in resolution and in the
       compression of elevation class input to non-elevation class
       output. (for CESM)

   * ``med_phases_prep_lnd``:

     * carry out land-ice to land mapping if land-ice is present (for CESM)
     * update the scalar data for the time of the next short
       wave calculation caried out by the atmosphere (this is needed to the
       land component to determine the zenith angle) (for CESM)

   * ``med_phases_prep_ocn``:

     * computation of net shortwave that is sent to the ocean.
     * apply precipitation fractor to scale rain and snow sent to ocean (for CESM)
     * carry out custom merges for NEMS coupling mode (for NEMS)

   * ``med_phases_prep_rof``:

     * reset the irrigation flux to the river model by pulling in
       irrigation out of the rof cells that are proportial to the
       river volumn in each cell (for CESM).

   * ``med_phases_prep_wav``:

     * currently there are no custom calculations.
