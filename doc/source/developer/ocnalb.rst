.. _ocean-albedo:

============
Ocean albedo
============

As summarized in :ref:`concepts`, the ocean albedo computation produces four
albedos, which are merged with the land and sea-ice albedos and sent to the
atmosphere. This page describes how they are computed, in
``med_phases_ocnalb_mod.F90``.

What is computed
================

Four ocean albedos are produced, spanning the two shortwave bands and the two
illumination types:

* ``So_avsdr`` — visible, direct
* ``So_avsdf`` — visible, diffuse
* ``So_anidr`` — near-infrared, direct
* ``So_anidf`` — near-infrared, diffuse

The albedos are computed **on the ocean grid** and held in the internal-state
field bundle ``FBMed_ocnalb_o``; they are then mapped to the atmosphere grid
(``FBMed_ocnalb_a``) for the ``prep_atm`` merge.

The run phase
=============

The albedo computations depend on the **solar zenith angle**, so they are
recomputed every coupling step for the current zenith angle. The zenith angle is
derived from the Earth's orbital parameters — eccentricity, obliquity and related
quantities — which are initialized by ``med_phases_ocnalb_orbital_init`` and
advanced over the run by ``med_phases_ocnalb_orbital_update``. The run sequence
drives the whole computation through ``med_phases_ocnalb_run``.

Options
=======

The albedo calculation is controlled by a few attributes and settings:

* ``ocean_albedo_scheme`` — the ocean albedo scheme:

  * ``0`` — Briegleb et al. (1986)
  * ``1`` — Taylor et al. (1996)

* ``flux_albav`` — if true, use **averaged** (zenith-angle-independent) diffuse
  and direct albedos instead of the zenith-angle-dependent calculation.
* ``use_min_albedo`` / ``min_albedo`` — apply a minimum floor to the direct
  visible and near-infrared albedos.
* ``albdir`` / ``albdif`` — the 60-degree reference albedos (direct and diffuse)
  used by the scheme.

Once computed and mapped to the atmosphere grid, the ocean albedos are merged
with the land and sea-ice albedos in ``prep_atm`` and sent to the atmosphere for
its shortwave radiation calculation (:ref:`concepts`).
