.. _aofluxes:

=======================
Atmosphere/ocean fluxes
=======================

The mediator computes the turbulent air-sea fluxes itself, from the atmosphere
and ocean states brought onto a common grid, using a bulk flux algorithm. The
**concept** —
why the mediator owns this calculation, which fields it produces, and how they
are delivered and averaged — is covered in :ref:`concepts`. This page describes
the **mechanism**: the interface in ``med_phases_aofluxes_mod.F90`` and the grid
on which the calculation is carried out.

Interface versus algorithm
==========================

``med_phases_aofluxes_mod.F90`` is an **interface**: it prepares the inputs,
selects the grid, invokes the flux algorithm, and distributes the outputs. The
bulk-flux **algorithm** itself lives outside the shared mediator, in
host-specific code — ``cesm/flux_atmocn`` for CESM/NorESM and ``ufs/`` for UFS
(see :ref:`concepts`). The mediator reaches it through
``med_aofluxes_init`` / ``med_aofluxes_update``.

The flux outputs are held in two internal-state field bundles,
``FBMed_aoflux_a`` and ``FBMed_aoflux_o`` — the flux fields on the atmosphere and
ocean grids.

The run phase
=============

``med_phases_aofluxes_run`` is the phase the run sequence calls
(``MED med_phases_aofluxes_run``). On the first call it initializes the flux
input/output field bundles (``med_phases_aofluxes_init_fldbuns``) and the flux
scheme; on every call it computes the fluxes with ``med_aofluxes_update``.

Which bulk-flux scheme is used, and its options, come from run-configuration
attributes:

* ``coupling_mode`` — the target modeling system that will compute the
  atmosphere-ocean fluxes. Valid values are ``cesm``, ``noresm``,
  ``ufs.(frac,nfrac).(aoflux)`` and ``hafs``.
* ``aoflux_grid`` — the grid the atmosphere-ocean flux computation is done on:
  ``agrid`` (atmosphere), ``ogrid`` (ocean), or ``xgrid`` (the ESMF exchange
  grid). This is described in detail in `The flux calculation grid`_ below.
* ``ocn_surface_flux_scheme`` — selects the ocean surface flux (bulk)
  algorithm:

  * ``0`` — Large and Pond
  * ``1`` — COARE algorithm
  * ``2`` — UA algorithm

  CESM and UFS use the Large and Pond scheme (``0``); NorESM uses the COARE
  scheme (``1``).
* ``aofluxes_use_shr_wv_sat`` — only applicable to CESM and NorESM. If
  ``true``, ``shr_wv_sat_mod`` is used to calculate the saturation specific
  humidity (``qsat``) for the atmosphere-ocean flux calculations; if ``false``,
  an older inline calculation of ``qsat`` is used, which uses a different
  formulation.
* ``add_gusts`` — only applicable to CESM and NorESM. Adds a wind gustiness
  factor. This should be ``false`` for ``ocn_surface_flux_scheme`` settings of
  ``1`` or ``2``.
* ``flux_max_iteration`` — the maximum number of iterations of the (iterative)
  flux solver.

The flux calculation grid
=========================

The fluxes may be computed on the atmosphere grid, the ocean grid, or a dedicated
**ESMF exchange grid**, selected by the ``aoflux_grid`` attribute; for cases with
**all prognostic components, the exchange grid is the default**. Because the
calculation needs both atmosphere and ocean state on a single grid, the interface
maps the inputs onto the chosen grid and maps the outputs back out. The three
cases:

``aoflux_grid = ogrid`` (ocean grid)
   Requires an atmosphere-to-ocean mapping: map the atmosphere inputs to the
   ocean grid, compute the fluxes there, and map the flux outputs back to the
   atmosphere grid.

``aoflux_grid = agrid`` (atmosphere grid)
   Requires an ocean-to-atmosphere mapping: map the ocean inputs to the
   atmosphere grid, compute the fluxes there, and map the flux outputs back to
   the ocean grid.

``aoflux_grid = xgrid`` (ESMF exchange grid)
   Map **both** the atmosphere and ocean inputs onto the exchange grid, compute
   the fluxes there, and map the flux outputs from the exchange grid to **both**
   the atmosphere and ocean grids.

The **ESMF exchange grid** is an ESMF construct formed from the overlap of the
atmosphere and ocean grids; the interface holds dedicated routehandles for it —
``rh_agrid2xgrid``, ``rh_ogrid2xgrid``, ``rh_xgrid2agrid`` and ``rh_xgrid2ogrid``
(with second-order-conservative and bilinear atmosphere-to-exchange variants).
Computing on the exchange grid avoids biasing the fluxes toward either
component's grid, at the cost of the extra mappings.

Whichever grid is chosen, the resulting fluxes end up in ``FBMed_aoflux_a`` and
``FBMed_aoflux_o`` for the ``prep_atm`` and ``prep_ocn`` phases to merge into the
atmosphere and ocean exports — on each component's grid and coupling interval, as
described in :ref:`concepts`.
