.. _phase-lifecycle:

===================
The phase lifecycle
===================

The mediator does its work through **phases** — routines that ``med.F90``
registers with NUOPC and that the :ref:`run sequence <run-sequence>` invokes in a
defined order. There are three kinds: initialization phases (run once),
run phases (run every coupling step), and finalization phases (run once at the
end). This page follows a run through those phases; it assumes the field-exchange
declarations from :ref:`field-exchange` and the internal-state structures from
:ref:`code-organization`.

Initialization phases
=====================

Initialization follows the NUOPC initialization protocol. The mediator registers
several phases in ``SetServices`` and NUOPC calls them in sequence:

Advertise
   ``AdvertiseFields`` calls ``esmFldsExchange_<host>`` with ``phase =
   'advertise'``, making all the ``addfld`` calls (see :ref:`field-exchange`).
   This advertises the superset of fields the mediator can exchange; nothing is
   allocated yet.

Realize
   After NUOPC negotiates which advertised fields are connected, the realize
   phases (``RealizeFieldsWithTransferProvided`` and
   ``RealizeFieldsWithTransferAccept``) create the ESMF fields for the connected
   fields. Because the mediator is mesh-based, this is also where component grids
   are **accepted and transferred onto ESMF meshes**, so the mediator has a mesh
   for each connected component.

Data initialize
   ``DataInitialize`` finishes setup and resolves data dependencies. It:

   * creates the mapping **routehandles** (``RH``) for every component pair and
     mapping type that the ``addmap`` declarations require — these are built once
     and reused every coupling step;
   * calls ``esmFldsExchange_<host>`` with ``phase = 'initialize'`` to populate
     the ``addmap`` / ``addmrg`` specifications;
   * initializes the surface **fractions** (``FBfrac``); and
   * ensures every connected field has an initial value before the run loop
     begins.

After initialization the mediator holds: a mesh per component, the realized
import/export states, the routehandle matrix, the fraction bundles, and the map
and merge specifications — everything the run phases need.

Run phases
==========

Every coupling step, the run sequence calls a series of ``MED`` phases. The bulk
of them are the per-component **post** and **prep** phases, which together move
data through the :ref:`FBImp matrix <code-organization>`:

* A **prep_<comp>** phase runs *before* the mediator exports fields to component
  ``<comp>``. It **merges** the already-mapped contributions from each source
  component (reading down the ``FBImp(:,comp)`` column) into ``FBExp(comp)``,
  using the ``addmrg`` declarations and the fraction weights. The connector then
  sends ``FBExp(comp)`` to the component.
* A **post_<comp>** phase runs *after* the mediator has imported fields from
  component ``<comp>``. It takes those fields (the diagonal ``FBImp(comp,comp)``,
  which points into ``NStateImp(comp)``) and **maps** them onto the grids of the
  components that will consume them — filling the off-diagonal
  ``FBImp(comp,k)``. In other words, ``post`` distributes a component's fresh
  output across the destination grids.

Neither kind of phase contains field-selection logic of its own — they execute
the maps and merges declared by :ref:`field-exchange`.

What each phase does
--------------------

The specific mappings and merges per component are:

``prep`` phases (build the export for a component by mapping+merging sources)
   * ``prep_atm`` — merges the land, ocean and ice surface fields into the
     atmosphere export.
   * ``prep_ocn`` — accumulates and merges the atmosphere/ocean fluxes and ice
     fields for the ocean. Split into an **accumulate** sub-phase (run every fast
     step) and an **average** sub-phase (run on the ocean coupling interval), it
     implements the accumulate-then-average pattern for the ocean.
   * ``prep_ice`` — merges atmosphere and ocean fields into the ice export.
   * ``prep_lnd`` — maps the atmosphere fields into the land export.
   * ``prep_rof`` — accumulates the land runoff fields, time-averages them, maps
     them to the river grid, and merges them into the river export.
   * ``prep_wav`` — merges atmosphere, ocean and ice fields into the wave export.
   * ``prep_glc`` — downscales the accumulated land fields onto the land-ice
     elevation classes (see :ref:`land-to-land-ice downscaling <concepts>`).

``post`` phases (distribute a component's imports onto destination grids)
   * ``post_atm`` — maps atmosphere fields to the ocean, ice and land (and wave).
   * ``post_ocn`` — maps ocean fields to the ice, and accumulates the ocean
     contribution destined for land-ice.
   * ``post_ice`` — updates the surface fractions and maps ice fields to the
     atmosphere and ocean.
   * ``post_lnd`` — maps/accumulates land fields destined for the river and
     land-ice.
   * ``post_rof`` — maps river fields to the land, ocean and ice.
   * ``post_glc`` — maps land-ice fields to the land and ocean.
   * ``post_wav`` — maps wave fields to their destinations.

Alongside these, the run phases also compute the atmosphere/ocean fluxes
(``med_phases_aofluxes``), the ocean albedos (``med_phases_ocnalb``), the
fractions, and the diagnostics — each covered in its own section.

Finalization phases
===================

At the end of the run the finalization phases are called once to write final
history and restart output and to release resources.
