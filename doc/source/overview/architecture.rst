.. _architecture:

=================================
Coupling architecture at a glance
=================================

This page sketches how a CMEPS-coupled run is organized: how the mediator code
is grouped, the *phases* it runs, and how those phases fit into a coupled time
loop. Each topic here is treated in depth in the Developer Guide.

How the code is organized
=========================

The mediator source lives in ``mediator/`` and falls into three broad groups.

**Generic infrastructure** — the reusable machinery that implements mediator
functionality independent of any particular host:

* The mediator component itself: initialization, the phase registrations, and
  the run loop.

  * ``med.F90``

* Mapping (interpolation) of fields between grids.

  * ``med_map_mod.F90``

* Merging of mapped source fields into a single destination field.

  * ``med_merge_mod.F90``

* Computing and maintaining the surface fractions.

  * ``med_fraction_mod.F90``

* The mediator's internal state and field bundles.

  * ``med_internalstate_mod.F90``

* I/O, helpers, time management, utilities, and kind/constant definitions.

  * ``med_io_mod.F90``, ``med_methods_mod.F90``, ``med_time_mod.F90``,
    ``med_utils_mod.F90``, ``med_constants_mod.F90``, ``med_kind_mod.F90``

* Water and energy budget diagnostics.

  * ``med_diag_mod.F90``

* History output, restarts, and timing.

  * ``med_phases_history_mod.F90``, ``med_phases_restart_mod.F90``,
    ``med_phases_profile_mod.F90``

**Application-specific exchange code** — determines *which* fields are
exchanged and *how* they are mapped and merged for a given host:

* The per-host field-exchange logic (which connections exist and how they map
  and merge).

  * ``esmFldsExchange_cesm_mod.F90`` (CESM / NorESM),
    ``esmFldsExchange_ufs_mod.F90`` (UFS),
    ``esmFldsExchange_hafs_mod.F90`` (HAFS)

* Shared definitions used by the exchange modules.

  * ``esmFlds.F90``

**Phase code** — the per-component work carried out each coupling step:

* *Prep* phases: map and merge fields from one or more source components into
  the fields the mediator exports to a destination component.

  * ``med_phases_prep_<comp>_mod.F90``

* *Post* phases: handle fields the mediator has just imported from a component —
  for example mapping them to the grids where they will later be used.

  * ``med_phases_post_<comp>_mod.F90``

* Atmosphere/ocean flux calculation.

  * ``med_phases_aofluxes_mod.F90``

* Ocean albedo calculation.

  * ``med_phases_ocnalb_mod.F90``

* Inline CDEPS functionality. This is used to fill unmapped regions
  that can arise in regional coupled configurations due to domain
  mismatch (i.e. between atmosphere-ocean)

  * ``med_phases_cdeps_mod.F90``

Here ``<comp>`` is one of ``atm``, ``ocn``, ``ice``, ``lnd``, ``rof``,
``wav`` or ``glc``.

.. note::

   The split is not perfectly clean: a few "generic" modules — notably
   ``med_phases_prep_ocn_mod.F90`` and ``med_fraction_mod.F90`` — also contain
   application-specific blocks.

Mediator phases
===============

The mediator does its work through **phases**: named routines that NUOPC calls
in a defined order. Broadly there are three kinds:

The **initialization phases** are run once and do the following:

* advertise the fields each component will exchange with the mediator — the
  host's field-exchange configuration determines which fields are connected;
* realize those fields on the mediator's meshes, accepting each component's grid
  and transferring it to an ESMF mesh;
* create the mapping routehandles needed to do both conservative and
  non-conservative mapping between source and destination components;
* compute the initial surface fractions; and
* resolve data dependencies so that every field has an initial value before the
  run loop begins.

The **run phases** are executed every coupling step and are dominated by the
``prep_*`` and ``post_*`` phases for each component, plus the flux, albedo and
fraction updates.

The **finalization phases** are run once at the end.

The *order* in which phases run — and the coupling intervals between components
— is not hard-coded in the mediator. It is set by the driver's **NUOPC run
sequence**, which is what makes different component combinations and coupling
frequencies possible without changing mediator code.

A coupled step, conceptually
============================

Within a single coupling interval, a typical sequence looks like:

#. Each component advances and returns its fields to the mediator.
#. A ``post_<comp>`` phase ingests the fields the mediator has just imported
   from that component (mapping as needed).
#. The mediator updates the surface fractions.
#. Optionally, the mediator can compute atmosphere/ocean fluxes and ocean albedos, depending
   on the component configuration.
#. A ``prep_<comp>`` phase maps and merges the required sources into the fields
   the mediator exports to each destination component, accumulating and
   averaging across coupling intervals as required.
#. The mediator writes history output and (periodically) restarts. For certain
   configurations it also writes water and energy budget diagnostics.

The details — which fields, which mappings, which merges — come from the
host's ``esmFldsExchange_*`` module and are the subject of the Developer Guide.

Host applications
=================

CMEPS serves several modeling systems from the same generic core:

* **CESM / NorESM** — the configuration this documentation treats as primary.
* **UFS** — NOAA's Unified Forecast System.
* **HAFS** — NOAA's Hurricane Analysis and Forecast System.

They share the generic infrastructure and phase machinery, and differ mainly in
their field dictionaries and ``esmFldsExchange_*`` exchange logic. Host-specific
behavior is flagged as such throughout the User and Developer guides.
