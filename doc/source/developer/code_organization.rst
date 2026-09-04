.. _code-organization:

=================
Code organization
=================

This page describes how the mediator code is laid out and introduces the
**internal state** — the data structure at the heart of the mediator. It builds
on the high-level grouping in :ref:`architecture`; read that first if you have
not. As noted in the :ref:`Developer Guide introduction <developer>`, this page
assumes familiarity with ESMF and the NUOPC Layer.

The ``mediator/`` directory
===========================

All of the shared mediator code lives in ``mediator/``. It falls into three
groups (described in the :ref:`architecture overview <architecture>`):

* **Generic infrastructure** — the reusable machinery:

  * the mediator component (``med.F90``)
  * the internal state (``med_internalstate_mod.F90``)
  * mapping (``med_map_mod.F90``)
  * merging (``med_merge_mod.F90``)
  * fractions (``med_fraction_mod.F90``)
  * I/O (``med_io_mod.F90``)
  * history output (``med_phases_history_mod.F90``)
  * restarts (``med_phases_restart_mod.F90``)
  * time management (``med_time_mod.F90``)
  * methods and helpers (``med_methods_mod.F90``)
  * shared utilities (``med_utils_mod.F90``)
  * constants (``med_constants_mod.F90``)
  * kind definitions (``med_kind_mod.F90``)
* **Application-specific exchange code** — decides *which* fields are exchanged
  and *how* they are mapped and merged:

  * the per-host exchange modules (``esmFldsExchange_<host>_mod.F90``)
  * shared exchange definitions (``esmFlds.F90``)
  * the per-host field dictionaries (``fd_<host>.yaml``)

* **Phase code** — the per-component and cross-cutting work run each coupling
  step:

  * per-component prep phases (``med_phases_prep_<comp>_mod.F90``)
  * per-component post phases (``med_phases_post_<comp>_mod.F90``)
  * the atmosphere/ocean flux interface (``med_phases_aofluxes_mod.F90``)
  * the ocean albedo phase (``med_phases_ocnalb_mod.F90``)
  * the diagnostics phase (``med_diag_mod.F90``)
  * the timing/profiling phase (``med_phases_profile_mod.F90``)

The host-specific driver (under ``cesm/``) and the UFS-specific flux code (under
``ufs/``) live outside ``mediator/`` and are not part of the shared mediator; see
:ref:`concepts` for the driver/mediator distinction.

The mediator component
======================

``med.F90`` is the mediator's NUOPC component — a ``NUOPC_Mediator`` (its
``SetServices`` specializes the generic ``NUOPC_Mediator``). It registers the
initialization and run phases with NUOPC (the phases the run sequence names),
and it owns the mediator's **internal state**. The registration and the
initialization sequence are covered on the phase-lifecycle page (next); this
page focuses on the data the phases operate on.

The internal state
==================

The mediator's working data is held in a single **internal state** derived type
(defined in ``med_internalstate_mod.F90`` and carried on the component as
``is_local%wrap``). Understanding its field bundles is the key to reading almost
any phase.

For a run with ``ncomps`` components, the central members are:

Per-component states
   * ``NStateImp(n)`` — the ESMF **import state** from component ``n``.
   * ``NStateExp(n)`` — the ESMF **export state** to component ``n``.

The import/export field bundles
   * ``FBImp(n,k)`` — a **matrix** of field bundles: the fields imported from
     component ``n``, interpolated to the grid of component ``k``. The diagonal
     ``FBImp(n,n)`` points into component ``n``'s fields on its own grid (the
     fields in ``NStateImp(n)``); the off-diagonal ``FBImp(n,k)`` holds those
     fields after mapping onto grid ``k``. This matrix is what lets a ``prep``
     phase pull any source component's fields already on the destination grid.
   * ``FBExp(n)`` — the field bundle of fields to be exported to component
     ``n``, on component ``n``'s grid; its fields point into ``NStateExp(n)``,
     the export state the connector sends to component ``n``.

Mediator-computed field bundles
   * ``FBMed_aoflux_a`` / ``FBMed_aoflux_o`` — the atmosphere/ocean flux outputs
     on the atmosphere and ocean grids.
   * ``FBMed_ocnalb_a`` / ``FBMed_ocnalb_o`` — the ocean albedos on the
     atmosphere and ocean grids.

Fractions and areas
   * ``FBfrac(n)`` — the surface fractions for component ``n`` on its grid.
   * ``FBArea(n)`` — component areas (used by the mediator history writes).

Mapping
   * ``RH(:,:,:)`` — the **routehandles** for each pair of components and each
     mapping type (conservative, bilinear, and so on). These are created once
     during initialization and reused every coupling step.

Accumulators
   * ``FBExpAccumOcn`` / ``FBExpAccumWav`` (with their counters
     ``ExpAccumOcnCnt`` / ``ExpAccumWavCnt``) — accumulators used to average an
     export over several fast coupling intervals before it is sent to a
     component coupled less frequently (the accumulate-then-average pattern seen
     in the :ref:`run sequence <run-sequence>` for the ocean).

The coupling matrix
-------------------

Which components exchange fields — and therefore which off-diagonal
``FBImp(n,k)`` entries are populated and which ``prep`` / ``post`` work runs — is
captured by a **coupling matrix** that the mediator prints at start-up. Rows are
the *source* (``from``) component and columns are the *destination* (``to``)
component; a ``T`` marks an active coupling from the row component to the column
component.

CMEPS prints two versions, which mirror the :ref:`advertise/realize
<field-exchange>` superset-vs-intersection idea.

The **allowed** coupling matrix is every coupling CMEPS supports (the superset):

.. code-block:: none

   from  to ->  med  atm  lnd  ocn  ice  rof  wav  glc1 glc2
   med           -    -    -    -    -    -    -    -    -
   atm           -    -    T    T    T    -    T    -    -
   lnd           -    T    -    -    -    T    -    T    T
   ocn           -    T    -    -    T    -    T    T    T
   ice           -    T    -    T    -    -    T    -    -
   rof           -    -    T    T    T    -    -    -    -
   wav           -    T    -    T    T    -    -    -    -
   glc1          -    -    T    -    T    T    -    -    -
   glc2          -    -    T    -    T    T    -    -    -

The **active** coupling matrix is the subset actually turned on for a given
configuration — mirroring how NUOPC realizes only the connected fields. Two
examples follow.

A fully coupled system without an active wave component (its row and column drop
out):

.. code-block:: none

   from  to ->  med  atm  lnd  ocn  ice  rof  wav  glc1 glc2
   med           -    -    -    -    -    -    -    -    -
   atm           -    -    T    T    T    -    -    -    -
   lnd           -    T    -    -    -    T    -    T    T
   ocn           -    T    -    -    T    -    -    T    T
   ice           -    T    -    T    -    -    -    -    -
   rof           -    -    T    T    T    -    -    -    -
   wav           -    -    -    -    -    -    -    -    -
   glc1          -    -    T    -    T    T    -    -    -
   glc2          -    -    T    -    T    T    -    -    -

A configuration with active **CAM**, **CLM** and **CICE** (atm, lnd, ice) and
**data** ocean (DOCN) and land-ice (DGLC):

.. code-block:: none

   from  to ->  med  atm  lnd  ocn  ice  rof  wav  glc1 glc2
   med           -    -    -    -    -    -    -    -    -
   atm           -    -    T    -    T    -    -    -    -
   lnd           -    T    -    -    -    T    -    T    T
   ocn           -    T    -    -    T    -    -    -    -
   ice           -    T    -    -    -    -    -    -    -
   rof           -    -    T    -    T    -    -    -    -
   wav           -    -    -    -    -    -    -    -    -
   glc1          -    -    T    -    T    T    -    -    -
   glc2          -    -    T    -    T    T    -    -    -

Notice how the **data** components change the picture. Nothing needs to *force*
the data ocean, so ``atm -> ocn``, ``ice -> ocn`` and ``rof -> ocn`` all go
inactive — but DOCN still provides its prescribed state *to* the atmosphere and
ice (``ocn -> atm``, ``ocn -> ice`` stay active). A data component is coupled
where it is a *source*, and decoupled where it would be a *destination* of
fields it does not consume.

Reading the matrix ties directly back to the internal state. Each active
coupling is stored as one off-diagonal ``FBImp`` entry:

.. code-block:: none

   an active coupling   atm -> ocn        (a "T" in the matrix above)
   is stored as         FBImp(atm, ocn)   =  atm fields on the ocn grid

and it ties back to the phases. For each active coupling ``X -> Y``:

* ``post_X`` maps ``X``'s fields onto ``Y``'s grid (filling ``FBImp(X,Y)``), and
* ``prep_Y`` merges those mapped fields — with the other active sources into
  ``Y`` — into ``FBExp(Y)``.

So a **row** shows everything a component feeds (its ``post`` targets), and a
**column** shows everything that feeds a component (its ``prep`` sources).
