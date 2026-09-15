.. _run-config:

==============================
NUOPC attributes (CESM/NorESM)
==============================

This page covers the CESM/NorESM run configuration and the attributes the
mediator obtains from it. For what a run configuration *is* and how it relates to
the driver and the mediator, see :ref:`concepts`.

For CESM/NorESM the run configuration is the ``nuopc.runconfig`` file, which the
driver ingests at start-up. As explained in the concepts, the mediator does not
read this file directly — the driver ingests it and *sets attributes* on the
components, and the mediator queries the attributes it needs.

Importantly, each component obtains its attributes **at initialization** through
calls to ``NUOPC_CompAttributeGet``. This is a general NUOPC mechanism, not
specific to the mediator: the driver makes the run-configuration attributes
available on every component, and each one — CMEPS **and** all of the model
components in the coupled configuration — obtains the attributes it uses this way.

The file is organized into the ``component_list`` (the active components, e.g.
``MED ATM LND ICE OCN ROF GLC WAV``), the :ref:`run sequence <run-sequence>`
``runSeq::`` block, and a set of **ESMF attribute groups**. Each group is written
as ``<name>::`` … ``::`` and holds the attributes for a particular scope:

* ``DRIVER_attributes`` — driver-level settings:

  * the model name
  * timing and log output paths
  * reproducible-sum options
  * the water-vapor saturation scheme
  * the coupler restart pointer (``drv_restart_pointer``)
* ``PELAYOUT_attributes`` — the processor layout:

  * the task and thread counts, root PE and stride for each component
    (``atm_ntasks``, ``atm_rootpe``, …)
  * the number of instances
  * asynchronous-PIO settings

* ``ALLCOMP_attributes`` — attributes **queried by every component**, not just
  the mediator:

  * the orbital parameters (``orb_*``)
  * the scalar-field descriptors (``ScalarField*``) that describe the scalar
    values exchanged between each component and the mediator
  * the per-component model names (``ATM_model``, ``OCN_model``, …)
  * the mesh mask, profiling, and other coupling-wide options

* ``MED_attributes`` — the mediator's own attributes, and by far the largest
  group. These are the attributes the mediator queries; they are detailed in
  **Attributes the mediator queries** below.

* ``CLOCK_attributes`` — time-management attributes, also **queried by every
  component**:

  * the per-component coupling timesteps (``atm_cpl_dt``, ``ocn_cpl_dt``, …)
  * the calendar
  * the start date and time
  * the stop and restart cadence (``stop_*``, ``read_restart``)
* ``ATM_attributes``, ``ICE_attributes``, ``GLC_attributes``,
  ``LND_attributes``, ``OCN_attributes``, ``ROF_attributes``,
  ``WAV_attributes`` — per-component attributes visible to the driver. Most hold
  little more than a verbosity level plus a few component-specific options (for
  example ``aqua_planet`` under ``ATM_attributes``).

Attributes the mediator queries
===============================

During its initialization phases the mediator retrieves attributes with
``NUOPC_CompAttributeGet``. Each lookup reports whether the attribute is present
and set, so optional attributes can fall back to defaults. The attributes CMEPS
consumes fall into a few groups:

* **Coupling mode** — ``coupling_mode`` selects which field exchange specification file is chosen
  (``esmFldsExchange_<host>_mod.F90``) and which custom ``prep`` calculations
  are used.
* **Atmosphere/ocean flux options** — for example ``aoflux_grid`` (the grid on
  which the air-sea fluxes are computed) and ``ocn_surface_flux_scheme``.
* **Diagnostics and budgets** — for example ``do_budgets`` (enable the water and
  energy budget diagnostics).
* **Auxiliary history** — the ``histaux_*`` settings for the mediator's
  auxiliary history streams.
* **Debugging** — for example ``dbug_flag`` and ``check_for_nans``.

Most of these attributes are set for you by the compset and case configuration.
Where you can change them yourself, it is through ``user_nl_cpl`` or
``xmlchange`` as described in :ref:`running-a-case`.

.. note::

   The complete list of attributes CMEPS requires or recognizes — including
   which are required, which are optional, and which are specific to a
   particular host — is given in the :ref:`Reference <reference>`. This page
   orients you to the attributes the mediator queries; the Reference is the
   authoritative catalogue.
