.. _fields:

===============================
Fields and the field dictionary
===============================

Components are connected by matching **standard field names**, not by position
or order (see :ref:`concepts`). This page describes the **field dictionary**
that defines those names, the **naming convention** they follow, and where field
names show up when you configure or debug a run.

The field dictionary
====================

For CESM/NorESM the field dictionary is the YAML file ``mediator/fd_cesm.yaml``.
It is a community-based dictionary of the shared coupling fields: every field
that can be exchanged has an entry giving its standard name, canonical units,
and (optionally) a longer alias and a description. For example:

.. code-block:: yaml

   - standard_name: Faox_evap
     alias: mean_evap_rate_atm_into_ocn
     canonical_units: kg m-2 s-1
     description: med export - atm/ocn evaporation water flux computed in mediator

Each entry has:

* **standard_name** — the short name used throughout CMEPS to make connections
  (for example ``Faox_evap``).
* **alias** — an optional longer, self-describing NUOPC-style name
  (``mean_evap_rate_atm_into_ocn``). Note that a given ``standard_name`` can have multiple aliases.
* **canonical_units** — the units the field is expected in.
* **description** — an optional human-readable description.

A connection between two components is made only when a field's standard name
appears in the export of one component and the import of another; the dictionary
is what guarantees both sides agree on the name and the units.

.. note::

   Each host application has its own field dictionary, and they are **not
   provided the same way**:

   * For **CESM/NorESM**, the dictionary is ``fd_cesm.yaml`` and it **ships with
     the CMEPS code** (``mediator/fd_cesm.yaml``).
   * For **UFS**, the dictionary is ``fd_ufs.yaml`` and it is **not** part of the
     CMEPS code — it is provided by the UFS application.

The naming convention
=====================

The CMEPS field-name convention is independent of the specific model components.
It distinguishes **state** fields from **flux** fields, and encodes which
components are involved using a one-letter designation:

.. code-block:: none

   a => atmosphere      o => ocean
   i => sea ice         r => river
   l => land            w => wave
   g => land-ice        x => mediator (merged, after mapping and merging)

State variables
---------------

A state variable has a three-character prefix ``S<c>_`` — where ``<c>`` is one of
the component letters above — followed by the field name. For example:

* ``Sa_t`` — temperature as a state from the atmosphere.
* ``Sx_t`` — the merged surface temperature (from land, sea ice and ocean) that
  the mediator sends to the atmosphere.

Flux variables
--------------

A flux variable has a five-character prefix ``F<c1><c2><n>_`` followed by the
flux name. Here ``<c1>`` and ``<c2>`` are the two components the flux is between,
and ``<n>`` is the component that computed it. For example:

**Imported to the mediator:**

.. code-block:: none

   Faxa_   atm flux computed by the atmosphere
   Fall_   lnd-atm flux computed by the land
   Fioi_   ice-ocn flux computed by the sea ice
   Faii_   ice-atm flux computed by the sea ice
   Flrr_   lnd-rof flux computed by the river
   Firr_   rof-ice flux computed by the river

**Exported from the mediator (merged):**

.. code-block:: none

   Faxx_   merged fluxes sent to the atmosphere
   Foxx_   merged fluxes sent to the ocean
   Fixx_   merged fluxes sent to the sea ice

So a name like ``Faxa_rainc`` tells you it is a flux (``F``) between the
atmosphere and the mediator, computed by the atmosphere (``Faxa_``) — here a
convective rain flux the atmosphere provides.

Where you encounter field names
===============================

Field names from the dictionary appear in several places you touch when
configuring or debugging a CESM/NorESM run:

* **Auxiliary history streams** — the ``histaux_*_flds`` settings in
  ``user_nl_cpl`` list fields by standard name (see :ref:`running-a-case`).
* **Mediator history and diagnostics output** — field names identify the coupled
  fields the mediator writes.
* **Debugging** — when a connection fails or a field is missing, the message
  refers to the standard name; matching it against ``fd_cesm.yaml`` is the first
  thing to check.

Adding a new field
==================

Adding a new coupled field starts with the field dictionary (a new entry in
``fd_cesm.yaml``) and then wiring it into the exchange logic that advertises,
maps and merges it. That process is a developer task, covered in the
:ref:`Developer Guide <developer>`.
