.. _what-is-cmeps:

==============
What is CMEPS?
==============

CMEPS, the Community Mediator for Earth Prediction Systems, is a
NUOPC-compliant *mediator* that uses `ESMF
<https://earthsystemmodeling.org/>`_ to couple Earth system model components
in a hub-and-spoke arrangement. The component models (atmosphere, ocean, sea
ice, land, river, wave, land-ice) are the spokes; CMEPS is the hub through
which they exchange information.

What a mediator does
====================

As a mediator, CMEPS is responsible for transferring field information from
one component to another. That transfer is rarely a simple copy: it usually
requires one or more of the following operations, each described later in this
documentation:

* **Mapping** (interpolation) of a field from a source component grid to a
  destination component grid.
* **Merging** of contributions from several source components into a single
  field for a destination component (for example, combining land, ocean and
  sea-ice surface fields into a single field seen by the atmosphere).
* **Time averaging and accumulation** of fields exchanged between components
  that are coupled at different frequencies.

The mediator also computes quantities that naturally live "between" components
— surface fractions, atmosphere/ocean fluxes, and ocean albedos — and manages
diagnostics, history output and restarts for the coupled fields.

Hub-and-spoke coupling
======================

Each component coupled through CMEPS is serviced by a NUOPC-compliant *cap* —
a thin software layer between the component's native code and the mediator.
Components exchange information with the mediator through **import** and
**export states**: ESMF containers that wrap native model data together with
metadata such as field names, the underlying grid or mesh, and the parallel
decomposition.

For every field the mediator connects between two components, the field is
placed in the appropriate import or export state within the component's cap.
The cap moves the data between those states and the component's own arrays.

.. important::

   Import and export states are named from the **owner's** point of view, so
   they are reversed between a component and the mediator: a component's
   *export* state is what the mediator *imports*, and the mediator's *export*
   state is what a component *imports*. Keep this in mind whenever "import" or
   "export" appears — see :ref:`concepts` for the full explanation.

.. note::

   CMEPS is itself a **mesh-based** mediator: internally all fields live on
   ESMF meshes. The components it couples, however, may be either grid-based
   or mesh-based — the caps present each component to the mediator as a mesh.

The driver
==========

The mediator and the component models do not run themselves — they are harnessed
by a NUOPC **driver**. The driver is the top-level application component: it
creates the mediator, the component models and the connectors between them, and
it runs the coupled time loop, invoking each component's run phase, each mediator
phase, and each field transfer in the order laid out by the **run sequence**
(described in the :ref:`User Guide <run-sequence>`).

The distinction between the driver and the mediator matters:

* The **mediator (CMEPS)** is shared across applications — the same mediator
  couples the components in CESM, NorESM, UFS and HAFS.
* The **driver** is *not* shared. Each application supplies its own driver (the
  driver under ``cesm/`` in this repository is used only by CESM and NorESM;
  UFS and HAFS use their own).

So CMEPS provides the coupling *machinery*, while the driver provides the
*harness* that assembles a particular application and steps it through time. See
:ref:`architecture` for how the pieces fit together at run time.

A community mediator
====================

Field connections between components rely on matching **standard field
names**. These names are defined in a *field dictionary*. Because CMEPS is a
community mediator serving several modeling systems, the field dictionary — and
some of the exchange logic — is specific to each application (CESM/NorESM, UFS,
HAFS). The concepts in this Overview apply to all of them; where behavior
differs by host, the User and Developer guides call it out explicitly.

Where to go next
================

* New to CMEPS? Continue with :ref:`concepts`, which defines the terms used
  throughout the rest of the documentation.
* Want the big picture of how a run is put together? See
  :ref:`architecture`.
* Configuring or running a coupled case? Jump to the :ref:`User Guide
  <users>`.
* Modifying the mediator? Start with the :ref:`Developer Guide <developer>`.
