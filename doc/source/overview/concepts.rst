.. _concepts:

=============
Core concepts
=============

This page defines the vocabulary that the User and Developer guides assume.
Read it once; later pages link back to these terms.

Components and caps
===================

A **component** is a model of one part of the Earth system — atmosphere
(``atm``), ocean (``ocn``), sea ice (``ice``), land (``lnd``), river runoff
(``rof``), waves (``wav``) or land-ice (``glc``). Each component is presented
to the mediator by a NUOPC **cap**, a small software layer that translates
between the component's native data and the ESMF states the mediator
understands.

Unlike a model component, CMEPS (the mediator) has no NUOPC cap of its own —
its data structures are already ESMF/NUOPC-compliant, so no translation layer
is needed.

The driver
==========

The **driver** is the top-level NUOPC application component that harnesses the
mediator and all of the component models. It creates them and the connectors
between them, and it runs the coupled time loop — invoking each component's run
phase, each mediator phase, and each field transfer in the order given by the
**run sequence**.

Crucially, the driver does **no science of its own** — it computes no physical
fields. Its role is purely **control flow and timekeeping**: creating the
components and connectors, advancing the clock, and invoking each run phase in
the prescribed order.

Whenever this documentation refers to "the driver," it means this host-specific
harness — not the mediator; the two are distinguished under **The run sequence**
below.

The run sequence
----------------

The **run sequence** is the ordered recipe that drives a coupled run. It lists,
in order, which component runs, which mediator phase executes, and which fields
are transferred between components — together with the **coupling intervals** at
which each of these happens. The driver (described above) ingests the run
sequence at start-up and uses it to step the coupled system through time.

The run sequence — not the mediator — determines the *order* of operations and
the coupling frequencies; changing which components are coupled, and how often,
is done by changing the run sequence.

.. important::

   Several things around the run sequence are **host-specific**, and confusing
   them with the mediator is a common source of error.

   **The driver — host-specific:**

   * The driver that ingests the run sequence is host-specific. The driver under
     ``cesm/`` in this repository is used only by **CESM and NorESM**; **UFS**
     uses its own driver.

   **The mediator (CMEPS) — shared across applications:**

   * The mediator is shared across applications such as **CESM**, **NorESM** and
     **UFS**; the same mediator code serves them all.
   * The mediator does not parse the run sequence; it simply registers the
     ``MED`` phases that the run sequence names.
   * The mediator obtains only the *attributes* the driver sets from the run
     configuration, not the host-specific configuration file itself.
   * The mediator is designed to be host-agnostic to the greatest extent
     possible — a small number of host-specific ``#ifdef`` blocks remain that
     could not be avoided, but the goal is to keep the mediator common across
     applications.

The run sequence itself is written in a **common NUOPC format** that is the same
whichever application's driver ingests it.

.. note::

   **A common misconception: the mediator is not "driving" the run.** A run
   sequence usually contains many ``MED`` phases — mapping, merging, the
   atmosphere/ocean flux calculation, accumulation, and so on — so it is easy to
   look at it and conclude that the mediator is orchestrating the coupled system.
   It is not. Each ``MED`` line is simply **the driver calling a mediator phase**
   at the point the run sequence specifies. The driver owns the time loop and
   decides what runs, in what order, and how often; the mediator only *provides*
   phases and executes them when it is called. The mediator never advances the
   clock, never chooses the order, and never invokes another component — that is
   entirely the driver's job. Put simply: the **driver is the conductor**,
   deciding who plays when, while the **mediator is the translator and combiner**
   that receives fields from each component, reshapes them (regridding, merging,
   unit conversion) and passes them on. The mediator is a *participant* in the
   run sequence, not its conductor.

The run configuration
---------------------

The **run configuration** is the second key input the driver ingests at
start-up (the run sequence being the first). Whereas the run sequence says *what
runs, in what order, and how often*, the run configuration supplies the
**attributes** that the driver and each component, including the mediator, can
query at run time in order to configure their target behavior.
For the driver, these attributes include the list of active components and the
processor layout.

The relationship to the mediator mirrors that of the run sequence:

* The mediator does **not** read the run-configuration file directly. The driver
  ingests the file and *sets attributes* on the components; the mediator then
  queries the attributes it needs (for example the coupling mode, the grid on
  which atmosphere/ocean fluxes are computed, or whether budget diagnostics are
  enabled, etc.).
* The run-configuration **file is host-specific** and is read by the
  driver. Because the mediator depends only on the *attributes*, not
  on the file that supplied them, the same mediator code works for
  every host.

The mediator
============

The mediator supplies all of the coupling machinery and also computes
atmosphere/ocean fluxes and ocean albedos where appropriate. The following is a
high level overview of the mediator functionality.


Import and export states
------------------------

Components and the mediator communicate through two ESMF states:

* an ESMF **import state**, holding fields flowing *into* a component, and
* an ESMF **export state**, holding fields flowing *out of* a component

An ESMF state wraps native model data as ESMF fields and also carries
metadata: the standard field names, the grid or mesh and its coordinates, and
the parallel decomposition.

.. important::

   **States are named from the component's point of view, so the
   mediator sees them reversed relative to a model component such as
   the atmosphere.** A component **sends out** fields in it's
   **export** state - and the mediator *receives* it, so from the
   mediator's point of view it arrives in the mediator's **import**
   state.  Likewise, a field the mediator sends to a component lives
   in the mediator's **export** state, and the component receives it
   into its **import** state.

   In short:

   * component **export** state  →  mediator **import** state
   * mediator **export** state   →  component **import** state

   Whenever "import" or "export" appears in this documentation, check whose
   state is meant: the direction is opposite for the mediator and the component
   exchanging the field.


Fields and the field dictionary
-------------------------------

Connections between components are made by matching **standard field names**
rather than by position or order. The full set of names, together with their
long names, units and metadata, is the `field dictionary
<https://earthsystemmodeling.org/docs/release/ESMF_8_4_2/NUOPC_refdoc/node3.html#SECTION00032000000000000000>`_.
CMEPS ships a
field dictionary per host application. Adding a new coupled field therefore
starts with the dictionary; this is covered in the Developer Guide.

Three mediator operations
-------------------------

Nearly everything the mediator does is built from three operations:

Mapping (interpolation)
   Moving a field from a source grid to a destination grid using precomputed
   ESMF *routehandles* (sparse-matrix weights). Weights can be conservative,
   bilinear, nearest-neighbor and so on. See the :ref:`mapping <mapping>` page
   in the Developer Guide.

Merging
   Combining several mapped source fields into one destination field — for
   example forming the field the atmosphere sees over a grid cell that is part
   land, part ocean and part sea ice. Merges are weighted by **surface
   fractions**. See the :ref:`merging <merging>` page in the Developer Guide.

Time averaging and accumulation
   Reconciling components coupled at different frequencies by accumulating a
   field over the fast-coupling periods and averaging it before it is passed to
   a component running on a slower coupling period. Accumulation and averaging
   are developed on the :ref:`conservation <conservation>` page in the Developer
   Guide.

Surface fractions
-----------------

Because a single atmosphere/land grid cell can overlap land, ocean and sea ice,
the mediator tracks the **fraction** of each surface type in every cell. For a
given atmosphere cell the land, ocean and ice fractions sum to one:

.. math::

   f_{al} + f_{ao} + f_{ai} = 1

Fractions are dynamic (the ice fraction changes as sea ice grows and melts) and
are themselves mapped between grids. They appear as weights in essentially
every merge, which is why a dedicated module maintains them. See the
:ref:`fractions <surface-fractions>` page in the Developer Guide for details.

Conservation
------------

A guiding requirement of the mediator is **conservation** of mass and energy:
a flux leaving one component must deposit the same integrated quantity in the
component that receives it, even after mapping and merging across mismatched
grids.

CMEPS's mapping weights are **pure area-overlap weights** — static geometry,
computed once and summing to one. Masking and normalization are deliberately
*not* folded into them, because those are **dynamic**: they use the ocean/ice
surface fraction, which changes in time as sea ice grows and melts, and a
time-varying fraction cannot be precomputed into static weights. So masking and
normalization are applied **at run time, when a field is mapped or merged**,
using the surface fractions — which enter in two distinct places:

* **Mapping** a field that covers only part of a cell — for example a
  sea-surface temperature, defined only where there is ocean — fraction-weights
  the field before mapping and normalizes by the mapped fraction, so undefined
  land values are not pulled in and the integral is preserved.
* **Merging** the land, ocean and ice contributions into the single field the
  atmosphere sees weights each surface by its fraction
  (``Fa = fal*Fal + fao*Fao + fai*Fai``).

These two fraction weightings look redundant, and in fact partly cancel; the
full derivation — why fraction-weighted normalized mapping conserves, how area
corrections are applied to fluxes, and how accumulation interacts with mapping —
is summarized in the :ref:`conservation <conservation>` page of the Developer
Guide. It is important background for anyone changing how fields are coupled, but
not required to run a coupled case.

.. note::

   Two rules of thumb that recur in the conservation discussion:

   * Accumulate the **fraction-weighted flux** (``f*F``), never the fraction
     and flux separately — ``sum(f*F)`` is not ``sum(f)*sum(F)``.
   * Area corrections apply to **fluxes only**, are computed once at
     initialization, and do not vary in time.

Beyond these three operations, the mediator also **computes** certain fields
itself rather than merely transferring them between components. These are not
simple map/merge/average steps; each is described in its own section below
(atmosphere/ocean fluxes, ocean albedos and land-to-land-ice downscaling) and
developed further in the :ref:`Developer Guide <developer>`.

Atmosphere/ocean fluxes
-----------------------

The turbulent fluxes at the air-sea interface — wind stress, and the sensible
heat, latent heat and evaporation (moisture) fluxes — are among the most
important quantities in a coupled simulation, and computing them is a central
responsibility of the mediator. Unlike most coupled fields, which are produced
by a component and simply transferred, the atmosphere/ocean fluxes belong to
neither the atmosphere nor the ocean: they describe the *interaction* between
the two. In many configurations they are therefore computed **inside the
mediator**, from the ocean state (for example sea-surface temperature and ocean
currents) and the mapped atmosphere state (for example near-surface winds,
temperature and humidity), using a bulk flux algorithm.

The calculation produces the turbulent air-sea fluxes:

* zonal surface wind stress
* meridional surface wind stress
* sensible heat flux
* latent heat flux
* evaporation (surface water flux)
* upward longwave radiation

Alongside these fluxes it also produces near-surface diagnostic states used by
other components — for example the 2 m reference temperature and humidity, the
10 m wind speed, and the friction velocity.

Computing the fluxes in one place, from a consistent set of states, keeps the
atmosphere and ocean seeing the same interface and makes conservation
tractable: the resulting fluxes are fraction-weighted, mapped conservatively
between grids, and merged with the ice/ocean fluxes so that each component
receives a consistent flux (see the discussion of **Conservation** above).

The fluxes are delivered to the two components on their own grids and at their
own coupling intervals: to the **atmosphere on the atmosphere grid at the
atmosphere coupling interval**, and to the **ocean on the ocean grid at the ocean
coupling interval**. Because the ocean coupling interval is often longer than the
atmosphere coupling interval, the atmosphere/ocean fluxes are **accumulated and
averaged over the ocean coupling interval** before being sent to the ocean — the
accumulate-then-average pattern visible in the run sequence.

Two aspects are configurable:

* **The grid on which the fluxes are computed.** The calculation may be done on
  the atmosphere grid, the ocean grid, or a dedicated *exchange grid* formed from
  the overlap of the atmosphere and ocean grids. The states needed are mapped to
  the chosen grid, and the resulting fluxes are mapped to the grids of the
  components that use them. This choice affects accuracy and cost.
* **The flux algorithm.** Different bulk-flux formulations are supported.

.. note::

   Although it is conceptually a mediator computation, the atmosphere/ocean flux
   calculation code currently lives outside the shared mediator tree, in a
   host-specific location: under ``cesm/`` (``cesm/flux_atmocn``) for CESM/NorESM
   and under ``ufs/`` (``ufs/flux_atmocn_mod.F90``, with a CCPP variant) for UFS.
   The mediator phase that drives it, ``med_phases_aofluxes``, is common; the
   full treatment is on the :ref:`atmosphere/ocean fluxes <aofluxes>` page of the
   Developer Guide.

Ocean albedos
-------------

The mediator computes the ocean surface albedos as a function of the solar
zenith angle. Four different ocean albedos are produced, spanning the two
shortwave bands and the two illumination types:

* direct-beam albedo, visible radiation
* diffuse albedo, visible radiation
* direct-beam albedo, near-infrared radiation
* diffuse albedo, near-infrared radiation

These ocean albedos are then merged with the corresponding albedos obtained from
the land and from the sea ice, and the merged albedos are sent to the atmosphere
for its shortwave radiation calculation. See the :ref:`ocean albedos
<ocean-albedo>` page in the Developer Guide.

Land to land-ice downscaling
----------------------------

The land component sends its fields to the mediator on a set of **elevation
classes**, but the land-ice (``glc``) component does not accept fields on
elevation classes. The mediator therefore **downscales** these fields — mapping
them from the land's elevation-class representation onto the land-ice grid,
using each land-ice cell's actual elevation to select and interpolate the
appropriate value. This is a vertical remapping driven by elevation, not a plain
horizontal interpolation.
