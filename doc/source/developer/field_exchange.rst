.. _field-exchange:

==============
Field exchange
==============

The mediator sets up all of its field exchange during initialization: which
fields it exchanges with each component, how each field is mapped between grids,
and how mapped fields are merged. To follow that setup, start with the NUOPC
two-phase field mechanism it is built on.

Advertise and realize
=====================

NUOPC (the coupling layer built on ESMF) splits component initialization into a
two-phase field setup, so components do not need to know each other's internals
up front:

* **Advertise** — each component lists the *names* of the fields it can
  potentially send or receive (for example ``sea_surface_temperature``,
  ``air_pressure``), without allocating any real memory or grid data behind
  them. Think of it as everyone posting a menu of "here is what I can offer /
  here is what I need," so the coupler (the mediator and driver) can figure out
  who talks to whom.
* **Realize** — after the coupler matches producers to consumers based on those
  advertised names, each component goes back and actually creates the real data
  structures for the fields that got matched: allocating arrays, attaching them
  to a grid or mesh, and connecting them to memory. Fields that were advertised
  but never matched to a partner are simply dropped instead of wastefully
  allocated.

Splitting setup this way keeps components decoupled from each other's grids and
internals during advertise, and only pays the memory and setup cost in realize,
once a field is known to be needed.

Within this scheme the mediator and the components advertise differently: the
mediator advertises **all possible fields** that could be exchanged with a
component — a set that is kept consistent with the ``fd_<host>.yaml`` field
dictionary — whereas each component advertises only the fields it actually
expects to send and receive.

.. important::

   Because the mediator advertises this superset while each component advertises
   only its own subset, NUOPC connects only the **intersection** — a field is
   realized only where both the mediator and a component advertised it. This is
   why the ``addfld`` calls below advertise so many fields: the mediator offers
   the entire menu, and the actual component set and configuration determine
   which offerings are taken up.

The field-exchange module
=========================

The mediator's side of all this lives in one place: the host's **field-exchange
module**, ``esmFldsExchange_<host>_mod.F90``. It declares *which* fields are
exchanged and *how* they are mapped and merged, using three routines provided by
the generic ``esmFlds.F90``:

* ``addfld`` — **advertise** a field the mediator can receive from, or send to, a
  component.
* ``addmap`` — declare **how** a source field is mapped from its grid to a
  destination grid.
* ``addmrg`` — declare **how** one or more mapped source fields are merged into a
  destination field.

Understanding these three routines is the key to the phases that follow: the
phases simply execute, every coupling step, the maps and merges declared here.

When ``esmFldsExchange`` runs
=============================

``esmFldsExchange_<host>`` is called twice during initialization, selected by a
``phase`` argument:

* ``phase == 'advertise'`` — makes all the ``addfld`` calls. This is the
  mediator's contribution to the NUOPC **advertise** phase above: it advertises
  the full superset of exchangeable fields, and NUOPC then determines which are
  connected and realized.
* ``phase == 'initialize'`` — runs after the connected fields have been realized,
  and makes the ``addmap`` and ``addmrg`` calls that set up mapping and merging
  for the fields that actually got connected.

Components are referred to by a **component index** — one of ``compatm``,
``compice``, ``compglc``, ``complnd``, ``compocn``, ``comprof`` or ``compwav``.

``addfld`` — advertise fields
=============================

The mediator advertises every field it *might* receive from or send to a
component. There are two directions:

.. code-block:: fortran

   call addfld_from(comp_index, 'field_name')   ! a field the mediator receives FROM comp_index
   call addfld_to(comp_index, 'field_name')     ! a field the mediator sends TO comp_index

``addfld_from`` populates the "from" field list for a component
(``fldListFr(comp_index)``) and ``addfld_to`` the "to" field list
(``fldListTo(comp_index)``). Advertising a field does not guarantee it is used —
only advertised fields that are connected on both sides are realized.

``addmap`` — map a source field to a destination
================================================

For each source field that must reach another component's grid, ``addmap``
declares the mapping:

.. code-block:: fortran

   call addmap_from(comp_index_src, 'field_name', comp_index_dst, maptype, mapnorm, mapfile)

A single source field may be mapped to several destinations, so it may appear in
more than one ``addmap`` call. The arguments are:

**maptype** — the interpolation method (from ESMF). Common values include:

* ``mapbilnr`` — bilinear
* ``mapconsf`` — first-order conservative, conservative-fraction normalization
* ``mapconsd`` — first-order conservative, destination-area normalization
* ``mappatch`` — patch recovery
* ``mapfcopy`` — redistribution (identical source and destination grids)
* ``mapnstod`` — nearest source-to-destination (and the
  ``mapnstod_consd`` / ``mapnstod_consf`` combinations)

**mapnorm** — the mapping normalization:

* ``unset`` — no normalization (only valid with ``mapfcopy``)
* ``none`` — no normalization (when ``maptype`` is not ``mapfcopy``)
* ``one`` — normalize by one
* ``lfrin`` / ``ifrac`` / ``ofrac`` — normalize by the land, ice or ocean
  fraction in ``FBFrac(...)`` (used when mapping land, ice or ocean fields to the
  atmosphere)
* ``custom`` — normalization is done by hand in the ``prep`` phase (used, for
  example, mapping glc to land)

.. note::

   When a fractional ``mapnorm`` is used, the field is **scaled by the relevant
   fraction before mapping and unscaled by the mapped fraction after mapping**.
   This fraction weighting is what keeps ocean/ice-to-atmosphere mappings
   accurate where surface fractions vary — see :ref:`surface fractions
   <concepts>`.

**mapfile** — whether to build the routehandle online or read weights from a file:

* ``unset`` — generate the routehandle at run time
* a pathname (supplied as a driver attribute) — read the mapping weights from
  that file

Example::

   call addmap_from(compice, 'Si_snowh', compatm, mapconsf, 'ifrac', 'unset')

This maps the sea-ice field ``Si_snowh`` conservatively to the atmosphere, using
ice-fraction normalization, with the routehandle generated at run time.

``addmrg`` — merge mapped fields into a destination
===================================================

Once source fields have been mapped onto a destination grid, ``addmrg`` declares
how they combine into the field the mediator exports to that component:

.. code-block:: fortran

   call addmrg_to(comp_index_dst, 'dst_field_name', &
                  mrg_from=comp_index_src, mrg_fld='src_field_name', &
                  mrg_type='...', mrg_fracname='...')

The arguments are:

* ``mrg_from`` — the source component index.
* ``mrg_fld`` — the source field name in that component's mapped field bundle.
* ``mrg_type`` — how the source contributes:

  * ``copy`` — copy the mapped source field into the destination.
  * ``copy_with_weights`` — copy, weighting by the source fraction on the
    destination grid (requires ``mrg_fracname``).
  * ``sum`` — cumulative sum of the mapped source fields.
  * ``sum_with_weights`` — cumulative sum, each source weighted by its fraction
    (requires ``mrg_fracname``).
  * ``merge`` — fraction-weighted merge of several sources into one destination.

* ``mrg_fracname`` — the fraction field in ``FBFrac(comp_index_dst)`` used as the
  weight for the weighted/merge types.

**Multi-source merges** are expressed as **one ``addmrg_to`` call per source**,
all targeting the same destination field. For example, the merged surface
temperature the mediator sends to the atmosphere is built from land, sea ice and
ocean::

   call addmrg_to(compatm, 'Sx_t', mrg_from=complnd, mrg_fld='Sl_t', &
                  mrg_type='merge', mrg_fracname=mrg_fracname_lnd2atm_state)
   call addmrg_to(compatm, 'Sx_t', mrg_from=compice, mrg_fld='Si_t', &
                  mrg_type='merge', mrg_fracname='ifrac')
   call addmrg_to(compatm, 'Sx_t', mrg_from=compocn, mrg_fld='So_t', &
                  mrg_type='merge', mrg_fracname='ofrac')

Each call adds one fraction-weighted contribution; together they form the
``Sx_t`` field over each atmosphere cell (the land/ice/ocean fractions sum to
one).

How it fits together
====================

For a single coupled field, the three routines line up like this:

#. ``addfld_from`` / ``addfld_to`` advertise the field on the source and
   destination sides.
#. ``addmap_from`` declares how the source field reaches each destination grid —
   these become the routehandles (``RH``) and populate the off-diagonal
   ``FBImp(n,k)`` bundles.
#. ``addmrg_to`` declares how the mapped contributions combine into ``FBExp`` for
   the destination component.

The run phases (next) do no field-selection logic of their own — they walk these
declarations and execute the maps and merges. See :ref:`code-organization` for
the ``FBImp`` / ``FBExp`` / ``FBFrac`` structures these routines populate, and
:ref:`fields` for the field-name convention.
