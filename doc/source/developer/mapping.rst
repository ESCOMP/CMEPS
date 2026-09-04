.. _mapping:

=======
Mapping
=======

**Mapping** is the interpolation of a field from a source grid to a destination
grid. It is declared per field by an :ref:`addmap <field-exchange>` call in
``esmFldsExchange_<host>_mod.F90`` and executed by the ``post`` and ``prep``
phases; the machinery lives in ``med_map_mod.F90``.

Routehandles
============

A mapping in ESMF is applied through a **routehandle** — a precomputed sparse
matrix of interpolation weights. CMEPS creates all the routehandles it needs
**once, during initialization** (``med_map_routehandles_init``, driven by the
``addmap`` declarations) and reuses them every coupling step. They are stored in
the internal-state array

.. code-block:: none

   RH(n, k, m)   ! a routehandle from grid n to grid k, map type m

indexed by the source grid ``n``, the destination grid ``k``, and the map type
``m``. ``med_map_rh_is_created`` tests whether a given routehandle exists.

A routehandle is created (in ``med_internalstate_mod.F90``) for every
**required** ``FBImp(n,k)`` mapping — that is, for each active coupling that
actually needs one — not for every possible pair of components.

Building weights once and reusing them is essential: weight generation is
expensive, while applying a routehandle (a sparse matrix–vector multiply) is
cheap.

Importantly, the routehandle (weight) creation is performed by **ESMF in
parallel at run time**, and is **scalable** both with the number of processors
and in the memory used. This is what makes it practical to generate the mapping
weights *online* — during initialization, on whatever processor layout the run
uses — rather than relying on precomputed offline weight files.

Map types
=========

Each ``maptype`` from ``addmap`` selects an ESMF regridding method:

============== =========================================================
``maptype``    ESMF regrid method
============== =========================================================
``mapbilnr``   ``ESMF_REGRIDMETHOD_BILINEAR``
``mapconsf``   ``ESMF_REGRIDMETHOD_CONSERVE`` (``NORMTYPE_FRACAREA``)
``mapconsd``   ``ESMF_REGRIDMETHOD_CONSERVE`` (``NORMTYPE_DSTAREA``)
``mappatch``   ``ESMF_REGRIDMETHOD_PATCH``
``mapnstod``   ``ESMF_REGRIDMETHOD_NEAREST_STOD``
``mapfcopy``   redistribution (identical source and destination grids)
============== =========================================================

The ``mapnstod_consd`` / ``mapnstod_consf`` variants apply a nearest
source-to-destination pass followed by a conservative pass, to fill destination
cells a purely conservative map would leave unmapped.

Fractional normalization
========================

State and flux fields that carry a surface fraction must be **fraction-weighted
and normalized** across the mapping to stay accurate and conservative (the
``mapnorm`` argument to ``addmap``). ``med_map_field_normalized`` implements this
as a scale–map–unscale sequence:

#. **scale** the source field by the relevant source fraction (for example
   ``FBfrac(compice)`` for an ice field);
#. **map** both the scaled field and the fraction to the destination grid; and
#. **unscale** the mapped field by the mapped fraction.

This is the "fraction-weighted normalized mapping" developed in
:ref:`concepts` — it is what makes ocean/ice-to-atmosphere mappings correct where
the surface fractions vary from cell to cell. Fields that do not need it are
mapped directly with ``med_map_field``.

Packed mapping
==============

Many fields typically share the same source, destination and map type — that is,
the same routehandle. Rather than apply the routehandle once per field, CMEPS
**packs** all such fields into a single field bundle with an extra undistributed
dimension and applies the routehandle **once** to the whole pack
(``med_map_packed_field_create`` and ``med_map_field_packed``, using the
internal-state ``packed_data(n,k,m)`` structure that mirrors ``RH``). Because
applying a routehandle is a sparse matrix multiply, packing turns many small
multiplies into one large one — a substantial performance win in a coupled run.

Vector mapping
==============

Vector fields (for example wind stress) cannot be interpolated component-by-
component on the sphere without introducing errors near the poles.
``med_map_uv_cart3d`` handles them by rotating the vector into 3-D Cartesian
components, mapping those, and rotating the result back on the destination grid.
