.. _merging:

=======
Merging
=======

**Merging** combines one or more mapped source fields into a single destination
field in a component's export bundle. It is declared per destination field by
:ref:`addmrg <field-exchange>` calls in ``esmFldsExchange_<host>_mod.F90`` and
executed by the ``prep`` phases; the machinery lives in ``med_merge_mod.F90``.
Where :ref:`mapping` fills the off-diagonal ``FBImp(n,k)`` entries, merging reads
down a destination column and writes ``FBExp``.

Automatic merging
=================

The ``prep`` phases call ``med_merge_auto``, which is *driven by the ``addmrg``
declarations* rather than by hand-written per-field logic. For each destination
field in ``fldListTo(comp)`` it looks up the source contributions that were
declared for that field — each contribution's source component
(``mrg_from``), source field (``mrg_fld``), merge type (``mrg_type``) and, where
needed, fraction name (``mrg_fracname``) — and applies them in turn to build the
output field. The mapped source fields are taken from the ``FBImp`` bundles
already on the destination grid; the weights are taken from ``FBfrac`` on the
destination grid.

The merge operations
====================

Every merge reduces to one of a small set of per-gridcell operations. Writing
``dst`` for the destination field, ``src`` for the mapped source field, and
``f`` for the fraction weight (``mrg_fracname`` in ``FBfrac`` of the destination):

.. list-table::
   :header-rows: 1
   :widths: 24 26 40

   * - ``mrg_type``
     - operation
     - meaning
   * - ``copy``
     - ``dst = src``
     - overwrite with the source
   * - ``copy_with_weights``
     - ``dst = src * f``
     - copy, fraction-weighted
   * - ``sum``
     - ``dst = dst + src``
     - cumulative sum (unweighted)
   * - ``sum_with_weights``
     - ``dst = dst + src*f``
     - cumulative sum, fraction-weighted
   * - ``merge``
     - ``dst = dst + src*f``
     - fraction-weighted accumulate

``merge`` and ``sum_with_weights`` perform the **same** arithmetic
(``dst = dst + src*f``); the different names document intent. The ``*_with_weights``
and ``merge`` types read their weight ``f`` from the destination's fraction
bundle; ``copy`` and ``sum`` use no weight.

Multi-source (surface) merges
=============================

A field that the mediator forms from several surfaces — for example the surface
temperature the atmosphere sees over a cell that is part land, part ocean and
part sea ice — is declared as **one ``addmrg_to`` call per source**, each with
``mrg_type = 'merge'`` and that surface's fraction. Because ``merge`` *accumulates*
``src*f``, the three calls together compute the fraction-weighted sum

.. code-block:: none

   Sx_t = f_lnd * Sl_t  +  f_ice * Si_t  +  f_ocn * So_t

Since the land, ice and ocean fractions sum to one over each atmosphere cell,
this is a proper fraction-weighted average — the gridcell-average surface
temperature. The first contribution behaves like an accumulate onto a
zero-initialized output, and each subsequent ``merge`` call adds its surface's
weighted contribution.

Merging and conservation
========================

The fraction-weighted merge here is one half of the conservation strategy
developed in :ref:`concepts`: the mediator maps fields with a **fraction-weighted,
normalized** map (:ref:`mapping`) and then merges them with **fraction weights**.
Applied consistently, the two keep coupled fluxes conserved and coupled states as
meaningful gridcell-average values across mismatched grids.
