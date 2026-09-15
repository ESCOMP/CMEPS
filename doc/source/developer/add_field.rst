.. _add-field:

======================
Adding a coupled field
======================

Adding a new field for the mediator to exchange between two components is one of
the most common CMEPS tasks, and it can be done on its own without touching the
rest of the mediator. This page is a practical recipe. It assumes the
:ref:`field-exchange` concepts (advertise/realize, ``addfld`` / ``addmap`` /
``addmrg``) and the :ref:`fields` naming convention.

What you will touch
===================

A new coupled field involves a handful of well-defined edits:

#. the **field dictionary** — declare the field's standard name and units;
#. the host **exchange module** — advertise, map and merge it in
   ``esmFldsExchange_<host>_mod.F90``;
#. the **component caps** — make sure each component actually advertises the
   field so NUOPC connects it.

The first two are inside CMEPS; the third is on the component side (and outside
CMEPS), but the connection will not form without it.

Step 1 — add the field to the dictionary
========================================

Add an entry for the field to the host field dictionary (``fd_cesm.yaml`` for
CESM/NorESM). Give it a standard name that follows the :ref:`naming convention
<fields>`, its canonical units, and optionally an alias and description:

.. code-block:: yaml

   - standard_name: Sa_foo
     canonical_units: K
     description: my new near-surface atmosphere state

Step 2 — advertise, map and merge it
====================================

All of this goes in the host exchange module, ``esmFldsExchange_<host>_mod.F90``.

**Advertise** the field on both sides, in the ``phase == 'advertise'`` block —
from the source component and to the destination component:

.. code-block:: fortran

   call addfld_from(compatm, 'Sa_foo')   ! atmosphere can send it
   call addfld_to(complnd, 'Sa_foo')     ! land can receive it

**Map** the field from the source grid to the destination grid, in the
``phase == 'initialize'`` block. Choose the map type and normalization
appropriate to the field (see :ref:`mapping`):

.. code-block:: fortran

   call addmap_from(compatm, 'Sa_foo', complnd, mapbilnr, 'none', 'unset')

**Merge** the mapped field into the destination's export, also in the
``initialize`` block. For a field that comes from a single source, a plain
``copy`` is enough; a field built from several surfaces uses one ``addmrg_to``
call per source with a fraction weight (see :ref:`merging`):

.. code-block:: fortran

   call addmrg_to(complnd, 'Sa_foo', mrg_from=compatm, mrg_fld='Sa_foo', &
                  mrg_type='copy')

Step 3 — advertise the field in the component caps
==================================================

The mediator advertises a superset of possible fields, but a connection is only
made where the component on the other side advertises the field too
(:ref:`field-exchange`). So the field must also be advertised by the component
caps — the atmosphere cap must advertise ``Sa_foo`` as an **export** and fill it,
and the land cap must advertise it as an **import** and use it. This is
component code, outside CMEPS, but the coupling will silently not form without
it.

Checklist
=========

* Field added to ``fd_<host>.yaml`` with correct name and units.
* ``addfld_from`` / ``addfld_to`` in the ``advertise`` phase.
* ``addmap_from`` with an appropriate map type and normalization.
* ``addmrg_to`` with the appropriate merge type (and fraction, for weighted
  merges).
* The source and destination **component caps** advertise the field.

If the field does not appear in the destination component, the first things to
check are that both caps advertised it (so NUOPC connected it) and that the name
matches the dictionary exactly.
