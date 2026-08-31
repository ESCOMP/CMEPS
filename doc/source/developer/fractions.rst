.. _surface-fractions:

=================
Surface fractions
=================

Because a single grid cell can be part land, part ocean and part sea ice, the
mediator tracks the **fraction** of each surface type on every component grid.
These fractions are the weights in essentially every :ref:`merge <merging>` and
in the fractional normalization of many :ref:`mappings <mapping>`, so getting
them right is central to conservation (:ref:`concepts`). The machinery lives in
``med_fraction_mod.F90``, and the fractions are stored per grid in the
internal-state fraction bundles ``FBfrac(comp)``.

Key assumptions
===============

The fraction bookkeeping rests on a few assumptions that are important to keep in
mind:

* **The ocean and sea ice must share the same grid** (and the same mask). Ice and
  ocean fractions are therefore computed together and mapped together.
* **The ocean grid has no partial land.** Every ocean gridcell is either entirely
  ocean or entirely land — the ocean mask is binary (0 or 1) per cell, never a
  fractional value.
* **The land fraction is derived from the ocean mask.** The fraction of land on
  the land grid is obtained by mapping the ocean mask onto the land grid (land is
  the complement of the mapped ocean fraction), rather than being prescribed
  independently.
* **The cross-grid maps have specific forms.** The ocean-to-atmosphere and
  ice-to-atmosphere maps are *masked* maps, while the land-to-atmosphere map is a
  *global* map.

The fraction fields
===================

The fractions for each component are held in an ESMF ``FieldBundle``. These are
the internal-state array ``FBfrac(:)`` — one bundle per component, allocated to
``ncomps`` in ``med_internalstate_mod.F90`` ("fraction data for various
components, on their grid"). They are initialized and subsequently updated in
``med_fraction_mod.F90``. Each such fraction bundle can carry the following
fields (not all on every grid):

* ``lfrac`` — fraction of land on the grid
* ``ofrac`` — fraction of ocean on the grid
* ``ifrac`` — fraction of sea ice on the grid
* ``ifrad`` / ``ofrad`` — the ice / ocean fractions as of the last radiation
  timestep

The set carried on each grid is fixed per component, for example:

.. code-block:: none

   atm : ifrac ofrac lfrac aofrac
   ocn : ifrac ofrac ifrad ofrad
   ice : ifrac ofrac
   lnd : lfrac
   glc : gfrac lfrac
   rof : rfrac lfrac

On the atmosphere grid the surface fractions are self-consistent and sum to one:

.. math::

   f_{lfrac} + f_{ofrac} + f_{ifrac} \approx 1

How the fractions are set
=========================

Relative versus total fractions
-------------------------------

Component fractions delivered at run time are **relative** fractions — the
fraction of the *active region* a surface covers, not the fraction of the whole
gridcell. For example, if a cell can be at most 50% ice and 50% land, an ice
model reporting "half covered" sends an ice fraction of 0.5 (a relative value
between 0 and 1). The mediator corrects these to **total** fractions so that, on
every grid,

.. math::

   f_{ifrac} + f_{ofrac} + f_{lfrac} = 1

The fractions are seeded from a few fields the surface components provide (read
from the diagonal ``FBImp`` bundles):

* the ocean mask ``So_omask`` (``FBImp(compocn,compocn)``)
* the land-grid land fraction ``Sl_lfrin`` (``FBImp(complnd,complnd)``), itself a
  map of the ocean mask onto the land grid
* the ice fraction ``Si_ifrac`` and ice mask ``Si_imask``
  (``FBImp(compice,compice)``)

In the sequences below, ``mapX2Y`` denotes a mapping from the ``X`` grid to the
``Y`` grid — for example ``mapo2a`` (ocean → atmosphere), ``mapl2a`` (land →
atmosphere), ``mapa2l`` (atmosphere → land), ``mapi2o`` (ice → ocean), ``mapi2a``
(ice → atmosphere), ``mapl2r`` (land → river) and ``mapl2g`` (land → land-ice).
These fraction and mask maps are all **first-order conservative** maps.

Initialization (``med_fraction_init``)
--------------------------------------

At initialization the ice fraction is taken to be zero everywhere; the ocean and
land-grid fractions are taken from the components, mapped to the atmosphere, and
the remaining fractions derived:

.. code-block:: none

   fractions_*(ifrac) = 0                            ! ice fraction initially zero
   fractions_o(ofrac) = So_omask                     ! ocean mask from the ocean
   fractions_l(lfrac) = Sl_lfrin                     ! land fraction (ocean-mask complement
                                                     !   on the land grid)

   fractions_a(ofrac) = mapo2a(fractions_o(ofrac))   ! ocean fraction -> atmosphere grid
   fractions_a(lfrac) = 1 - fractions_a(ofrac)       ! land is the complement
                                                     !   (truncated to 0 below 0.001)

   fractions_r(lfrac) = mapl2r(fractions_l(lfrac))   ! land fraction -> river grid
   fractions_g(lfrac) = mapl2g(fractions_l(lfrac))   ! land fraction -> land-ice grid

The land-fraction truncation (values below 0.001 set to zero) helps preserve
non-land gridcells as fully non-land.

Run time (``med_fraction_set``)
-------------------------------

Every coupling step the time-evolving fractions are recomputed from the current
ice fraction — this is what the ``post_ice`` phase triggers when it "updates the
surface fractions":

.. code-block:: none

   fractions_i(ifrac) = Si_ifrac                     ! ice fraction from ice model
   fractions_i(ofrac) = 1 - fractions_i(ifrac)       ! ocean = complement of ice
                                                     !   (relative -> total corrected)
   fractions_o(ifrac) = mapi2o(fractions_i(ifrac))   ! -> ocean grid
   fractions_o(ofrac) = mapi2o(fractions_i(ofrac))
   fractions_a(ifrac) = mapi2a(fractions_i(ifrac))   ! -> atmosphere grid
   fractions_a(ofrac) = mapi2a(fractions_i(ofrac))

Because the land fraction is static, only the ice and ocean fractions are
updated at run time.

Where the fractions are used
============================

The fractions feed the two operations of the previous pages:

* **Merges.** A merge *to* the atmosphere is weighted by ``fractions_a`` (``lfrac``,
  ``ofrac``, ``ifrac``); a merge *to* the ocean is weighted by ``fractions_o``
  (``ofrac``, ``ifrac``) normalized to one.
* **Mapping normalization.** Fraction-weighted mappings scale by a fraction
  before mapping and unscale by the mapped fraction after (:ref:`mapping`) — for
  example ``mapo2a`` scales by ``fractions_o(ofrac)`` and divides by
  ``fractions_a(ofrac)``; ``mapi2a`` uses ``ifrac``; ``mapl2a`` uses ``lfrac``.
