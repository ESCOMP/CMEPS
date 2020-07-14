.. _fractions:

Fractions
=========

The component fractions on their corresponding meshes are defined and
updated in ``med_fractions_mod.F90.`` We now describe how these fractions are defined and computed.

CMEPS component fractions are defined as follows:

* Below, ``comp_index`` is the component index and can be one of ``compatm``, ``compice``, ``complnd``, ``compglc``, ``compocn``, ``comprof`` or ``compwav``.

* An array of field bundles, ``Frac(:)`` is created, where the size of
  ``Frac`` corresponds to the number of active components.

* For each active component, a fraction field bundle is created, ``Frac(comp_index)``, where the fields in the field bundle are unique.
  Below, ``Frac(comp_index)[fieldname]`` refers to the field in the ``Frac(comp_index)`` field bundle that has the name ``fieldname``.

* The following are the field names for each component of FBFrac::

    Frac(compatm) = afrac,ifrac,ofrac,lfrac,lfrin
    Frac(compocn) = afrac,ifrac,ofrac,ifrad,ofrad
    Frac(compice) = afrac,ifrac,ofrac
    Frac(complnd) = afrac,lfrac,lfrin
    Frac(compglc) = gfrac,lfrac
    Frac(comprof) = lfrac,rfrac
    Frac(compwav) = wfrac

    where

    afrac = fraction of atm on a grid
    lfrac = fraction of lnd on a grid
    ifrac = fraction of ice on a grid
    ofrac = fraction of ocn on a grid
    lfrin = land fraction defined by the land model
    rfrac = fraction of rof on a grid
    wfrac = fraction of wav on a grid
    ifrad = fraction of ocn on a grid at last radiation time
    ofrad = fraction of ice on a grid at last radiation time

    As an example, ``Frac(compatm)[lfrac]`` is the land fraction on
    the atmosphere mesh.

* ``lfrin`` and ``lfrac`` can be different from ``lfrac`` when the
  atmosphere and land meshes are different.  ``lfrac`` is the land
  fraction consistent with the ocean mask where ``lfrin`` is the land
  fraction in the land component.

* ``ifrad`` and ``ofrad`` are fractions at the last radiation
  timestep.  These fractions preserve conservation of heat in the net
  shortwave calculation because the net shortwave calculation is one
  timestep behind the ice fraction evolution in the system.

The following assumptions are made regarding fractions:

* The ocean and ice are on the same meshes with same masks
* The ice fraction can evolve in time
* The land fraction does not evolve in time
* the ocean fraction is just the complement of the ice fraction over the region
  of the ocean/ice mask.
* The component fractions are always the relative fraction covered.
  For example, if an ice cell can be up to 50% covered in
  ice and 50% land, then the ice domain should have a fraction
  value of 0.5 at that grid cell. At run time though, the ice
  fraction will be between 0.0 and 1.0 meaning that grid cells
  is covered with between 0.0 and 0.5 by ice.  The "relative" fractions
  sent at run-time are corrected by the model to be total fractions
  such that in general, on every mesh cell:

  * ``Frac(:)[afrac]`` = 1.0
  * ``Frac(:)[ifrac]`` + ``Frac(:)[ofrac]`` + ``Frac(:)[lfrac]`` = 1.0

Initialization of the fractions occurs as follows (note that all mapping is first order conservative):

* ``Frac(compatm)[afrac]`` = 1.0

* ``Frac(compocn)[afrac]`` = map atm -> ocn ``Frac(compatm)[afrac]``

* ``Frac(compice)[afrac]`` = map atm -> ice ``Frac(compatm)[afrac]``

* ``Frac(complnd)[afrac]`` = map atm -> lnd ``Frac(compatm)[afrac]``

* ``FBfrac(:)[ifrac]``     = 0.0

* ``Frac(compocn)[ofrac]`` = ocean mask provided by ocean

* ``Frac(complnd)[lfrin]`` = land fraction provided by land

* ``Frac(compatm)[ofrac]`` = map ocn -> atm ``Frac(compocn)[ofrac]``

* ``Frac(compatm)[lfrin]`` = map lnd -> atm ``Frac(complnd)[lfrin]``

* ``Frac(compatm)[lfrac]`` = 1.0 - ``Frac(compatm)[ofrac]``
  (this is truncated to zero for very small values (< 0.001) to attempt to preserve non-land gridcells.)

* ``Frac(complnd)[lfrac]`` = map atm -> lnd ``Frac(compatm)[lfrac]``

* ``Frac(comprof)[lfrac]`` = map lnd -> rof ``Frac(complnd)[lfrac]``

* ``Frac(compglc)[lfrac]`` = map lnd -> glc ``Frac(complnd)[lfrac]``

Run time calculation of fractions is as follows:

* ``Frac(compice)[ofrac]`` = 1.0 - ``Frac(compice)[ifrac]``
  (Note: the relative fractions are corrected to total fractions)

* ``Frac(compocn)[ifrac]`` = map ice -> ocn ``Frac(compice)[ifrac]``

* ``Frac(compocn)[ofrac]`` = map ice -> ocn ``Frac(compice)[ofrac]``

* ``Frac(compatm)[ifrac]`` = map ice -> atm ``Frac(compice)[ifrac]``

* ``Frac(compatm)[ofrac]`` = map ice -> atm ``Frac(compice)[ofrac]``

* ``Frac(compatm)[lfrac]`` + ``Frac(compatm)[ofrac]`` + ``Frac(compatm)[ifrac]`` ~ 1.0
  (0.0-eps < Frac(:)[*] < 1.0+eps)
