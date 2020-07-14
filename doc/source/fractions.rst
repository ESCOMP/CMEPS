.. _fractions:

Fractions
=========

Sets fractions on all component grids.
the fractions fields are now ``afrac, ifrac, ofrac, lfrac, and lfrin``::

 afrac = fraction of atm on a grid
 lfrac = fraction of lnd on a grid
 ifrac = fraction of ice on a grid
 ofrac = fraction of ocn on a grid
 lfrin = land fraction defined by the land model
 ifrad = fraction of ocn on a grid at last radiation time
 ofrad = fraction of ice on a grid at last radiation time

      afrac, lfrac, ifrac, and ofrac:
         are the self-consistent values in the system
      lfrin:
         is the fraction on the land grid and is allowed to
         vary from the self-consistent value as descibed below.
      ifrad and ofrad:
         are needed for the swnet calculation.

    the fractions fields are defined for each grid in the fraction bundles as
      needed as follows.
      character(*),parameter :: fraclist_a = 'afrac:ifrac:ofrac:lfrac:lfrin'
      character(*),parameter :: fraclist_o = 'afrac:ifrac:ofrac:ifrad:ofrad'
      character(*),parameter :: fraclist_i = 'afrac:ifrac:ofrac'
      character(*),parameter :: fraclist_l = 'afrac:lfrac:lfrin'
      character(*),parameter :: fraclist_g = 'gfrac:lfrac'
      character(*),parameter :: fraclist_r = 'lfrac:rfrac'

    we assume ocean and ice are on the same grids, same masks
    we assume ocn2atm and ice2atm are masked maps
    we assume lnd2atm is a global map
    we assume that the ice fraction evolves in time but that
      the land model fraction does not.  the ocean fraction then
      is just the complement of the ice fraction over the region
      of the ocean/ice mask.
    we assume that component fractions sent at runtime
      are always the relative fraction covered.
      for example, if an ice cell can be up to 50% covered in
      ice and 50% land, then the ice domain should have a fraction
      value of 0.5 at that grid cell.  at run time though, the ice
      fraction will be between 0.0 and 1.0 meaning that grid cells
      is covered with between 0.0 and 0.5 by ice.  the "relative" fractions
      sent at run-time are corrected by the model to be total fractions
      such that in general, on every grid,
         fractions_*(afrac) = 1.0
         fractions_*(ifrac) + fractions_*(ofrac) + fractions_*(lfrac) = 1.0
    where fractions_* are a bundle of fractions on a particular grid and
      *frac (ie afrac) is the fraction of a particular component in the bundle.

    the fractions are computed fundamentally as follows (although the
      detailed implementation might be slightly different)

    initialization:
      afrac is set on all grids
        fractions_a(afrac) = 1.0
        fractions_o(afrac) = mapa2o(fractions_a(afrac))
        fractions_i(afrac) = mapa2i(fractions_a(afrac))
        fractions_l(afrac) = mapa2l(fractions_a(afrac))
      initially assume ifrac on all grids is zero
        fractions_*(ifrac) = 0.0
      fractions/masks provided by surface components
        fractions_o(ofrac) = ocean "mask" provided by ocean
        fractions_l(lfrin) = land "fraction  provided by land
      then mapped to the atm model
        fractions_a(ofrac) = mapo2a(fractions_o(ofrac))
        fractions_a(lfrin) = mapl2a(fractions_l(lfrin))
      and a few things are then derived
        fractions_a(lfrac) = 1.0 - fractions_a(ofrac)
             this is truncated to zero for very small values (< 0.001)
             to attempt to preserve non-land gridcells.
        fractions_l(lfrac) = mapa2l(fractions_a(lfrac))
        fractions_r(lfrac) = mapl2r(fractions_l(lfrac))
        fractions_g(lfrac) = mapl2g(fractions_l(lfrac))

    run-time (frac_set):
      update fractions on ice grid
        fractions_i(ifrac) = i2x_i(Si_ifrac)   ice frac from ice model
        fractions_i(ofrac) = 1.0 - fractions_i(ifrac)
          note: the relative fractions are corrected to total fractions
        fractions_o(ifrac) = mapi2o(fractions_i(ifrac))
        fractions_o(ofrac) = mapi2o(fractions_i(ofrac))
        fractions_a(ifrac) = mapi2a(fractions_i(ifrac))
        fractions_a(ofrac) = mapi2a(fractions_i(ofrac))

    fractions used in merging are as follows
    merge to atm   uses fractions_a(lfrac,ofrac,ifrac)
    merge to ocean uses fractions_o(ofrac,ifrac) normalized to one

    fraction corrections in mapping are as follows
      mapo2a uses *fractions_o(ofrac) and /fractions_a(ofrac)
      mapi2a uses *fractions_i(ifrac) and /fractions_a(ifrac)
      mapl2a uses *fractions_l(lfrin) and /fractions_a(lfrin)
      mapl2g weights by fractions_l(lfrac) with normalization and multiplies by fractions_g(lfrac)
      mapa2* should use *fractions_a(afrac) and /fractions_*(afrac) but this
        has been defered since the ratio always close to 1.0

    run time:
        fractions_a(lfrac) + fractions_a(ofrac) + fractions_a(ifrac) ~ 1.0
        0.0-eps < fractions_*(*) < 1.0+eps

   Note that the following FBImp field names are current hard-wired below
   TODO: this needs to be generalized - these names should be set dynamically at run time in the
   source component
      is_local%wrap%FBImp(compglc,compglc) => 'frac'
      is_local%wrap%FBImp(complnd,complnd) => 'Sl_lfrin'
      is_local%wrap%FBImp(compice,compice) => 'Si_imask'
      is_local%wrap%FBImp(compocn,compocn) => 'So_omask'
      is_local%wrap%FBImp(compice,compice) => 'Si_ifrac' (runtime)
