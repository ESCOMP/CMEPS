.. _api-for-esmflds:

================================
 CMEPS application specfic code
================================

CMEPS contains two coupled model specific files that determine:

* the allowed field names in the mediator and aliases for those names that the components might have
* the fields that are exchanged between components
* how source fields are mapped to destination fields
* how source fields are merged after mapping to destination fields

This occurs via the following files:

* **esmFldsExchange_cesm_mod.F90**
* **esmFldsExchange_nems_mod.F90**
* **esmFldsExchange_hafs_mod.F90**

CMEPS advertises **all possible fields** that can be imported to and
exported by the mediator for the target coupled system. Not all of
these fields will be connected to the various components. The
connections will be determined by what the components advertise in
their respective advertise phase.

The mediator variable names can be seen in the application specific
YAML field dictionary. Currently, three field dictionaries are
supported for the target coupled model applications:

* **fd_cesm.yaml** for CESM
* **fd_nems.yaml** for UFS-S2S
* **fd_hafs.yaml** for UFS-HAFS

Details of the naming conventions and API's of this file can be found
in the description of the :ref:`exchange of fields in
CMEPS<api-for-esmflds>`.

Field Naming Convention
-----------------------

The CMEPS field name convention in these YAML files is independent of the model components.
The convention differentiates between variables that are state fields versus flux fields.

The CMEPS naming convention assumes the following one letter designation for the various components as
well as the mediator. Furthermore, the naming convention distinguishes between state and flux variables.

**import to mediator**::

  a => atmosphere
  i => sea-ice
  l => land
  g => land-ice
  o => ocean
  r => river
  w => wave

**export from mediator (after  mapping and merging)**::

  x => mediator

**State Variables** ::

  State variables have a 3 character prefix followed by the state
  name.  The prefix has the form S[a,i,l,g,o,r,w]_ and is followed by
  the field name.  As an example, ``Sx_t`` is the merged surface
  temperature from land, ice and ocean sent to the atmopshere for CESM.

**Flux variables** ::

  Specify both source and destination components and have a 5 character prefix::
  flux-prefix
    first 5 characters: Flmn_
    lm => between components l and m
    n  => computed by component n
  flux-name
     what follows the flux-prefix
  mediator import flux prefix:
    Faxa_, atm flux computed by atm
    Fall_, lnd-atm flux computed by lnd
    Fioi_, ice-ocn flux computed by ice
    Faii_, ice_atm flux computed by ice
    Flrr_, lnd-rof flux computed by rof
    Firr_, rof-ice flux computed by rof
  mediator export flux prefix:
    Faxx_, mediator merged fluxes sent to the atm
    Foxx_, mediator merged fluxes sent to the ocn
    Fixx_, mediator merged fluxes sent to the ice

Exchange of fields
------------------

The application specific module, ``esmFldsExchange_xxx.F90`` contains
all of the information that determines how the mediator performs the
exchange of fields between components. In particular, this module uses the subroutines
``addfld``, ``addmap`` and ``addmrg`` to do the following:

* ``addfld`` advertises all possible fields that the mediator can send
  to and receive from each component that is part of the target
  application

* ``addmap`` determines how each source field is mapped from its
  source mesh to a target destinatinon mesh. Note that a given source
  field may be mapped to more than one destination meshes and so there
  can be more than one call to ``addmap`` for that source field.

* ``addmrg`` subsequent to the field mapping, how a collection of source fields
  is merged to the target destination field.

This section describes the API for the calls that determine the above
information. All of the API's discussed below use the code in the
module ``esmFlds.F90``.

`addfld`
~~~~~~~~~~
CMEPS advertises all possible fields that it can receive from a component or send to any component via a call to ``addfld``.
The API for this call is:

.. code-block:: Fortran

   call addfld(fldListFr(comp_index)%flds, 'field_name')

   where

   comp_index = component index, can be any of [compatm, compice, compglc, complnd, compocn, comprof, compwav]
   field_name = the field name that will be advertised

`addmap`
~~~~~~~~~~
CMEPS determines how to map each source field from its source mesh to a target destination mesh via a call to ``addmap``.
The API for this call is:

.. code-block:: Fortran

   call addmap(FldListFr(comp_index_src)%flds, 'field_name', comp_index_dst, maptype, mapnorm, mapfile)

where

* ``comp_index_src`` is the  source component index and can be one of [compatm, compice, compglc, complnd, compocn, comprof, compwav]

* ``comp_index_dst`` is the  destination component index and can be one of [compatm, compice, compglc, complnd, compocn, comprof, compwav]

* **maptype** determines the mapping type and can have values of:

  * ``mapbilnr``: bilinear mapping

  * ``mapconsf``: first order conservative mapping with normalization type of conservative fraction.

  * ``mapconsd``: first order conservative mapping with normalization type of conservative fraction.

  * ``mappatch``: patch mapping

  * ``mapfcopy``: redist mapping

  * ``mapnstod``: nearest source to desintation mapping

  * ``mapnstod_consd``: nearest source to destination followed by conservative destination

  * ``mapnstod_consf``: nearest source to destination followed by conservative fraction

* **mapnorm** determines the  mapping normalization and can have values of:

  * ``unset`` : no normalization is set, should only be used only if maptype is 'mapfcopy'

  * ``none``  : no normalization is done, should only be used if maptype is not 'mapfcopy'

  * ``one``   : normalize by 1. (see description below for normalization)

  * ``lfrin`` : normalize by the ``lfrin`` field in FBFrac(complnd) (i.e. ``Ffrac(complnd)[lfrin]``).
    Used to map lnd->atm, (see description of :ref:`fractions<fractions>`).
    Scale the field by ``FBfrac(complnd)[`lfrin`]`` before mapping and unscale it by the mapped ``FBFrac(complnd)[lfrin]`` after mapping.

  * ``ifrac`` : normalize by the 'ifrac' field in FBFrac(compice). Used to map ice->atm, (see description of :ref:`fractions<fractions>`).
    Used to map lnd->atm, (see description of :ref:`fractions<fractions>`).
    Scale the field by ``FBfrac(compice)[`ifrac`]`` before mapping and unscale it by the mapped ``FBFrac(compice)[ifrac]`` after mapping.

  * ``ofrac`` : normalize by the 'ofrac' field in FBFrac(compocn). Used to map ice->atm, (see description of :ref:`fractions<fractions>`).
    Used to map lnd->atm, (see description of :ref:`fractions<fractions>`).
    Scale the field by ``FBfrac(compocn)[`ofrac`]`` before mapping and unscale it by the mapped ``FBFrac(compice)[ofrac]`` after mapping.

  * ``custom`` : custom mapping and normalization will be done in the prep phase for the corresponding field (used to map glc->lnd).

* **mapfile**  determines if a mapping file will be read in or the route handle will be generated at run time:

  * ``unset``  : online route handles will be generated

  * ``<filename>``: read in corresponding full pathname

Fractional normalization is needed to improve the accuracy field
exchanges between ice and ocean and atmosphere.  Consider now the case
where one cell has an ice fraction of 0.3 and the other has a fraction
of 0.5.  Mapping the ice fraction to the atmospheric cell results in a
value of 0.4.  If the same temperatures are mapped in the same way, a
temperature of -1.5 results which is reasonable, but not entirely
accurate.  Because of the relative ice fractions, the weight of the
second cell should be greater than the weight of the first cell.
Taking this into account properly results in a fraction weighted ice
temperature of -1.625 in this example.  This is the fraction
correction that is carried out whenever ocean and ice fields are
mapped to the atmosphere grid.

Time varying fraction corrections are not required in other mappings
to improve accuracy because their relative fractions remain static.

An example of addmap for CESM is shown below:

.. code-block:: Fortran

   call addmap(fldListFr(compice)%flds, 'Si_snowh', compatm, mapconsf, 'ifrac', ice2atm_fmap)

This will create an entry in ``fldListFr(compatm)`` specifying that the
``Si_snowh`` field from the ice should be mapped conservatively to the
atmosphere using fractional normalization where the ice fraction is
obtained from ``Frac(compice)[snowh]``. The ``ice2atm_fmap`` is a
character string obtained as an attribute from the driver and in
general is set to ``unset`` specifying that an run time route handle
for this mapping needs to be created.

`addmrg`
~~~~~~~~~~
CMEPS determines how to map a set of one or more mapped source fields to create the target destination field in the export state.
The API for this call is:

.. code-block:: Fortran

   call addmrg(fldListTo(comp_index_dst)%flds, dst_fieldname, &
               mrg_from1, mrg_fld1, mrg_type1, mrg_fracname1, &
               mrg_from2, mrg_fld2, mrg_type2, mrg_fracname2, &
               mrg_from3, mrg_fld3, mrg_type3, mrg_fracname3, &
               mrg_from4, mrg_fld4, mrg_type4, mrg_fracname4)

where

* ``mrg_fromN``, ``mrgfldN``, ``mrgtypeN`` and ``mrg_fracnameN``, where ``N=[1,2,3,4]``, are optional arguments.
  ``mrgfrom1`` is corresponds to the first source component index (e.g. ``compatm``).

* **mrg_fromN**: is an integer corresponding to the source component index

* **mrg_fldN** : is a character string corresponding to the field name in the mapped field bundle of the source component with index ``mrg_fromN``

* **mrg_typeN**: the type of merging that will be carried out for component with index ``mrg_fromN``. The allowed values are:

  * ``copy``: simply copy the source mapped field into the destination field bundle

  * ``copy_with_weights``: weight the mapped source field by its fraction on the destination mesh.

  * ``sum_with_weights``: do a cumulative sum of all the mapped source fields where each field is weighed by by its fraction on the destination mesh.

  * ``sum_with_weights``: do a cumulative sum of all the mapped source fields.

For ``copy_with_weights`` and ``sum_with_weights``, the mapped source field is weighted by ``mrg_fracnameN`` in ``FBFrac(comp_index_dst)``. If copy_with_weights is chose as the ``mrg_typeN`` value then ``mrg_fracnameN`` is also required as an argument. If sum_with_weights is chose as the ``mrg_typeN`` value then ``mrg_fracnameN`` is also required as an argument.
