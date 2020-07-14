.. _api-for-esmflds:

Exchange of fields in CMEPS
===========================

The application specific module, ``esmFldsExchange_xxx.F90`` contains
all of the information that determines how the mediator performs the
exchange of fields between components. In particular, this module uses the subroutines
``addfld``, ``addmap`` and ``addmrg`` to do the following:

* ``addfld`` advertises all possible fields that the mediator can send
  to and receive from each component that is part of the target
  application

* ``addmap`` determines how each source field is mapped from its
  source mesh to a target destinatinon mesh. note that a given source
  field may be mapped to more than one destination meshes and so there
  will be more than one call to ``addmap`` for that source field.

* subsequent to the field mapping, how a collection of source fields
  is merged to the target destination field.

This section describes the API for the calls that determine this
information. All of the API's discussed below use the code in the
module ``esmFlds.F90``.

``addfld``
----------
CMEPS advertises all possible fields that it can receive from a component or send to any component via calling ``addfld``.
The API for this call is:

.. code-block:: Fortran

   call addfld(fldListFr(comp_index)%flds, 'field_name')

      where

   comp_index = component index, can be any of [compatm, compice, compglc, complnd, compocn, comprof, compwav]
   field_name = the field name that will be advertised

``addmap``
----------
CMEPS determines how to map each source field from its source mesh to a target destination mesh via calling ``addmap``.
The API for this call is:

.. code-block:: Fortran

   call addmap(FldListFr(comp_index_src)%flds, 'field_name', comp_index_dst, maptype, mapnorm, mapfile)

where

*   comp_index_src: source component index
      can be one of [compatm, compice, compglc, complnd, compocn, comprof, compwav]
*   comp_index_dst: destination component index
      can be one of [compatm, compice, compglc, complnd, compocn, comprof, compwav]
*   maptype: mapping type, can be one of
      * mapbilnr: bilinear mapping
      * mapconsf: first order conservative mapping with normalization type of conservative fraction.
      * mapconsd: first order conservative mapping with normalization type of conservative fraction.
      * mappatch: patch mapping
      * mapfcopy: redist mapping
      * mapnstod: nearest source to desintation mapping
      * mapnstod_consd: nearest source to destination followed by conservative destination
      * mapnstod_consf: nearest source to destination followed by conservative fraction
*   mapnorm: mapping normalization, can be one of
      * unset : no normalization is set, should only be used only if maptype is 'mapfcopy'
      * none  : no normalization is done, should only be used if maptype is not 'mapfcopy'
      * one   : normalize by 1. (see description below for normalization)
      * lfrin : normalize by the 'lfrin' field in FBFrac(complnd). Used to map lnd->atm, (see description of :ref:`fractions<fractions>`).
      * ifrac : normalize by the 'ifrac' field in FBFrac(compice). Used to map ice->atm, (see description of :ref:`fractions<fractions>`).
      * ofrac : normalize by the 'ofrac' field in FBFrac(compocn). Used to map ocn->atm, (see description of :ref:`fractions<fractions>`).
      * custom: custom mapping and normalization will be done in the prep phase for the corresponding field (used to map glc->lnd).
*   mapfile: determines if a mapping file will be read in or the route handle will be generated at run time
      * unset     : online route handles will be generated
      * <filename>: read in corresponding full pathname

``addmrg``
----------
CMEPS determines how to map a set of one or more mapped source fields to create the target destination field in the export state.
The API for this call is:

.. code-block:: Fortran

   call addmrg(fldListTo(comp_index_dst)%flds, dst_fieldname, &
               mrg_from1, mrg_fld1, mrg_type1, mrg_fracname1, &
               mrg_from2, mrg_fld2, mrg_type2, mrg_fracname2, &
               mrg_from3, mrg_fld3, mrg_type3, mrg_fracname3, &
               mrg_from4, mrg_fld4, mrg_type4, mrg_fracname4)

where

* the arguments ``mrg_fromN``, ``mrgfldN``, ``mrgtypeN`` and ``mrg_fracnameN``, where ``N=[1,2,3,4]``, are optional arguments.
* ``mrg_fromN``: is an integer corresponding to the source component index
* ``mrg_fldN`` : is a character string corresponding to the field name in the mapped field bundle of the source component with index ``mrg_fromN``
* ``mrg_type1``: the type of merging that will be carried out for component with index ``mrg_fromN``. The allowed values are:
   * copy: simply copy the source mapped field into the destination field bundle
   * copy_with_weights: weight the mapped source field by its fraction on the destination mesh.
     This is given by the field ``mrg_fracnameN`` in ``FBFrac(comp_index_dst)``.
     If copy_with_weights is chose as the ``mrg_typeN`` value then ``mrg_fracnameN`` is also required as an argument.
   * sum_with_weights: do a cumulative sum of all the mapped source fields where each field is weighed by by its fraction on the destination mesh.
     As mentioned above, this is given by the field ``mrg_fracnameN`` in ``FBFrac(comp_index_dst)``.
     If sum_with_weights is chose as the ``mrg_typeN`` value then ``mrg_fracnameN`` is also required as an argument.
   * sum_with_weights: do a cumulative sum of all the mapped source fields.
