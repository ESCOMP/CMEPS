.. _generic_modules:

=========================
 CMEPS `generic` modules
=========================

The following describes in some detail the CMEPS modules that are not
application specific and provide general functionality.

**med.F90**

  This module is initializes the CMEPS mediator functionality by performing the following functions:

  * adding a namespace (i.e. nested state) for each import and export
    component state in the mediator's InternalState

  * initializing the mediator component specific fields via a call to
    ``esmFldsExchange_xxx_`` (where currently xxx can be ``cesm``, ``nems`` or ``hafs``).

  * determining which components are present

  * advertising the import and export mediator fields

  * creating import (``FBImp``), export (``FBExp``) and accumulation (``FBExpAccum``) field bundles

  * initializing the mediatory route handles and field bundles needed for normalization

  * initializing component ``FBFrac`` field bundles

  * reading mediator restarts

  * optionally carrying out initializations for atmosphere/ocean flux
    calculations and ocean albedo calculations (these are needed by CESM)

  * carrying out the NUOPC data initialization via the ``DataInitialize`` routine.

    .. note:: After the first DataInitialize() of CMEPS returns,
	          NUOPC will note that its InitializeDataComplete is not yet true. The
	          NUOPC Driver will then execute the Run() phase of all of the Connectors that
	          fit the xxx-TO-MED pattern. After that, it will call CMEPS
	          DataInitialize() again. Note that the time stamps are only set
	          when the Run() phase of all the connectors are run and the
	          Connectors Run() phase is called before the second call of the
	          CMEPS DataInitialize phase. As a result, CMEPS will see the
	          correct timestamps, which also indicates that the actual data has
	          been transferred reliably, and CMEPS can safely use it.

**med_map_mod.F90**

  This module creates the required route handles that are needed for
  the model run.  The route handles are stored in the multi-dimensional array
  ``RH(ncomps,ncomps,nmappers)`` in the module ``med_internal_state_mod.F90``.

  ``nmappers`` is the total number of mapping types that CMEPS  supports (currently 8).
  These are described in :ref:`mapping types<addmap>`.

  ``ncomps,ncomps`` corresponds to the source and destination component indices.

  As an example ``RH(compatm,compocn,mapbilnr)`` is the atm->ocn bilinear route handle.

  **med_map_mod.F90** also initializes additional field bundles that
  are needed for mapping fractional normalization (see the
  :ref:`mapping normalization <normalization>`). Normalization is
  normally done using the relevant field from ``FBFrac(:)``.

  The default call to carry out mediator mapping is done in the
  :ref:`prep_modules<prep_modules>` by calling
  ``med_map_FB_Regrid_Norm``.  Mapping is done by using the
  ``fldListFr(:)`` data that was initialized in the
  ``esmFldsExchange_xxxx_mod.F90`` calls to ``addmap``.

**med_merge_mod.F90**

   This module carries out merging of one or more mapped source fields
   to the target destination field (see :ref:`merging
   types<addmrg>`). The :ref:`prep_modules<prep_modules>` carry out
   merging via the call to ``med_merge_auto`` Merging is done by using
   the ``fldListTo(:)`` data that was initialized in the
   ``esmFldsExchange_xxx_mod.F90`` calls to ``addmrg``.

**med_io_mod.F90**

   CMEPS uses the PIO2 parallel library to carry out all IO. PIO
   provides a netCDF-like API, and allows users to designate some
   subset of processors to perform IO. Computational code calls
   netCDF-like functions to read and write data, and PIO uses the IO
   processors to perform all necessary IO. This module contains
   wrapper layers to PIO for writing and reading mediator restart
   files and for writing mediator history files.

.. _history_writes:

**med_phases_history_mod.F90**

   This module writes mediator history files. The freqency of CMEPS
   history writes is controlled via the NUOPC attributes
   ``history_option`` and ``history_n``. These attributes control
   instantaneous mediator history output as follows:

   ==============  =============================================================
   history_option          description
   ==============  =============================================================
   none		   do not write any history files
   never	   do not write any history files
   nsteps	   write files every ``history_n`` mediator coupling intervals
   nseconds	   write files every ``history_n`` seconds
   nminutes	   write files every ``history_n`` minutes
   nhours	   write files every ``history_n`` hours
   ndays	   write files every ``history_n`` days
   nmonths	   write files every ``history_n`` months
   nyears	   write files every ``history_n`` years
   monthly	   write files on the month boundary
   yearly	   write files on the year boundary
   ==============  =============================================================

   .. note:: It is assumed that the NUOPC attributes ``history_option`` and ``history_n``
	     are obtained by the model driver and passed down to the mediator.

.. _restart_writes:

**med_phases_restart_mod.F90**

   This module reads and writes mediator restart files. The freqency of CMEPS
   restart writes is controlled via the NUOPC attributes
   ``restart_option`` and ``restart_n``. These attributes control
   instantaneous mediator history output as follows:

   ==============  =============================================================
   restart_option          description
   ==============  =============================================================
   none		   do not write any restart files
   never	   do not write any restart files
   nsteps	   write files every ``restart_n`` mediator coupling intervals
   nseconds	   write files every ``restart_n`` seconds
   nminutes	   write files every ``restart_n`` minutes
   nhours	   write files every ``restart_n`` hours
   ndays	   write files every ``restart_n`` days
   nmonths	   write files every ``restart_n`` months
   nyears	   write files every ``restart_n`` years
   monthly	   write files on the month boundary
   yearly	   write files on the year boundary
   ==============  =============================================================

   .. note:: It is assumed that the NUOPC attributes ``restart_option`` and ``restart_n``
	     are obtained by the model driver and passed down to the mediator.
