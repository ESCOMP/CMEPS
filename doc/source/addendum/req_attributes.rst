.. _attributes:

==========================================
 CMEPS Application Independent Attributes
==========================================

The following attributes are obtained from the respective driver and
available to all components that the driver uses.  In the case of
UFS, the UFS driver ingests these attributes from the
``ufs.configure`` file.  In the case of CESM, the CESM driver ingests
these attributes from the ``nuopc.runconfig`` file.  The list of
attributes below are separated into application independent attributes
and at this time additional attributes required by CESM. There are no
UFS-specific attributes required by the UFS application.


General
-------

**coupling_mode** (required)

  The coupling_mode attribute determines which
  ``esmFlds_exchange_xxx_mod.F90`` and ``fd_xxx.yaml`` is used by
  CMEPS and is also leveraged in some of the custom calculations in
  the ``prep`` modules.

  The currently supported values for ``coupling_mode`` are ``cesm``, ``ufs.(frac,nfrac).(aoflux)``, and ``hafs``.

Scalar attributes
-----------------

**ScalarFieldCount**
  The maximum number of scalars that are going to be communicated
  between the mediator and a component.  Currently scalar values are
  put into a field bundle that only contains an undistributed
  dimension equal to the size of ``ScalarFieldCount`` and communicated
  between the component and the mediator on the `main task` of each
  component.

**ScalarFieldName** (required)
  This is the name of the scalar field bundle. By default it is ``cpl_scalars``.

**ScalarFieldIdxGridNX**, **ScalarFieldIdxGridNY** (required)
  The global number of longitude and latitude points. For unstructured grids::

    ScalarFieldIdxGridNY = 1
    ScalarFieldIdxGridNX = global size of mesh

  For cases where ``ScalarFieldIdxGridNY`` is not 1, this scalar data
  is needed by the mediator for the history output.

**ScalarFieldIdxNextSwCday** (optional)
 Send by the atmosphere component to specify the calendar day of its
 next short wave computation.  This is subsequently used by other
 components (e.g. cesm-land and sea-ice) in determining the zenith
 angle for its albedo calculation. It is also used in the mediator
 routine ``med_phases_ocnalb_mod.F90`` to determine the zenith angle
 in the ocean albedo calculation.

Mediator history and restart attributes
---------------------------------------

**history_option**, **history_n** (required)
  Determines the write frequency for a mediator history file (see :ref:`mediator history writes<history_writes>`).
**restart_option**, **restart_n** (required)
  Determines the write frequency for a mediator restart file (see :ref:`mediator restart writes<restart_writes>`).
**read_restart** (required)
  Determines if a mediator restart file is read in.
