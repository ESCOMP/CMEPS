.. _glossary:

========
Glossary
========

.. glossary::

   cap
      A NUOPC-compliant software layer between a component's native code and the
      mediator, presenting the component to the coupled system through import and
      export states. See :ref:`concepts`.

   CDEPS
      The Community Data Models for Earth Prediction Systems — data components
      (DATM, DOCN, …) that can replace an active component, for example by
      replaying saved coupler forcing (see the ``histaux`` discussion in
      :ref:`history`).

   component
      A model of one part of the Earth system — atmosphere (``atm``), ocean
      (``ocn``), sea ice (``ice``), land (``lnd``), river (``rof``), wave
      (``wav``) or land-ice (``glc``).

   coupling interval
      The period between successive couplings of a component, set by the
      ``<COMP>_NCPL`` settings and reflected in the run sequence. See
      :ref:`running-a-case`.

   driver
      The top-level NUOPC application component that harnesses the mediator and
      the component models and runs the coupled time loop. It is host-specific,
      unlike the shared mediator. See :ref:`concepts`.

   exchange grid
      An ESMF grid formed from the overlap of two component grids, on which
      atmosphere/ocean fluxes may be computed. See :ref:`aofluxes`.

   ``FBExp``
      The internal-state field bundle holding the fields the mediator exports to
      a component. See :ref:`code-organization`.

   ``FBImp``
      The internal-state matrix of field bundles ``FBImp(n,k)`` — component
      ``n``'s fields on grid ``k``. See :ref:`code-organization`.

   field dictionary
      The per-host YAML file (``fd_<host>.yaml``) defining the standard names,
      units and metadata of the coupled fields. See :ref:`fields`.

   fraction
      The fraction of a surface type (land, ocean, sea ice) in a grid cell, used
      as a weight in merges and mappings. See :ref:`surface-fractions`.

   ``histaux``
      Auxiliary coupler history streams that record the fields a source component
      sends into the mediator, for reuse as forcing via CDEPS. See
      :ref:`history`.

   mapping
      Interpolation of a field from a source grid to a destination grid through
      an ESMF routehandle. See :ref:`mapping`.

   mediator
      The CMEPS component that couples the component models — mapping, merging,
      computing fluxes and albedos, and maintaining fractions. Shared across host
      applications.

   merging
      Combining one or more mapped source fields into a single destination
      field. See :ref:`merging`.

   phase
      A named routine the mediator registers with NUOPC and that the run sequence
      invokes (for example ``med_phases_prep_ocn``). See :ref:`phase-lifecycle`.

   routehandle
      A precomputed sparse-matrix set of interpolation weights used to apply a
      mapping. See :ref:`mapping`.

   run configuration
      The host-specific file (``nuopc.runconfig`` for CESM/NorESM,
      ``ufs.configure`` for UFS) supplying the attributes the driver and mediator
      read. See :ref:`run-config`.

   run sequence
      The ordered recipe of component runs, mediator phases and field transfers,
      with their coupling intervals, that the driver ingests. See
      :ref:`run-sequence`.
