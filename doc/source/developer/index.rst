.. _developer:

###############
Developer Guide
###############

.. important::

   **Everything in the Developer Guide assumes familiarity with ESMF and the
   NUOPC Layer.** The mediator is built entirely on ESMF/NUOPC data types and
   conventions — in particular ESMF ``States``, ``Fields``, ``FieldBundles`` and
   ``RouteHandles``, and the NUOPC concepts of components, caps, connectors and
   the run sequence. If these are unfamiliar, start with the `ESMF reference
   manual
   <https://earthsystemmodeling.org/docs/nightly/develop/ESMF_refdoc/node2.html#SECTION02065000000000000000>`_
   and the `NUOPC reference manual
   <https://earthsystemmodeling.org/docs/nightly/develop/NUOPC_refdoc/>`_
   before continuing.

Understanding this part is a prerequisite for any CMEPS customization or
modification. It covers the code organization in depth, the phase lifecycle,
how fields are declared and
exchanged for each host (``addfld`` / ``addmap`` / ``addmrg``), the mapping,
merging and fraction machinery, atmosphere/ocean flux and albedo calculation,
diagnostics/history/restart, and a walkthrough for adding a new coupled field.
It also develops the full conservation, area-correction and accumulation
derivations summarized in :ref:`concepts`.

Unlike the User Guide, the Developer Guide is **not** host-specific: it describes
the shared mediator code that serves all host applications. Where behavior
differs by host, it is in the ``esmFldsExchange_<host>`` exchange modules and the
``fd_<host>.yaml`` dictionaries, which are called out where relevant.

.. toctree::
   :maxdepth: 2

   code_organization
   field_exchange
   conservation
   phase_lifecycle
   mapping
   merging
   fractions
   aofluxes
   ocnalb
   diagnostics
   history
   restart
   add_field
