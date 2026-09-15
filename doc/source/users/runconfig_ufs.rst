.. _run-config-ufs:

======================
NUOPC attributes (UFS)
======================

.. note::

   **This section is to be filled in.**

This page will describe the UFS run configuration, the counterpart to the
CESM/NorESM :ref:`nuopc.runconfig <run-config>` page.

For what a run configuration *is* and how it relates to the driver and the
mediator, see :ref:`concepts`. As on every host, each component obtains its
attributes **at initialization** through calls to ``NUOPC_CompAttributeGet`` —
the driver makes the run-configuration attributes available on every component,
and the mediator does not read the configuration file directly.

For UFS the run configuration is the ``ufs.configure`` file (the counterpart to
CESM/NorESM's ``nuopc.runconfig``). *The contents and attribute groups of
``ufs.configure`` remain to be documented here.*
