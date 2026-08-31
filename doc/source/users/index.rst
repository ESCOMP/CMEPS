.. _users:

##########
User Guide
##########

.. important::

   **The driver, the run-configuration file, and the run sequence file are
   host-specific.** Where applicable, differences between CESM/NorESM and UFS
   are called out in host-specific sections.

   The mediator itself is shared across applications, so the :ref:`Overview
   <overview>` and :ref:`Developer Guide <developer>` are not host-specific.

A coupled run is driven by **two key inputs that the NUOPC driver ingests** at
start-up:

* the **run sequence** — the ordered recipe of component runs, mediator phases
  and field transfers, together with the coupling intervals (:ref:`run-sequence`); and
* the **run configuration** — the ESMF attributes the driver reads at
  initialization and makes available to the mediator and the other components
  (:ref:`run-config`).

The pages below describe each of these two inputs in turn, and then show how
CESM/NorESM and UFS cases are configured to produce them.

.. toctree::
   :maxdepth: 2

   run_sequence
   runconfig
   runconfig_ufs
   running
   running_ufs
   fields
