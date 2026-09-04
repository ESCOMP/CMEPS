.. _index:

###################
CMEPS documentation
###################

The Community Mediator for Earth Prediction Systems (CMEPS) is a
NUOPC-compliant mediator used to couple the components of an Earth
system model. It is used by, among others, NCAR's Community Earth System
Model (CESM), NorESM, and NOAA's UFS and HAFS coupled systems.

This documentation is organized into four parts:

* **Overview & Concepts** introduces what the mediator does and the ideas
  the rest of the documentation relies on. Start here.
* **User Guide** is for people configuring and running a coupled system
  that uses CMEPS. It is specific to **CESM and NorESM**; other applications
  (e.g. UFS, HAFS) provide their own user documentation.
* **Developer Guide** is for people modifying the mediator: its code
  structure, coupling phases, and how fields are exchanged, mapped and merged.
* **Reference** collects field-name and attribute tables and a glossary.

.. toctree::
   :maxdepth: 2
   :numbered:

   overview/index
   users/index
   developer/index
   reference/index
