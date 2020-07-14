Introduction
============

CMEPS is a mediator for use in hub and spoke system.

When you check out the code you will get the following mediator files
that you obtain can be separated into several groups:

1. application specific code that determines what fields are exchanged
between components and how they are merged and mapped

2 "prep" phase modules that carry out the mapping and merging from one
or more source components to the destination component (the rules
for these operations are determined in 1. above)

3. a "fraction" module that determine the fractions of different
source model components on every source desintation mesh

4. totally generic components that carry out the mediator
functionality such as mapping, merging, restarts and history writes.

.. todo::

   add description of what nuopc is, including advertise, etc.

CMEPS application specfic code
----------------------------------------

CMEPS is a community mediator in that it supports multiple coupling systems.
CMEPS contains two coupled model specific files that determine:
- the fields that are exchanged between components
- how source fields are mapped to destination fields
- how source fields are merged after mapping to destination fields

This occurs via the following files:

- ``esmFldsExchange_cesm_mod.F90``
- ``esmFldsExchange_nems_mod.F90``
- ``esmFldsExchange_hafs_mod.F90``

CMEPS advertises **all possible fields** that can be imported to and
exported by the mediator for the target coupled system. Not all of
these fields will be connected to the various components. The
connections will be determined by what the components advertise in
their respective advtise phase.

The mediator variable names can be seen in the application specific
YAML field dictionary. Currently, three field dictionaries are
supported for the target coupled model applications:

- ``fd_cesm.yaml`` for CESM
- ``fd_nems.yaml`` for UFS-S2S
- ``fd_hafs.yaml`` for UFS-HAFS

Details of the naming conventions and API's of this file can be found in the description of the :ref:`exchange of fields in CMEPS<api-for-esmflds>`.

CMEPS generic code
------------------

The following modules comprise the generic code::

    ESMFConvenienceMacros.h
    ESMFVersionDefine.h
    Makefile
    esmFlds.F90
    med_constants_mod.F90
    med_phases_history_mod.F90
    med_io_mod.F90
    med_internalstate_mod.F90
    med_kind_mod.F90
    med_map_mod.F90
    med_merge_mod.F90
    med_methods_mod.F90
    med_phases_ocnalb_mod.F90
    med_phases_aofluxes_mod.F90
    med_phases_profile_mod.F90
    med_time_mod.F90
    med_utils_mod.F90
    med_phases_restart_mod.F90

CMEPS generic code
------------------

    med.F90
    med_fraction_mod.F90
    med_phases_prep_atm_mod.F90
    med_phases_prep_glc_mod.F90
    med_phases_prep_ice_mod.F90
    med_phases_prep_lnd_mod.F90
    med_phases_prep_ocn_mod.F90
    med_phases_prep_rof_mod.F90
    med_phases_prep_wav_mod.F90
