Introduction
============

CMEPS is a NUOPC-compliant mediator which uses ESMF to couple earth grid components in a hub and spoke system.

As a mediator, CMEPS is responsible  for transferring field information from one
model component to another. This transfer  can require one or more operations on
the transferred  fields, including  mapping of  fields between  component grids,
merging  of fields  between different  components and  time-averaging of  fields
over varying coupling periods.



Components share information via import  and export states, which are containers
for  ESMF data  types  that wrap  native  model data.  The  states also  contain
metadata, which  includes physical  field names,  the underlying  grid structure
and coordinates,  and information on  the parallel decomposition of  the fields.
Note that while CMEPS itself is  a mesh based mediator, component models coupled
by the CMEPS mediator can be either grid or mesh based.



Each component model  using the CMEPS mediator is serviced  by a NUOPC-compliant
cap. The NUOPC cap  is a small software layer between  the underlying model code
and  the mediator.  Fields  for  which the  mediator  has  created a  connection
between model components are placed in either  the import or export state of the
component within  the NUOPC cap.  The information contained within  these states
is then passed into  native model arrays or structures for  use by the component
model.



Field  connections  made  by  the  CMEPS mediator  between  components  rely  on
matching of  standard field names. These  standard names are defined  in a field
dictionary.  Because CMEPS  is a  community mediator,  these standard  names are
specific to each application.


 .. todo::

   add description of what nuopc is, including advertise, etc.

Organization of the CMEPS mediator code
#######################################


When you check out the code you  will files, which can be separated into several
groups:

* totally generic components that carry  out the mediator functionality such as mapping, 
  merging, restarts and history writes. Included here is a a  "fraction" module that 
  determine  the fractions of different  source model components on every source 
  desinatation mesh.

* application specific  code that determines what fields  are exchanged between 
  components and how they are merged and mapped.

* prep phase modules  that carry out the mapping and merging  from one or more 
  source components to  the destination component.

=========================== ============================ ===========================
  Generic Code               Application Specific Code   Prep Phase Code
=========================== ============================ ===========================
med.F90                     esmFldsExchange_cesm_mod.F90 med_phases_prep_atm_mod.F90
esmFlds.F90                 esmFldsExchange_nems_mod.F90 med_phases_prep_ice_mod.F90
med_map_mod.F90             esmFldsExchange_hafs_mod.F90 med_phases_prep_ocn_mod.F90
med_merge_mod.F90           fd_cesm.yaml                 med_phases_prep_glc_mod.F90
med_frac_mod.F90            fd_nems.yaml                 med_phases_prep_lnd_mod.F90                          
med_internalstate_mod.F90   fd_hafs.yaml                 med_phases_prep_rof_mod.F90               
med_methods_mod.F90.                         
med_phases_aofluxes_mod.F90 
med_phases_ocnalb_mod.F90
med_phases_history_mod.F90
med_phases_restart_mod.F90
med_phases_profile_mod.F90
med_io_mod.F90
med_constants_mod.F90
med_kind_mod.F90
med_time_mod.F90
med_utils_mod.F90
=========================== ============================ ===========================

.. note:: Some modules, such as med_phases_prep_ocn.F90 and med_frac_mod.F90 contain application specific-code blocks.
