.. _field_naming_convention:

Application Specific Field Exchange Specification
=================================================

CMEPS contains two component specific files that determine:
 - the fields that are exchanged between components
 - how source fields are mapped to destination fields
 - how source fields are merged after mapping to destination fields


Field Naming Convention
=======================

The mediator variable names can be seen in the application specific YAML field dictionary. Currently, three
field dictionaries are supported::

  fd_cesm.yaml
  fd_nems.yaml

The CMEPS field name convention in these YAML files is independent of the model components.
The convention differentiates between variables that are state fields versus flux fields.

State variables have a prefix that always start with an ``S`` followned by a two character string::

  state-prefix
    first 3 characters: Sx_, Sa_, Si_, Sl_, So_
    one letter indices: x,a,l,i,o,g,r
    x => mediator (mapping, merging)
    a => atmosphere
    l => land
    i => sea-ice
    o => ocean
    g => land-ice
    r => river
    w => wave

  state-name
    what follows state prefix

As an example, ``Sx_t`` is the merged surface temperature from land, ice and ocean sent to the atmopshere for CESM.

Flux variables that specifies both source and destination components and have a 5 character prefix::

  flux-prefix
    first 5 characters: Flmn_
    lm => between components l and m
    n  => computed by component n
    example: Fioi => ice-ocn flux computed by ice
    example: Fall => atm-lnd flux computed by lnd

   flux-name
     what follows flux-prefix
