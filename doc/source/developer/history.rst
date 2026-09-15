.. _history:

==============
History output
==============

``med_phases_history_mod.F90`` writes the mediator's history files — snapshots of
the coupled fields, on their grids, at a requested frequency. It is the coupler's
own diagnostic output, separate from each component's history.

Instantaneous and time-averaged streams
========================================

History output comes in two kinds, maintained as separate streams:

* **instantaneous** — the field values at the output time; and
* **time-averaged** — the fields averaged over the output interval (the mediator
  accumulates them between writes and divides at write time).

Each stream has its own output frequency, set by a ``freq_option`` / ``freq_n``
pair (for example "every ``nmonths``", "every ``nsteps``"). These come from the
``history_option_*`` / ``history_n_*`` settings (see :ref:`running-a-case`).

The write phases
================

The run sequence calls the appropriate write phase:

* ``med_phases_history_write`` — instantaneous output of all variables;
* ``med_phases_history_write_comp`` — instantaneous, averaged and auxiliary
  output for a component (``..._comp_inst`` and ``..._comp_avg`` do the actual
  writes);
* ``med_phases_history_write_med`` — the mediator-computed fields (the
  atmosphere/ocean fluxes and ocean albedos);
* ``med_phases_history_write_data2glc`` — the (normally yearly) average of the
  fields sent to the land-ice component.

``med_phases_history_init`` sets up the streams and their alarms.

File names
==========

The mediator history files are named after the case with a coupler (``.cpl``)
tag and the output timestamp, ending in ``.nc`` — for example
``<case>.cpl[.<inst_tag>].h*.<yyyy-mm-dd-sssss>.nc``.

Auxiliary history (``histaux``) files
=====================================

A ``histaux`` stream records, at the frequency they are exchanged, the fields
that a source component sends into the mediator. For example,
``histaux_atm2med_*`` saves the atmosphere-to-mediator fields each coupling
period.

Because these files capture the exact forcing a component provided, they can be
**replayed later**: a CDEPS data model reads them back as an input stream and
uses them to drive a *target* component, without rerunning the original source
component. This is a cheap way to **spin up** a slow component.

A common example is the CESM/CTSM land biogeochemistry, which needs many years of
atmospheric forcing to equilibrate. Rather than rerunning CAM each time, one run
with an active CAM saves the atmosphere-to-mediator fields as ``histaux`` files.
CDEPS/DATM (a data atmosphere) then repeatedly cycles over those saved years — for
instance a 20-year block — reading the ``histaux`` files as its forcing to drive
CTSM on its own.

Auxiliary history files carry an ``.hx.`` tag in their name. For example, the
land-ice files written by ``med_phases_history_write_data2glc`` name the source
of the mapping: ``<case>.cpl.hx.lnd2glc.<time>.nc`` and
``<case>.cpl.hx.ocn2glc.<time>.nc``.

Per-slot settings
-----------------

Each ``histaux`` stream is defined by one or more numbered *file slots*
(``histaux_<stream>_file<N>_*``). Each slot has these settings:

* ``enabled`` — whether the slot is written.
* ``flds`` — the ``:``-separated list of fields to write (or ``all``).
* ``history_option`` / ``history_n`` — the output frequency.
* ``doavg`` — instantaneous (``.false.``) or time-averaged (``.true.``).
* ``auxname`` — a tag placed in the auxiliary file name.
* ``ntperfile`` — the number of time samples per file.

Default histaux streams
-----------------------

These settings are available as **NUOPC attributes**, obtained during
initialization in ``med_phases_history_mod.F90``. The default file slots are
summarized below; all can be overridden in ``user_nl_cpl``.

.. list-table::
   :header-rows: 1
   :widths: 12 16 10 8 8 46

   * - Import type
     - ``auxname``
     - Frequency
     - Type
     - Samples / file
     - Default fields written out for NorESM/CESM
   * - ``atm2med``
     - ``atm.1h.inst``
     - 1 h
     - inst
     - 24
     - ``Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf``
   * - ``atm2med``
     - ``atm.1h.avrg``
     - 1 h
     - avg
     - 24
     - ``Sa_u:Sa_v``
   * - ``atm2med``
     - ``atm.3hprec.avrg``
     - 3 h
     - avg
     - 8
     - ``Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl``
   * - ``atm2med``
     - ``atm.3h.avrg``
     - 3 h
     - avg
     - 8
     - ``Sa_z:Sa_topo:Sa_u:Sa_v:Sa_tbot:Sa_ptem:Sa_shum:Sa_dens:Sa_pbot:Sa_pslv:Faxa_lwdn:Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl:Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf:Sa_co2prog``
   * - ``atm2med``
     - ``atm.24h.avrg``
     - 1 day
     - avg
     - 1
     - ``Faxa_bcph:Faxa_ocph:Faxa_dstwet:Faxa_dstdry:Faxa_ndep:Sa_co2diag``
   * - ``lnd2med``
     - ``lnd.ncpl.inst``
     - coupling step
     - inst
     - 1
     - ``all``
   * - ``lnd2med``
     - ``lnd2rof.24hr.avg``
     - 1 day
     - avg
     - 365
     - ``Flrl_rofsur:Flrl_rofi:Flrl_rofgwl:Flrl_rofsub``
   * - ``ocn2med``
     - ``ocn.24h.avg``
     - 1 day
     - avg
     - 30
     - ``So_bldepth:So_t:So_u:So_v``
   * - ``rof2med``
     - ``rof.24h.avrg``
     - 3 h
     - avg
     - 2
     - ``all``
   * - ``wav2med``
     - ``wav.24h.avg``
     - 1 day
     - inst
     - 30
     - *(none set)*

The ``histaux_*`` streams are configured through ``user_nl_cpl`` (see
:ref:`running-a-case`).
