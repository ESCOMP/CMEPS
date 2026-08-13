.. _restart:

========
Restarts
========

``med_phases_restart_mod.F90`` writes and reads the mediator restart files that
let a run continue exactly where it stopped. ``med_phases_restart_write`` writes
the mediator's internal state — the field bundles, accumulators and fractions it
must carry across a restart — at the restart cadence, and
``med_phases_restart_read`` restores that state when a run is continued
(controlled by the ``read_restart`` attribute).

Timing profile
==============

``med_phases_profile_mod.F90`` writes the mediator's timing profile to the log
(``med_phases_profile``, with ``med_phases_profile_finalize`` at the end of the
run), giving a per-phase picture of where the mediator spends its time.
