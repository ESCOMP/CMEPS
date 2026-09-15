.. _running-a-case:

============================================
Configuring and running a case (CESM/NorESM)
============================================

.. note::

   This page describes the **CESM/NorESM** workflow, which uses the CIME
   case-control system. It covers only the parts that affect the CMEPS mediator;
   for the general CIME workflow (``create_newcase``, ``case.setup``,
   ``case.build``, ``case.submit``) refer to the CIME project at
   https://github.com/ESMCI/cime. The UFS workflow is covered separately in
   :ref:`Configuring and running a case (UFS) <running-a-case-ufs>`.

From a CMEPS point of view, two kinds of configuration matter: the **coupling
intervals** (which shape the :ref:`run sequence <run-sequence>`) and the
**mediator attributes and namelist settings** (which reach the mediator through
the :ref:`run configuration <run-config>`). Both are set in the case directory —
some through xml variables changed with ``xmlchange``, others through the driver
namelist file ``user_nl_cpl``.

Where settings live: xml versus ``user_nl_cpl``
===============================================

CESM/NorESM configuration is split between two mechanisms, and it matters which
one a given setting uses:

* **xml variables**, changed with ``./xmlchange VAR=value``. These control
  case-level structure — components, processor layout, coupling intervals, run
  length — and are used to *generate* the run sequence and parts of the run
  configuration.
* **the driver namelist**, edited in ``user_nl_cpl``. These set individual
  mediator namelist values in ``drv_in`` (and in turn the attributes CMEPS
  reads). Some driver values are defined by xml variables and **must** be
  changed with ``xmlchange`` instead of ``user_nl_cpl`` — ``user_nl_cpl`` itself
  documents which.

Setting the coupling intervals
==============================

Coupling intervals are set per component with the ``<COMP>_NCPL`` xml variables,
relative to a base period:

* ``NCPL_BASE_PERIOD`` — the base period the couplings are counted against (for
  example one day).
* ``ATM_NCPL``, ``OCN_NCPL``, ``ICE_NCPL``, ``LND_NCPL``, ``ROF_NCPL``,
  ``GLC_NCPL``, ``WAV_NCPL`` — the number of times each component is coupled per
  base period.

A component's coupling interval is ``NCPL_BASE_PERIOD / <COMP>_NCPL``. For
example, with a base period of one day, ``ATM_NCPL=48`` gives a 30-minute
atmosphere coupling interval and ``OCN_NCPL=24`` gives a one-hour ocean coupling
interval::

  ./xmlchange NCPL_BASE_PERIOD=day
  ./xmlchange ATM_NCPL=48
  ./xmlchange OCN_NCPL=24

These intervals are exactly what determine the nesting of the ``@<dt>`` time
loops in the generated run sequence — the faster-coupled components end up in the
inner loops. See :ref:`run-sequence` for how to read the result.

Controlling the mediator through ``user_nl_cpl``
================================================

Mediator behavior is adjusted by adding driver namelist settings to
``user_nl_cpl``. The most commonly used groups are:

* **Water and energy budgets** — ``do_budgets`` turns the budget diagnostics on,
  and the ``budget_*`` settings (``budget_inst``, ``budget_daily``,
  ``budget_month``, ``budget_ann``, ``budget_ltann``, ``budget_ltend``, …)
  select which averaging periods are printed.
* **Atmosphere/ocean fluxes** — ``aoflux_grid`` selects the grid on which the
  mediator computes air-sea fluxes.
* **History and auxiliary history** — the ``history_option_*`` / ``history_n_*``
  settings control mediator history output, and the ``histaux_*`` settings
  configure auxiliary history streams (for example ``histaux_atm2med_*``).
* **Debugging** — ``dbug_flag`` sets the mediator's diagnostic verbosity and
  ``check_for_nans`` enables NaN checking on coupled fields.
* **Special couplings** — for example ``ocn2glc_coupling`` and
  ``wav_coupling_to_cice`` enable particular exchanges when the relevant
  components are active.

.. note::

   ``coupling_mode`` — which selects the host exchange logic — is normally set
   by the compset rather than changed by hand. See :ref:`run-config`.

How it fits together
====================

Putting the pieces together for a CESM/NorESM case:

#. ``xmlchange`` sets the components, coupling intervals and run length.
#. CIME **generates** the run sequence from those intervals and writes the run
   configuration (``nuopc.runconfig``), folding in the driver and mediator
   attributes — including any values you set in ``user_nl_cpl``.
#. At run time the driver ingests the run configuration and run sequence, and
   the mediator obtains the attributes it needs and executes the ``MED`` phases in
   the order the run sequence specifies.

For the details of the two files the driver ingests, see :ref:`run-sequence` and
:ref:`run-config`.

Branch runs and the coupler restart pointer
===========================================

A **branch run** continues exactly from another case's restarts, so the mediator
reads its coupler restart. In addition to the usual branch settings, you must set
the **coupler restart pointer** explicitly to the dated ``rpointer.cpl`` file of
the reference case:

.. code-block:: bash

   ./xmlchange RUN_STARTDATE=1990-01-01
   ./xmlchange RUN_REFDATE=1990-01-01
   ./xmlchange DRV_RESTART_POINTER=rpointer.cpl.1990-01-01-00000

The coupler restart pointer files are **dated** — for example
``rpointer.cpl.1850-01-06-00000`` for a restart written at that timestamp.
Setting ``DRV_RESTART_POINTER`` explicitly for a branch run is **by design**:

* On a **continue** run (``CONTINUE_RUN=TRUE``), CIME already knows the exact
  restart it just wrote and automatically sets ``DRV_RESTART_POINTER`` to that
  dated pointer, so you never touch it.
* On a **branch** run, you are instead choosing *which* restart of a reference
  case to resume from — information CIME does not have on its own
  (``RUN_REFDATE`` gives the date, but the reference case may hold several
  restart snapshots). So you name the pointer yourself, consistent with the
  ``RUN_REFCASE`` / ``RUN_REFDATE`` you already set.

.. note::

   A **hybrid** run does not read the coupler restart at all, so the pointer does
   not apply there.
