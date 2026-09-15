.. _diagnostics:

========================
Water and energy budgets
========================

One of the mediator's roles is to check that the coupled system conserves water
and energy. ``med_diag_mod.F90`` computes spatial and time averages of the
fluxed quantities needed for these **budget diagnostics**.

Sign convention
===============

Fluxes are accounted with a **positive-downward** sign convention and a fixed
component hierarchy — atm, glc, lnd, rof, ice, ocn. With this convention a flux
that one component gains is the same flux another loses, so the totals across the
coupled system should close. A budget that does not close points to a
conservation error somewhere in the coupling, which makes these diagnostics a
primary tool for validating conservation.

The budget phases
=================

The budget work is spread across several phases that the run sequence calls:

* per-component accumulation phases — ``med_phases_diag_atm``,
  ``med_phases_diag_lnd``, ``med_phases_diag_ocn``, and so on — accumulate each
  component's fluxed quantities as they pass through the mediator;
* ``med_phases_diag_accum`` accumulates those contributions over the averaging
  periods; and
* ``med_phases_diag_print`` prints the budget tables.

``med_diag_init`` and ``med_diag_zero`` set up and reset the accumulators.

The budget tables printed by ``med_phases_diag_print`` are written to files of
the form ``diags.log.<runid>.<yyyymmdd-hhmmss>``.

Reading a budget table
======================

For each averaging period the mediator prints a set of budgets — a **net area**
budget (m²/m²), a **net heat** budget (W/m²) and a **net water** budget
(kg/m²/s, scaled by 1e6). In each table the **columns** are the per-component
contributions (``atm``, ``lnd``, ``rof``, ``ocn``, ``ice nh``, ``ice sh``,
``glc``) plus a ``*SUM*`` column, and the **rows** are the individual flux terms
plus a ``*SUM*`` row.

Conservation shows up directly in the sums: because a flux one component gains is
another's loss, the ``*SUM*`` column (across components) for each term, and the
``*SUM*`` row, should be **essentially zero**. A term whose row sum is far from
zero points at a non-conserving exchange. For example, from a NorESM run
(``diags.log`` excerpt), the annual **heat** budget (W/m²):

.. code-block:: none

   (med_diag_print_summary) NET HEAT BUDGET (W/m2): period =   annual: date =  11060101     0
                              atm            lnd            rof            ocn         ice nh         ice sh            glc        *SUM*
       hfreeze          0.00000000     0.00000000     0.00000000     0.03292080    -0.00913038    -0.02379042     0.00000000     0.00000000
       hmelt            0.00000000     0.00000000     0.00000000    -0.89090241     0.32600589     0.56489652     0.00000000     0.00000000
       hnetsw        -164.81223433    42.73795062     0.00000000   120.82313101     0.62392301     0.63004087     0.00000000     0.00281117
       hlwdn         -336.69371959    87.02240716     0.00000000   239.28275576     4.67536186     5.71298315     0.00000000    -0.00021166
       hlwup          394.51485762  -107.98891375     0.00000000  -274.37211846    -5.48118805    -6.67262512     0.00000000     0.00001223
       hlatvap         85.44674862   -12.15707541     0.00000000   -73.10714372    -0.05853463    -0.12425713     0.00000000    -0.00026227
       hlatfus          0.82811998    -0.28266769     0.00000000    -0.36550691    -0.04546068    -0.13447928     0.00000000     0.00000542
       hiroff           0.00000000     0.00024877     0.00000000    -0.05026752     0.00000000     0.00000000     0.05001961     0.00000087
       hsen            20.86523625    -9.39739411     0.00000000   -11.46215568    -0.03417290     0.02776099     0.00000000    -0.00072545
              *SUM*     0.14900854    -0.06544441     0.00000000    -0.10928714    -0.00319588    -0.01947042     0.05001961     0.00163030

and the annual **water** budget (kg/m²/s, scaled by 1e6):

.. code-block:: none

   (med_diag_print_summary) NET WATER BUDGET (kg/m2s*1e6): period =   annual: date =  11060101     0
                              atm            lnd            rof            ocn         ice nh         ice sh            glc        *SUM*
       wfreeze          0.00000000     0.00000000     0.00000000    -0.08728167     0.02420702     0.06307465     0.00000000     0.00000000
       wmelt            0.00000000     0.00000000     0.00000000     0.55264305    -0.16932198    -0.38332107     0.00000000    -0.00000000
       wrain          -31.68859545     7.10225360     0.00000000    24.51612745     0.04130430     0.02895740     0.00000000     0.00004730
       wsnow           -2.48163015     0.84707129     0.00000000     1.09531588     0.13623219     0.40299455     0.00000000    -0.00001623
       wevap           34.15015152    -4.85457580     0.00000000   -29.23116502    -0.02073919    -0.04377629     0.00000000    -0.00010478
       weqsaltf         0.00000000     0.00000000     0.00000000     0.00072424     0.00102805    -0.00175229     0.00000000     0.00000000
       wrunoff          0.00000000    -2.84704267     0.00760572     2.83936982     0.00000000     0.00000000     0.00000000    -0.00006712
       wfrzrof          0.00000000    -0.00074550     0.00000000     0.15063685     0.00000000     0.00000000    -0.14989395    -0.00000260
              *SUM*    -0.02007408     0.24696091     0.00760572    -0.16362939     0.01271039     0.06617694    -0.14989395    -0.00014344

The rows are the individual flux terms — net shortwave, up/down longwave, latent
and sensible heat, freeze/melt and runoff for heat; rain, snow, evaporation,
runoff, freeze/melt for water. The small ``*SUM*`` values confirm the coupled
budgets close to round-off. The area budget has the same layout.

Enabling the budgets
====================

The budgets are enabled by the ``do_budgets`` attribute, and the ``budget_*``
settings select which averaging periods (instantaneous, daily, monthly, annual
and long-term) are printed. These are configured in :ref:`running-a-case`.
