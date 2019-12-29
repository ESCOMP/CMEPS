#!/usr/bin/env python

import os, shutil, sys
from CIME.utils import expect
from gen_runseq import RunSeq

_CIMEROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir,os.pardir,os.pardir,os.pardir))
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *

#pylint:disable=undefined-variable
logger = logging.getLogger(__name__)

def runseq(case, coupling_times):

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")
    comp_rof  = case.get_value("COMP_ROF")
    comp_glc  = case.get_value("COMP_GLC")

    rof_cpl_dt = coupling_times["rof_cpl_dt"]
    atm_cpl_dt = coupling_times["atm_cpl_dt"]
    glc_cpl_dt = coupling_times["glc_cpl_dt"]

    # NOTE: in the mediator the glc_avg_period will be set as an alarm
    # on the mediator clock - when this alarm rings - the averaging will be done AND an
    # attribute will be set on the glc export state from the mediator
    # saying that the data coming to glc is valid
    glc_couplint_time = glc_cpl_dt

    # ASSUME that will not run I compset with data runoff - either stub or prognostic
    if (comp_rof == 'srof'):
        couple_rof = False
    else:
        couple_rof = True

    cism_evolve = False
    if (comp_glc == 'sglc'):
        glc_to_med = False
        med_to_glc = False
    elif (comp_glc == 'cism'):
        cism_evolve = case.get_value("CISM_EVOLVE")
        glc_to_med = True
        if cism_evolve:
            med_to_glc = True
        else:
            med_to_glc = False
    else:
        expect(False, "Invalid comp_glc %s " %compglc)
        
    # If CISM is not evolving only get data back from cism at the initial time
    # However will still need to call the exchange at the end if the stop_option
    # is nsteps or days - or otherwise just every ndays
    # Note that nsteps is the minimum component coupling time
    if (comp_glc == 'cism'):
        glc_coupling_time = glc_cpl_dt
        if not cism_evolve:
            stop_option = case.get_value('STOP_OPTION')
            stop_n = case.get_value('STOP_N')
            if stop_option == 'nsteps':
                glc_coupling_time = stop_n * atm_cpl_dt
            elif stop_option == 'ndays':
                glc_coupling_time = stop_n * 86400
            else:
                glc_coupling_time = 86400
    else:
        glc_coupling_time = 0

    # ----------------------------------
    # Write nuopc.runseq
    # ----------------------------------

    with RunSeq(os.path.join(caseroot, "CaseDocs", "nuopc.runseq")) as runseq:

        runseq.enter_time_loop(glc_coupling_time, active=med_to_glc, if_newtime=(glc_to_med or med_to_glc)) 

        runseq.add_action("MED med_phases_prep_glc_avg"    , med_to_glc)
        runseq.add_action("MED -> GLC :remapMethod=redist" , med_to_glc)
        runseq.add_action("GLC"                            , med_to_glc)
        runseq.add_action("GLC -> MED :remapMethod=redist" , glc_to_med)

        runseq.enter_time_loop(rof_cpl_dt, if_newtime=couple_rof) 

        runseq.add_action("MED med_phases_prep_rof_avg"    , couple_rof)
        runseq.add_action("MED -> ROF :remapMethod=redist" , couple_rof)
        runseq.add_action("ROF"                            , couple_rof)

        runseq.enter_time_loop(atm_cpl_dt, if_newtime=((not couple_rof) or (atm_cpl_dt < rof_cpl_dt)))

        runseq.add_action("MED med_phases_prep_lnd"            , True)
        runseq.add_action("MED -> LND :remapMethod=redist"     , True)
        runseq.add_action("LND"                                , True)
        runseq.add_action("LND -> MED :remapMethod=redist"     , True)
        runseq.add_action("MED med_phases_prep_rof_accum"      , couple_rof )
        runseq.add_action("MED med_phases_prep_glc_accum"      , med_to_glc )
        runseq.add_action("ATM"                                , True)
        runseq.add_action("ATM -> MED :remapMethod=redist"     , True)
       #runseq.add_action("MED med_phases_history_write"       , atm_cpl_dt < rof_cpl_dt)

        runseq.leave_time_loop(couple_rof and (atm_cpl_dt < rof_cpl_dt))

        runseq.add_action("ROF -> MED :remapMethod=redist"     , couple_rof)
        runseq.add_action("MED med_phases_history_write"       , couple_rof)

        runseq.leave_time_loop(couple_rof and (atm_cpl_dt < rof_cpl_dt))

    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
    

