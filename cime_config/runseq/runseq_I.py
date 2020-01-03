#!/usr/bin/env python

import os, shutil, sys
from CIME.utils import expect
from gen_runseq import RunSeq
from driver_config import DriverConfig

_CIMEROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir,os.pardir,os.pardir,os.pardir))
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *

#pylint:disable=undefined-variable
logger = logging.getLogger(__name__)

def runseq(case, coupling_times):

    rundir   = case.get_value("RUNDIR")
    caseroot = case.get_value("CASEROOT")

    driver_config = DriverConfig(case, coupling_times)
    atm_to_med, med_to_atm, atm_coupling_time = driver_config['atm']
    ice_to_med, med_to_ice, ice_coupling_time = driver_config['ice']
    glc_to_med, med_to_glc, glc_coupling_time = driver_config['glc']
    lnd_to_med, med_to_lnd, lnd_coupling_time = driver_config['lnd']
    ocn_to_med, med_to_ocn, ocn_coupling_time = driver_config['ocn']
    rof_to_med, med_to_rof, rof_coupling_time = driver_config['rof']
    wav_to_med, med_to_wav, wav_coupling_time = driver_config['wav']

    with RunSeq(os.path.join(caseroot, "CaseDocs", "nuopc.runseq")) as runseq:

        runseq.enter_time_loop(glc_coupling_time, active=med_to_glc, newtime=(glc_to_med or med_to_glc)) 

        runseq.add_action("MED med_phases_prep_glc_avg"    , med_to_glc)
        runseq.add_action("MED -> GLC :remapMethod=redist" , med_to_glc)
        runseq.add_action("GLC"                            , med_to_glc)
        runseq.add_action("GLC -> MED :remapMethod=redist" , glc_to_med)

        runseq.enter_time_loop(rof_coupling_time, newtime=rof_to_med) 

        runseq.add_action("MED med_phases_prep_rof_avg"    , rof_to_med)
        runseq.add_action("MED -> ROF :remapMethod=redist" , rof_to_med)
        runseq.add_action("ROF"                            , rof_to_med)

        runseq.enter_time_loop(atm_coupling_time, newtime=((not rof_to_med) or (atm_coupling_time < rof_coupling_time)))

        runseq.add_action("MED med_phases_prep_lnd"        , med_to_lnd)
        runseq.add_action("MED -> LND :remapMethod=redist" , med_to_lnd)
        runseq.add_action("LND"                            , lnd_to_med)
        runseq.add_action("LND -> MED :remapMethod=redist" , lnd_to_med)
        runseq.add_action("MED med_phases_prep_rof_accum"  , rof_to_med)
        runseq.add_action("MED med_phases_prep_glc_accum"  , med_to_glc)
        runseq.add_action("ATM"                            , atm_to_med)
        runseq.add_action("ATM -> MED :remapMethod=redist" , atm_to_med)
       #runseq.add_action("MED med_phases_history_write"   , atm_coupling_time < rof_coupling_time)

        runseq.leave_time_loop(rof_to_med and (atm_coupling_time < rof_coupling_time))

        runseq.add_action("ROF -> MED :remapMethod=redist" , rof_to_med)
        runseq.add_action("MED med_phases_history_write"   , rof_to_med)

        runseq.leave_time_loop(rof_to_med and (atm_coupling_time < rof_coupling_time))

    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
    

