#!/usr/bin/env python3

import os, shutil, sys
from CIME.utils import expect
from gen_runseq import RunSeq
from driver_config import DriverConfig

_CIMEROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir,os.pardir,os.pardir,os.pardir))
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *

#pylint:disable=undefined-variable
logger = logging.getLogger(__name__)

def gen_runseq(case, coupling_times):

    rundir         = case.get_value("RUNDIR")
    caseroot       = case.get_value("CASEROOT")

    driver_config = DriverConfig(case, coupling_times)
    run_glc, med_to_glc, glc_cpl_time = driver_config['glc']
    run_lnd,  _        , lnd_cpl_time = driver_config['lnd']

    if lnd_cpl_time != glc_cpl_time:
        expect(False,"for TG compset require that lnd_cpl_time equal glc_cpl_time")

    with RunSeq(os.path.join(caseroot, "CaseDocs", "nuopc.runseq")) as runseq:

        runseq.enter_time_loop(glc_cpl_time, True)

        runseq.add_action ("LND"                            , run_lnd)
        runseq.add_action ("LND -> MED :remapMethod=redist" , run_lnd)
        runseq.add_action ("MED med_phases_post_lnd"        , run_lnd)
        runseq.add_action ("MED med_phases_prep_glc"        , med_to_glc)
        runseq.add_action ("MED -> GLC :remapMethod=redist" , med_to_glc)
        runseq.add_action ("GLC"                            , run_glc)
        # Need to do GLC -> MED even if not running GLC; otherwise, we get a
        # failure in InitializeRealize ("Object being used before creation")
        runseq.add_action ("GLC -> MED :remapMethod=redist" , med_to_glc)
        runseq.add_action ("MED med_phases_history_write"   , True)

        runseq.leave_time_loop(True)

    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
