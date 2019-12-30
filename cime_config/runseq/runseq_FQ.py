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

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")

    driver_config = DriverConfig(case, coupling_times)
    atm_to_med, med_to_atm, atm_coupling_time = driver_config['atm']
    ice_to_med, med_to_ice, ice_coupling_time = driver_config['ice']
    glc_to_med, med_to_glc, glc_coupling_time = driver_config['glc']
    lnd_to_med, med_to_lnd, lnd_coupling_time = driver_config['lnd']
    ocn_to_med, med_to_ocn, ocn_coupling_time = driver_config['ocn']
    rof_to_med, med_to_rof, rof_coupling_time = driver_config['rof']
    wav_to_med, med_to_wav, wav_coupling_time = driver_config['wav']

    with RunSeq(os.path.join(caseroot, "CaseDocs", "nuopc.runseq")) as runseq:

        runseq.enter_time_loop(glc_coupling_time, active=med_to_glc, if_newtime=(glc_to_med or med_to_glc))

        runseq.add_action("MED med_phases_prep_glc_avg"    , med_to_glc)
        runseq.add_action("MED -> GLC :remapMethod=redist" , med_to_glc)
        runseq.add_action("GLC"                            , med_to_glc)
        runseq.add_action("GLC -> MED :remapMethod=redist" , glc_to_med)

        runseq.enter_time_loop(rof_coupling_time, if_newtime=rof_to_med)

        runseq.add_action("MED med_phases_prep_rof_avg"    , rof_to_med)
        runseq.add_action("MED -> ROF :remapMethod=redist" , rof_to_med)
        runseq.add_action("ROF"                            , rof_to_med)

        runseq.enter_time_loop(atm_coupling_time, if_newtime=((not rof_to_med) or (atm_coupling_time < rof_coupling_time)))

        runseq.add_action("MED med_phases_prep_ocn_avg"        , med_to_ocn)
        runseq.add_action("MED -> OCN :remapMethod=redist"     , med_to_ocn)

        runseq.add_action("MED med_phases_prep_lnd"            , lnd_to_med)
        runseq.add_action("MED -> LND :remapMethod=redist"     , lnd_to_med)

        runseq.add_action("MED med_phases_prep_ice"            , ice_to_med)
        runseq.add_action("MED -> ICE :remapMethod=redist"     , ice_to_med)

        runseq.add_action("ICE"                                , ice_to_med)
        runseq.add_action("LND"                                , lnd_to_med)
        runseq.add_action("OCN"                                , ocn_to_med)

        runseq.add_action("OCN -> MED :remapMethod=redist"     , ocn_to_med)

        runseq.add_action("ICE -> MED :remapMethod=redist"     , ice_to_med)
        runseq.add_action("MED med_fraction_set"               , ice_to_med or ocn_to_med)
        runseq.add_action("MED med_phases_aofluxes_run"        , ocn_to_med)
        runseq.add_action("MED med_phases_ocnalb_run"          , ocn_to_med)

        runseq.add_action("LND -> MED :remapMethod=redist"     , lnd_to_med)
        runseq.add_action("MED med_phases_prep_rof_accum"      , med_to_rof)
        runseq.add_action("MED med_phases_prep_glc_accum"      , med_to_glc)

        runseq.add_action("MED med_phases_prep_atm"            , med_to_atm)
        runseq.add_action("MED -> ATM :remapMethod=redist"     , med_to_atm)
        runseq.add_action("ATM"                                , atm_to_med)
        runseq.add_action("ATM -> MED :remapMethod=redist"     , atm_to_med)

        runseq.add_action("MED med_phases_prep_ocn_map"        , med_to_ocn)
        runseq.add_action("MED med_phases_prep_ocn_merge"      , med_to_ocn)
        runseq.add_action("MED med_phases_prep_ocn_accum_fast" , med_to_ocn)
        runseq.add_action("MED med_phases_history_write"       , not rof_to_med)

        runseq.leave_time_loop(rof_to_med and (atm_coupling_time < rof_coupling_time))

        runseq.add_action("ROF -> MED :remapMethod=redist"     , rof_to_med)
        runseq.add_action("MED med_phases_history_write"       , rof_to_med)

        runseq.leave_time_loop(rof_to_med and (atm_coupling_time < rof_coupling_time))

    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
