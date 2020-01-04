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

# The following is for RASM_OPTION1

def gen_runseq(case, coupling_times):

    rundir         = case.get_value("RUNDIR")
    caseroot       = case.get_value("CASEROOT")
    cpl_seq_option = case.get_value('CPL_SEQ_OPTION')

    driver_config = DriverConfig(case, coupling_times)
    run_atm, atm_to_med, med_to_atm, atm_cpl_time = driver_config['atm']
    run_ice, ice_to_med, med_to_ice, ice_cpl_time = driver_config['ice']
    run_ocn, ocn_to_med, med_to_ocn, ocn_cpl_time = driver_config['ocn']
    run_rof, rof_to_med, _         , rof_cpl_time = driver_config['rof']
    run_wav, wav_to_med, med_to_wav, wav_cpl_time = driver_config['wav']

    # error checking
    if cpl_seq_option != 'RASM':
        expect(False,"for C or G compset require that cpl_seq_option is RASM")
    if wav_cpl_time != atm_cpl_time:
        expect(False,"for C or G compset require that wav_cpl_time equal atm_cpl_time")
    if ice_cpl_time != atm_cpl_time:
        expect(False,"for C or G compset require that ice_cpl_time equals atm_cpl_time")
    if rof_cpl_time != atm_cpl_time:
        expect(False,"for C or G compset require that rof_cpl_time equals atm_cpl_time")
    if atm_cpl_time > ocn_cpl_time:
        expect(False,"for C or G compset require that ocn_cpl_time >= atm_cpl_time")

    with RunSeq(os.path.join(caseroot, "CaseDocs", "nuopc.runseq")) as runseq:

        runseq.enter_time_loop(ocn_cpl_time)

        runseq.add_action ("MED med_phases_prep_ocn_accum_avg" , med_to_ocn and atm_cpl_time < ocn_cpl_time)
        runseq.add_action ("MED -> OCN :remapMethod=redist"    , med_to_ocn and atm_cpl_time < ocn_cpl_time)

        runseq.enter_time_loop(atm_cpl_time, newtime=(atm_cpl_time < ocn_cpl_time))

        runseq.add_action ("MED med_phases_prep_ocn_map"        , med_to_ocn)
        runseq.add_action ("MED med_phases_aofluxes_run"        , med_to_ocn)
        runseq.add_action ("MED med_phases_prep_ocn_merge"      , med_to_ocn)
        runseq.add_action ("MED med_phases_prep_ocn_accum_fast" , med_to_ocn)
        runseq.add_action ("MED med_phases_ocnalb_run"          , med_to_ocn)
        runseq.add_action ("MED med_phases_prep_ocn_accum_avg"  , med_to_ocn and (atm_cpl_time == ocn_cpl_time))
        runseq.add_action ("MED -> OCN :remapMethod=redist"     , med_to_ocn and (atm_cpl_time == ocn_cpl_time))
        runseq.add_action ("MED med_phases_prep_ice"            , med_to_ice)
        runseq.add_action ("MED -> ICE :remapMethod=redist"     , med_to_ice)
        runseq.add_action ("MED med_phases_prep_wav"            , med_to_wav)
        runseq.add_action ("MED -> WAV :remapMethod=redist"     , med_to_wav)
        runseq.add_action ("ICE"                                , run_ice)
        runseq.add_action ("ROF"                                , run_rof)
        runseq.add_action ("ATM"                                , run_atm)
        runseq.add_action ("WAV"                                , run_wav)
        runseq.add_action ("ICE -> MED :remapMethod=redist"     , ice_to_med)
        runseq.add_action ("MED med_fraction_set"               , ice_to_med)
        runseq.add_action ("OCN"                                , run_ocn    and ocn_cpl_time == atm_cpl_time)
        runseq.add_action ("OCN -> MED :remapMethod=redist"     , ocn_to_med and ocn_cpl_time == atm_cpl_time)
        runseq.add_action ("ROF -> MED :remapMethod=redist"     , rof_to_med)
        runseq.add_action ("ATM -> MED :remapMethod=redist"     , atm_to_med)
        runseq.add_action ("WAV -> MED :remapMethod=redist"     , wav_to_med)
        runseq.add_action ("MED med_phases_history_write"       , atm_cpl_time == ocn_cpl_time)

        runseq.leave_time_loop(atm_cpl_time == rof_cpl_time)

        runseq.add_action ("OCN"                            , run_ocn    and ocn_cpl_time > atm_cpl_time)
        runseq.add_action ("OCN -> MED :remapMethod=redist" , ocn_to_med and ocn_cpl_time > atm_cpl_time)
        runseq.add_action ("MED med_phases_history_write"   , ocn_cpl_time > atm_cpl_time)

        runseq.leave_time_loop(True)

    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
