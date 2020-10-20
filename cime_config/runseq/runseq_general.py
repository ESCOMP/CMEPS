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

def gen_runseq(case, coupling_times):

    rundir         = case.get_value("RUNDIR")
    caseroot       = case.get_value("CASEROOT")
    cpl_seq_option = case.get_value('CPL_SEQ_OPTION')
    cpl_add_aoflux = case.get_value('ADD_AOFLUX_TO_RUNSEQ')
    coupling_mode  = case.get_value('COUPLING_MODE')
    diag_mode      = case.get_value('BUDGETS')

    # It is assumed that if a component will be run it will send information to the mediator
    # so the flags run_xxx and xxx_to_med are redundant

    driver_config = DriverConfig(case, coupling_times)
    run_atm, med_to_atm, atm_cpl_time = driver_config['atm']
    run_ice, med_to_ice, ice_cpl_time = driver_config['ice']
    run_glc, med_to_glc, glc_cpl_time = driver_config['glc']
    run_lnd, med_to_lnd, lnd_cpl_time = driver_config['lnd']
    run_ocn, med_to_ocn, ocn_cpl_time = driver_config['ocn']
    run_rof, med_to_rof, rof_cpl_time = driver_config['rof']
    run_wav, med_to_wav, wav_cpl_time = driver_config['wav']

    # Note: assume that atm_cpl_dt, lnd_cpl_dt, ice_cpl_dt and wav_cpl_dt are the same

    if lnd_cpl_time != atm_cpl_time:
        expect(False, "assume that lnd_cpl_time is equal to atm_cpl_time")
    if ice_cpl_time != atm_cpl_time:
        expect(False, "assume that ice_cpl_time is equal to atm_cpl_time")
    if wav_cpl_time != atm_cpl_time:
        expect(False, "assume that wav_cpl_time is equal to atm_cpl_time")
    if rof_cpl_time < ocn_cpl_time:
        expect(False, "assume that rof_cpl_time is always greater than or equal to ocn_cpl_time")

    rof_outer_loop = run_rof and rof_cpl_time > atm_cpl_time
    ocn_outer_loop = run_ocn and ocn_cpl_time > atm_cpl_time

    inner_loop = ((atm_cpl_time < ocn_cpl_time) or
                  (atm_cpl_time < rof_cpl_time) or
                  (run_glc and atm_cpl_time < glc_cpl_time) or
                  atm_cpl_time == ocn_cpl_time)

    with RunSeq(os.path.join(caseroot, "CaseDocs", "nuopc.runseq")) as runseq:

        #------------------
        runseq.enter_time_loop(glc_cpl_time, newtime=run_glc, active=med_to_glc)
        #------------------

        runseq.add_action("MED med_phases_prep_glc_avg"    , med_to_glc)
        runseq.add_action("MED -> GLC :remapMethod=redist" , med_to_glc)
        runseq.add_action("GLC"                            , run_glc and med_to_glc)
        runseq.add_action("GLC -> MED :remapMethod=redist" , run_glc)

        #------------------
        runseq.enter_time_loop(rof_cpl_time, newtime=rof_outer_loop)
        #------------------

        runseq.add_action("MED med_phases_prep_rof_avg"    , med_to_rof and rof_outer_loop)
        runseq.add_action("MED -> ROF :remapMethod=redist" , med_to_rof and rof_outer_loop)
        runseq.add_action("ROF"                            , rof_outer_loop)

        #------------------
        runseq.enter_time_loop(ocn_cpl_time, newtime=ocn_outer_loop)
        #------------------

        runseq.add_action("MED med_phases_prep_ocn_accum_avg"  , med_to_ocn and ocn_outer_loop)
        runseq.add_action("MED -> OCN :remapMethod=redist"     , med_to_ocn and ocn_outer_loop)

        #------------------
        runseq.enter_time_loop(atm_cpl_time, newtime=inner_loop)
        #------------------

        if (cpl_seq_option == 'RASM'):
            runseq.add_action("MED med_phases_prep_ocn_map"        , med_to_ocn)
            if cpl_add_aoflux:
                runseq.add_action("MED med_phases_aofluxes_run"        , run_ocn and run_atm and (med_to_ocn or med_to_atm))
            runseq.add_action("MED med_phases_prep_ocn_merge"      , med_to_ocn)
            runseq.add_action("MED med_phases_prep_ocn_accum_fast" , med_to_ocn)
            runseq.add_action("MED med_phases_ocnalb_run"          , run_ocn and run_atm and (med_to_ocn or med_to_atm))
        runseq.add_action("MED med_phases_prep_lnd"                , med_to_lnd)
        runseq.add_action("MED -> LND :remapMethod=redist"         , med_to_lnd)
        runseq.add_action("MED med_phases_prep_ice"                , med_to_ice)
        runseq.add_action("MED -> ICE :remapMethod=redist"         , med_to_ice)
        runseq.add_action("MED med_phases_diag_ice_med2ice"        , run_ice and diag_mode)
        runseq.add_action("MED med_phases_prep_wav"                , med_to_wav)
        runseq.add_action("MED -> WAV :remapMethod=redist"         , med_to_wav)
        runseq.add_action("MED med_phases_prep_rof_avg"            , med_to_rof and not rof_outer_loop)
        runseq.add_action("MED -> ROF :remapMethod=redist"         , med_to_rof and not rof_outer_loop)
        runseq.add_action("MED med_phases_prep_ocn_accum_avg"      , med_to_ocn and not ocn_outer_loop)
        runseq.add_action("MED -> OCN :remapMethod=redist"         , med_to_ocn and not ocn_outer_loop)
        runseq.add_action("ICE"                                    , run_ice)
        runseq.add_action("LND"                                    , run_lnd)
        runseq.add_action("ROF"                                    , run_rof and not rof_outer_loop)
        runseq.add_action("WAV"                                    , run_wav)
        runseq.add_action("OCN"                                    , run_ocn and not ocn_outer_loop)
        if coupling_mode == 'hafs':
            runseq.add_action("OCN -> MED :remapMethod=redist:ignoreUnmatchedIndices=true"         , run_ocn and not ocn_outer_loop)
        else:
            runseq.add_action("OCN -> MED :remapMethod=redist"         , run_ocn and not ocn_outer_loop)
        if (cpl_seq_option == 'TIGHT'):
            runseq.add_action("MED med_phases_prep_ocn_map"        , med_to_ocn)
            if cpl_add_aoflux:
                runseq.add_action("MED med_phases_aofluxes_run"        , run_ocn and run_atm)
            runseq.add_action("MED med_phases_prep_ocn_merge"      , med_to_ocn)
            runseq.add_action("MED med_phases_prep_ocn_accum_fast" , med_to_ocn)
            runseq.add_action("MED med_phases_ocnalb_run"          , run_ocn and run_atm)
        runseq.add_action("MED med_phases_diag_ocn"                , run_ocn and diag_mode and not ocn_outer_loop)
        runseq.add_action("LND -> MED :remapMethod=redist"         , run_lnd)
        runseq.add_action("ICE -> MED :remapMethod=redist"         , run_ice)
        runseq.add_action("MED med_phases_diag_ice_ice2med"        , run_ice and diag_mode)
        runseq.add_action("MED med_fraction_set"                   , run_ice)

        runseq.add_action("MED med_phases_prep_rof_accum"          , med_to_rof)
        runseq.add_action("MED med_phases_prep_glc_accum"          , med_to_glc)
        runseq.add_action("MED med_phases_prep_atm"                , med_to_atm)
        runseq.add_action("MED -> ATM :remapMethod=redist"         , med_to_atm)
        runseq.add_action("ATM"                                    , run_atm)
        runseq.add_action("ATM -> MED :remapMethod=redist"         , run_atm)
        runseq.add_action("WAV -> MED :remapMethod=redist"         , run_wav)
        runseq.add_action("ROF -> MED :remapMethod=redist"         , run_rof and not rof_outer_loop)
        runseq.add_action("MED med_phases_diag_atm"                , run_atm and diag_mode)
        runseq.add_action("MED med_phases_diag_lnd"                , run_lnd and diag_mode)
        runseq.add_action("MED med_phases_diag_rof"                , run_rof and diag_mode)
        runseq.add_action("MED med_phases_diag_glc"                , run_glc and diag_mode)
        runseq.add_action("MED med_phases_diag_accum"              , diag_mode)
        runseq.add_action("MED med_phases_diag_print"              , diag_mode)
        #------------------
        runseq.leave_time_loop(inner_loop)
        #------------------

        runseq.add_action("OCN"                                    , run_ocn and ocn_outer_loop)

        if coupling_mode == 'hafs':
            runseq.add_action("OCN -> MED :remapMethod=redist:ignoreUnmatchedIndices=true"         , run_ocn and ocn_outer_loop)
        else:
            runseq.add_action("OCN -> MED :remapMethod=redist"         , run_ocn and ocn_outer_loop)


        #------------------
        runseq.leave_time_loop(ocn_outer_loop)
        #------------------

        runseq.add_action("ROF -> MED :remapMethod=redist"         , run_rof and rof_outer_loop)

        #------------------
        runseq.leave_time_loop(rof_outer_loop)
        #------------------

    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
