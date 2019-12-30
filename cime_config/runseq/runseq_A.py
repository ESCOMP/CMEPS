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
    comp_atm  = case.get_value("COMP_ATM")
    comp_ice  = case.get_value("COMP_ICE")
    comp_lnd  = case.get_value("COMP_LND")
    comp_ocn  = case.get_value("COMP_OCN")
    comp_rof  = case.get_value("COMP_ROF")
    comp_wav  = case.get_value("COMP_WAV")

    atm_cpl_dt = coupling_times["atm_cpl_dt"]
    ocn_cpl_dt = coupling_times["ocn_cpl_dt"]
    wav_cpl_dt = coupling_times["wav_cpl_dt"]
    lnd_cpl_dt = coupling_times["lnd_cpl_dt"]

    #force_prognostic = case.get_value("FORCE_PROGNOSTIC_FOR_DATAMODELS")
    force_prognostic = False

    if (comp_atm == 'datm'):
        couple_atm = True
        med_to_atm = force_prognostic
    else:
        couple_atm = False
        med_to_atm = False

    if (comp_lnd == 'dlnd'):
        couple_lnd = True
        med_to_lnd = force_prognostic
    else:
        couple_lnd = False
        med_to_lnd = False

    if (comp_ice == 'dice'):
        couple_ice = True
        med_to_ice = force_prognostic
    else:
        couple_ice = False
        med_to_ice = False

    if (comp_ocn == 'docn'):
        couple_ocn = True
        docn_mode = case.get_value("DOCN_MODE")
        med_to_ocn = ('som' in docn_mode or 'interannual' in docn_mode or force_prognostic)
    else:
        couple_ocn = False
        med_to_ocn = False 

    if (comp_rof == 'drof'):
        couple_rof = True
        med_to_rof = force_prognostic
    else:
        couple_rof = False
        med_to_rof = False

    if (comp_wav == 'dwav'):
        couple_wav = True
        med_to_wav = force_prognostic
    else:
        couple_wav = False
        med_to_wav = False

    # set coupling time
    if couple_wav:
        coupling_time = wav_cpl_dt
    elif couple_lnd:
        coupling_time = lnd_cpl_dt
    else:
        coupling_time = ocn_cpl_dt

    with RunSeq(os.path.join(caseroot, "CaseDocs", "nuopc.runseq")) as runseq:

        runseq.enter_time_loop(coupling_time)

        runseq.add_action ("MED med_phases_prep_ocn_accum_avg"  , med_to_ocn)
        runseq.add_action ("MED -> OCN :remapMethod=redist"     , med_to_ocn)

        runseq.add_action ("MED med_phases_prep_ocn_map"        , med_to_ocn)
        runseq.add_action ("MED med_phases_aofluxes_run"        , med_to_atm or med_to_ocn)
        runseq.add_action ("MED med_phases_prep_ocn_merge"      , med_to_ocn)
        runseq.add_action ("MED med_phases_prep_ocn_accum_fast" , med_to_ocn)
        runseq.add_action ("MED med_phases_ocnalb_run"          , med_to_atm or med_to_ocn)

        runseq.add_action ("MED med_phases_prep_ice"            , couple_ice)
        runseq.add_action ("MED -> ICE :remapMethod=redist"     , couple_ice)
        runseq.add_action ("MED med_phases_prep_rof_accum"      , med_to_rof)
        runseq.add_action ("MED med_phases_prep_rof_avg"        , med_to_rof)
        runseq.add_action ("MED -> ROF :remapMethod=redist"     , med_to_rof)
        runseq.add_action ("ICE"                                , couple_ice)
        runseq.add_action ("ROF"                                , couple_ice)
        runseq.add_action ("ICE -> MED :remapMethod=redist"     , couple_ice)
        runseq.add_action ("MED med_fraction_set"               , couple_ice)
        runseq.add_action ("ROF -> MED :remapMethod=redist"     , couple_rof)
        runseq.add_action ("MED med_phases_prep_atm"            , med_to_atm)
        runseq.add_action ("MED -> ATM :remapMethod=redist"     , med_to_atm)
        runseq.add_action ("ATM"                                , couple_atm)
        runseq.add_action ("ATM -> MED :remapMethod=redist"     , couple_atm)
        runseq.add_action ("OCN"                                , couple_ocn)
        runseq.add_action ("OCN -> MED :remapMethod=redist"     , couple_ocn)
        runseq.add_action ("WAV"                                , couple_wav)
        runseq.add_action ("WAV -> MED :remapMethod=redist"     , couple_wav)
        runseq.add_action ("LND"                                , couple_lnd)
        runseq.add_action ("LND -> MED :remapMethod=redist"     , couple_lnd)
        runseq.add_action ("MED med_fraction_set"               , couple_lnd)
        runseq.add_action ("MED med_phases_history_write"       , True)

        runseq.leave_time_loop(True, if_write_hist_rest=True)

    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
