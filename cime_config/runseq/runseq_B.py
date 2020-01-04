#!/usr/bin/env python

import os, shutil, sys

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *
#pylint:disable=undefined-variable
logger = logging.getLogger(__name__)

def gen_runseq(case, coupling_times):

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")

    ocn_cpl_dt = coupling_times["ocn_cpl_dt"]
    atm_cpl_dt = coupling_times["atm_cpl_dt"]

    outfile   = open(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), "w")

    outfile.write ("runSeq::                                   \n")
    outfile.write ("@" + str(ocn_cpl_dt) + "                   \n")
    outfile.write ("  MED med_phases_prep_ocn_accum_avg        \n")
    outfile.write ("  MED -> OCN :remapMethod=redist           \n")
    if atm_cpl_dt < ocn_cpl_dt:
        outfile.write ("@" + str(atm_cpl_dt) + "               \n")
    outfile.write ("  MED med_phases_prep_ocn_map              \n")
    outfile.write ("  MED med_phases_aofluxes_run              \n")
    outfile.write ("  MED med_phases_prep_ocn_merge            \n")
    outfile.write ("  MED med_phases_prep_ocn_accum_fast       \n")
    outfile.write ("  MED med_phases_ocnalb_run                \n")
    outfile.write ("  MED med_phases_prep_lnd                  \n")
    outfile.write ("  MED -> LND :remapMethod=redist           \n")
    outfile.write ("  MED med_phases_prep_ice                  \n")
    outfile.write ("  MED -> ICE :remapMethod=redist           \n")
    outfile.write ("  ICE                                      \n")
    outfile.write ("  LND                                      \n")
    outfile.write ("  ICE -> MED :remapMethod=redist           \n")
    outfile.write ("  MED med_fraction_set                     \n")
    outfile.write ("  LND -> MED :remapMethod=redist           \n")
    outfile.write ("  MED med_phases_prep_atm                  \n")
    outfile.write ("  MED -> ATM :remapMethod=redist           \n")
    outfile.write ("  ATM                                      \n")
    outfile.write ("  ATM -> MED :remapMethod=redist           \n")
    outfile.write ("  MED med_phases_history_write             \n")
    outfile.write ("  MED med_phases_profile                   \n")
    if atm_cpl_dt < ocn_cpl_dt:
        outfile.write ("@                                      \n")
    outfile.write ("  OCN                                      \n")
    outfile.write ("  OCN -> MED :remapMethod=redist           \n")
    outfile.write ("  MED med_phases_restart_write             \n")
    outfile.write ("@                                          \n")
    outfile.write ("::                                         \n")

    outfile.close()
    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
