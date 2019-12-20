#!/usr/bin/env python

import os, shutil, sys
from CIME.utils import expect

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *
#pylint:disable=undefined-variable
logger = logging.getLogger(__name__)

def runseq(case, coupling_times):

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")
    comp_rof  = case.get_value("COMP_ROF")
    comp_glc  = case.get_value("COMP_GLC")

    outfile   = open(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), "w")

    rof_cpl_dt = coupling_times["rof_cpl_dt"]
    atm_cpl_dt = coupling_times["atm_cpl_dt"]
    glc_cpl_dt = coupling_times["glc_cpl_dt"]

    glc_couplint_time = glc_cpl_dt

    # NOTE: in the mediator the glc_avg_period will be set as an alarm
    # on the mediator clock - when this alarm rings - the averaging will be done AND an
    # attribute will be set on the glc export state from the mediator
    # saying that the data coming to glc is valid

    if (comp_rof == 'srof'):
        med_to_rof = False
        rof_to_med = False
    elif (comp_rof == 'drof'):
        med_to_rof = False
        rof_to_med = True
    else:
        med_to_rof = True
        rof_to_med = True

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

    outfile.write ("runSeq::                                  \n" )

    if glc_to_med or med_to_glc:
        outfile.write ("@" + str(glc_coupling_time) + "       \n" )
    if med_to_glc:
        outfile.write ("  MED med_phases_prep_glc_accum       \n" )
        outfile.write ("  MED med_phases_prep_glc_avg         \n" )
        outfile.write ("  MED -> GLC :remapMethod=redist      \n" )
        outfile.write ("  GLC                                 \n" )
    if glc_to_med:
        outfile.write ("  GLC -> MED :remapMethod=redist      \n" )

    if med_to_rof:
        outfile.write ("@" + str(rof_cpl_dt) + "              \n" )
        outfile.write ("  MED med_phases_prep_rof_avg         \n" )
        outfile.write ("  MED -> ROF :remapMethod=redist      \n" )
        outfile.write ("  ROF                                 \n" )
        if atm_cpl_dt < rof_cpl_dt:
            outfile.write ("@" + str(atm_cpl_dt) + "          \n" )
    else:
        outfile.write ("@" + str(atm_cpl_dt) + "              \n" )

    outfile.write ("  MED med_phases_prep_lnd                 \n" )
    outfile.write ("  MED -> LND :remapMethod=redist          \n" )
    outfile.write ("  LND                                     \n" )
    outfile.write ("  LND -> MED :remapMethod=redist          \n" )
    outfile.write ("  MED med_phases_prep_rof_accum           \n" )
    outfile.write ("  ATM                                     \n" )
    outfile.write ("  ATM -> MED :remapMethod=redist          \n" )
    outfile.write ("  MED med_phases_profile                  \n" )

    if rof_to_med:
        if atm_cpl_dt < rof_cpl_dt:
            outfile.write ("@                                 \n" )
        outfile.write ("  ROF -> MED :remapMethod=redist      \n" )

    outfile.write ("  MED med_phases_history_write            \n" )
    outfile.write ("  MED med_phases_restart_write            \n" )
    outfile.write ("  MED med_phases_profile                  \n" )

    if med_to_rof:
        outfile.write ("@                                     \n" )
        
    if glc_to_med or med_to_glc:
        outfile.write ("@                                     \n" )

    outfile.write ("::                                        \n" )

    outfile.close()
    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
