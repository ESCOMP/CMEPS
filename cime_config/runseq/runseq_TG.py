#!/usr/bin/env python

import os, shutil, sys

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *

logger = logging.getLogger(__name__)

def runseq(case, coupling_times):

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")

    outfile   = open(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), "w")

    glc_cpl_dt = coupling_times["glc_cpl_dt"]

    outfile.write ("runSeq::                         \n" )
    outfile.write ("@" + str(glc_cpl_dt) + "         \n" )
    outfile.write ("  LND                            \n" )
    outfile.write ("  LND -> MED :remapMethod=redist \n" )
    outfile.write ("  MED med_fraction_set           \n" )
    outfile.write ("  MED med_phases_prep_glc_accum  \n" )
    outfile.write ("  MED med_phases_prep_glc_avg    \n" )
    outfile.write ("  MED -> GLC :remapMethod=redist \n" )
    outfile.write ("  GLC                            \n" )
    outfile.write ("  GLC -> MED :remapMethod=redist \n" )
    outfile.write ("  MED med_phases_history_write   \n" )
    outfile.write ("  MED med_phases_profile         \n" )
    outfile.write ("  MED med_phases_restart_write   \n" )
    outfile.write ("@                                \n" )
    outfile.write ("::                               \n" )

    outfile.close()
    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)




