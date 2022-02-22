#!/usr/bin/env python

###############################################################################
# Used modules                                                                #
###############################################################################

import os

###############################################################################
# Query required information/s                                                 #
###############################################################################

fv3_path = os.environ['FV3_PATH']

###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "CMEPS"

# Add all files with metadata tables on the host model side and in CCPP,
# relative to basedir = top-level directory of host model. This includes
# kind and type definitions used in CCPP physics. Also add any internal
# dependencies of these files to the list.
VARIABLE_DEFINITION_FILES = [
    # actual variable definition files
    '{}/ccpp/framework/src/ccpp_types.F90'.format(fv3_path),
    '{}/ccpp/physics/physics/machine.F'.format(fv3_path),
    'CMEPS/ufs/ccpp/data/MED_typedefs.F90',
    'CMEPS/ufs/ccpp/data/MED_data.F90'
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_t' : 'cdata',
        'ccpp_types' : '',
        },
    'machine' : {
        'machine' : '',
        },
    'MED_typedefs' : {
        'MED_init_type' : 'physics%init',
        'MED_statein_type' : 'physics%Statein',
        'MED_interstitial_type' : 'physics%Interstitial',
        'MED_control_type' : 'physics%Model',
        'MED_coupling_type' : 'physics%Coupling',
        'MED_grid_type' : 'physics%Grid',
        'MED_sfcprop_type' : 'physics%Sfcprop',
        'MED_diag_type' : 'physics%Diag',
        'MED_typedefs' : '',
        },
    'MED_data' : {
        'MED_data' : '',
        'physics_type' : 'physics',
        }
    }

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    '{}/ccpp/physics/physics/sfc_ocean.F'.format(fv3_path),
    '{}/ccpp/physics/physics/sfc_diff.f'.format(fv3_path),
    '{}/ccpp/physics/physics/GFS_surface_loop_control.F90'.format(fv3_path),
    '{}/ccpp/physics/physics/GFS_surface_composites.F90'.format(fv3_path)
    ]

# Default build dir, relative to current working directory,
# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'CMEPS'

# Auto-generated makefile/cmakefile snippets that contain all type definitions
TYPEDEFS_MAKEFILE   = '{build_dir}/physics/CCPP_TYPEDEFS.mk'
TYPEDEFS_CMAKEFILE  = '{build_dir}/physics/CCPP_TYPEDEFS.cmake'
TYPEDEFS_SOURCEFILE = '{build_dir}/physics/CCPP_TYPEDEFS.sh'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE   = '{build_dir}/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE  = '{build_dir}/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE = '{build_dir}/physics/CCPP_SCHEMES.sh'

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE   = '{build_dir}/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE  = '{build_dir}/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE = '{build_dir}/physics/CCPP_CAPS.sh'

# Directory where to put all auto-generated physics caps
CAPS_DIR = '{build_dir}/physics'

# Directory where the suite definition files are stored
SUITES_DIR = 'CMEPS/ufs/ccpp/suites'

# Directory where to write static API to
STATIC_API_DIR = '{build_dir}/physics'
STATIC_API_CMAKEFILE  = '{build_dir}/physics/CCPP_STATIC_API.cmake'
STATIC_API_SOURCEFILE = '{build_dir}/physics/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
METADATA_HTML_OUTPUT_DIR = '{build_dir}/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = '{build_dir}/physics/CCPP_VARIABLES_CMEPS.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = '{build_dir}/framework/doc/DevelopersGuide/CCPP_VARIABLES_CMEPS.tex'
