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
    'CMEPS/ufs/ccpp/data/GFS_typedefs.F90',
    'CMEPS/ufs/ccpp/data/med_typedefs.F90'
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_t' : 'cdata',
        'ccpp_types' : '',
        },
    'machine' : {
        'machine' : '',
        },
    'GFS_typedefs' : {
        'GFS_statein_type' : 'physics%Statein',
        'GFS_typedefs' : '',
        },
    'med_typedefs' : {
        'med_typedefs' : '',
        'physics_type' : 'physics',
        }
    }

# Add all physics scheme files relative to basedir
SCHEME_FILES = ['{}/ccpp/physics/physics/sfc_ocean.F'.format(fv3_path)]
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    #'{}/ccpp/physics/physics/GFS_DCNV_generic.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_GWD_generic.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_MP_generic.F90'.format(fv3_pathmt(fv3_path),
    #'{}/ccpp/physics/physics/GFS_PBL_generic.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_SCNV_generic.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_debug.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_phys_time_vary.fv3.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rad_time_vary.fv3.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_radiation_surface.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmg_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmg_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmg_setup.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_stochastics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_suite_interstitial.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_surface_generic.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_surface_composites.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_surface_loop_control.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_time_vary_pre.fv3.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cires_ugwp.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cires_ugwp_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/unified_ugwp.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/unified_ugwp_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/ugwpv1_gsldrag.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/ugwpv1_gsldrag_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cnvc90.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/cs_conv.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cs_conv_aw_adj.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cu_ntiedtke_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cu_ntiedtke.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cu_ntiedtke_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/dcyc2.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/drag_suite.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/gcm_shoc.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/get_prs_fv3.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/gfdl_cloud_microphys.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/gfdl_fv_sat_adj.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/gfdl_sfc_layer.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/gscond.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/gwdc.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/gwdps.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/h2ophys.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/samfdeepcnv.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/samfshalcnv.f',
    #'{}/ccpp/physics/physics/sascnvn.F'.format(fv3_path),
    #'{}/ccpp/physics/physics/shalcnv.F'.format(fv3_path),
    #'{}/ccpp/physics/physics/maximum_hourly_diagnostics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/m_micro.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/m_micro_interstitial.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cu_gf_driver_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cu_gf_driver.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/cu_gf_driver_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/moninedmf.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/moninshoc.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/satmedmfvdif.F'.format(fv3_path),
    #'{}/ccpp/physics/physics/satmedmfvdifq.F'.format(fv3_path),
    #'{}/ccpp/physics/physics/shinhongvdif.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/ysuvdif.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/module_MYNNPBL_wrapper.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/module_MYNNSFC_wrapper.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/module_SGSCloud_RadPre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/module_SGSCloud_RadPost.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/module_MYJSFC_wrapper.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/module_MYJPBL_wrapper.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/mp_thompson_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/mp_thompson.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/mp_thompson_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/ozphys.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/ozphys_2015.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/precpd.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/phys_tend.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/radlw_main.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/radsw_main.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rascnv.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rayleigh_damp.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmg_lw_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmg_lw_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmg_sw_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmg_sw_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_diag.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_diag_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_drv_ruc.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_cice.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_diff.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_drv.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_noahmp_drv.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/flake_driver.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_nst.f'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_ocean.F'.format(fv3_path),
    #'{}/ccpp/physics/physics/sfc_sice.f'.format(fv3_path),
    ## HAFS FER_HIRES
    #'{}/ccpp/physics/physics/mp_fer_hires.F90'.format(fv3_path),
    ## RRTMGP
    #'{}/ccpp/physics/physics/rrtmgp_lw_gas_optics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_lw_cloud_optics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_sw_gas_optics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_sw_cloud_optics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_sw_aerosol_optics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_lw_rte.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_sw_rte.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_lw_aerosol_optics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_setup.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_lw_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_sw_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_lw_post.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_lw_cloud_sampling.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/rrtmgp_sw_cloud_sampling.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_cloud_diagnostics.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_thompsonmp_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_gfdlmp_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_zhaocarr_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_cloud_overlap_pre.F90'.format(fv3_path),
    #'{}/ccpp/physics/physics/GFS_rrtmgp_sw_post.F90'.format(fv3_path)
    #]

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
SUITES_DIR = '{}/ccpp/suites'.format(fv3_path)

# Directory where to write static API to
STATIC_API_DIR = '{build_dir}/physics'
STATIC_API_SRCFILE = '{build_dir}/physics/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
METADATA_HTML_OUTPUT_DIR = '{build_dir}/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = '{build_dir}/physics/CCPP_VARIABLES_CMEPS.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = '{build_dir}/framework/doc/DevelopersGuide/CCPP_VARIABLES_CMEPS.tex'
