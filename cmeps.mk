#-----------------------------------------------
# NUOPC/ESMF self-describing build dependency
# makefile fragment for CMEPS
#-----------------------------------------------

# component module name
#MED_ESMF_DEP_FRONT := MED
#MED_ESMF_DEP_SHRD_PATH := /opt/PIO/lib /tmp/CMEPS/mediator /tmp/CMEPS/shared /tmp/CMEPS/ufs /opt/libFMS/intel/32bit /opt/nceplibs/lib
#MED_ESMF_DEP_SHRD_LIBS := pioc piof cmeps shared cmeps_share

# component module name
MED_ESMF_DEP_FRONT := MED
# component module path
MED_ESMF_DEP_INCPATH := /opt/PIO/include
# component module objects
MED_ESMF_DEP_CMPL_OBJS := /tmp/CMEPS/mediator/*.o /tmp/CMEPS/shared/*.o /tmp/CMEPS/ufs/*.o
# component object/archive list
MED_ESMF_DEP_LINK_OBJS := /tmp/CMEPS/mediator/libcmeps.a /tmp/CMEPS/shared/libshared.a /tmp/CMEPS/ufs/libcmeps_share.a
MED_ESMF_DEP_SHRD_PATH := /opt/PIO/lib  /opt/libFMS/intel/32bit /opt/nceplibs/lib /tmp/CMEPS/mediator /tmp/CMEPS/shared /tmp/CMEPS/ufs
MED_ESMF_DEP_SHRD_LIBS := pioc piof cmeps shared cmeps_share
