# BASE_DIR points to root of CMEPS clone
BASE_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)

ifndef FC
$(error FC not defined)
endif

ifndef CC
$(error CC not defined)
endif

ifndef CXX
$(error CXX not defined)
endif

MEDIATOR_DIR := $(BASE_DIR)/mediator
LIBRARY_MEDIATOR := $(MEDIATOR_DIR)/libcmeps.a
LIBRARY_UTIL := $(BASE_DIR)/util/libcmeps_util.a
PIO_INSTALL_DIR := $(BASE_DIR)/lib/ParallelIO/install
PIO_INSTALL_LIBS := $(PIO_INSTALL_DIR)/lib/libpiof.a
PIO_INCLUDE_DIR := $(PIO_INSTALL_DIR)/include

all default: install

install: $(LIBRARY_MEDIATOR)
ifndef INSTALLDIR
	$(error INSTALLDIR not defined for CMEPS installation location)
else
	rm -f cmeps.mk.install
	@echo "# ESMF self-describing build dependency makefile fragment" > cmeps.mk.install
	@echo "# src location: $(PWD)" >> cmeps.mk.install
	@echo  >> cmeps.mk.install
	@echo "ESMF_DEP_FRONT     = MED" >> cmeps.mk.install
	@echo "ESMF_DEP_INCPATH   = $(INSTALLDIR)/include" >> cmeps.mk.install
	@echo "ESMF_DEP_CMPL_OBJS = " >> cmeps.mk.install
	@echo "ESMF_DEP_LINK_OBJS = $(INSTALLDIR)/libcmeps.a $(INSTALLDIR)/libcmeps_util.a $(INSTALLDIR)/libgptl.a $(INSTALLDIR)/libpiof.a $(INSTALLDIR)/libpioc.a $(INSTALLDIR)/libgptl.a $(PNETCDF_LD_OPTS)" >> cmeps.mk.install
	mkdir -p $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/include
	cp -f $(PIO_INSTALL_DIR)/lib/*.a $(INSTALLDIR)
	cp -f $(PIO_INSTALL_DIR)/include/* $(INSTALLDIR)/include
	cp -f $(LIBRARY_UTIL) $(INSTALLDIR)
	cp -f $(LIBRARY_MEDIATOR) $(INSTALLDIR)
	cp -f mediator/*.mod $(INSTALLDIR)/include
	cp -f util/*.mod $(INSTALLDIR)/include
	cp -f cmeps.mk.install $(INSTALLDIR)/cmeps.mk
endif

$(LIBRARY_MEDIATOR): $(LIBRARY_UTIL) .FORCE
	cd mediator ;\
	exec $(MAKE) PIO_INCLUDE_DIR=$(PIO_INCLUDE_DIR)

$(LIBRARY_UTIL): $(PIO_INSTALL_LIBS) .FORCE
	cd util ;\
	exec $(MAKE) PIO_INCLUDE_DIR=$(PIO_INCLUDE_DIR) 

$(PIO_INSTALL_LIBS):
	cd lib ;\
	exec $(MAKE) install FC="$(FC)" CC="$(CC)" CXX="$(CXX)" PIO_INSTALL_DIR=$(PIO_INSTALL_DIR)

.FORCE:

clean:
	cd mediator; \
	exec $(MAKE) clean PIO_INCLUDE_DIR=$(PIO_INCLUDE_DIR)
	cd util; \
	exec $(MAKE) clean PIO_INCLUDE_DIR=$(PIO_INCLUDE_DIR)
	cd lib; \
	exec $(MAKE) clean PIO_INSTALL_DIR=$(PIO_INSTALL_DIR)

