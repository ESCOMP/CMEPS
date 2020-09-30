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
LIBRARY_UTIL := $(BASE_DIR)/share/libcmeps_util.a

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
	@echo "ESMF_DEP_LINK_OBJS = $(INSTALLDIR)/libcmeps.a $(INSTALLDIR)/libcmeps_util.a $(INSTALLDIR)/libpiof.a $(INSTALLDIR)/libpioc.a $(PNETCDF_LD_OPTS)" >> cmeps.mk.install
	mkdir -p $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/include
	cp -f $(LIBRARY_UTIL) $(INSTALLDIR)
	cp -f $(LIBRARY_MEDIATOR) $(INSTALLDIR)
	cp -f mediator/*.mod $(INSTALLDIR)/include
	cp -f share/*.mod $(INSTALLDIR)/include
	cp -f cmeps.mk.install $(INSTALLDIR)/cmeps.mk
endif

$(LIBRARY_MEDIATOR): $(LIBRARY_UTIL) .FORCE
	cd mediator ;\
	exec $(MAKE) PIO_INCLUDE_DIR=$(PIO_INCLUDE_DIR) INTERNAL_PIO_INIT=1

$(LIBRARY_UTIL): $(PIO_INSTALL_LIBS) .FORCE
	cd share ;\
	exec $(MAKE) PIO_INCLUDE_DIR=$(PIO_INCLUDE_DIR)

.FORCE:

clean:
	cd mediator; \
	exec $(MAKE) clean PIO_INCLUDE_DIR=$(PIO_INCLUDE_DIR)
	cd share; \
	exec $(MAKE) clean PIO_INCLUDE_DIR=$(PIO_INCLUDE_DIR)
