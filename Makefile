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

ifndef PIO_LIBDIR
$(error PIO_LIBDIR should point to PIO library directory.)
endif

ifndef PIO_INC
$(error PIO_INC should point to PIO include directory.)
endif

ifndef INTERNAL_PIO_INIT
INTERNAL_PIO_INIT := 1
endif
$(info INTERNAL_PIO_INIT is set to $(INTERNAL_PIO_INIT))

MEDIATOR_DIR := $(BASE_DIR)/mediator
LIBRARY_MEDIATOR := $(MEDIATOR_DIR)/libcmeps.a
LIBRARY_UTIL := $(BASE_DIR)/util/libcmeps_util.a

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
	@echo "ESMF_DEP_LINK_OBJS = $(INSTALLDIR)/libcmeps.a $(INSTALLDIR)/libcmeps_util.a $(PIO_LIBDIR)/libpiof.a $(PIO_LIBDIR)/libpioc.a $(PNETCDF_LD_OPTS)" >> cmeps.mk.install
	mkdir -p $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/include
	cp -f $(LIBRARY_UTIL) $(INSTALLDIR)
	cp -f $(LIBRARY_MEDIATOR) $(INSTALLDIR)
	cp -f mediator/*.mod $(INSTALLDIR)/include
	cp -f util/*.mod $(INSTALLDIR)/include
	cp -f cmeps.mk.install $(INSTALLDIR)/cmeps.mk
endif

$(LIBRARY_MEDIATOR): $(LIBRARY_UTIL) .FORCE
	cd mediator ;\
	exec $(MAKE) PIO_INC=$(PIO_INC) INTERNAL_PIO_INIT=$(INTERNAL_PIO_INIT)

$(LIBRARY_UTIL): .FORCE
	cd util ;\
	exec $(MAKE) PIO_INC=$(PIO_INC)

.FORCE:

clean:
	cd mediator; \
	exec $(MAKE) clean
	cd util; \
	exec $(MAKE) clean

