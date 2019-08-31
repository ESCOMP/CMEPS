ifndef NETCDF
ifndef PNETCDF
$(error NETCDF or PNETCDF must be set to build PIO)
endif
endif

ifndef PIO_INSTALL_DIR
$(error PIO_INSTALL_DIR must be set)
endif

BASE_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PIO_SRC_DIR := "$(BASE_DIR)/ParallelIO"
PIO_BUILD_DIR := "$(BASE_DIR)/ParallelIO/build"

CMAKE_ENV_VARS := CC=$(CC)   \
                  CXX=$(CXX) \
                  FC=$(FC)

#CMAKE_OPTS += -D CMAKE_Fortran_FLAGS:STRING="$(FFLAGS) $(EXTRA_PIO_FPPDEFS) $(INCLDIR)" \
#              -D CMAKE_C_FLAGS:STRING="$(CFLAGS) $(EXTRA_PIO_CPPDEFS) $(INCLDIR)" \
#              -D CMAKE_CXX_FLAGS:STRING="$(CXXFLAGS) $(EXTRA_PIO_CPPDEFS) $(INCLDIR)" \
#              -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
#              -D GPTL_PATH:STRING=$(INSTALL_SHAREDPATH) \
#              -D PIO_ENABLE_TESTS:BOOL=OFF \
#              -D PIO_USE_MALLOC:BOOL=OFF \
#              -D USER_CMAKE_MODULE_PATH:LIST="$(CIMEROOT)/src/CMake;$(CIMEROOT)/src/externals/pio2/cmake" \

CMAKE_OPTS := -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
              -D PIO_ENABLE_TESTS:BOOL=OFF \
              -D PIO_USE_MALLOC:BOOL=OFF \
              -D CMAKE_INSTALL_PREFIX=$(PIO_INSTALL_DIR)  

CMAKE_OPTS += -D NetCDF_PATH:PATH=$(NETCDF)
CMAKE_OPTS += -D PnetCDF_PATH:STRING=$(PNETCDF)

.PHONY: install build clean

install: build

build: $(PIO_BUILD_DIR)/Makefile
	cd $(PIO_BUILD_DIR); \
	exec $(MAKE) install

$(PIO_BUILD_DIR)/Makefile:
	mkdir -p $(PIO_BUILD_DIR); \
	cd $(PIO_BUILD_DIR); \
	$(CMAKE_ENV_VARS) cmake $(CMAKE_OPTS) ../

clean:
	rm -rf $(PIO_BUILD_DIR)
	rm -rf $(PIO_INSTALL_DIR)
