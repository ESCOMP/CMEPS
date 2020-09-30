# - Try to find PIO
#
# This can be controled by setting PIO_PATH or PIO_<lang>_PATH Cmake variables,
# where <lang> is the COMPONENT language one needs.
#
# Once done, this will define:
#
#  PIO_<lang>_FOUND   (BOOL) - system has PIO
#  PIO_<lang>_IS_SHARED (BOOL) - whether the library is shared/dynamic
#  PIO_<lang>_INCLUDE_DIR (PATH) - Location of the header files and modules
#  PIO_<lang>_LIBRARY     (File) - Path to the <lang> library files
#  PIO_<lang>_LIBRARIES   (List) - link these to use PIO
#
# Available COMPONENTS are: C Fortran
# If no components are specified only C is assumed
include (LibFind)
include (LibCheck)

# Define PIO C Component
define_package_component(PIO DEFAULT
                         COMPONENT C
                         INCLUDE_NAMES pio.h
                         LIBRARY_NAMES pioc)

# Define PIO Fortran Component
define_package_component(PIO
                         COMPONENT Fortran
                         INCLUDE_NAMES pio.mod pio.inc
                         LIBRARY_NAMES piof)

# Search for list of valid components requested
find_valid_components(PIO)

#==============================================================================
# SEARCH FOR VALIDATED COMPONENTS
foreach (pcomp IN LISTS PIO_FIND_VALID_COMPONENTS)

    # If not found already, search...
    if (NOT PIO_${pcomp}_FOUND)

        # Manually add the MPI include and library dirs to search paths
        # and search for the package component
        if (MPI_${pcomp}_FOUND)
            initialize_paths (PIO_${pcomp}_PATHS
                              INCLUDE_DIRECTORIES ${MPI_${pcomp}_INCLUDE_PATH}
                              LIBRARIES ${MPI_${pcomp}_LIBRARIES})
            find_package_component(PIO COMPONENT ${pcomp}
                                   PATHS ${PIO_${pcomp}_PATHS})
        else ()
            find_package_component(PIO COMPONENT ${pcomp} HINT PIO_${pcomp}_PATH=${PIO_PATH})
        endif ()

        # Continue only if component found
        if (PIO_${pcomp}_FOUND)

            # Checks
            if (pcomp STREQUAL C)

                # Check version
                check_version (PIO
                               NAME "pio_meta.h"
                               HINTS ${PIO_C_INCLUDE_DIRS}
                               MACRO_REGEX "PIO_VERSION_")

            endif ()

            # Dependencies
            if (pcomp STREQUAL C AND NOT PIO_C_IS_SHARED)

                # DEPENDENCY: PnetCDF (if PnetCDF enabled)
                check_macro (PIO_HAS_PNETCDF
                             NAME TryPIO_PNETCDF.c
                             HINTS ${CMAKE_MODULE_PATH}
                             DEFINITIONS -I${PIO_C_INCLUDE_DIR}
                             COMMENT "whether PIO has PnetCDF support")
                if (PIO_HAS_PNETCDF)
                    find_package (PnetCDF COMPONENTS C)
                endif ()


            elseif (pcomp STREQUAL Fortran AND NOT PIO_Fortran_IS_SHARED)

                # DEPENDENCY: PIO
                set (orig_comp ${pcomp})
                set (orig_comps ${PIO_FIND_VALID_COMPONENTS})
                find_package (PIO COMPONENTS C)
                set (PIO_FIND_VALID_COMPONENTS ${orig_comps})
                set (pcomp ${orig_comp})
                if (PIO_C_FOUND)
                    list (APPEND PIO_Fortran_INCLUDE_DIRS ${PIO_C_INCLUDE_DIRS})
                    list (APPEND PIO_Fortran_LIBRARIES ${PIO_C_LIBRARIES})
                endif ()

            endif ()

        endif ()

    endif ()

endforeach ()
message("PIO_C_FOUND ${PIO_C_FOUND}")
message("PIO_Fortran_FOUND ${PIO_Fortran_FOUND}")
message("PIO_Fortran_INCLUDE_DIR ${PIO_Fortran_INCLUDE_DIR}")
