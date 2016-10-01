# - Locate GFAL library
# Defines:
#
#  GFAL_FOUND
#  GFAL_INCLUDE_DIR
#  GFAL_INCLUDE_DIRS (not cached)
#  GFAL_LIBRARIES (not cached)

find_path(GFAL_INCLUDE_DIR NAMES gfal_api.h
          PATH_SUFFIXES . gfal gfal2
          HINTS ${GFAL_DIR}/include $ENV{GFAL_DIR}/include)
find_library(GFAL_LIBRARY NAMES gfal gfal2
             HINTS ${GFAL_DIR}/lib $ENV{GFAL_DIR}/lib)
find_path(SRM_IFCE_INCLUDE_DIR  gfal_srm_ifce_types.h 
          HINTS ${SRM_IFCE_DIR}/include $ENV{SRM_IFCE_DIR}/include)

set(GFAL_LIBRARIES ${GFAL_LIBRARY})
set(GFAL_INCLUDE_DIRS ${GFAL_INCLUDE_DIR} ${SRM_IFCE_INCLUDE_DIR})

if(GFAL_LIBRARY MATCHES gfal2)
  # use pkg-config to get the directories for glib and then use these values
  find_package(PkgConfig)
  pkg_check_modules(GLIB2 glib-2.0)
  list(APPEND GFAL_INCLUDE_DIRS ${GLIB2_INCLUDE_DIRS})
endif()

# handle the QUIETLY and REQUIRED arguments and set GFAL_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GFAL DEFAULT_MSG GFAL_INCLUDE_DIR SRM_IFCE_INCLUDE_DIR GFAL_LIBRARY)

mark_as_advanced(GFAL_FOUND GFAL_INCLUDE_DIR GFAL_LIBRARIES SRM_IFCE_INCLUDE_DIR GLIB_INCLUDE_DIR)
