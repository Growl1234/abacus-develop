# Find ELPA, preferring pkg-config unless deprecated manual inputs are supplied.
# Both discovery paths provide ELPA::ELPA with complete usage requirements.

# Migrate caches created when FindPkgConfig used the public ELPA prefix.
if(DEFINED CACHE{ELPA_LIBRARIES})
  get_property(
    _elpa_libraries_cache_type
    CACHE ELPA_LIBRARIES
    PROPERTY TYPE)
  if(_elpa_libraries_cache_type STREQUAL "INTERNAL")
    unset(ELPA_LIBRARIES CACHE)
    unset(ELPA_INCLUDE_DIRS CACHE)
    unset(ELPA_LINK_LIBRARIES CACHE)
  endif()
endif()
unset(_elpa_libraries_cache_type)

if(DEFINED ELPA_DIR)
  message(
    WARNING "ELPA_DIR is deprecated and will be removed in a future release.")
endif()

set(_elpa_use_pkgconfig FALSE)
if(NOT ELPA_DIR
   AND NOT ELPA_INCLUDE_DIR
   AND NOT ELPA_INCLUDE_DIRS
   AND NOT ELPA_LIBRARIES
   AND NOT ELPA_LINK_LIBRARIES)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    if(ENABLE_OPENMP)
      pkg_search_module(ELPA_PKG QUIET IMPORTED_TARGET GLOBAL elpa_openmp elpa)
    else()
      pkg_search_module(ELPA_PKG QUIET IMPORTED_TARGET GLOBAL elpa)
    endif()
    set(_elpa_use_pkgconfig ${ELPA_PKG_FOUND})
  endif()
endif()

if(_elpa_use_pkgconfig)
  set(ELPA_INCLUDE_DIRS ${ELPA_PKG_INCLUDE_DIRS})
  set(ELPA_LINK_LIBRARIES ${ELPA_PKG_LINK_LIBRARIES})
  set(ELPA_LIBRARIES ${ELPA_PKG_LIBRARIES})
  set(ELPA_VERSION ${ELPA_PKG_VERSION})
  set(_elpa_required_vars ELPA_PKG_FOUND)
else()
  if(ELPA_INCLUDE_DIR)
    set(ELPA_INCLUDE_DIRS ${ELPA_INCLUDE_DIR})
  endif()
  if(ELPA_LIBRARIES)
    set(ELPA_LINK_LIBRARIES ${ELPA_LIBRARIES})
  endif()
  find_path(
    ELPA_INCLUDE_DIRS elpa/elpa.h
    HINTS ${ELPA_DIR}
    PATH_SUFFIXES include include/elpa)
  if(ENABLE_OPENMP)
    set(_elpa_names elpa_openmp elpa)
  else()
    set(_elpa_names elpa)
  endif()
  find_library(
    ELPA_LINK_LIBRARIES
    NAMES ${_elpa_names}
    HINTS ${ELPA_DIR}
    PATH_SUFFIXES lib lib64)
  set(_elpa_required_vars ELPA_LINK_LIBRARIES ELPA_INCLUDE_DIRS)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  ELPA
  REQUIRED_VARS ${_elpa_required_vars}
  VERSION_VAR ELPA_VERSION)

if(ELPA_FOUND)
  list(GET ELPA_LINK_LIBRARIES 0 ELPA_LIBRARY)
  set(ELPA_INCLUDE_DIR ${ELPA_INCLUDE_DIRS})
  if(NOT TARGET ELPA::ELPA)
    if(_elpa_use_pkgconfig)
      add_library(ELPA::ELPA ALIAS PkgConfig::ELPA_PKG)
    else()
      add_library(ELPA::ELPA INTERFACE IMPORTED)
      set_target_properties(
        ELPA::ELPA
        PROPERTIES INTERFACE_LINK_LIBRARIES "${ELPA_LINK_LIBRARIES}"
                   INTERFACE_INCLUDE_DIRECTORIES "${ELPA_INCLUDE_DIRS}")
    endif()
  endif()

  include(CheckCXXSourceCompiles)
  include(CMakePushCheckState)
  cmake_push_check_state(RESET)
  # Older CMake versions cannot export alias targets into try_compile projects.
  get_target_property(_elpa_check_target ELPA::ELPA ALIASED_TARGET)
  if(NOT _elpa_check_target)
    set(_elpa_check_target ELPA::ELPA)
  endif()
  set(CMAKE_REQUIRED_LIBRARIES ${_elpa_check_target})
  unset(ELPA_VERSION_SATISFIES CACHE)
  check_cxx_source_compiles(
    "#include <elpa/elpa_version.h>
     #if ELPA_API_VERSION < 20210430
     #error ELPA version is too old.
     #endif
     int main() {}"
    ELPA_VERSION_SATISFIES)
  cmake_pop_check_state()
  if(NOT ELPA_VERSION_SATISFIES)
    message(FATAL_ERROR "ELPA must provide API version 20210430 or newer.")
  endif()
endif()

unset(_elpa_use_pkgconfig)
unset(_elpa_required_vars)
unset(_elpa_names)
unset(_elpa_check_target)
mark_as_advanced(ELPA_INCLUDE_DIR ELPA_INCLUDE_DIRS ELPA_LIBRARY
                 ELPA_LINK_LIBRARIES)
