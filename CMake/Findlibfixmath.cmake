# - Try to find libfixmath
# Once done this will define
#  LIBFIXMATH_FOUND - System has libfixmath
#  LIBFIXMATH_INCLUDE_DIRS - The libfixmath include directories
#  LIBFIXMATH_DEFINITIONS - Compiler switches required for using libfixmath

find_package(PkgConfig)
pkg_check_modules(PC_LIBFIXMATH QUIET libfixmath)
set(LIBFIXMATH_DEFINITIONS ${PC_LIBFIXMATH_CFLAGS_OTHER})

find_path(LIBFIXMATH_INCLUDE_DIR
        NAMES fix16.h fix16.hpp
        HINTS ${PC_LIBFIXMATH_INCLUDEDIR} ${PC_LIBFIXMATH_INCLUDE_DIRS}
        PATH_SUFFIXES libfixmath)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBFIXMATH_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(libfixmath  DEFAULT_MSG
        LIBFIXMATH_INCLUDE_DIR)

mark_as_advanced(LIBFIXMATH_INCLUDE_DIR)

set(LIBFIXMATH_INCLUDE_DIRS ${LIBFIXMATH_INCLUDE_DIR} )