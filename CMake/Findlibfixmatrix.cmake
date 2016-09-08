# - Try to find libfixmatrix
# Once done this will define
#  LIBFIXMATRIX_FOUND - System has libfixmatrix
#  LIBFIXMATRIX_INCLUDE_DIRS - The libfixmatrix include directories
#  LIBFIXMATRIX_DEFINITIONS - Compiler switches required for using libfixmatrix

find_package(PkgConfig)
pkg_check_modules(PC_LIBFIXMATRIX QUIET libfixmatrix)
set(LIBFIXMATRIX_DEFINITIONS ${PC_LIBFIXMATRIX_CFLAGS_OTHER})

find_path(LIBFIXMATRIX_INCLUDE_DIR fixmatrix.h
        HINTS ${PC_LIBFIXMATRIX_INCLUDEDIR} ${PC_LIBFIXMATRIX_INCLUDE_DIRS}
        PATH_SUFFIXES libfixmatrix )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBFIXMATRIX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(libfixmatrix  DEFAULT_MSG
        LIBFIXMATRIX_INCLUDE_DIR)

mark_as_advanced(LIBFIXMATRIX_INCLUDE_DIR)

set(LIBFIXMATRIX_INCLUDE_DIRS ${LIBFIXMATRIX_INCLUDE_DIR} )