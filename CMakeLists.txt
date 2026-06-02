cmake_minimum_required(VERSION 3.14)  # FetchContent_Populate / *_SOURCE_DIR
project(libfixkalman C)

# ---------------------------------------------------------------------------
# Build-time configuration. These map directly onto the library's compile
# macros (see settings.h / fixkalman.h) and are applied PUBLIC-ly so that the
# library and every consumer compile with the *same* configuration. Mixing
# them would link inline (header) code against differently-compiled object
# code - silent breakage.
# ---------------------------------------------------------------------------
set(FIXMATRIX_MAX_SIZE "8" CACHE STRING
    "Maximum matrix dimension; must be >= max(#states, #inputs, #measurements)")

option(KALMAN_JOSEPH_FORM    "Use the Joseph-form covariance update"           OFF)
option(KALMAN_TIME_VARYING   "Covariance prediction uses B*Q*B' (control)"     OFF)
option(KALMAN_DISABLE_C      "Drop filter functions WITH control inputs"       OFF)
option(KALMAN_DISABLE_UC     "Drop filter functions WITHOUT control inputs"    OFF)
option(KALMAN_DISABLE_LAMBDA "Drop certainty (lambda) tuning"                  OFF)

option(KALMAN_BUILD_EXAMPLES "Build the gravity example and package test"      ON)

# ---------------------------------------------------------------------------
# Dependencies. libfixmath / libfixmatrix are not on Conan Center or vcpkg
# (the old Bintray packages are gone), so we fetch the upstream sources and
# compile the parts we need directly. libfixmatrix ships no CMake build, so
# FetchContent only *populates* the sources here.
#
# Commits are pinned for reproducible CI; bump deliberately.
# ---------------------------------------------------------------------------
include(FetchContent)

FetchContent_Declare(libfixmath
    GIT_REPOSITORY https://github.com/PetteriAimonen/libfixmath.git
    GIT_TAG        9318ccd9428114d318eb8703a7e40c152b9059c7)
FetchContent_Declare(libfixmatrix
    GIT_REPOSITORY https://github.com/PetteriAimonen/libfixmatrix.git
    GIT_TAG        a8d583e4fe1fa27ba2e0ec55970f73a1b1e6928b)

FetchContent_GetProperties(libfixmath)
if(NOT libfixmath_POPULATED)
    FetchContent_Populate(libfixmath)
endif()
FetchContent_GetProperties(libfixmatrix)
if(NOT libfixmatrix_POPULATED)
    FetchContent_Populate(libfixmatrix)
endif()

set(LIBFIXMATH_DIR   "${libfixmath_SOURCE_DIR}/libfixmath")
set(LIBFIXMATRIX_DIR "${libfixmatrix_SOURCE_DIR}")

# ---------------------------------------------------------------------------
# fixkalman static library (bundles the dependency sources it depends on)
# ---------------------------------------------------------------------------
add_library(fixkalman STATIC
    fixkalman.c
    ${LIBFIXMATRIX_DIR}/fixmatrix.c
    ${LIBFIXMATRIX_DIR}/fixarray.c
    ${LIBFIXMATH_DIR}/fix16.c
    ${LIBFIXMATH_DIR}/fix16_sqrt.c
    ${LIBFIXMATH_DIR}/fix16_str.c
    ${LIBFIXMATH_DIR}/fix16_trig.c
    ${LIBFIXMATH_DIR}/uint32.c
    ${LIBFIXMATH_DIR}/fract32.c)

target_include_directories(fixkalman PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${LIBFIXMATRIX_DIR}
    ${LIBFIXMATH_DIR})

target_compile_definitions(fixkalman PUBLIC
    FIXMATRIX_MAX_SIZE=${FIXMATRIX_MAX_SIZE}
    FIXMATH_NO_CACHE)

foreach(_opt KALMAN_JOSEPH_FORM KALMAN_TIME_VARYING
             KALMAN_DISABLE_C KALMAN_DISABLE_UC KALMAN_DISABLE_LAMBDA)
    if(${_opt})
        target_compile_definitions(fixkalman PUBLIC ${_opt})
    endif()
endforeach()

target_link_libraries(fixkalman PUBLIC m)
set_target_properties(fixkalman PROPERTIES OUTPUT_NAME fixkalman DEBUG_POSTFIX d)

# ---------------------------------------------------------------------------
# Examples + package test. Both use the control-input-free (UC) API, so they
# cannot be built when that API is disabled.
# ---------------------------------------------------------------------------
if(KALMAN_BUILD_EXAMPLES AND NOT KALMAN_DISABLE_UC)
    enable_testing()

    add_executable(gravity_example example_gravity.c)
    target_link_libraries(gravity_example PRIVATE fixkalman)

    add_executable(package_test test_package/main.c)
    target_link_libraries(package_test PRIVATE fixkalman)

    # package_test reports failure via its exit code (works under NDEBUG too);
    # gravity_example asserts internally and returns 0 on success.
    add_test(NAME gravity_example COMMAND gravity_example)
    add_test(NAME package_test    COMMAND package_test)
endif()
