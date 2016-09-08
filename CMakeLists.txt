cmake_minimum_required(VERSION 3.6)
project(libfixkalman)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    compiler.h
    example_gravity.c
    fixkalman.c
    fixkalman.h)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

find_package (libfixmath REQUIRED)
find_package (libfixmatrix REQUIRED)

if (LIBFIXMATH_FOUND)
    include_directories(${LIBFIXMATH_INCLUDE_DIRS})
    set(SOURCE_FILES
            ${SOURCE_FILES}
            ${LIBFIXMATH_INCLUDE_DIR}/fix16
            ${LIBFIXMATH_INCLUDE_DIR}/fix16_sqrt)
endif (LIBFIXMATH_FOUND)

if (LIBFIXMATRIX_FOUND)
    include_directories(${LIBFIXMATRIX_INCLUDE_DIRS})
    set(SOURCE_FILES
            ${SOURCE_FILES}
            ${LIBFIXMATRIX_INCLUDE_DIR}/fixarray
            ${LIBFIXMATRIX_INCLUDE_DIR}/fixmatrix)
endif (LIBFIXMATRIX_FOUND)

add_executable(libfixkalman ${SOURCE_FILES})
