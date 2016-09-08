cmake_minimum_required(VERSION 3.6)
project(libfixkalman)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_DEBUG_POSTFIX d)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY         ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY         ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

set(SOURCE_FILES
    fixkalman.c)

set(HEADER_FILES
    compiler.h
    settings.h
    fixkalman.h)

set(GRAVITY_EXAMPLE_FILES
    example_gravity.c)

add_library(libfixkalman STATIC ${HEADER_FILES} ${SOURCE_FILES})
add_executable(gravity_example ${HEADER_FILES} ${GRAVITY_EXAMPLE_FILES})

target_link_libraries(libfixkalman ${CONAN_LIBS})
target_link_libraries(gravity_example ${CONAN_LIBS} libfixkalman)