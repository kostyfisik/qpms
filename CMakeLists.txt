cmake_minimum_required(VERSION 3.0.2)
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake-scripts")
include(CMakeAddFortranSubdirectory)
include(version.cmake)
include(GNUInstallDirs)

project (QPMS)

set(CMAKE_BUILD_TYPE Debug)

macro(use_c99)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    endif ()
  else ()
    set (CMAKE_C_STANDARD 99)
  endif ()
endmacro(use_c99)

# We need the amos libraries to be compiled with -fPIC
# but at the same time, we do not want to make a .so file,
# so distutils package it with the rest of qpms c/cython lib.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set (QPMS_VERSION_MAJOR 0)
#set (QPMS_VERSION_MINOR 3)
cmake_add_fortran_subdirectory (amos
	PROJECT amos
	LIBRARIES amos
	NO_EXTERNAL_INSTALL)
add_subdirectory (qpms)


enable_testing()
add_subdirectory (tests EXCLUDE_FROM_ALL)

#add_subdirectory (apps/transop-ewald)
