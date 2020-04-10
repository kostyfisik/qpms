cmake_minimum_required(VERSION 3.0.2)

option(QPMS_USE_FORTRAN_AMOS "Use the original AMOS Fortran libraries instead of the C ones" OFF)

if (QPMS_USE_FORTRAN_AMOS)
include(CMakeAddFortranSubdirectory)
endif (QPMS_USE_FORTRAN_AMOS)
include(version.cmake)
include(GNUInstallDirs)

project (QPMS)



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

if (QPMS_USE_FORTRAN_AMOS)
	cmake_add_fortran_subdirectory (amos
		PROJECT amos
		LIBRARIES amos
		NO_EXTERNAL_INSTALL)
	set(QPMS_AMOSLIB amos)
else (QPMS_USE_FORTRAN_AMOS)
	set(CAMOS_BUILD_STATIC ON)
	add_subdirectory (camos)
	set(QPMS_AMOSLIB camos)
endif (QPMS_USE_FORTRAN_AMOS)

add_subdirectory (qpms)


enable_testing()
add_subdirectory (tests EXCLUDE_FROM_ALL)

#add_subdirectory (apps/transop-ewald)
