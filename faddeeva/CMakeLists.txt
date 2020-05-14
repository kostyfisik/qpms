cmake_minimum_required(VERSION 3.0)
include(GNUInstallDirs)

project(Faddeeva VERSION 1.0 LANGUAGES C)

option(FADDEEVA_BUILD_STATIC "Build Faddeeva as static library" OFF)

if (FADDEEVA_BUILD_STATIC)
	add_library(Faddeeva STATIC Faddeeva.h Faddeeva.c Faddeeva.cc)
else (FADDEEVA_BUILD_STATIC)
	add_library(Faddeeva SHARED Faddeeva.c)
	set_target_properties(Faddeeva PROPERTIES PUBLIC_HEADER "Faddeeva.h")
	install(TARGETS Faddeeva
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif (FADDEEVA_BUILD_STATIC)
target_include_directories(Faddeeva PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	

