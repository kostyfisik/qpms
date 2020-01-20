#packages
find_package(GSL 2.0 REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# disable an annoying warning that gives false positives probably due to a bug in gcc
# and other not very relevant warnings
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-int-in-bool-context -Wno-comment")

#includes
set (DIRS ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
include_directories(${DIRS})

add_library (qpms SHARED translations.c tmatrices.c vecprint.c vswf.c wigner.c ewald.c
	ewaldsf.c pointgroups.c latticegens.c
	lattices2d.c gaunt.c error.c legendre.c symmetries.c vecprint.c
	bessel.c own_zgemm.c parsing.c scatsystem.c materials.c drudeparam_data.c
	lll.c beyn.c
	)
use_c99()

set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

target_link_libraries (qpms
	gsl
	lapack
	blas
	amos
	)
target_include_directories (qpms PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(qpms PRIVATE -Wall -Wno-return-type -Wno-unused-variable -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-label)
target_compile_definitions(qpms PRIVATE LATTICESUMS32 QPMS_VECTORS_NICE_TRANSFORMATIONS
	EWALD_AUTO_CUTOFF
	)

install(TARGETS qpms
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
