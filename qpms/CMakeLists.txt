#packages
find_package(GSL 2.0 REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GMP REQUIRED)


#includes
set (DIRS ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
include_directories(${DIRS})

add_library (qpms SHARED translations.c tmatrices.c vecprint.c vswf.c wigner.c ewald.c
	ewaldsf.c pointgroups.c latticegens.c
	lattices2d.c gaunt.c error.c legendre.c symmetries.c vecprint.c
	bessel.c own_zgemm.c parsing.c scatsystem.c materials.c drudeparam_data.c
	lll.c polynomials.c
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

target_compile_options(qpms PRIVATE -Wall -Wno-return-type)
target_compile_definitions(qpms PRIVATE LATTICESUMS32 QPMS_VECTORS_NICE_TRANSFORMATIONS
	EWALD_AUTO_CUTOFF
	)

install(TARGETS qpms
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
