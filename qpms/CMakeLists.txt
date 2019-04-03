#packages
find_package(GSL 2.0 REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_definitions(-DLATTICESUMS32)
add_definitions(-DQPMS_VECTORS_NICE_TRANSFORMATIONS)

#includes
set (DIRS ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
include_directories(${DIRS})

add_library (qpms translations.c tmatrices.c vecprint.c vswf.c wigner.c
	lattices2d.c gaunt.c error.c legendre.c symmetries.c vecprint.c
	bessel.c own_zgemm.c)
use_c99()

set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

target_link_libraries (qpms
	gsl
	lapack
	blas
	amos
	)
target_include_directories (qpms PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
