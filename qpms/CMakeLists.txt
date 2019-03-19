add_library (qpms translations.c tmatrices.c vecprint.c vswf.c wigner.c
	lattices2d.c gaunt.c error.c legendre.c symmetries.c vecprint.c)
use_c99()
target_include_directories (qpms PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
