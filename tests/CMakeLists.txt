
add_executable( test_vswf_translations_array test_vswf_translations_array.c )
target_link_libraries( test_vswf_translations_array qpms gsl blas lapacke amos m )
target_include_directories ( test_vswf_translations_array PRIVATE .. )

add_executable( test_vswf_translations test_vswf_translations.c )
target_link_libraries( test_vswf_translations qpms gsl blas lapacke amos m )
target_include_directories ( test_vswf_translations PRIVATE .. )
