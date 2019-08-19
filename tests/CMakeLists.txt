
add_executable( test_vswf_translations_array test_vswf_translations_array.c )
target_link_libraries( test_vswf_translations_array qpms gsl blas lapacke amos m )
target_include_directories ( test_vswf_translations_array PRIVATE .. )

add_executable( test_vswf_translations test_vswf_translations.c )
target_link_libraries( test_vswf_translations qpms gsl blas lapacke amos m )
target_include_directories ( test_vswf_translations PRIVATE .. )

add_executable( test_single_translations_vs_calc single_translations_vs_calc.c )
target_link_libraries( test_single_translations_vs_calc qpms gsl lapacke amos m )
target_include_directories( test_single_translations_vs_calc PRIVATE .. )

add_custom_target( mytests DEPENDS test_single_translations_vs_calc test_vswf_translations test_vswf_translations_array )

add_test( NAME single_vs_array_translation_coeffs COMMAND test_single_translations_vs_calc )
