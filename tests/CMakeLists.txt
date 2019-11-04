
add_executable( test_vswf_translations_array test_vswf_translations_array.c )
target_link_libraries( test_vswf_translations_array qpms gsl blas lapacke amos m )
target_include_directories ( test_vswf_translations_array PRIVATE .. )

add_executable( test_vswf_translations test_vswf_translations.c )
target_link_libraries( test_vswf_translations qpms gsl blas lapacke amos m )
target_include_directories ( test_vswf_translations PRIVATE .. )

add_executable( test_single_translations_vs_calc single_translations_vs_calc.c )
target_link_libraries( test_single_translations_vs_calc qpms gsl lapacke amos m )
target_include_directories( test_single_translations_vs_calc PRIVATE .. )

add_executable( test_scalar_ewald32 test_scalar_ewald32.c )
target_link_libraries( test_scalar_ewald32 qpms gsl lapacke amos m )
target_include_directories( test_scalar_ewald32 PRIVATE .. )

add_executable( kidneycontour kidneycontour.c )
target_link_libraries( kidneycontour qpms gsl lapacke amos m )
target_include_directories( kidneycontour PRIVATE .. )

add_executable( tbeyn tbeyn.c )
target_link_libraries( tbeyn qpms gsl lapacke amos m )
target_include_directories( tbeyn PRIVATE .. )

add_executable( tbeyn2 tbeyn2.c )
target_link_libraries( tbeyn2 qpms gsl lapacke amos m )
target_include_directories( tbeyn2 PRIVATE .. )

add_executable( tbeyn3a tbeyn3.c )
target_link_libraries( tbeyn3a qpms gsl lapacke amos m )
target_include_directories( tbeyn3a PRIVATE .. )
target_compile_definitions( tbeyn3a PRIVATE VARIANTA )

add_executable( tbeyn3a_implus tbeyn3.c )
target_link_libraries( tbeyn3a_implus qpms gsl lapacke amos m )
target_include_directories( tbeyn3a_implus PRIVATE .. )
target_compile_definitions( tbeyn3a_implus PRIVATE VARIANTA IMPLUS )

add_executable( tbeyn3a_kidney tbeyn3.c )
target_link_libraries( tbeyn3a_kidney qpms gsl lapacke amos m )
target_include_directories( tbeyn3a_kidney PRIVATE .. )
target_compile_definitions( tbeyn3a_kidney PRIVATE VARIANTA IMPLUS_KIDNEY )

add_executable( tbeyn3b tbeyn3.c )
target_link_libraries( tbeyn3b qpms gsl lapacke amos m )
target_include_directories( tbeyn3b PRIVATE .. )
target_compile_definitions( tbeyn3b PRIVATE VARIANTB RXSMALL )

add_executable( tbeyn3bfail tbeyn3.c )
target_link_libraries( tbeyn3bfail qpms gsl lapacke amos m )
target_include_directories( tbeyn3bfail PRIVATE .. )
target_compile_definitions( tbeyn3bfail PRIVATE VARIANTB )

add_executable( tbeyn3c tbeyn3.c )
target_link_libraries( tbeyn3c qpms gsl lapacke amos m )
target_include_directories( tbeyn3c PRIVATE .. )
target_compile_definitions( tbeyn3c PRIVATE VARIANTC RXSMALL )

add_executable( tbeyn3d tbeyn3.c )
target_link_libraries( tbeyn3d qpms gsl lapacke amos m )
target_include_directories( tbeyn3d PRIVATE .. )
target_compile_definitions( tbeyn3d PRIVATE VARIANTD RXSMALL )

add_executable( tbeyn3e tbeyn3.c )
target_link_libraries( tbeyn3e qpms gsl lapacke amos m )
target_include_directories( tbeyn3e PRIVATE .. )
target_compile_definitions( tbeyn3e PRIVATE VARIANTE RXSMALL )

add_executable( tbeyn3f tbeyn3.c )
target_link_libraries( tbeyn3f qpms gsl lapacke amos m )
target_include_directories( tbeyn3f PRIVATE .. )
target_compile_definitions( tbeyn3f PRIVATE VARIANTF )

add_executable( tbeyn_gsl tbeyn_gsl.c )
target_link_libraries( tbeyn_gsl qpms gsl lapacke amos m )
target_include_directories( tbeyn_gsl PRIVATE .. )

add_custom_target( mytests DEPENDS test_single_translations_vs_calc test_vswf_translations test_vswf_translations_array tbeyn )

add_test( NAME single_vs_array_translation_coeffs COMMAND test_single_translations_vs_calc )

add_test( NAME scalar_ewald32_square_realk1 COMMAND test_scalar_ewald32 
#	lMax b1.x b1.y b2.x b2.y wavenum.real wavenum.imag k.x k.y particle_shift.x particle_shift.y csphase rtol atol maxR maxK eta1 ...
	3    1    0    0    1    2.3          0            2.7   1              0.5           0.1325      -1 1e-8 1e-10  20  160  1.5   1.6 2.5 2.6
	)

add_test( NAME scalar_ewald32_square_realk2 COMMAND test_scalar_ewald32 
#	lMax b1.x b1.y b2.x b2.y wavenum.real wavenum.imag k.x k.y particle_shift.x particle_shift.y csphase rtol atol maxR maxK eta1 ...
	3    1    0    0    1    2.3          0            2.7   1              0.0           0.0      -1 1e-8 1e-10  20  160  1.5   1.6 2.5 2.6
	)

add_test( NAME scalar_ewald32_square_cplxk_positive_Im1 COMMAND test_scalar_ewald32 
#	lMax b1.x b1.y b2.x b2.y wavenum.real wavenum.imag k.x k.y particle_shift.x particle_shift.y csphase rtol atol maxR maxK eta1 ...
	3    1    0    0    1    2.3          0.1          2.7   1              0.5           0.1325      -1 1e-8 1e-10  20  160  1.5   1.6 2.5 2.6
	)

add_test( NAME scalar_ewald32_square_cplxk_negative_Im1 COMMAND test_scalar_ewald32 
#	lMax b1.x b1.y b2.x b2.y wavenum.real wavenum.imag k.x k.y particle_shift.x particle_shift.y csphase rtol atol maxR maxK eta1 ...
	3    1    0    0    1    2.3          -0.1          2.7   1              0.5           0.1325      -1 1e-8 1e-10  20  160  1.5   1.6 2.5 2.6
	)

add_test( NAME scalar_ewald32_square_cplxk_positive_Im2 COMMAND test_scalar_ewald32 
#	lMax b1.x b1.y b2.x b2.y wavenum.real wavenum.imag k.x k.y particle_shift.x particle_shift.y csphase rtol atol maxR maxK eta1 ...
	3    1    0    0    1    2.3          0.1          2.7   1              0.0           0.0         -1 1e-8 1e-10  20  160  1.5   1.6 2.5 2.6
	)

add_test( NAME scalar_ewald32_square_cplxk_negative_Im2 COMMAND test_scalar_ewald32 
#	lMax b1.x b1.y b2.x b2.y wavenum.real wavenum.imag k.x k.y particle_shift.x particle_shift.y csphase rtol atol maxR maxK eta1 ...
	3    1    0    0    1    2.3          -0.1          2.7   1              0.0           0.0        -1 1e-8 1e-10  20  160  1.5   1.6 2.5 2.6
	)

add_test( NAME scalar_ewald32_square_zerok COMMAND test_scalar_ewald32 
#	lMax b1.x b1.y b2.x b2.y wavenum.real wavenum.imag k.x k.y particle_shift.x particle_shift.y csphase rtol atol maxR maxK eta1 ...
	3    1    0    0    1    2.3          0            0   0                0.5           0.1325      -1 1e-8 1e-10  40  160  1.5   1.6 2.5 2.6
	)
