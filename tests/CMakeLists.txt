
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

add_executable( tbeyn tbeyn.c )
target_link_libraries( tbeyn qpms gsl lapacke amos m )
target_include_directories( tbeyn PRIVATE .. )

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
