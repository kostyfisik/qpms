# Kde co je
## Staré věci
hexpoints.py
legacy.py
qpms_p.py (až na změny souřadnic???)

## Nové věci
lattices2d.py
scripts_common.py
timetrack.py
tmatrices.py
types.py
svwf.c

## Smíšené / v přepisu
scattering.py
qpms_c.pyx

## ???
hexpoints_c.pyx

# hexpoints.py
Asi hlavně starý kód pro vytváření trojúhelníkových a hexagonálních mřížek
a počítání (a ukládání) interakčních matic

## funkce
generate_trianglepoints
generate_trianglepoints_hexcomplement
hexlattice_precalc_AB_save
hexlattice_precalc_AB_loadunwrap
hexlattice_get_AB
hexlattice_zsym_getSVD
	 
# hexpoints_c.pyx
Obsahuje pouze jedinou funkci (která je i v hexpoints.py). 
Používá se tohle vůbec někde?
## funkce
hexlattice_zsym_getSVD

# lattices2d.py
Nový kód, manipulace s basemi, vytváření mřížek atd.

## třídy
LatticeType(Enum)

## funkce
reduceBasisSingle
shortestBase3
shortestbase46
is_obtuse
classifyLatticeSingle
range2D
generateLattice
generateLatticeDisk
cellCornersWS
cutWS
filledWS
filledWS2
change_basis

# legacy.py
Stařičký kód

## funkce
q_max
a_q
Ã
B̃
G_Mie_scat_precalc_cart
G_Mie_scat_cart
scatter_plane_wave
scatter_plane_wave_rectarray
scatter_constmultipole_rectarray
hexlattice_precalc_AB_save2
hexlattice_precalc_AB_save_purepy
hexlattice_precalc_AB_loadunwrap
hexlattice_get_AB

# qpms_p.py
## funkce
cart2sph
sph2cart
sph_loccart2cart
sph_loccart_basis
nelem2lMax
lpy
lpy1
vswf_yr
_sph_zn_1
_sph_zn_2
_sph_zn_3
_sph_zn_4
zJn
π̃_pilim
τ̃_zerolim
τ̃_pilim
get_π̃τ̃_y1
vswf_yr1
zplane_pq_y
plane_pq_y
ε_drude
mie_coefficients
G_Mie_scat_precalc_cart_new
Grr_Delga
Grr_Delga
G0_dip_1
_P
_Q
G0_analytical
G0L_analytical
G0T_analytical
G0_sum_1_slow

# scattering.py
## třídy
Scattering
LatticeScattering (neimplementováno nic, asi zrovna rozepsáno)
Scattering_2D_zsym

# scripts_common.py
## funkce
make_action_sharedlist
add_argparse_k_output_options
add_argparse_unitcell_definitions
add_argparse_infinite_lattice_options
add_argparse_output_options
add_argparse_common_options
arg_preprocess_particles

# timetrack.py

# tmatrices.py
## funkce
WignerD_mm
WignerD_mm_fromvector
WignerD_yy
WignerD_yy_fromvector
xflip_yy
xflip_tyy
xflip_tyty
yflip_yy
yflip_tyy
yflip_tyty
zflip_yy
zflip_tyy
zflip_tyty
parity_yy
loadScuffTMatrices
apply_matrix_left
apply_ndmatrix_left
symz_indexarrays
get_TMatrix_fromspec
perform_tmspect

## třídy
TMatrix

# types.py
## třídy
NormalizationT
BesselT
TMatrixOp
TMatrixSpec
ParticleSpec
LatticeSpec

#qpms_c.pyx
## funkce
get_mn_y
get_mn_y_unsigned
q_max
loop_D_iiiidddii_As_D_lllldddbl a jiné pomocné funkce (pro ufunc)

## třídy
trans_calculator
