Known bugs
===========

Wrong factor on B coefficient
-----------------------------
(Probably fixed in the "calculator object" versions!)
Under Kristensson normalisation (with CS = -1), my code gives
B(1,0,n,n)/B(1,0,n,-n) == -(2n)! at (x,y,z) = (x,0,0)
(expected plus or minus 1).
A-coefficients seem to behave correctly.

Xu's antinormalisation
----------------------
"Xu's antinormalisation" is broken (most likely in legendre.c and maybe
also in qpms_types.h) – the plane wave test fails and the spherical wave
reconstruction as well (but the translation coefficients match the 
Xu's tables).

Translation coefficients inconsistent
-------------------------------------
The translation coefficients currently do not work correctly except for certain
combinations of additional i-factors (defined by the [AB][NMF][123] macros
in translations.c) and only for certain normalisations.
QPMS_NORMALISATION_KRISTENSSON_CS does not work at all
QPMS_NORMALISATION_NONE_CS does not work at all
QPMS_NORMALISATION_TAYLOR_CS works for the following macros defined:
  AN1 AM0 BN1 BM0 BF0
  AN1 AM0 BN3 BM2 BF2
  AN3 AM2 BN1 BM0 BF0
  AN3 AM2 BN3 BM2 BF2
QPMS_NORMALISATION_TAYLOR works for the following macros defined:
  AN1 AM2 BN1 BM3 BF0
  AN1 AM2 BN3 BM0 BF2
  AN3 AM0 BN1 BM2 BF0
  AN3 AM0 BN3 BM0 BF2

The default behaviour is now that the QPMS_NORMALISATION_TAYLOR_CS works.

Longitudinal waves
------------------
Plane wave decompositions gives wrong value on the longitudinal part.
The implementation of the L coefficients OR the longitudinal waves
is thus probably wrong.

Scattering result asymmetry
---------------------------
The lattice scattering code (such as finitesqlatzsym-scattery.py) produces
asymmetric results where one should not get them, due to the system symmetry.

It seems that the asymmetry appears mostly in the y-direction (i.e.
for example the scattering/absorption cross section at k = (kx, ky, kz)
is not exactly the same as k = (kx, -ky, kz).

What has been checked (hopefully):
 - The flip operators for electric waves
 - Some weird kind of rounding or other numerical error depending on
   the position ordering of the matrix (randomized indices give
   the same asymmetry).

What has not been checked (so well):
 - The x,y,z-flip operators for magnetic waves (i.e. are they really 
   supposet to bejust the
   same as the operators for electric waves, just with opposite sign?) 
 - zplane_pq_y
 - the translation operators


Singular value asymmetry
------------------------
Similar as the scattering result asymmetry, although not necessarily 
only in the y-direction?

Overflows etc.
--------------
Assertion failed in gaunt_xu for test_vswf_translations.c and high values of LMAX
(LMAX=25)


