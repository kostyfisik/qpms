Known bugs
===========

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

