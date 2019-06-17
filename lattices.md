Using QPMS library for finding modes of 2D-periodic systems
===========================================================

Calculating modes of infinite 2D arrays is now done 
in several steps (assuming the T-matrices have already
been obtained using `scuff-tmatrix` or can be obtained
from Lorenz-Mie solution (spherical particles)):

 1. Sampling the *k*, *ω* space.
 2. Pre-calculating the
    Ewald-summed translation operators.
 3. For each *k*, *ω* pair, build the LHS operator
    for the scattering problem (TODO reference), optionally decomposed
    into suitable irreducible representation subspaces.
 4. Evaluating the singular values and finding their minima.

The steps above may (and will) change as more user-friendly interface
will be developed.


Preparation: compile the `ew_gen_kin` utility
---------------------------------------------

This will change, but at this point, the lattice-summed
translation operators are computed using the `ew_gen_kin`
utility located in the `qpms/apps` directory. It has to be built
manually like this:

```bash
  cd qpms/apps
  c99 -o ew_gen_kin -Wall -I ../.. -I ../../amos/ -O2 -ggdb -DQPMS_VECTORS_NICE_TRANSFORMATIONS -DLATTICESUMS32 2dlattice_ewald.c ../translations.c ../ewald.c ../ewaldsf.c ../gaunt.c ../lattices2d.c ../latticegens.c ../bessel.c -lgsl -lm -lblas ../../amos/libamos.a -lgfortran ../error.c
``` 

Step 1: Sampling the *k*, *ω* space
--------------------------------------

`ew_gen_kin` expects a list of (*k_x*, *k_y*)
pairs on standard input (separated by whitespaces),
the rest is specified via command line arguments.

So if we want to examine the line between the Г point and the point
\f$ k = (0, 10^5\,\mathrm{m}^{-1}) \f$, we can generate an input
running
```bash
  for ky in $(seq 0 1e3 1e5); do
    echo 0 $ky >> klist
  done
```

It also make sense to pre-generate the list of *ω* values,
e.g. 
```bash
  seq 6.900 0.002 7.3 | sed -e 's/,/./g' > omegalist
```


Step 2: Pre-calculating the translation operators
-------------------------------------------------

`ew_gen_kin` currently uses command-line arguments in
an atrocious way with a hard-coded order:
```
  ew_gen_kin outfile b1.x b1.y b2.x b2.y lMax scuffomega refindex npart part0.x part0.y [part1.x part1.y [...]]
```
where `outfile` specifies the path to the output, `b1` and `b2` are the
direct lattice vectors, `lMax` is the multipole degree cutoff,
`scuffomega` is the frequency in the units used by `scuff-tmatrix`
(TODO specify), `refindex` is the refractive index of the background
medium, `npart` number of particles in the unit cell, and `partN` are 
the positions of these particles inside the unit cell.

Assuming we have the `ew_gen_kin` binary in our `${PATH}`, we can
now run e.g.
```bash
  for omega in $(cat omegalist); do
    ew_gen_kin $omega 621e-9 0 0 571e-9 3 w_$omega 1.52 1 0 0 < klist
  done
```
This pre-calculates the translation operators for a simple (one particle per unit cell)
621 nm × 571 nm rectangular lattice inside a medium with refractive index 1.52, 
up to the octupole (`lMax` = 3) order, yielding one file per frequency. 
This can take some time and
it makes sense to run a parallelised `for`-loop instead; this is a stupid but working
way to do it in bash:
```bash
  N=4   # number of parallel processes
  for omega in $(cat omegalist); do
    ((i=i%N)); ((i++==0)) && wait
    ew_gen_kin $omega 621e-9 0 0 571e-9 3 w_$omega 1.52 1 0 0 < klist
    echo $omega   # optional, to follow progress
  done
```

When this is done, we convert all the text output files into 
numpy's binary format in order to speed up loading in the following steps.
This is done using the processWfiles_sortnames.py script located in the 
`misc` directory. Its usage pattern is
```
  processWfiles_sortnames.py npart dest src1 [src2 ...]
```
where `npart` is the number of particles in the unit cell, `dest`
is the destination path for the converted data (this will be 
a directory), and the remaining arguments are paths to the
files generated by `ew_gen_kin`. In the case above, one could use
```
  processWfiles_sortnames.py 1 all w_*
```
which would create a directory named `all` containing several
.npy files.


Steps 3, 4
----------

TODO. For the time being, see e.g. the `SaraRect/dispersions.ipynb` jupyter notebook
from the `qpms_ipynotebooks` repository
for the remaining steps.
