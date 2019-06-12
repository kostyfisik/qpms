QPMS README
===========

QPMS is a toolkit for frequency-domain simulations of photonic systems
consisting of compact objects (particles) inside a homogeneous medium. Scattering
properties of the individual particles are described by their T-matrices
(which can be obtained e.g. with the `scuff-tmatrix` tool from 
the [SCUFF-EM] suite).

QPMS handles the multiple scattering of electromagnetic radiation between 
the particles. The system can consist either of a finite number of particles
or an infinite number of periodically arranged lattices (with finite number
of particles in a single unit cell).

Features
========

Finite systems
--------------
 * Computing multipole excitations *and fields (TODO)* scattered from nanoparticle
   clusters illuminated by plane, spherical or *cylindrical (TODO)* waves.
 * Finding eigenmodes.
 * *Calculating cross sections (TODO).*
 * Reducing numerical complexity of the computations by exploiting
   symmetries of the cluster (decomposition to irreducible representations).

Infinite systems (lattices)
---------------------------
 * 2D-periodic systems supported. (TODO 1D and 3D.)
 * *Calculation of transmission and reflection properties (TODO).*
 * Finding eigenmodes and calculating dispersion relations.
 * *Calculation of far-field radiation patterns of an excited array (TODO).*
 * Reducing numerical complexity of the computations by exploiting
   symmetries of the lattice (decomposition to irreducible representations).


Installation
============
The package depends on several python modules and GSL (>= 2.0).
The python module dependencies should be installed automatically when running
the installation script. If you have a recent enough OS,
you can get GSL easily from the repositories; on Debian and derivatives,
just run `apt-get install libgsl-dev` under root. Alternatively,
you can [get the source and compile it yourself][GSL].

You also need a fresh enough version of [cmake][].

After GSL is installed, you can install qpms to your local python library using::

```
  cmake .
  make amos
  python3 setup.py install --user
```

If GSL is not installed the standard library path on your system, you might 
need to pass it to the installation script using the
`LIBRARY_PATH` and `LD_LIBRARY_PATH` environment
variables.

Special care has often be taken when installing QPMS in cluster environments.
Specific installation instructions for Aalto University's Triton cluster
can be found in a [separate document][TRITON-README].

Documentation
=============

Documentation of QPMS is a work in progress. Most of the newer code
is documented using doxygen comments. To build the documentation, just run
`doxygen`
in the root directory; the documentation will then be found in 
`docs/html/index.html`.

Of course, the prerequisite of this is having doxygen installed.
If you don't, you will probably find it easily in your OS's
repositories. On Debian and derivatives, simply run `apt-get install doxygen`
under root.


Tutorials
---------

  * [Infinite system (lattice) tutorial][tutorial-infinite]
  * [Finite system tutorial][tutorial-finite]

[SCUFF-EM]: https://homerreid.github.io/scuff-em-documentation/
[GSL]: https://www.gnu.org/software/gsl/
[cmake]: https://cmake.org
[TRITON-README]: README.Triton.md
[tutorial-finite]: finite_systems.md
[tutorial-infinite]: lattices.md

