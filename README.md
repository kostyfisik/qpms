QPMS README
===========

QPMS (standing for QPMS Photonic Multiple Scattering) 
is a toolkit for frequency-domain simulations of photonic systems
consisting of compact objects (particles) inside a homogeneous medium. Scattering
properties of the individual particles are described by their T-matrices
(which can be obtained using one of the built-in generators or
 e.g. with the `scuff-tmatrix` tool from 
the [SCUFF-EM] suite).

QPMS handles the multiple scattering of electromagnetic radiation between 
the particles. The system can consist either of a finite number of particles
or an infinite number of periodically arranged lattices (with finite number
of particles in a single unit cell).


Features
========


Finite systems
--------------
 * Computing multipole excitations and fields scattered from nanoparticle
   clusters illuminated by plane, spherical or *cylindrical (TODO)* waves.
 * Finding eigenmodes (optical resonances).
 * Calculating cross sections.
 * Reducing numerical complexity of the computations by exploiting
   symmetries of the cluster (decomposition to irreducible representations).


Infinite systems (lattices)
---------------------------
 * 2D-periodic systems with arbitrary unit cell geometry supported. (TODO 1D and 3D.)
 * Computing multipole excitations and fields scattered from nanoparticle
 * Finding eigenmodes and calculating dispersion relations.
 * Calculation of the scattered fields.
 * *Calculation of total transmission and reflection properties (TODO).*
 * *Reducing numerical complexity of the computations by exploiting
   symmetries of the lattice (decomposition to irreducible representations) (in development).* 


Getting the code
================

The main upstream public repository is located at <https://repo.or.cz/qpms.git>.
Just clone the repository with `git` and proceed to the installation instructions
below.


Installation
============
The package depends on several python modules, a BLAS/LAPACK library with 
the respective C bindings (incl. the `lapacke.h` and `cblas.h` headers;
[OpenBLAS][OpenBLAS] does have it all and is recommended) and GSL (>= 2.0).
The python module dependencies should be installed automatically when running
the installation script. If you have a recent enough OS,
you can get GSL easily from the repositories; on Debian and derivatives,
just run `apt-get install libgsl-dev` under root. Alternatively,
you can [get the source and compile it yourself][GSL].

You also need a fresh enough version of [cmake][].

QPMS uses a C version of the Amos library for calculating Bessel function
from a submodule. Before proceeding with running `cmake`, the submodules
need to be downloaded first (in the QPMS source root directory):

```{.sh}
  git submodule init
  git submodule update
```

After GSL is installed and submodules updated, you can install qpms to your local python library using

```{.sh}
  cmake -DCMAKE_INSTALL_PREFIX=${YOUR_PREFIX} .
  make install
  python3 setup.py install --user
```
Above, replace `${YOUR_PREFIX}` with the path to where you want to install the shared library;
you will also need to make sure that the linker can find it;
on Linux, this means the path `${YOUR_PREFIX}/lib` is included in your
`LIBRARY_PATH` and `LD_LIBRARY_PATH` environment variables. The same applies
to the GSL and OpenBLAS dependencies: they must be installed where the
installation scripts and linker can find them (setting the `C_INCLUDE_PATH` environment
variable might be necessary as well).

Special care might need to be taken when installing QPMS in cluster environments.
Specific installation instructions for Aalto University's Triton cluster
can be found in a [separate document][TRITON-README].


Documentation
=============

Documentation of QPMS is a work in progress. Most of the newer code
is documented using [doxygen][] comments. To build the documentation, just run
`doxygen`
in the root directory; the documentation will then be found in 
`docs/html/index.html`.

Of course, the prerequisite of this is having doxygen installed.
If you don't, you will probably find it easily in your OS's
repositories. On Debian and derivatives, simply run `apt-get install doxygen`
under root.


Tutorials
---------

  * [Finite system tutorial][tutorial-finite]

See also the examples directory in the source repository.

Acknowledgments
================

This software has been developed in the [Quantum Dynamics research group][QD],
Aalto University, Finland. If you use the code in your work, please cite 
**M. Nečada and P. Törmä, Multiple-scattering T-matrix simulations for nanophotonics: symmetries and periodic lattices, [arXiv: 2006.12968][lepaper] (2020)**
in your publications, presentations, and similar. 

Please also have a look at other publications by the group 
(google scholar Päivi Törmä), they may be useful for your work as well. 


Bug reports
===========

If you believe that some parts of QPMS behave incorrectly, please mail
a bug report to <marek@necada.org>. To ensure that your message is not
considered spam, please start the subject line with `QPMS`.

If you were able to fix a bug yourself, please include the patch as well,
see below.


Contributions
=============

Contributions to QPMS are welcome, be it bug fixes, improvements to the
documentation, code quality, or new features.

You can send patches prepared using the 
[`git format-patch`](https://git-scm.com/docs/git-format-patch) tool
to <marek@necada.org>. 

If you plan to contribute with major changes to the codebase, it is 
recommended to discuss that first (see the contact information below).


Contact & discussion
====================

You can contact the main author e.g. via [e-mail](marek@necada.org) 
or [Telegram](https://t.me/necadam).

You are also warmly welcome to the [QPMS user chat](https://t.me/QPMScattering)
in Telegram!



[SCUFF-EM]: https://homerreid.github.io/scuff-em-documentation/
[OpenBLAS]: https://www.openblas.net/
[GSL]: https://www.gnu.org/software/gsl/
[cmake]: https://cmake.org
[tRITON-README]: README.Triton.md
[tutorial-finite]: finite_systems.md
[tutorial-infinite]: lattices.md
[doxygen]: http://doxygen.nl/
[QD]: https://www.aalto.fi/en/department-of-applied-physics/quantum-dynamics-qd
[lepaper]: https://arxiv.org/abs/2006.12968
