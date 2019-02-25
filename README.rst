Quantum photonic multiple scattering
====================================

TODO description

Installation
============
The package depends on numpy, scipy, cython and GSL (>= 2.0).
The first three can be obtained by pip. If you have a recent enough OS,
you can get GSL easily from the repositories; on Debian and derivatives,
just run `apt-get install libgsl-dev` under root. Alternatively,
you can `get the source 
<https://www.gnu.org/software/gsl/>`_ get the source and compile it yourself.

After all dependencies are installed, install qpms to your local python library using::

  python3 setup.py install --user


Easiest installation ever 
=========================
(Just skip those you have already installed.)

::

  pip3 install --user numpy
  pip3 install --user scipy
  pip3 install --user cython
  pip3 install --user git+https://github.com/moble/quaternion.git
  pip3 install --user git+https://github.com/moble/spherical_functions.git
  python3 setup.py install --user


  
Documentation
=============

Documentation of QPMS is a work in progress. Most of the newer code
is documented using doxygen comments. To build the documentation, just run::

  doxygen

in the root directory; the documentation will then be found in 
``docs/index.html``.
