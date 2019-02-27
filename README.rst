Quantum photonic multiple scattering
====================================

TODO description

Installation
============
The package depends on several python modules and GSL (>= 2.0).
The python module dependencies should be installed automatically when running
the installation script. If you have a recent enough OS,
you can get GSL easily from the repositories; on Debian and derivatives,
just run ``apt-get install libgsl-dev`` under root. Alternatively,
you can `get the source 
<https://www.gnu.org/software/gsl/>`_ get the source and compile it yourself.

After GSL is installed, you can install qpms to your local python library using::

  python3 setup.py install --user

If GSL is not installed the standard library path on your system, you might 
need to pass it to the installation script using the LD_LIBRARY_PATH environment
variable.

Documentation
=============

Documentation of QPMS is a work in progress. Most of the newer code
is documented using doxygen comments. To build the documentation, just run::

  doxygen

in the root directory; the documentation will then be found in 
``docs/index.html``.
