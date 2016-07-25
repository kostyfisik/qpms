Quantum photonic multiple scattering
====================================

TODO description

Installation
============
The package depends on numpy, scipy, cython and customized version of py_gmm.
The first three can be obtained by pip, the last one can be obtained from github:

git clone --branch standalone_mie  https://github.com/texnokrates/py_gmm.git 

After all dependencies are installed, install qpms to your local python library using

python3 setup.py install --user


Easiest installation ever 
=========================
(Just skip those you have already installed.)

pip3 install --user numpy
pip3 install --user scipy
pip3 install --user cython
pip3 install --user git+https://github.com/moble/quaternion.git
pip3 install --user git+https://github.com/moble/spherical_functions.git
pip3 install --user git+https://github.com/texnokrates/py_gmm.git@standalone_mie
python3 setup.py install --user
