from setuptools import setup#, Extension
from Cython.Distutils import build_ext
from distutils.extension import Extension
# setuptools DWIM monkey-patch madness
# http://mail.python.org/pipermail/distutils-sig/2007-September/thread.html#8204
#import sys
#if 'setuptools.extension' in sys.modules:
#        m = sys.modules['setuptools.extension']
#        m.Extension.__dict__ = m._Extension.__dict__

qpms_c = Extension('qpms_c',
        sources = ['qpms/qpms_c.pyx'])

setup(name='qpms',
        version = "0.1.8",
        packages=['qpms'],
#        setup_requires=['setuptools_cython'],
        install_requires=['cython>=0.21','quaternion','spherical_functions','py_gmm'],
        dependency_links=['https://github.com/texnokrates/py_gmm','https://github.com/moble/quaternion','https://github.com/moble/spherical_functions'],
        ext_modules=[qpms_c],
        cmdclass = {'build_ext': build_ext},
        )
        

