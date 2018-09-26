from setuptools import setup#, Extension
from Cython.Distutils import build_ext
from distutils.extension import Extension
# setuptools DWIM monkey-patch madness
# http://mail.python.org/pipermail/distutils-sig/2007-September/thread.html#8204
#import sys
#if 'setuptools.extension' in sys.modules:
#        m = sys.modules['setuptools.extension']
#        m.Extension.__dict__ = m._Extension.__dict__

# TODO CHECK THIS OUT http://stackoverflow.com/questions/4056657/what-is-the-easiest-way-to-make-an-optional-c-extension-for-a-python-package
# also this: https://docs.python.org/2/extending/building.html
import os

print("You might want to add additional library path to LD_LIBRARY_PATH (especially if you are not using"
        " GNU GSL in your system library path) and if import fails. ")
if("LD_LIBRARY_PATH" in os.environ):
    print(os.environ['LD_LIBRARY_PATH'].split(':'))

qpms_c = Extension('qpms_c',
        sources = ['qpms/qpms_c.pyx', #'qpms/hexpoints_c.pyx',
            'qpms/gaunt.c',#'qpms/gaunt.h','qpms/vectors.h','qpms/translations.h',
            # FIXME http://stackoverflow.com/questions/4259170/python-setup-script-extensions-how-do-you-include-a-h-file
            'qpms/translations.c'],
        extra_compile_args=['-std=c99','-ggdb', '-O0',
            '-DQPMS_COMPILE_PYTHON_EXTENSIONS', # this is required
            #'-DQPMS_USE_OMP',
            '-DDISABLE_NDEBUG', # uncomment to enable assertions in the modules
            #'-fopenmp',
            ],
        libraries=['gsl', 'blas', 'gslcblas', #'omp'
            # TODO resolve the problem with openblas (missing gotoblas symbol) and preferable use other blas library
	],
        runtime_library_dirs=os.environ['LD_LIBRARY_PATH'].split(':') if 'LD_LIBRARY_PATH' in os.environ else []
        )

setup(name='qpms',
        version = "0.2.992",
        packages=['qpms'],
#        setup_requires=['setuptools_cython'],
        install_requires=['cython>=0.21','quaternion','spherical_functions','py_gmm'],
        dependency_links=['https://github.com/texnokrates/py_gmm','https://github.com/moble/quaternion','https://github.com/moble/spherical_functions'],
        ext_modules=[qpms_c],
        cmdclass = {'build_ext': build_ext},
        )
        

