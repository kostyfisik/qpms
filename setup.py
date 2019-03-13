from setuptools import setup#, Extension
from Cython.Build import cythonize, build_ext
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
            'qpms/translations.c',
            'qpms/symmetries.c',
            'qpms/wigner.c',
            'qpms/scatsystem.c',
            'qpms/vswf.c', # FIXME many things from vswf.c are not required by this module, but they have many dependencies (following in this list); maybe I want to move all the  "basespec stuff"
            'qpms/legendre.c',
            'qpms/tmatrix_io.c',
            'qpms/error.c'],
        extra_compile_args=['-std=c99','-ggdb', '-O0',
            '-DQPMS_COMPILE_PYTHON_EXTENSIONS', # this is required
            #'-DQPMS_USE_OMP',
            '-DDISABLE_NDEBUG', # uncomment to enable assertions in the modules
            #'-fopenmp',
            ],
        libraries=['gsl', 'lapacke', 'blas', 'gslcblas', #'omp'
            # TODO resolve the problem with openblas (missing gotoblas symbol) and preferable use other blas library
	],
        runtime_library_dirs=os.environ['LD_LIBRARY_PATH'].split(':') if 'LD_LIBRARY_PATH' in os.environ else []
        )

setup(name='qpms',
        version = "0.2.995",
        packages=['qpms'],
        setup_requires=['cython>=0.28',],
        install_requires=['cython>=0.28','quaternion','spherical_functions','scipy>=0.18.0', 'sympy>=1.2'],
        dependency_links=['https://github.com/moble/quaternion/archive/v2.0.tar.gz','https://github.com/moble/spherical_functions/archive/master.zip'],
        ext_modules=cythonize([qpms_c], include_path=['qpms'], gdb_debug=True),
        cmdclass = {'build_ext': build_ext},
        )
        

