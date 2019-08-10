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
# AND THIS https://groups.google.com/forum/#!topic/cython-users/GAAPYb2X304
# also this: https://docs.python.org/2/extending/building.html
import os

#print("You might want to add additional library path to LD_LIBRARY_PATH (especially if you are not using"
#        " GNU GSL in your system library path) and if import fails. ")
#if("LD_LIBRARY_PATH" in os.environ):
#    print(os.environ['LD_LIBRARY_PATH'].split(':'))



'''
amos_sources = [
        'amos/d1mach.f',
        'amos/dgamln.f',
        'amos/i1mach.f',
        'amos/xerror.f',
        'amos/zabs.f',
        'amos/zacai.f',
        'amos/zacon.f',
        'amos/zairy.f',
        'amos/zasyi.f',
        'amos/zbesh.f',
        'amos/zbesi.f',
        'amos/zbesj.f',
        'amos/zbesk.f',
        'amos/zbesy.f',
        'amos/zbinu.f',
        'amos/zbiry.f',
        'amos/zbknu.f',
        'amos/zbuni.f',
        'amos/zbunk.f',
        'amos/zdiv.f',
        'amos/zexp.f',
        'amos/zkscl.f',
        'amos/zlog.f',
        'amos/zmlri.f',
        'amos/zmlt.f',
        'amos/zrati.f',
        'amos/zs1s2.f',
        'amos/zseri.f',
        'amos/zshch.f',
        'amos/zsqrt.f',
        'amos/zuchk.f',
        'amos/zunhj.f',
        'amos/zuni1.f',
        'amos/zuni2.f',
        'amos/zunik.f',
        'amos/zunk1.f',
        'amos/zunk2.f',
        'amos/zuoik.f',
        'amos/zwrsk.f',
        ]

libqpms_sources = [
            #'qpms/hexpoints_c.pyx',
            'qpms/gaunt.c',#'qpms/gaunt.h','qpms/vectors.h','qpms/translations.h',
            # FIXME http://stackoverflow.com/questions/4259170/python-setup-script-extensions-how-do-you-include-a-h-file
            'qpms/translations.c',
            'qpms/symmetries.c',
            'qpms/wigner.c',
            'qpms/scatsystem.c',
            'qpms/vswf.c', # FIXME many things from vswf.c are not required by this module, but they have many dependencies (following in this list); maybe I want to move all the  "basespec stuff"
            'qpms/legendre.c',
            'qpms/tmatrices.c',
            'qpms/materials.c',
            'qpms/error.c',
            'qpms/bessel.c',
            'qpms/own_zgemm.c',
            'qpms/pointgroups.c',
            ]
'''

cycommon = Extension('qpms.cycommon',
        sources = ['qpms/cycommon.pyx'],
        extra_link_args=['qpms/libqpms.a'],
        libraries=['gsl', 'lapacke', 'blas', 'gslcblas', 'pthread',]
        )
cytranslations = Extension('qpms.cytranslations',
        sources = ['qpms/cytranslations.pyx',
            'qpms/translations_python.c',
            ],
        extra_compile_args=['-std=c99',
            '-DQPMS_COMPILE_PYTHON_EXTENSIONS', # This is needed to enable it in translations.h
            ],
        extra_link_args=['qpms/libqpms.a', 'amos/libamos.a'],
        libraries=['gsl', 'lapacke', 'blas', 'gslcblas', 'pthread',]
        )
cybspec = Extension('qpms.cybspec',
        sources = ['qpms/cybspec.pyx'],
        extra_link_args=['qpms/libqpms.a'],
        libraries=['gsl', 'lapacke', 'blas', 'gslcblas', 'pthread',]
        )
cymaterials = Extension('qpms.cymaterials',
        sources = ['qpms/cymaterials.pyx'],
        extra_link_args=['qpms/libqpms.a'],
        libraries=['gsl', 'lapacke', 'blas', 'gslcblas', 'pthread',]
        )
cyquaternions = Extension('qpms.cyquaternions',
        sources = ['qpms/cyquaternions.pyx'],
        extra_link_args=['amos/libamos.a', 'qpms/libqpms.a'],
        libraries=['gsl', 'lapacke', 'blas', 'gslcblas', 'pthread',]
        )

qpms_c = Extension('qpms.qpms_c',
        sources = [
            'qpms/qpms_c.pyx',
        ],
        libraries=['gsl', 'lapacke', 'blas', 'gslcblas', 'pthread', #'omp'
            #('amos', dict(sources=amos_sources) ),
	],
        include_dirs=['amos', 'qpms'],
        extra_link_args=[ 'qpms/libqpms.a','amos/libamos.a', ],
        #runtime_library_dirs=os.environ['LD_LIBRARY_PATH'].split(':') if 'LD_LIBRARY_PATH' in os.environ else [],
        #extra_objects = ['amos/libamos.a'], # FIXME apparently, I would like to eliminate the need to cmake/make first
        )

setup(name='qpms',
        version = "0.2.996",
        packages=['qpms'],
#        libraries = [('amos', {'sources': amos_sources} )],
        setup_requires=['cython>=0.28',],
        install_requires=['cython>=0.28',
            #'quaternion','spherical_functions',
            'scipy>=0.18.0', 'sympy>=1.2'],
        #dependency_links=['https://github.com/moble/quaternion/archive/v2.0.tar.gz','https://github.com/moble/spherical_functions/archive/master.zip'],
        ext_modules=cythonize([qpms_c, cytranslations, cycommon, cyquaternions, cybspec, cymaterials], include_path=['qpms', 'amos'], gdb_debug=True),
        cmdclass = {'build_ext': build_ext},
        zip_safe=False
        )
        

