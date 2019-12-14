from pkg_resources import get_distribution
__version__ = get_distribution('qpms').version

import os as __os
from sys import platform as __platform

import warnings as __warnings

try:
    from .qpms_c import PointGroup, FinitePointGroup, FinitePointGroupElement, Particle, scatsystem_set_nthreads, ScatteringSystem, ScatteringMatrix, pitau, set_gsl_pythonic_error_handling, pgsl_ignore_error, gamma_inc, lll_reduce
except ImportError as ex:
    if __platform == "linux" or __platform == "linux2":
        if 'LD_LIBRARY_PATH' not in __os.environ.keys():
            __warnings.warn("Environment variable LD_LIBRARY_PATH has not been set. Make it point to a directory where you installed libqpms and run python again")
        else:
            __warnings.warn("Does your LD_LIBRARY_PATH include a directory where you installed libqpms? Check and run python again."
                '\nCurrently, I see LD_LIBRARY_PATH="%s"' % __os.environ['LD_LIBRARY_PATH'])        
    raise ex
from .cyquaternions import CQuat, IRot3
from .cybspec import VSWFNorm, BaseSpec, default_bspec
from .cytmatrices import CTMatrix, TMatrixInterpolator, TMatrixGenerator
from .cytranslations import trans_calculator
from .cymaterials import MaterialInterpolator, EpsMu, LorentzDrudeModel, lorentz_drude, EpsMuGenerator
from .cycommon import dbgmsg_enable, dbgmsg_disable, dbgmsg_active, BesselType, VSWFType
from .cywaves import vswf_single

from .qpms_p import * # maybe don't import automatically in the future (adds around 0.5 s delay)
from .constants import *

# legacy code which brutally slows down the whole package init:
#from .lattices2d import *
#from .hexpoints import *
#from .tmatrices import *

