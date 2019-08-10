from pkg_resources import get_distribution
__version__ = get_distribution('qpms').version

from .qpms_c import *
from .qpms_p import *
from .cyquaternions import CQuat, IRot3
from .cybspec import VSWFNorm, BaseSpec
from .lattices2d import *
from .hexpoints import *
from .tmatrices import *

