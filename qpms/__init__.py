from pkg_resources import get_distribution
__version__ = get_distribution('qpms').version

from qpms_c import *
from .qpms_p import *
