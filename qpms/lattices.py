'''
Object oriented approach for the classical multiple scattering problem.
'''

import numpy as np
from qpms_c import * # TODO be explicit about what is imported
from .qpms_p import * # TODO be explicit about what is imported

class Scatterers(object):
    '''
    This is the most general class for a system of scatterers
    in a non-lossy homogeneous background
    to be solved with the multiple_scattering method. The scatterers,
    as long as they comply with the disjoint circumscribed sphere
    hypothesis, can each have any position in the 3D space and
    any T-matrix.

    Note that this object describes the scattering problem only for
    a single given frequency, as the T-matrices and wavelenght
    otherwise differ and all the computationally demanding
    parts have to be done for each frequency. However,
    the object can be recycled for many incident field shapes
    at the given frequency.


    '''

