# Cythonized parts of QPMS here
# -----------------------------

import numpy as np
import cmath
from qpms_cdefs cimport *
from cybspec cimport *
from cycommon import *
from cycommon cimport make_c_string
cimport cython
import enum
import warnings
import os
from libc.stdlib cimport malloc, free, calloc, abort


cdef class MaterialInterpolator:
    '''
    Wrapper over the qpms_permittivity_interpolator_t structure.
    '''

    def __cinit__(self, filename, *args, **kwargs):
        '''Creates a permittivity interpolator.'''
        cdef char *cpath = make_c_string(filename)
        self.interp = qpms_permittivity_interpolator_from_yml(cpath, gsl_interp_cspline)
        if not self.interp: 
            raise IOError("Could not load permittivity data from %s" % filename)
        self.omegamin = qpms_permittivity_interpolator_omega_min(self.interp)
        self.omegamax = qpms_permittivity_interpolator_omega_max(self.interp)

    def __dealloc__(self):
        qpms_permittivity_interpolator_free(self.interp)

    def __call__(self, double freq):
        '''Returns interpolated permittivity, corresponding to a given angular frequency.'''
        if freq < self.omegamin or freq > self.omegamax:
            raise ValueError("Input frequency %g is outside the interpolator domain (%g, %g)."
                    % (freq, self.minomega, self.freqs[self.maxomega]))
        return qpms_permittivity_interpolator_eps_at_omega(self.interp, freq)

    property freq_interval:
        def __get__(self):
            return [self.omegamin, self.omegamax]

