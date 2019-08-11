# Cythonized parts of QPMS here
# -----------------------------

import numpy as np
import cmath
from .qpms_cdefs cimport qpms_permittivity_interpolator_from_yml, qpms_permittivity_interpolator_free, qpms_permittivity_interpolator_omega_min, qpms_permittivity_interpolator_omega_max, gsl_interp_type, qpms_permittivity_interpolator_t, gsl_interp_cspline, qpms_permittivity_interpolator_eps_at_omega, qpms_epsmu_const_g, qpms_permittivity_interpolator_epsmu_g, qpms_epsmu_const_g, qpms_lorentzdrude_epsmu_g, qpms_ldparams_triple_t, qpms_lorentzdrude_eps
from .cycommon cimport make_c_string
cimport cython
import enum
import warnings
import os
from scipy.constants import e as eV, hbar
from libc.stdlib cimport malloc, free, calloc, abort

class EpsMuGeneratorType(enum.IntEnum):
    CONSTANT = 1
    PERMITTIVITY_INTERPOLATOR = 2
    LORENTZ_DRUDE = 3

cdef class EpsMu:
    def __init__(self, *args ,**kwargs):
        self.em.eps = 1
        self.em.mu = 1
        if(len(args)>=1):
            self.em.eps = args[0]
        if(len(args)>=2):
            self.em.mu = args[1]
        if 'eps' in kwargs.keys():
            self.em.eps = kwargs['eps']
        if 'mu' in kwargs.keys():
            self.em.mu = kwargs['mu']
        return

    def __repr__(self):
        return 'EpsMu(' + repr(self.em.eps) + ', ' + repr(self.em.mu) + ')'

cdef class LorentzDrudeModel:
    def __cinit__(self, eps_inf, omega_p, f_arr, omega_arr, gamma_arr):
        cdef size_t n = len(omega_arr)
        if (n != len(gamma_arr) or n != len(f_arr)):
            raise ValueError("omega_arr, gamma_arr and f_arr must have equal lengths!")
        cdef qpms_ldparams_t *p
        p = <qpms_ldparams_t *>malloc(sizeof(qpms_ldparams_t) + sizeof(qpms_ldparams_triple_t) * n)
        p[0].eps_inf = eps_inf
        p[0].omega_p = omega_p
        p[0].n = n
        cdef size_t i
        for i in range(0,n):
            p[0].data[i].f = f_arr[i]
            p[0].data[i].omega = omega_arr[i]
            p[0].data[i].gamma = gamma_arr[i]
        self.params = p

    def __dealloc__(self):
        free(self.params)
        self.params = NULL

    def __call__(self, omega):
        return qpms_lorentzdrude_eps(omega, self.params)

cdef double eh = eV/hbar

# Some basic Lorentz-Drude parameters
lorentz_drude = {
        'Au' : 
            LorentzDrudeModel(1, 9.03*eh,
                (0.76, 0.024, 0.01, 0.071, 0.601, 4.384),
                (0, 0.415*eh, 0.83*eh, 2.969*eh, 4.304*eh, 13.32*eh),
                (0.053*eh, 0.241*eh, 0.345*eh, 0.87*eh, 2.494*eh, 2.214*eh)),
        'Ag' : 
            LorentzDrudeModel(1, 9.01*eh,
                (0.84, 0.065,0.124, 0.111, 0.840, 5.646),
                (0, 0.816*eh,4.481*eh, 8.185*eh, 9.083*eh, 20.29*eh),
                    (0.053*eh, 3.886*eh, 0.452*eh,0.065*eh, 0.916*eh, 2.419*eh)),
}

cdef class EpsMuGenerator:
    def __init__(self, what):
        if isinstance(what, EpsMu):
            self.holder = what
            self.g.function = qpms_epsmu_const_g
            self.g.params = (<EpsMu?>self.holder).rawpointer()
        elif isinstance(what, LorentzDrudeModel):
            self.holder = what
            self.g.function = qpms_lorentzdrude_epsmu_g
            self.g.params = (<LorentzDrudeModel?>self.holder).rawpointer()
        elif isinstance(what, MaterialInterpolator):
            self.holder = what
            self.g.function = qpms_permittivity_interpolator_epsmu_g
            self.g.params = (<MaterialInterpolator?>self.holder).rawpointer()
        else:
            raise ValueError("Must be constructed from EpsMu, LorentzDrudeModel or MaterialInterpolator")

    property typ:
        def __get__(self):
            if(self.g.function == qpms_epsmu_const_g):
                return EpsMuGeneratorType.CONSTANT
            elif(self.g.function == qpms_lorentzdrude_epsmu_g):
                return EpsMuGeneratorType.LORENTZ_DRUDE
            elif(self.g.function == qpms_permittivity_interpolator_epsmu_g):
                return EpsMuGeneratorType.PERMITTIVITY_INTERPOLATOR
            else:
                raise ValueError("Damn, this is a bug.")

    def __call__(self, omega):
        cdef qpms_epsmu_t em
        if self.g.function == qpms_permittivity_interpolator_epsmu_g:
            i = self.holder.freq_interval
            if(omega < i[0] or omega > i[1]):
                raise ValueError("Input frequency %g is outside the interpolator domain (%g, %g)."
                    % (omega, i[0], i[1]))
        em = self.g.function(omega, self.g.params)
        return EpsMu(em.eps, em.mu)

    cdef qpms_epsmu_generator_t raw(self):
        return self.g


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

