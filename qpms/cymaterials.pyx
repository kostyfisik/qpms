# Cythonized parts of QPMS here
# -----------------------------

import numpy as np
import cmath
from .qpms_cdefs cimport qpms_permittivity_interpolator_from_yml, qpms_permittivity_interpolator_free, qpms_permittivity_interpolator_omega_min, qpms_permittivity_interpolator_omega_max, gsl_interp_type, qpms_permittivity_interpolator_t, gsl_interp_cspline, qpms_permittivity_interpolator_eps_at_omega, qpms_epsmu_const_g, qpms_permittivity_interpolator_epsmu_g, qpms_epsmu_const_g, qpms_lorentzdrude_epsmu_g, qpms_ldparams_triple_t, qpms_lorentzdrude_eps, cdouble
from .cycommon cimport make_c_string
cimport cython
import enum
import warnings
import os
from scipy.constants import e as eV, hbar, c
from libc.stdlib cimport malloc, free, calloc, abort

class EpsMuGeneratorType(enum.Enum):
    CONSTANT = 1
    PERMITTIVITY_INTERPOLATOR = 2
    LORENTZ_DRUDE = 3
    PYTHON_CALLABLE = 4

cdef class EpsMu:
    """Permittivity and permeability of an isotropic material.

    This wraps the C structure qpms_epsmu_t.

    See Also
    --------
    EpsMuGenerator : generates EpsMu objects as a function of frequency.
    """
    def __init__(self, *args ,**kwargs):
        """EpsMu object constructor

        Parameters
        ----------
        eps : complex
            Relative electric permittivity.
        mu : complex
            Relative magnetic permeability.
        """
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

    def k(self, omega):
        """Wavenumber of the material at frequency omega

        Parameters
        ----------
        omega : complex
            Angular frequency in 1/s.

        Returns
        -------
        out : complex
            Wavenumber of the material at omega assuming permittivity and permeability
            from self.
        """
        return self.n * omega / c

    property n:
        """Refractive index of a material specified by this permittivity and permeability."""
        def __get__(self):
            return (self.em.eps * self.em.mu)**.5
    property Z:
        """Wave impedance of a material specified by this permittivity and permeability."""
        def __get__(self):
            return (self.em.mu / self.em.eps)**.5

cdef class LorentzDrudeModel:
    """Lorentz-Drude model of material permittivity.

    This wraps the C structure qpms_ldparams_t.

    Some materials are available in the `lorentz_drude` dictionary.
    """

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
        """Evaluates the permittivity at a given frequency

        Parameters
        ----------
        omega : complex
            Angular frequency in 1/s.

        Returns
        -------
        eps : complex
            Relative permittivity from the Lorentz-Drude model
        """ 
        return qpms_lorentzdrude_eps(omega, self.params)

cdef class _CLorentzDrudeModel:
    """Lorentz-Drude model of material permittivity.

    This is an auxilliary class making the pre-compiled C Lorentz-Drude models accessible.
    For defining own Lorentz-Drude models from python, use LorentzDrudeModel class instead.
    """
    def __cinit__(self):
        "Do not use directly. Do not use from Python. Use the link() method instead."
        pass

    @staticmethod
    cdef link(const qpms_ldparams_t *params):
        self = _CLorentzDrudeModel()
        self.params = params
        return self

    def __call__(self, omega):
        """Evaluates the permittivity at a given frequency

        Parameters
        ----------
        omega : complex
            Angular frequency in 1/s.

        Returns
        -------
        eps : complex
            Relative permittivity from the Lorentz-Drude model
        """
        return qpms_lorentzdrude_eps(omega, self.params)

cdef double eh = eV/hbar

# Some basic Lorentz-Drude parameters
lorentz_drude = {
        'Au_py' :  # This should give the same results as 'Au'; to be removed.
            LorentzDrudeModel(1, 9.03*eh,
                (0.76, 0.024, 0.01, 0.071, 0.601, 4.384),
                (0, 0.415*eh, 0.83*eh, 2.969*eh, 4.304*eh, 13.32*eh),
                (0.053*eh, 0.241*eh, 0.345*eh, 0.87*eh, 2.494*eh, 2.214*eh)),
        'Ag_py' :  # This should give the same results as 'Ag'; to be removed.
            LorentzDrudeModel(1, 9.01*eh,
                (0.84, 0.065,0.124, 0.111, 0.840, 5.646),
                (0, 0.816*eh,4.481*eh, 8.185*eh, 9.083*eh, 20.29*eh),
                    (0.053*eh, 3.886*eh, 0.452*eh,0.065*eh, 0.916*eh, 2.419*eh)),
        'Au' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_AU),
        'Ag' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_AG),
        'Cu' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_CU),
        'Al' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_AL),
        'Cr' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_CR),
        'Ti' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_TI),
        'Be' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_BE),
        'Ni' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_NI),
        'Pd' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_PD),
        'Pt' : _CLorentzDrudeModel.link(QPMS_LDPARAMS_PT),
        'W'  : _CLorentzDrudeModel.link(QPMS_LDPARAMS_W),

}

cdef qpms_epsmu_t python_epsmu_generator(cdouble omega, const void *params):
    cdef object fun = <object> params
    cdef qpms_epsmu_t em
    em.eps, em.mu = fun(omega)
    return em

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
        elif isinstance(what, _CLorentzDrudeModel):
            self.holder = what
            self.g.function = qpms_lorentzdrude_epsmu_g
            self.g.params = (<_CLorentzDrudeModel?>self.holder).rawpointer()
        elif isinstance(what, MaterialInterpolator):
            self.holder = what
            self.g.function = qpms_permittivity_interpolator_epsmu_g
            self.g.params = (<MaterialInterpolator?>self.holder).rawpointer()
        elif isinstance(what, EpsMuGenerator): # Copy constructor
            self.holder = (<EpsMuGenerator?>what).holder
            self.g = (<EpsMuGenerator?>what).g
        elif callable(what):
            warnings.warn("Custom python (eps,mu) generators are an experimental feature")
            self.holder = what
            self.g.function = python_epsmu_generator
            self.g.params = <const void *> what
        else:
            raise ValueError("Must be constructed from EpsMu, LorentzDrudeModel or MaterialInterpolator, or a python callable object that returns an (epsilon, mu) tuple.")

    property typ:
        def __get__(self):
            if(self.g.function == qpms_epsmu_const_g):
                return EpsMuGeneratorType.CONSTANT
            elif(self.g.function == qpms_lorentzdrude_epsmu_g):
                return EpsMuGeneratorType.LORENTZ_DRUDE
            elif(self.g.function == qpms_permittivity_interpolator_epsmu_g):
                return EpsMuGeneratorType.PERMITTIVITY_INTERPOLATOR
            elif(self.g.function == python_epsmu_generator):
                return EpsMuGeneratorType.PYTHON_CALLABLE
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

    def n(self, omega):
        return self(omega).n
    def Z(self, omega):
        return self(omega).Z
    def k(self, omega):
        return self(omega).k(omega)

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

