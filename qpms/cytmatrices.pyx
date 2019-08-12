import numpy as np
from .qpms_cdefs cimport *
from .cybspec cimport BaseSpec
from .cycommon import *
from .cycommon cimport make_c_string
from .cymaterials cimport EpsMuGenerator
from .qpms_c cimport FinitePointGroup
import warnings
import os
from libc.stdlib cimport free

cdef class TMatrixInterpolator:
    '''
    Wrapper over the qpms_tmatrix_interpolator_t structure.
    '''
    def __cinit__(self, filename, BaseSpec bspec,  *args, **kwargs):
        '''Creates a T-matrix interpolator object from a scuff-tmatrix output'''
        global qpms_load_scuff_tmatrix_crash_on_failure
        qpms_load_scuff_tmatrix_crash_on_failure = False
        self.spec = bspec
        cdef char * cpath = make_c_string(filename)
        retval = qpms_load_scuff_tmatrix(cpath, self.spec.rawpointer(),
                &(self.nfreqs), &(self.freqs), &(self.freqs_su),
                &(self.tmatrices_array), &(self.tmdata))
        if (retval != QPMS_SUCCESS):
            raise IOError("Could not read T-matrix from %s: %s" % (filename, os.strerror(retval)))
        if 'symmetrise' in kwargs:
            sym = kwargs['symmetrise']
            if isinstance(sym, FinitePointGroup):
                if QPMS_SUCCESS != qpms_symmetrise_tmdata_finite_group(
                        self.tmdata, self.nfreqs, self.spec.rawpointer(),
                        (<FinitePointGroup?>sym).rawpointer()):
                    raise Exception("This should not happen.")
                atol = kwargs['atol'] if 'atol' in kwargs else 1e-16
                qpms_czero_roundoff_clean(self.tmdata, self.nfreqs * len(bspec)**2, atol)
            else:
                warnings.warn('symmetrise argument type not supported; ignoring.')
        self.interp = qpms_tmatrix_interpolator_create(self.nfreqs,
                self.freqs, self.tmatrices_array, gsl_interp_cspline)
        if not self.interp: raise Exception("Unexpected NULL at interpolator creation.")
    def __call__(self, double freq):
        '''Returns a TMatrix instance, corresponding to a given frequency.'''
        if freq < self.freqs[0] or freq > self.freqs[self.nfreqs-1]:# FIXME here I assume that the input is already sorted
            raise ValueError("input frequency %g is outside the interpolator domain (%g, %g)"
                    % (freq, self.freqs[0], self.freqs[self.nfreqs-1]))
        # This is a bit stupid, I should rethink the CTMatrix constuctors
        cdef qpms_tmatrix_t *t = qpms_tmatrix_interpolator_eval(self.interp, freq)
        cdef CTMatrix res = CTMatrix(self.spec, <cdouble[:len(self.spec),:len(self.spec)]>(t[0].m))
        qpms_tmatrix_free(t)
        return res
    def __dealloc__(self):
        qpms_tmatrix_interpolator_free(self.interp)
        free(self.tmatrices_array)
        free(self.tmdata)
        free(self.freqs_su)
        free(self.freqs)
    property freq_interval:
        def __get__(self):
            return [self.freqs[0], self.freqs[self.nfreqs-1]]

cdef class CTMatrix: # N.B. there is another type called TMatrix in tmatrices.py!
    '''
    Wrapper over the C qpms_tmatrix_t stucture. 
    '''

    def __cinit__(CTMatrix self, BaseSpec spec, matrix):
        self.spec = spec
        self.t.spec = self.spec.rawpointer();
        if (matrix is None) or not np.any(matrix):
            self.m = np.zeros((len(spec),len(spec)), dtype=complex, order='C')
        else:
            # The following will raise an exception if shape is wrong
            self.m = np.array(matrix, dtype=complex, copy=True, order='C').reshape((len(spec), len(spec)))
        #self.m.setflags(write=False) # checkme
        cdef cdouble[:,:] m_memview = self.m
        self.t.m = &(m_memview[0,0])
        self.t.owns_m = False # Memory in self.t.m is "owned" by self.m, not by self.t...

    property rawpointer:
        def __get__(self):
            return <uintptr_t> &(self.t)

    # Transparent access to the T-matrix elements.
    def __getitem__(self, key):
        return self.m[key]
    def __setitem__(self, key, value):
        self.m[key] = value

    def as_ndarray(CTMatrix self):
        ''' Returns a copy of the T-matrix as a numpy array.'''
        # Maybe not totally needed after all, as np.array(T[...]) should be equivalent and not longer
        return np.array(self.m, copy=True)

    def spherical_fill(CTMatrix self, double radius, cdouble k_int,
            cdouble k_ext, cdouble mu_int = 1, cdouble mu_ext = 1):
        ''' Replaces the contents of the T-matrix with those of a spherical particle.'''
        qpms_tmatrix_spherical_fill(&self.t, radius, k_int, k_ext, mu_int, mu_ext)

    def spherical_perm_fill(CTMatrix self, double radius, double freq, cdouble epsilon_int,
            cdouble epsilon_ext):
        '''Replaces the contenst of the T-matrix with those of a spherical particle.'''
        qpms_tmatrix_spherical_mu0_fill(&self.t, radius, freq, epsilon_int, epsilon_ext)
        
    @staticmethod
    def spherical(BaseSpec spec, double radius, cdouble k_int, cdouble k_ext, 
            cdouble mu_int = 1, cdouble mu_ext = 1):
        ''' Creates a T-matrix of a spherical nanoparticle. '''
        tm = CTMatrix(spec, 0)
        tm.spherical_fill(radius, k_int, k_ext, mu_int, mu_ext)
        return tm
    
    @staticmethod
    def spherical_perm(BaseSpec spec, double radius, double freq, cdouble epsilon_int, cdouble epsilon_ext):
        '''Creates a T-matrix of a spherical nanoparticle.'''
        tm = CTMatrix(spec, 0)
        tm.spherical_perm_fill(radius, freq, epsilon_int, epsilon_ext)
        return tm

cdef class __MieParams:
    # Not directly callable right now, serves just to be used by TMatrixGenerator.
    cdef qpms_tmatrix_generator_sphere_param_t cparam
    cdef EpsMuGenerator outside
    cdef EpsMuGenerator inside
    cdef inline  void *rawpointer(self):
        return <void *>&(self.cparam)
    
    def __init__(self, outside, inside, r):
        self.inside = inside
        self.outside = outside
        self.cparam.inside = self.inside.raw();
        self.cparam.outside = self.outside.raw();
        self.cparam.radius = r

    property r:
        def __get__(self):
            return self.cparam.radius
        def __set__(self, val):
            self.cparam.radius = val

cdef class __ArcCylinder:
    cdef qpms_arc_cylinder_params_t p
    cdef inline void *rawpointer(self):
        return <void *> &(self.p)
    def __init__(self, R, h):
        self.p.R = R
        self.p.h = h

cdef class __ArcSphere:
    cdef double r
    cdef inline void *rawpointer(self):
        return <void *> &(self.r)
    def __init__(self, r):
        self.r = r

cdef qpms_arc_function_retval_t userarc(double theta, const void *params):
    cdef object fun = <object> params
    cdef qpms_arc_function_retval_t retval
    retval.r, retval.beta = fun(theta)
    return retval


cdef class ArcFunction:
    cdef qpms_arc_function_t g
    cdef object holder
    def __init__(self, what):
        if isinstance(what, __ArcCylinder):
            self.holder = what
            self.g.function = qpms_arc_cylinder
            self.g.params = (<__ArcCylinder?>self.holder).rawpointer()
        elif isinstance(what, __ArcSphere):
            self.holder = what
            self.g.function = qpms_arc_sphere
            self.g.params = (<__ArcSphere?>self.holder).rawpointer()
        elif callable(what):
            warnings.warn("Custom python (r, beta) arc functions are an experimental feature. Also expect it to be slow.")
            self.holder = what
            self.g.function = userarc
            self.g.params = <const void *> self.holder
        elif isinstance(what, ArcFunction): #Copy constructor
            self.holder = what.holder
            self.g = (<ArcFunction?>what).g
            self.holder.rawpointer()

cdef class __AxialSymParams:
    cdef qpms_tmatrix_generator_axialsym_param_t p
    cdef EpsMuGenerator outside
    cdef EpsMuGenerator inside
    cdef ArcFunction shape
    cdef void * rawpointer(self):
        return <void *> &(self.p)
    property lMax_extend:
        def __get__(self):
            return self.p.lMax_extend
        def __set__(self, val):
            self.p.lMax_extend = val
    def __init__(self, outside, inside, shape, *args, **kwargs):
        self.outside = outside
        self.p.outside = self.outside.g
        self.inside = inside
        self.p.inside = self.inside.g
        self.shape = shape
        self.p.shape = self.shape.g
        if len(args)>0:
            self.lMax_extend = args[0]
        if 'lMax_extend' in kwargs.keys():
            self.lMax_extend = kwargs['lMax_extend']

cdef class TMatrixGenerator:
    cdef qpms_tmatrix_generator_t g
    cdef object holder
    cdef qpms_tmatrix_generator_t raw(self):
        return self.g
    def __init__(self, what):
        if isinstance(what, __MieParams):
            self.holder = what
            self.g.function = qpms_tmatrix_generator_sphere
            self.g.params = (<__MieParams?>self.holder).rawpointer()
        elif isinstance(what,__AxialSymParams):
            self.holder = what
            self.g.function = qpms_tmatrix_generator_axialsym
            self.g.params = (<__AxialSymParams?>self.holder).rawpointer()
        # TODO INTERPOLATOR
        else:
            raise ValueError("Can't construct TMatrixGenerator from that")

    def __call__(self, arg, cdouble omega):
        cdef CTMatrix tm
        if isinstance(arg, CTMatrix): # fill the matrix
            tm = arg
            if self.g.function(tm.rawpointer(), omega, self.g.params) != 0:
                raise ValueError("Something went wrong")
            return
        elif isinstance(arg, BaseSpec): # Make a new CTMatrix
            tm = CTMatrix(arg, None)
            if self.g.function(tm.rawpointer(), omega, self.g.params) != 0:
                raise ValueError("Something went wrong")
            return tm
        else:
            raise ValueError("Must specify CTMatrix or BaseSpec")

    # Better "constructors":
    @staticmethod
    def sphere(outside, inside, r):
        return TMatrixGenerator(__MieParams(EpsMuGenerator(outside),
                    EpsMuGenerator(inside), r))
    @staticmethod
    def sphere_asarc(outside, inside, r, *args, **kwargs):
        return TMatrixGenerator(__AxialSymParams(
            EpsMuGenerator(outside), EpsMuGenerator(inside),
            ArcFunction(__ArcSphere(r)), *args, **kwargs))
    @staticmethod
    def cylinder(outside, inside, r, h, *args, **kwargs):
        return TMatrixGenerator(__AxialSymParams(
            EpsMuGenerator(outside), EpsMuGenerator(inside),
            ArcFunction(__ArcCylinder(r, h)), *args, **kwargs))

