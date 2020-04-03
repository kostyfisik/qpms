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
    property omega_table:
        def __get__(self):
            cdef size_t i
            omegas = np.empty((self.nfreqs,), dtype=float)
            cdef double[:] omegas_view = omegas
            for i in range(self.nfreqs):
                omegas_view[i] = self.freqs[i]
            return omegas


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
        if self.lMax_extend == 0:
            self.lMax_extend = 1
    def Q_transposed(self, cdouble omega, norm):
        cdef size_t n = 2*(self.p.lMax_extend*(self.p.lMax_extend +2))
        cdef np.ndarray[np.complex_t, ndim=2] arr = np.empty((n,n), dtype=complex, order='C')
        cdef cdouble[:,::1] arrview = arr
        qpms_tmatrix_generator_axialsym_RQ_transposed_fill(&arrview[0][0], omega, &self.p, norm, QPMS_HANKEL_PLUS)
        return arr
    def R_transposed(self, cdouble omega, norm):
        cdef size_t n = 2*(self.p.lMax_extend*(self.p.lMax_extend +2))
        cdef np.ndarray[np.complex_t, ndim=2] arr = np.empty((n,n), dtype=complex, order='C')
        cdef cdouble[:,::1] arrview = arr
        qpms_tmatrix_generator_axialsym_RQ_transposed_fill(&arrview[0][0], omega, &self.p, norm, QPMS_BESSEL_REGULAR)
        return arr

cdef class TMatrixFunction:
    '''
    Wrapper over qpms_tmatrix_function_t. The main functional difference between this
    and TMatrixGenerator is that this remembers a specific BaseSpec
    and its __call__ method takes only one mandatory argument (in addition to self).
    '''
    def __init__(self, TMatrixGenerator tmg, BaseSpec spec):
        self.generator = tmg
        self.spec = spec
        self.f.gen = self.generator.rawpointer()
        self.f.spec = self.spec.rawpointer()

    def __call__(self, cdouble omega, fill = None):
        cdef CTMatrix tm
        if fill is None: # make a new CTMatrix
            tm = CTMatrix(self.spec, None)
        else: # TODO check whether fill has the same bspec as self?
            tm = fill
        if self.f.gen.function(tm.rawpointer(), omega, self.f.gen.params) != 0:
            raise ValueError("Something went wrong")
        else:
            return tm


cdef class TMatrixGenerator:
    def __init__(self, what):
        if isinstance(what, __MieParams):
            self.holder = what
            self.g.function = qpms_tmatrix_generator_sphere
            self.g.params = (<__MieParams?>self.holder).rawpointer()
        elif isinstance(what,__AxialSymParams):
            self.holder = what
            self.g.function = qpms_tmatrix_generator_axialsym
            self.g.params = (<__AxialSymParams?>self.holder).rawpointer()
        elif isinstance(what, CTMatrix):
            self.holder = what
            self.g.function = qpms_tmatrix_generator_constant
            self.g.params = <void*>(<CTMatrix?>self.holder).rawpointer()
        elif isinstance(what, TMatrixInterpolator):
            self.holder = what
            self.g.function = qpms_tmatrix_generator_interpolator
            self.g.params = <void*>(<TMatrixInterpolator?>self.holder).rawpointer()
        else:
            raise TypeError("Can't construct TMatrixGenerator from that")

    def __call__(self, arg, cdouble omega):
        """Produces a T-matrix at a given frequency.

        Parameters
        ----------
        arg : CTMatrix or BaseSpec
            If arg is a CTMatrix, its contents will be replaced.
            If arg is a BaseSpec, a new CTMatrix instance will be created.
        omega : complex
            Angular frequency at which the T-matrix shall be evaluated.

        Returns
        -------
        t : CTMatrix
        """
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

    def Q_transposed(self, cdouble omega, norm):
        if self.g.function != qpms_tmatrix_generator_axialsym:
            raise TypeError("Applicable only for axialsym generators")
        return self.holder.Q_transposed(omega, norm)
    def R_transposed(self, cdouble omega, norm):
        if self.g.function != qpms_tmatrix_generator_axialsym:
            raise TypeError("Applicable only for axialsym generators")
        return self.holder.R_transposed(omega, norm)

    # Better "constructors":
    @staticmethod
    def sphere(outside, inside, r):
        """Creates a T-matrix generator for a spherical scatterer.

        This method uses analytical Mie-Lorentz formulae.

        Parameters:
        -----------
        outside : EpsMuGenerator or EpsMu
            Optical properties of the surrounding medium.
        inside : EpsMuGenerator or EpsMu
            Optical properties of the material inside the sphere.
        r : double
            Sphere radius
        """
        return TMatrixGenerator(__MieParams(EpsMuGenerator(outside),
                    EpsMuGenerator(inside), r))

    @staticmethod
    def sphere_asarc(outside, inside, r, *args, **kwargs):
        """Creates a T-matrix generator for a spherical scatterer.

        This method uses numerical evaluation for generic axially-symmetric
        scatterer, and is intended for testing and benchmarking only.
        For regular use, see TMatrigGenerator.sphere() instead.

        Parameters
        ----------
        outside : EpsMuGenerator or EpsMu
            Optical properties of the surrounding medium.
        inside : EpsMuGenerator or EpsMu
            Optical properties of the material inside the sphere.
        r : double
            Sphere radius

        Returns
        -------
        tmgen : TMatrixGenerator

        See Also
        --------
        TMatrigGenerator.sphere : Faster and more precise method.
        """
        return TMatrixGenerator(__AxialSymParams(
            EpsMuGenerator(outside), EpsMuGenerator(inside),
            ArcFunction(__ArcSphere(r)), *args, **kwargs))

    @staticmethod
    def cylinder(outside, inside, r, h, *args, **kwargs):
        """Creates a T-matrix generator for a right circular cylinder.

        Parameters:
        -----------
        outside : EpsMuGenerator or EpsMu
            Optical properties of the surrounding medium.
        inside : EpsMuGenerator or EpsMu
            Optical properties of the material inside the cylinder.
        r : double
            Cylinder base radius.
        h : double
            Cylinder height.

        Returns
        -------
        tmgen : TMatrixGenerator
        """
        return TMatrixGenerator(__AxialSymParams(
            EpsMuGenerator(outside), EpsMuGenerator(inside),
            ArcFunction(__ArcCylinder(r, h)), *args, **kwargs))

