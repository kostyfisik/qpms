import numpy as np
from qpms_cdefs cimport *
from cybspec cimport BaseSpec
from cycommon import *
from cycommon cimport make_c_string
from qpms_c cimport FinitePointGroup
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
