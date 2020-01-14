cimport numpy as np
from .qpms_cdefs cimport qpms_tmatrix_t, cdouble, qpms_tmatrix_interpolator_t, qpms_tmatrix_generator_t
from .cybspec cimport BaseSpec

cdef class TMatrixInterpolator:
    #cdef readonly np.ndarray m # Numpy array holding the matrix data
    cdef readonly BaseSpec spec # Here we hold the base spec for the correct reference counting; TODO check if it gets copied
    cdef qpms_tmatrix_t *tmatrices_array
    cdef cdouble *tmdata
    cdef double *freqs
    cdef double *freqs_su
    cdef size_t nfreqs
    cdef qpms_tmatrix_interpolator_t *interp
    cdef inline qpms_tmatrix_interpolator_t *rawpointer(self):
        return self.interp

cdef class TMatrixGenerator:
    cdef qpms_tmatrix_generator_t g
    cdef object holder
    cdef inline qpms_tmatrix_generator_t raw(self):
        return self.g
    cdef inline qpms_tmatrix_generator_t *rawpointer(self):
        return &(self.g)

cdef class TMatrixGeneratorTransformed:
    pass

cdef class CTMatrix: # N.B. there is another type called TMatrix in tmatrices.py!
    cdef readonly np.ndarray m # Numpy array holding the matrix data
    cdef readonly BaseSpec spec # Here we hold the base spec for the correct reference counting; TODO check if it gets copied
    cdef qpms_tmatrix_t t


    cdef inline qpms_tmatrix_t *rawpointer(CTMatrix self):
        '''Pointer to the qpms_tmatrix_t structure.
        Don't forget to reference the BaseSpec object itself when storing the pointer anywhere!!!
        '''
        return &(self.t)




