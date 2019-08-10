from .qpms_cdefs cimport qpms_permittivity_interpolator_t

cdef class MaterialInterpolator:
    cdef qpms_permittivity_interpolator_t *interp
    cdef readonly double omegamin
    cdef readonly double omegamax


