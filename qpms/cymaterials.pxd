from .qpms_cdefs cimport qpms_permittivity_interpolator_t, qpms_epsmu_generator_t, qpms_epsmu_t, qpms_ldparams_t

cdef class MaterialInterpolator:
    cdef qpms_permittivity_interpolator_t *interp
    cdef readonly double omegamin
    cdef readonly double omegamax
    cdef inline void *rawpointer(self):
        return <void *>&(self.interp)

cdef class EpsMu:
    cdef public qpms_epsmu_t em
    cdef inline void *rawpointer(self):
        return <void *>&(self.em)

cdef class LorentzDrudeModel:
    cdef const qpms_ldparams_t *params
    #cdef bint owns_params
    cdef inline void *rawpointer(self):
        return <void *>&(self.params)

cdef class EpsMuGenerator:
    cdef qpms_epsmu_generator_t g
    cdef object holder
    cdef qpms_epsmu_generator_t raw(self)

