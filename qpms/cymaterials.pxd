from .qpms_cdefs cimport qpms_permittivity_interpolator_t, qpms_epsmu_generator_t, qpms_epsmu_t, qpms_ldparams_t

cdef extern from "materials.h":
    const qpms_ldparams_t *const QPMS_LDPARAMS_AG
    const qpms_ldparams_t *const QPMS_LDPARAMS_AU
    const qpms_ldparams_t *const QPMS_LDPARAMS_CU
    const qpms_ldparams_t *const QPMS_LDPARAMS_AL
    const qpms_ldparams_t *const QPMS_LDPARAMS_CR
    const qpms_ldparams_t *const QPMS_LDPARAMS_TI
    const qpms_ldparams_t *const QPMS_LDPARAMS_BE
    const qpms_ldparams_t *const QPMS_LDPARAMS_NI
    const qpms_ldparams_t *const QPMS_LDPARAMS_PD
    const qpms_ldparams_t *const QPMS_LDPARAMS_PT
    const qpms_ldparams_t *const QPMS_LDPARAMS_W

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
        return <void *>(self.params)

cdef class _CLorentzDrudeModel:
    ''' Drude-Lorentz parameters initialised from raw C structure. Private, do not use. '''
    cdef const qpms_ldparams_t *params
    cdef inline void *rawpointer(self):
        return <void *>(self.params)
    @staticmethod 
    cdef link(const qpms_ldparams_t *params) # The actual constructor

cdef class EpsMuGenerator:
    cdef qpms_epsmu_generator_t g
    cdef object holder
    cdef qpms_epsmu_generator_t raw(self)

