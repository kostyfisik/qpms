from .qpms_cdefs cimport *

cimport numpy as np

cdef class BaseSpec:
    cdef qpms_vswf_set_spec_t s
    cdef np.ndarray __ilist

    cdef qpms_vswf_set_spec_t *rawpointer(BaseSpec self)
