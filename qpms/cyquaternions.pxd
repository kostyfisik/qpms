from .qpms_cdefs cimport qpms_quat_t, qpms_irot3_t

cdef class CQuat:
    cdef readonly qpms_quat_t q

cdef class IRot3:
    cdef readonly qpms_irot3_t qd
    cdef void cset(self, qpms_irot3_t qd)
