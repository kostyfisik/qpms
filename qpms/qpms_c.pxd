from .qpms_cdefs cimport qpms_finite_group_t

cdef class FinitePointGroup:
    cdef readonly bint owns_data
    cdef qpms_finite_group_t *G

    cdef inline qpms_finite_group_t *rawpointer(self):
        return self.G
