from .qpms_cdefs cimport *
from libc.stdlib cimport malloc, free, calloc
import numpy as np

cdef extern from "ewald.h":
    void ewald3_2_sigma_long_Delta(qpms_csf_result *target, int maxn, cdouble x, cdouble z)
    int complex_gamma_inc_e(double a, cdouble x, int m, qpms_csf_result *result)

def e32_Delta(int maxn, cdouble x, cdouble z):
    cdef qpms_csf_result *target = <qpms_csf_result *>malloc((maxn+1)*sizeof(qpms_csf_result))
    cdef np.ndarray[cdouble, ndim=1] target_np = np.empty((maxn+1,), dtype=complex, order='C')
    ewald3_2_sigma_long_Delta(target, maxn, x, z)
    cdef int i
    for i in range(maxn+1):
        target_np[i] = target[i].val
    free(target)
    return target_np

def gamma_inc(double a, cdouble x, int m=0):
    cdef qpms_csf_result res
    complex_gamma_inc_e(a, x, m, &res)
    return res.val


