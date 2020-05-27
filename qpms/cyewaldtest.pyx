from .qpms_cdefs cimport *
from libc.stdlib cimport malloc, free, calloc
import numpy as np

cdef extern from "ewald.h":
    void ewald3_2_sigma_long_Delta(cdouble *target, double *err,  int maxn, cdouble x, cdouble z)
    void ewald3_2_sigma_long_Delta_series(cdouble *target, double *err,  int maxn, cdouble x, cdouble z)
    void ewald3_2_sigma_long_Delta_recurrent(cdouble *target, double *err,  int maxn, cdouble x, cdouble z)
    int complex_gamma_inc_e(double a, cdouble x, int m, qpms_csf_result *result)

def e32_Delta(int maxn, cdouble x, cdouble z, get_err=True, method='auto'):
    cdef np.ndarray[double, ndim=1] err_np
    cdef double[::1] err_view
    cdef np.ndarray[np.complex_t, ndim=1] target_np = np.empty((maxn+1,), dtype=complex, order='C')
    cdef cdouble[::1] target_view = target_np
    if get_err:
        err_np = np.empty((maxn+1,), order='C')
        err_view = err_np
    if method == 'recurrent':
        ewald3_2_sigma_long_Delta_recurrent(&target_view[0], &err_view[0] if get_err else NULL, maxn, x, z)
    elif method == 'series':
        ewald3_2_sigma_long_Delta_series(&target_view[0], &err_view[0] if get_err else NULL, maxn, x, z)
    else:
        ewald3_2_sigma_long_Delta(&target_view[0], &err_view[0] if get_err else NULL, maxn, x, z)
    if get_err:
        return target_np, err_np
    else:
        return target_np

def gamma_inc(double a, cdouble x, int m=0):
    cdef qpms_csf_result res
    complex_gamma_inc_e(a, x, m, &res)
    return res.val


