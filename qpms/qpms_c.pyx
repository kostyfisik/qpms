# Cythonized parts of QPMS here
# -----------------------------

import numpy as np
cimport numpy as np
cimport cython
from cython.parallel cimport parallel, prange
#cimport openmp
#openmp.omp_set_dynamic(1)

## Auxillary function for retrieving the "meshgrid-like" indices; inc. nmax
@cython.boundscheck(False)
def get_mn_y(int nmax):
    """
    Auxillary function for retreiving the 'meshgrid-like' indices from the flat indexing; 
    inc. nmax.
    ('y to mn' conversion)
    
    Parameters
    ----------

    nmax : int
        The maximum order to which the VSWFs / Legendre functions etc. will be evaluated.
        
    Returns
    -------
    
    output : (m, n)
        Tuple of two arrays of type np.array(shape=(nmax*nmax + 2*nmax), dtype=np.int),
        where [(m[y],n[y]) for y in range(nmax*nmax + 2*nma)] covers all possible 
        integer pairs n >= 1, -n <= m <= n.
    """
    cdef Py_ssize_t nelems = nmax * nmax + 2 * nmax
    cdef np.ndarray[np.int_t,ndim=1] m_arr = np.empty([nelems], dtype=np.int)
    cdef np.ndarray[np.int_t,ndim=1] n_arr = np.empty([nelems], dtype=np.int)
    cdef Py_ssize_t i = 0
    cdef np.int_t n, m
    for n in range(1,nmax+1):
        for m in range(-n,n+1):
            m_arr[i] = m
            n_arr[i] = n
            i = i + 1
    return (m_arr, n_arr)

def get_y_mn_unsigned(int nmax): 
    """
    Auxillary function for mapping 'unsigned m', n indices to the flat y-indexing.
    For use with functions as scipy.special.lpmn, which have to be evaluated separately
    for positive and negative m.
    
    Parameters
    ----------

    nmax : int
        The maximum order to which the VSWFs / Legendre functions etc. will be evaluated.
        
    output : (ymn_plus, ymn_minus)
        Tuple of two arrays of shape (nmax+1,nmax+1), containing the flat y-indices corresponding
        to the respective (m,n) and (-m,n). The elements for which |m| > n are set to -1.
        (Therefore, the caller must not use those elements equal to -1.)
    """
    cdef np.ndarray[np.intp_t, ndim=2] ymn_plus = np.full((nmax+1,nmax+1),-1, dtype=np.intp)
    cdef np.ndarray[np.intp_t, ndim=2] ymn_minus = np.full((nmax+1,nmax+1),-1, dtype=np.intp)
    cdef Py_ssize_t i = 0
    cdef np.int_t n, m
    for n in range(1,nmax+1):
        for m in range(-n,0):
            ymn_minus[-m,n] = i
            i = i + 1
        for m in range(0,n+1):
            ymn_plus[m,n] = i
            i = i + 1
    return(ymn_plus, ymn_minus)

cdef int q_max(int m, int n, int mu, int nu):
    return min(n,nu,(n+nu-abs(m+mu)//2))

"""
Now we generate our own universal functions to be used with numpy.

Good way to see how this is done is to look at scipy/scipy/special/generate_ufuncs.py
and scipy/scipy/special/generate_ufuncs.py

In simple words, it works like this:
- Let's have a single element function. This can be function which returns or a "subroutine".
- Then we need a loop function; this is a wrapper that gets bunch of pointers from numpy and
  has to properly call the single element function.
- From those two, we build a python object using PyUFunc_FromFuncAndData.
  * If the ufunc is supposed to work on different kinds of input/output types,
    then a pair of single-element and loop functions is o be provided for
    each combination of types. However, the single-element function can be reused if
    the corresponding loop functions do the proper casting.
"""

## as in scipy/special/_ufuncs_cxx.pyx
##-------------------------------------
#cdef extern from "numpy/ufuncobject.h":
#    int PyUFunc_getfperr() nogil
#
#cdef public int wrap_PyUFunc_getfperr() nogil:
#    """
#    Call PyUFunc_getfperr in a context where PyUFunc_API array is initialized;
#
#    """
#    return PyUFunc_getfperr()
#
#cimport sf_error
#-------------------------------------


ctypedef double complex cdouble

cdef void loop_D_iiiidddii_As_D_lllldddbl(char **args, np.npy_intp *dims, np.npy_intp *steps, void *data) nogil:
    cdef np.npy_intp i, n = dims[0]
    cdef void *func = (<void**>data)#[0]
    #cdef char *func_name= <char*>(<void**>data)[1] # i am not using this, nor have I saved func_name to data
    cdef char *ip0 = args[0]
    cdef char *ip1 = args[1]
    cdef char *ip2 = args[2]
    cdef char *ip3 = args[3]
    cdef char *ip4 = args[4]
    cdef char *ip5 = args[5]
    cdef char *ip6 = args[6]
    cdef char *ip7 = args[7]
    cdef char *ip8 = args[8]
    cdef char *op0 = args[9]
    cdef cdouble ov0
    for i in range(n): # iterating over dimensions
        ov0 = (<double complex(*)(int, int, int, int, double, double, double, int, int) nogil>func)(
            <int>(<np.npy_long*>ip0)[0],
            <int>(<np.npy_long*>ip1)[0],
            <int>(<np.npy_long*>ip2)[0],
            <int>(<np.npy_long*>ip3)[0],
            <double>(<np.npy_double*>ip4)[0],
            <double>(<np.npy_double*>ip5)[0],
            <double>(<np.npy_double*>ip6)[0],
            <int>(<np.npy_bool*>ip7)[0],
            <int>(<np.npy_long*>ip8)[0],
        )
        (<cdouble *>op0)[0] = <cdouble>ov0
        ip0 += steps[0]
        ip1 += steps[1]
        ip2 += steps[2]
        ip3 += steps[3]
        ip4 += steps[4]
        ip5 += steps[5]
        ip6 += steps[6]
        ip7 += steps[7]
        ip8 += steps[8]
        op0 += steps[9]
# FIXME ERROR HANDLING!!! requires correct import and different data passed (see scipy's generated ufuncs)
#    sf_error.check_fpe(func_name)


#cdef extern from "numpy/arrayobject.h":
#    cdef enum NPY_TYPES:
#        NPY_DOUBLE
#        NPY_CDOUBLE # complex double
#        NPY_LONG # int
#    ctypedef int npy_intp


cdef extern from "translations.h":
    cdouble qpms_trans_single_A_Taylor_ext(int m, int n, int mu, int nu,
        double r, double th, double ph, int r_ge_d, int J) nogil
    cdouble qpms_trans_single_B_Taylor_ext(int m, int n, int mu, int nu,
        double r, double th, double ph, int r_ge_d, int J) nogil
    struct qpms_trans_calculator:
        pass
    enum qpms_normalization_t:
        pass
    qpms_trans_calculator* qpms_trans_calculator_init(int lMax, int nt) # should be qpms_normalization_t
    void qpms_trans_calculator_free(qpms_trans_calculator* c)
    cdouble qpms_trans_calculator_get_A_ext(const qpms_trans_calculator* c,
            int m, int n, int mu, int nu, double kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J)
    cdouble qpms_trans_calculator_get_B_ext(const qpms_trans_calculator* c,
            int m, int n, int mu, int nu, double kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J)
    int qpms_trans_calculator_get_AB_p_ext(const qpms_trans_calculator* c,
            cdouble *Adest, cdouble *Bdest,
            int m, int n, int mu, int nu, double kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J)




# Module initialisation
# ---------------------

np.import_array()  # not sure whether this is really needed
np.import_ufunc()

# Arrays passed to PyUFunc_FromFuncAndData()
# ------------------------------------------

# BTW, aren't there anonymous arrays in cython?

cdef np.PyUFuncGenericFunction trans_X_taylor_loop_func[1]
cdef void *trans_A_taylor_elementwise_funcs[1]
cdef void *trans_B_taylor_elementwise_funcs[1]

trans_X_taylor_loop_func[0] = loop_D_iiiidddii_As_D_lllldddbl

# types to be used for all of the single-type translation 
# coefficient retrieval ufuncs called like
# coeff = func(m, n, mu, nu, r, theta, phi, r_ge_d, J)
# currently supported signatures: (D_lllldddbl)
cdef char ufunc__get_either_trans_coeff_types[10]
ufunc__get_either_trans_coeff_types[0] = np.NPY_LONG
ufunc__get_either_trans_coeff_types[1] = np.NPY_LONG
ufunc__get_either_trans_coeff_types[2] = np.NPY_LONG
ufunc__get_either_trans_coeff_types[3] = np.NPY_LONG
ufunc__get_either_trans_coeff_types[4] = np.NPY_DOUBLE
ufunc__get_either_trans_coeff_types[5] = np.NPY_DOUBLE
ufunc__get_either_trans_coeff_types[6] = np.NPY_DOUBLE
ufunc__get_either_trans_coeff_types[7] = np.NPY_BOOL
ufunc__get_either_trans_coeff_types[8] = np.NPY_LONG
ufunc__get_either_trans_coeff_types[9] = np.NPY_CDOUBLE

# types to be used for all of the both-type translation 
# coefficient retrieval ufuncs called like
# errval = func(m, n, mu, nu, r, theta, phi, r_ge_d, J, &A, &B)
# currently supported signatures: (lllldddbl_DD)
cdef char ufunc__get_both_coeff_types[11]
ufunc__get_both_coeff_types[0] = np.NPY_LONG
ufunc__get_both_coeff_types[1] = np.NPY_LONG
ufunc__get_both_coeff_types[2] = np.NPY_LONG
ufunc__get_both_coeff_types[3] = np.NPY_LONG
ufunc__get_both_coeff_types[4] = np.NPY_DOUBLE
ufunc__get_both_coeff_types[5] = np.NPY_DOUBLE
ufunc__get_both_coeff_types[6] = np.NPY_DOUBLE
ufunc__get_both_coeff_types[7] = np.NPY_BOOL
ufunc__get_both_coeff_types[8] = np.NPY_LONG
ufunc__get_both_coeff_types[9] = np.NPY_CDOUBLE
ufunc__get_both_coeff_types[10] = np.NPY_CDOUBLE


trans_A_taylor_elementwise_funcs[0] = <void*> qpms_trans_single_A_Taylor_ext
trans_B_taylor_elementwise_funcs[0] = <void*> qpms_trans_single_B_Taylor_ext

trans_A_Taylor = np.PyUFunc_FromFuncAndData(
        trans_X_taylor_loop_func, # func
        trans_A_taylor_elementwise_funcs, #data 
        ufunc__get_either_trans_coeff_types, # types
        1, # ntypes: number of supported input types
        9, # nin: number of input args
        1, # nout: number of output args
        0, # identity element, unused
        "trans_A_Taylor", # name
        """
        TODO computes the E-E or M-M translation coefficient in Taylor's normalisation
        """, # doc
        0 # unused, for backward compatibility of numpy c api
        )

trans_B_Taylor = np.PyUFunc_FromFuncAndData(
        trans_X_taylor_loop_func,
        trans_B_taylor_elementwise_funcs,
        ufunc__get_either_trans_coeff_types,
        1, # number of supported input types
        9, # number of input args
        1, # number of output args
        0, # identity element, unused
        "trans_B_Taylor",
        """
        TODO computes the E-E or M-M translation coefficient in Taylor's normalisation
        """,
        0 # unused
        )

# ---------------------------------------------
# Wrapper for the qpms_trans_calculator "class"
# ---------------------------------------------
ctypedef struct trans_calculator_get_X_data_t:
    qpms_trans_calculator* c
    void* cmethod

cdef void trans_calculator_loop_D_Ciiiidddii_As_D_lllldddbl(char **args, np.npy_intp *dims, np.npy_intp *steps, void *data) nogil:
    cdef np.npy_intp i, n = dims[0]
    cdef void *func = (<trans_calculator_get_X_data_t*>data)[0].cmethod
    #cdef cdouble (*func)(qpms_trans_calculator*, int, int, int, int, double, double, double, int, int) nogil = (<trans_calculator_get_X_data_t*>data)[0].cmethod
    cdef qpms_trans_calculator* c = (<trans_calculator_get_X_data_t*>data)[0].c
    #cdef char *func_name= <char*>(<void**>data)[1] # i am not using this, nor have I saved func_name to data
    cdef char *ip0 = args[0]
    cdef char *ip1 = args[1]
    cdef char *ip2 = args[2]
    cdef char *ip3 = args[3]
    cdef char *ip4 = args[4]
    cdef char *ip5 = args[5]
    cdef char *ip6 = args[6]
    cdef char *ip7 = args[7]
    cdef char *ip8 = args[8]
    cdef char *op0 = args[9]
    cdef cdouble ov0
    for i in range(n): # iterating over dimensions
        #ov0 = func( 
        ov0 = (<double complex(*)(qpms_trans_calculator*, int, int, int, int, double, double, double, int, int) nogil>func)(
            c,
            <int>(<np.npy_long*>ip0)[0],
            <int>(<np.npy_long*>ip1)[0],
            <int>(<np.npy_long*>ip2)[0],
            <int>(<np.npy_long*>ip3)[0],
            <double>(<np.npy_double*>ip4)[0],
            <double>(<np.npy_double*>ip5)[0],
            <double>(<np.npy_double*>ip6)[0],
            <int>(<np.npy_bool*>ip7)[0],
            <int>(<np.npy_long*>ip8)[0],
        )
        (<cdouble *>op0)[0] = <cdouble>ov0
        ip0 += steps[0]
        ip1 += steps[1]
        ip2 += steps[2]
        ip3 += steps[3]
        ip4 += steps[4]
        ip5 += steps[5]
        ip6 += steps[6]
        ip7 += steps[7]
        ip8 += steps[8]
        op0 += steps[9]
# FIXME ERROR HANDLING!!! requires correct import and different data passed (see scipy's generated ufuncs)
#    sf_error.check_fpe(func_name)


cdef void trans_calculator_loop_E_C_DD_iiiidddii_As_lllldddbl_DD(char **args, np.npy_intp *dims, np.npy_intp *steps, void *data) nogil:
    # E stands for error value (int), C for qpms_trans_calculator*
    cdef np.npy_intp i, n = dims[0]
    cdef void *func = (<trans_calculator_get_X_data_t*>data)[0].cmethod
    #cdef complex double (*func)(qpms_trans_calculator*, double complex *, double complex *, int, int, int, int, double, double, double, int, int) nogil = (<trans_calculator_get_X_data_t*>data)[0].cmethod
    cdef qpms_trans_calculator* c = (<trans_calculator_get_X_data_t*>data)[0].c
    #cdef char *func_name= <char*>(<void**>data)[1] # i am not using this, nor have I saved func_name to data
    cdef char *ip0 = args[0]
    cdef char *ip1 = args[1]
    cdef char *ip2 = args[2]
    cdef char *ip3 = args[3]
    cdef char *ip4 = args[4]
    cdef char *ip5 = args[5]
    cdef char *ip6 = args[6]
    cdef char *ip7 = args[7]
    cdef char *ip8 = args[8]
    cdef char *op0 = args[9]
    cdef char *op1 = args[10]
    cdef cdouble ov0
    cdef int errval
    for i in range(n): # iterating over dimensions
        #errval = func( 
        errval = (<int(*)(qpms_trans_calculator*, double complex *, double complex *, int, int, int, int, double, double, double, int, int) nogil>func)(
            c,
            <cdouble *> op0,
            <cdouble *> op1,
            <int>(<np.npy_long*>ip0)[0],
            <int>(<np.npy_long*>ip1)[0],
            <int>(<np.npy_long*>ip2)[0],
            <int>(<np.npy_long*>ip3)[0],
            <double>(<np.npy_double*>ip4)[0],
            <double>(<np.npy_double*>ip5)[0],
            <double>(<np.npy_double*>ip6)[0],
            <int>(<np.npy_bool*>ip7)[0],
            <int>(<np.npy_long*>ip8)[0],
        )
        ip0 += steps[0]
        ip1 += steps[1]
        ip2 += steps[2]
        ip3 += steps[3]
        ip4 += steps[4]
        ip5 += steps[5]
        ip6 += steps[6]
        ip7 += steps[7]
        ip8 += steps[8]
        op0 += steps[9]
        op1 += steps[10]
        # TODO if (errval != 0): ...
# FIXME ERROR HANDLING!!! requires correct import and different data passed (see scipy's generated ufuncs)
#    sf_error.check_fpe(func_name)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void trans_calculator_parallel_loop_E_C_DD_iiiidddii_As_lllldddbl_DD(char **args, np.npy_intp *dims, np.npy_intp *steps, void *data) nogil:
    # E stands for error value (int), C for qpms_trans_calculator*
    cdef np.npy_intp i, n = dims[0]
    cdef void *func = (<trans_calculator_get_X_data_t*>data)[0].cmethod
    #cdef complex double (*func)(qpms_trans_calculator*, double complex *, double complex *, int, int, int, int, double, double, double, int, int) nogil = (<trans_calculator_get_X_data_t*>data)[0].cmethod
    cdef qpms_trans_calculator* c = (<trans_calculator_get_X_data_t*>data)[0].c
    #cdef char *func_name= <char*>(<void**>data)[1] # i am not using this, nor have I saved func_name to data

    cdef char *ip0
    cdef char *ip1
    cdef char *ip2
    cdef char *ip3
    cdef char *ip4
    cdef char *ip5
    cdef char *ip6
    cdef char *ip7
    cdef char *ip8
    cdef char *op0
    cdef char *op1
    cdef int errval
    for i in prange(n): # iterating over dimensions
        ip0 = args[0] + i * steps[0]
        ip1 = args[1] + i * steps[1]
        ip2 = args[2] + i * steps[2]
        ip3 = args[3] + i * steps[3]
        ip4 = args[4] + i * steps[4]
        ip5 = args[5] + i * steps[5]
        ip6 = args[6] + i * steps[6]
        ip7 = args[7] + i * steps[7]
        ip8 = args[8] + i * steps[8]
        op0 = args[9] + i * steps[9]
        op1 = args[10] + i * steps[10]
        #errval = func( 
        errval = (<int(*)(qpms_trans_calculator*, double complex *, double complex *, int, int, int, int, double, double, double, int, int) nogil>func)(
            c,
            <cdouble *> op0,
            <cdouble *> op1,
            <int>(<np.npy_long*>ip0)[0],
            <int>(<np.npy_long*>ip1)[0],
            <int>(<np.npy_long*>ip2)[0],
            <int>(<np.npy_long*>ip3)[0],
            <double>(<np.npy_double*>ip4)[0],
            <double>(<np.npy_double*>ip5)[0],
            <double>(<np.npy_double*>ip6)[0],
            <int>(<np.npy_bool*>ip7)[0],
            <int>(<np.npy_long*>ip8)[0],
        )
       # TODO if (errval != 0): ...
# FIXME ERROR HANDLING!!! requires correct import and different data passed (see scipy's generated ufuncs)
#    sf_error.check_fpe(func_name)


cdef np.PyUFuncGenericFunction trans_calculator_get_X_loop_funcs[1]
trans_calculator_get_X_loop_funcs[0] = trans_calculator_loop_D_Ciiiidddii_As_D_lllldddbl

cdef np.PyUFuncGenericFunction trans_calculator_get_AB_loop_funcs[1]
#trans_calculator_get_AB_loop_funcs[0] = trans_calculator_parallel_loop_E_C_DD_iiiidddii_As_lllldddbl_DD
trans_calculator_get_AB_loop_funcs[0] = trans_calculator_loop_E_C_DD_iiiidddii_As_lllldddbl_DD
cdef void *trans_calculator_get_AB_elementwise_funcs[1]
trans_calculator_get_AB_elementwise_funcs[0] = <void *>qpms_trans_calculator_get_AB_p_ext

cdef class trans_calculator:
    cdef qpms_trans_calculator* c
    cdef trans_calculator_get_X_data_t get_A_data[1]
    cdef trans_calculator_get_X_data_t* get_A_data_p[1]

    cdef trans_calculator_get_X_data_t get_B_data[1]
    cdef trans_calculator_get_X_data_t* get_B_data_p[1]

    cdef trans_calculator_get_X_data_t get_AB_data[1]
    cdef trans_calculator_get_X_data_t* get_AB_data_p[1]
    cdef public: # TODO CHECK FOR CORRECT REFERENCE COUNTING AND LEAKS
        # have to be cdef public in order that __init__ can set these attributes
        object get_A, get_B, get_AB

    def __cinit__(self, int lMax, int normalization = 1):
        self.c = qpms_trans_calculator_init(lMax, normalization)

    def __init__(self, int lMax, int normalization = 1):
        self.get_A_data[0].c = self.c
        self.get_A_data[0].cmethod = <void *>qpms_trans_calculator_get_A_ext
        self.get_A_data_p[0] = &(self.get_A_data[0])
        self.get_A = <object>np.PyUFunc_FromFuncAndData(# TODO CHECK FOR CORRECT REFERENCE COUNTING AND LEAKS
                trans_calculator_get_X_loop_funcs, # func
                <void **>self.get_A_data_p, #data
                ufunc__get_either_trans_coeff_types, #types
                1, # ntypes: number of supported input types
                9, # nin: number of input args
                1, # nout: number of output args
                0, # identity element, unused
                "get_A", #name
                """
                TODO doc
                """, # doc
                0 # unused
                )
        self.get_B_data[0].c = self.c
        self.get_B_data[0].cmethod = <void *>qpms_trans_calculator_get_B_ext
        self.get_B_data_p[0] = &(self.get_B_data[0])
        self.get_B = <object>np.PyUFunc_FromFuncAndData(# TODO CHECK FOR CORRECT REFERENCE COUNTING AND LEAKS
                trans_calculator_get_X_loop_funcs, # func
                <void **>self.get_B_data_p, #data
                ufunc__get_either_trans_coeff_types, #types
                1, # ntypes: number of supported input types
                9, # nin: number of input args
                1, # nout: number of output args
                0, # identity element, unused
                "get_B", #name
                """
                TODO doc
                """, # doc
                0 # unused
                )
        self.get_AB_data[0].c = self.c
        self.get_AB_data[0].cmethod = <void *>qpms_trans_calculator_get_AB_p_ext
        self.get_AB_data_p[0] = &(self.get_AB_data[0])
        self.get_AB = <object>np.PyUFunc_FromFuncAndData(# TODO CHECK FOR CORRECT REFERENCE COUNTING AND LEAKS
                trans_calculator_get_AB_loop_funcs, # func
                <void **>self.get_AB_data_p, #data
                ufunc__get_both_coeff_types, #types
                1, # ntypes: number of supported input types
                9, # nin: number of input args
                2, # nout: number of output args
                0, # identity element, unused
                "get_AB", #name
                """
                TODO doc
                """, # doc
                0 # unused
                )
    def __dealloc__(self):
        qpms_trans_calculator_free(self.c)
        # TODO Reference counts to get_A, get_B, get_AB?


    # TODO make possible to access the attributes (to show normalization etc)

    

