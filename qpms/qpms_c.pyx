# Cythonized parts of QPMS here
# -----------------------------

import numpy as np
import cmath
from qpms_cdefs cimport *
cimport cython
from cython.parallel cimport parallel, prange

import math # for copysign in crep methods
#import re # TODO for crep methods?

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

def get_nelem(unsigned int lMax):
    return lMax * (lMax + 2)

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

'''
cdef extern from "numpy/ndarrayobject.h":
    struct PyArrayInterface:
        int itemsize
        np.npy_uintp *shape
        np.npy_uintp *strides
        void *data
'''


from libc.stdlib cimport malloc, free, calloc, abort



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
        if (lMax <= 0):
            raise ValueError('lMax has to be greater than 0.')
        self.c = qpms_trans_calculator_init(lMax, normalization)
        if self.c is NULL:
            raise MemoryError

    def __init__(self, int lMax, int normalization = 1):
        if self.c is NULL:
            raise MemoryError()
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
        if self.c is not NULL:
            qpms_trans_calculator_free(self.c)
        # TODO Reference counts to get_A, get_B, get_AB?

    def lMax(self):
        return self.c[0].lMax

    def nelem(self):
        return self.c[0].nelem

    def get_AB_arrays(self, r, theta, phi, r_ge_d, int J, 
            destaxis=None, srcaxis=None, expand=True):
        """
        Returns arrays of translation coefficients, inserting two new nelem-sized axes
        (corresponding to the destination and source axes of the translation matrix,
        respectively).

        By default (expand==True), it inserts the new axes. or it can be provided with
        the resulting shape (with the corresponding axes dimensions equal to one).
        The provided axes positions are for the resulting array.

        If none axis positions are provided, destaxis and srcaxis will be the second-to-last
        and last, respectively.
        """
        # TODO CHECK (and try to cast) INPUT ARRAY TYPES (now is done)
        # BIG FIXME: make skalars valid arguments, now r, theta, phi, r_ge_d have to be ndarrays
        cdef:
            int daxis, saxis, smallaxis, bigaxis, resnd, i, j, d, ax, errval
            np.npy_intp sstride, dstride, longi
            int *local_indices
            char *r_p
            char *theta_p
            char *phi_p
            char *r_ge_d_p
            char *a_p
            char *b_p
        # Process the array shapes
        baseshape = np.broadcast(r,theta,phi,r_ge_d).shape # nope, does not work as needed
        '''
        cdef int r_orignd = r.ndim if hasattr(r, "ndim") else 0
        cdef int theta_orignd = theta.ndim if hasattr(theta, "ndim") else 0
        cdef int phi_orignd = phi.ndim if hasattr(phi, "ndim") else 0
        cdef int r_ge_d_orignd = r_ge_d.ndim if hasattr(r_ge_d, "__len__") else 0
        cdef int basend = max(r_orignd, theta_orignd, phi_orignd, r_ge_d_orignd)
        baseshape = list()
        for d in range(basend):
            baseshape.append(max(
                    r.shape[d+r_orignd-basend] if d+r_orignd-basend >= 0 else 1,
                    theta.shape[d+theta_orignd-basend] if d+theta_orignd-basend >= 0 else 1,
                    phi.shape[d+phi_orignd-basend] if d+phi_orignd-basend >= 0 else 1,
                    r_ge_d.shape[d+r_ge_d_orignd-basend] if d+r_ge_d_orignd-basend >= 0 else 1,
                    ))
        baseshape = tuple(baseshape)
        '''
        if not expand:
            resnd = len(baseshape)
            if resnd < 2:
                raise ValueError('Translation matrix arrays must have at least 2 dimensions!')
            daxis = (resnd-2) if destaxis is None else destaxis
            saxis = (resnd-1) if srcaxis is None else srcaxis
            if daxis < 0:
                daxis = resnd + daxis
            if saxis < 0:
                saxis = resnd + saxis
            if daxis < 0 or saxis < 0 or daxis >= resnd or saxis >= resnd or daxis == saxis:
                raise ValueError('invalid axes provided (destaxis = %d, srcaxis = %d, # of axes: %d'
                        % (daxis, saxis, resnd)) 
            if baseshape[daxis] != 1 or baseshape[saxis] != 1:
                raise ValueError('dimension mismatch (input argument dimensions have to be 1 both at'
                        'destaxis (==%d) and srcaxis (==%d) but are %d and %d' % 
                        (daxis, saxis, baseshape[daxis], baseshape[saxis]))
            resultshape = list(baseshape)
        else:
            resnd = len(baseshape)+2
            daxis = (resnd-2) if destaxis is None else destaxis
            saxis = (resnd-1) if srcaxis is None else srcaxis
            if daxis < 0:
                daxis = resnd + daxis
            if saxis < 0:
                saxis = resnd + saxis
            if daxis < 0 or saxis < 0 or daxis >= resnd or saxis >= resnd or daxis == saxis:
                raise ValueError('invalid axes provided') # TODO better error formulation
            resultshape = list(baseshape)
            if daxis > saxis:
                smallaxis = saxis
                bigaxis = daxis
            else:
                smallaxis = daxis
                bigaxis = saxis
            resultshape.insert(smallaxis,1)
            resultshape.insert(bigaxis,1)
            r = np.expand_dims(np.expand_dims(r.astype(np.float_, copy=False), smallaxis), bigaxis)
            theta = np.expand_dims(np.expand_dims(theta.astype(np.float_, copy=False), smallaxis), bigaxis)
            phi = np.expand_dims(np.expand_dims(phi.astype(np.float_, copy=False), smallaxis), bigaxis)
            r_ge_d = np.expand_dims(np.expand_dims(r_ge_d.astype(np.bool_, copy=False), smallaxis), bigaxis)

        resultshape[daxis] = self.c[0].nelem
        resultshape[saxis] = self.c[0].nelem
        cdef np.ndarray r_c = np.broadcast_to(r,resultshape)
        cdef np.ndarray theta_c = np.broadcast_to(theta,resultshape)
        cdef np.ndarray phi_c = np.broadcast_to(phi,resultshape)
        cdef np.ndarray r_ge_d_c = np.broadcast_to(r_ge_d, resultshape)
        cdef np.ndarray a = np.empty(resultshape, dtype=complex)
        cdef np.ndarray b = np.empty(resultshape, dtype=complex)
        dstride = a.strides[daxis]
        sstride = a.strides[saxis]
        with nogil: 
            errval = qpms_cython_trans_calculator_get_AB_arrays_loop(
                    self.c, J, resnd,
                    daxis, saxis,
                    a.data, a.shape, a.strides,
                    b.data, b.shape, b.strides,
                    r_c.data, r_c.shape, r_c.strides,
                    theta_c.data, theta_c.shape, theta_c.strides,
                    phi_c.data, phi_c.shape, phi_c.strides,
                    r_ge_d_c.data, r_ge_d_c.shape, r_ge_d_c.strides
                    )
        return a, b

    # TODO make possible to access the attributes (to show normalization etc)


def complex_crep(complex c, parentheses = False, shortI = True, has_Imaginary = False):
    '''
    Return a C-code compatible string representation of a (python) complex number.
    '''
    return ( ('(' if parentheses else '')
            + repr(c.real)
            + ('+' if math.copysign(1, c.imag) >= 0 else '')
            + repr(c.imag)
            + ('*I' if shortI else '*_Imaginary_I' if has_Imaginary else '*_Complex_I')
            + (')' if parentheses else '')
        )

# Quaternions from wigner.h 
# (mainly for testing; use moble's quaternions in python)

cdef class cquat:
    '''
    Wrapper of the qpms_quat_t object, with the functionality
    to evaluate Wigner D-matrix elements.
    '''
    cdef readonly qpms_quat_t q

    def __cinit__(self, double w, double x, double y, double z):
        cdef qpms_quat4d_t p
        p.c1 = w
        p.ci = x
        p.cj = y
        p.ck = z
        self.q = qpms_quat_2c_from_4d(p)

    def copy(self):
        res = cquat(0,0,0,0)
        res.q = self.q
        return res

    def __repr__(self): # TODO make this look like a quaternion with i,j,k
        return repr(self.r)

    def __add__(cquat self, cquat other):
        # TODO add real numbers
        res = cquat(0,0,0,0)
        res.q = qpms_quat_add(self.q, other.q)
        return res

    def __mul__(cquat self, cquat other):
        res = cquat(0,0,0,0)
        res.q = qpms_quat_mult(self.q, other.q)
        return res

    def __neg__(cquat self):
        res = cquat(0,0,0,0)
        res.q = qpms_quat_rscale(-1, self.q)
        return res

    def __sub__(cquat self, cquat other):
        res = cquat(0,0,0,0)
        res.q = qpms_quat_add(self.q, qpms_quat_rscale(-1,other.q))
        return res

    def __abs__(self):
        return qpms_quat_norm(self.q)

    def norm(self):
        return qpms_quat_norm(self.q)

    def imnorm(self):
        return qpms_quat_imnorm(self.q)

    def exp(self):
        res = cquat(0,0,0,0)
        res.q = qpms_quat_exp(self.q)
        return res

    def log(self):
        res = cquat(0,0,0,0)
        res.q = qpms_quat_exp(self.q)
        return res

    def __pow__(cquat self, double other, _):
        res = cquat(0,0,0,0)
        res.q = qpms_quat_pow(self.q, other)
        return res

    def normalise(self):
        res = cquat(0,0,0,0)
        res.q = qpms_quat_normalise(self.q)
        return res

    def isclose(cquat self, cquat other, rtol=1e-5, atol=1e-8):
        '''
        Checks whether two quaternions are "almost equal".
        '''
        return abs(self - other) <= (atol + rtol * abs(other))

    property c:
        '''
        Quaternion representation as two complex numbers
        '''
        def __get__(self):
            return (self.q.a, self.q.b)
        def __set__(self, RaRb):
            self.q.a = RaRb[0]
            self.q.b = RaRb[1]

    property r:
        '''
        Quaternion representation as four real numbers
        '''
        def __get__(self):
            cdef qpms_quat4d_t p
            p = qpms_quat_4d_from_2c(self.q)
            return (p.c1, p.ci, p.cj, p.ck)
        def __set__(self, wxyz):
            cdef qpms_quat4d_t p
            p.c1 = wxyz[0]
            p.ci = wxyz[1]
            p.cj = wxyz[2]
            p.ck = wxyz[3]
            self.q = qpms_quat_2c_from_4d(p)

    def crepr(self):
        '''
        Returns a string that can be used in C code to initialise a qpms_irot3_t
        '''
        return '{' + complex_crep(self.q.a) + ', ' + complex_crep(self.q.b)  + '}'

    def wignerDelem(self, qpms_l_t l, qpms_m_t mp, qpms_m_t m):
        '''
        Returns an element of a bosonic Wigner matrix.
        '''
        # don't crash on bad l, m here
        if (abs(m) > l or abs(mp) > l):
            return 0
        return qpms_wignerD_elem(self.q, l, mp, m)

cdef class irot3:
    '''
    Wrapper over the C type qpms_irot3_t.
    '''
    cdef qpms_irot3_t qd

    def __cinit__(self, cquat q, short det): 
        # TODO implement a constructor with
        #  - no arguments (returns identity)
        #  - irot3 as argument (makes a copy)
        #  - cquat as argument (returns a corresponding proper rotation)
        #  - tuple as argument ...?
        if (det != 1 and det != -1):
            raise ValueError("Improper rotation determinant has to be 1 or -1")
        self.qd.rot = q.normalise().q
        self.qd.det = det

    def copy(self):
        res = irot3(cquat(1,0,0,0),1)
        res.qd = self.qd
        return res

    property rot:
        '''
        The proper rotation part of the irot3 type.
        '''
        def __get__(self):
            res = cquat(0,0,0,0)
            res.q = self.qd.rot
            return res
        def __set__(self, cquat r):
            # TODO check for non-zeroness and throw an exception if norm is zero
            self.qd.rot = r.normalise().q

    property det:
        '''
        The determinant of the improper rotation.
        '''
        def __get__(self):
            return self.qd.det
        def __set__(self, d):
            d = int(d)
            if (d != 1 and d != -1):
                raise ValueError("Improper rotation determinant has to be 1 or -1")
            self.qd.det = d

    def __repr__(self): # TODO make this look like a quaternion with i,j,k
        return '(' + repr(self.rot) + ', ' + repr(self.det) + ')'

    def crepr(self):
        '''
        Returns a string that can be used in C code to initialise a qpms_irot3_t
        '''
        return '{' + self.rot.crepr() + ', ' + repr(self.det) + '}'

    def __mul__(irot3 self, irot3 other):
        res = irot3(cquat(1,0,0,0), 1) 
        res.qd = qpms_irot3_mult(self.qd, other.qd)
        return res

    def __pow__(irot3 self, n, _):
        cdef int nint
        if (n % 1 == 0):
            nint = n
        else:
            raise ValueError("The exponent of an irot3 has to have an integer value.")
        res = irot3(cquat(1,0,0,0), 1)
        res.qd = qpms_irot3_pow(self.qd, n)
        return res

    def isclose(irot3 self, irot3 other, rtol=1e-5, atol=1e-8):
        '''
        Checks whether two (improper) rotations are "almost equal".
        Returns always False if the determinants are different.
        '''
        if self.det != other.det: 
            return False
        return self.rot.isclose(other.rot, rtol=rtol, atol=atol)

