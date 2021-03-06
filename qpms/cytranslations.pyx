import numpy as np
import cmath
from .qpms_cdefs cimport *
from .cycommon import *
from .cybspec cimport *
cimport cython
from cython.parallel cimport parallel, prange
from libc.stdlib cimport malloc, free, calloc, abort

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


# This one is probably not used anymore an can perhaps be removed:
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

    def __cinit__(self, int lMax, int normalization = QPMS_NORMALISATION_DEFAULT):
        if (lMax <= 0):
            raise ValueError('lMax has to be greater than 0.')
        self.c = qpms_trans_calculator_init(lMax, normalization)
        if self.c is NULL:
            raise MemoryError

    def __init__(self, int lMax, int normalization = QPMS_NORMALISATION_DEFAULT):
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

    # THIS FUNCTION MIGHT BE OBSOLETE; NOT SURE WHETHER IT'S WORKING ANY MORE
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
    
    def get_trans_array_bspec_sph(self, BaseSpec destspec, BaseSpec srcspec,
            kdlj, qpms_bessel_t J = QPMS_HANKEL_PLUS):
        kdlj = np.array(kdlj)
        if kdlj.shape != (3,):
            raise ValueError("Array of shape (3,) with spherical coordinates of the translation expected")
        cdef size_t destn = len(destspec)
        cdef size_t srcn = len(srcspec)
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (destn, srcn), dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        cdef csph_t kdlj_sph
        kdlj_sph.r = kdlj[0]
        kdlj_sph.theta = kdlj[1].real
        kdlj_sph.phi = kdlj[2].real
        qpms_trans_calculator_get_trans_array(self.c, &target_view[0][0], 
                destspec.rawpointer(), srcn, srcspec.rawpointer(), 1, 
                kdlj_sph, False, J)
        return target

    def get_trans_array_bspec_c3pos(self, BaseSpec destspec, BaseSpec srcspec,
            cdouble k, destpos, srcpos, qpms_bessel_t J = QPMS_HANKEL_PLUS):
        destpos = np.array(destpos)
        srcpos = np.array(srcpos)
        if destpos.shape != (3,) or srcpos.shape != (3,):
            raise ValueError("Array of shape (3,) with cartesian coordinates of the particle position expected")
        cdef size_t destn = len(destspec)
        cdef size_t srcn = len(srcspec)
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (destn, srcn), dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        cdef cart3_t srcp, destp
        srcp.x = srcpos[0]
        srcp.y = srcpos[1]
        srcp.z = srcpos[2]
        destp.x = destpos[0]
        destp.y = destpos[1]
        destp.z = destpos[2]
        qpms_trans_calculator_get_trans_array_lc3p(self.c, &target_view[0][0], 
                destspec.rawpointer(), srcn, srcspec.rawpointer(), 1, k,
                destp, srcp, J)
        return target


    # TODO make possible to access the attributes (to show normalization etc)

