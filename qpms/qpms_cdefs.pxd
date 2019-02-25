cimport numpy as np

cdef extern from "qpms_types.h":
    cdef struct cart3_t:
        double x
        double y
        double z
    cdef struct cart2_t:
        double x
        double y
    cdef struct sph_t:
        double r
        double theta
        double phi
    cdef struct pol_t:
        double r
        double phi
    cdef union anycoord_point_t:
        double z
        cart3_t cart3
        cart2_t cart2
        pol_t pol
    ctypedef enum qpms_normalisation_t:
        QPMS_NORMALISATION_XU
        QPMS_NORMALISATION_XU_CS
        QPMS_NORMALISATION_NONE
        QPMS_NORMALISATION_NONE_CS
        QPMS_NORMALISATION_KRISTENSSON
        QPMS_NORMALISATION_KRISTENSSON_CS
        QPMS_NORMALISATION_POWER
        QPMS_NORMALISATION_POWER_CS
        QPMS_NORMALISATION_TAYLOR
        QPMS_NORMALISATION_TAYLOR_CS
        QPMS_NORMALISATION_SPHARM
        QPMS_NORMALISATION_SPHARM_CS
        QPMS_NORMALISATION_UNDEF
    ctypedef int qpms_lm_t
    ctypedef int qpms_l_t
    ctypedef int qpms_m_t
    # maybe more if needed

# Point generators from lattices.h
cdef extern from "lattices.h":
    ctypedef enum PGenPointFlags:
        pass
    struct PGenReturnData:
        pass
    struct PGenZReturnData:
        pass
    struct PGenPolReturnData:
        pass
    struct PGenSphReturnData:
        pass
    struct PGenCart2ReturnData:
        pass
    struct PGenCart3ReturnData:
        pass
    struct PGenClassInfo: # maybe important
        pass
    struct PGen: # probably important
        PGenClassInfo* c
        void *statedata
    void PGen_destroy(PGen *g)

    # now the individual PGen implementations:
    # FIXME Is bint always guaranteed to be equivalent to _Bool? (I dont't think so.)
    PGen PGen_xyWeb_new(cart2_t b1, cart2_t b2, double rtol, cart2_t offset,
            double minR, bint inc_minR, double maxR, bint inc_maxR)
    ctypedef enum PGen_1D_incrementDirection:
        PGEN_1D_INC_FROM_ORIGIN
        PGEN_1D_INC_TOWARDS_ORIGIN
    PGen PGen_1D_new_minMaxR(double period, double offset, double minR, bint inc_minR,
            double maxR, bint inc_maxR, PGen_1D_incrementDirection incdir)

ctypedef double complex cdouble

cdef extern from "wigner.h":
    struct qpms_quat_t:
        cdouble a
        cdouble b
    struct qpms_quat4d_t:
        double c1
        double ci
        double cj
        double ck
    qpms_quat_t qpms_quat_2c_from_4d(qpms_quat4d_t q)
    qpms_quat4d_t qpms_quat_4d_from_2c(qpms_quat_t q)
    qpms_quat_t qpms_quat_mult(qpms_quat_t p, qpms_quat_t q)
    qpms_quat_t qpms_quat_add(qpms_quat_t p, qpms_quat_t q)
    qpms_quat_t qpms_quat_rscale(double s, qpms_quat_t q) 
    qpms_quat_t qpms_quat_conj(qpms_quat_t q) 
    double qpms_quat_norm(qpms_quat_t q) 
    double qpms_quat_imnorm(qpms_quat_t q)
    qpms_quat_t qpms_quat_normalise(qpms_quat_t q) 
    qpms_quat_t qpms_quat_log(qpms_quat_t q)
    qpms_quat_t qpms_quat_exp(qpms_quat_t q)
    qpms_quat_t qpms_quat_pow(qpms_quat_t q, double exponent)
    cdouble qpms_wignerD_elem(qpms_quat_t q, qpms_l_t l,
                           qpms_m_t mp, qpms_m_t m)


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
        int lMax
        size_t nelem
        cdouble** A_multipliers
        cdouble** B_multipliers
    enum qpms_normalization_t:
        pass
    qpms_trans_calculator* qpms_trans_calculator_init(int lMax, int nt) # should be qpms_normalization_t
    void qpms_trans_calculator_free(qpms_trans_calculator* c)
    cdouble qpms_trans_calculator_get_A_ext(const qpms_trans_calculator* c,
            int m, int n, int mu, int nu, double kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J) nogil
    cdouble qpms_trans_calculator_get_B_ext(const qpms_trans_calculator* c,
            int m, int n, int mu, int nu, double kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J) nogil
    int qpms_trans_calculator_get_AB_p_ext(const qpms_trans_calculator* c,
            cdouble *Adest, cdouble *Bdest,
            int m, int n, int mu, int nu, double kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J) nogil
    int qpms_trans_calculator_get_AB_arrays_ext(const qpms_trans_calculator *c,
            cdouble *Adest, cdouble *Bdest,
            size_t deststride, size_t srcstride,
            double kdlj_r, double kdlj_theta, double kdlj_phi,
            int r_ge_d, int J) nogil
    int qpms_cython_trans_calculator_get_AB_arrays_loop(qpms_trans_calculator *c,
            int J, int resnd,
            int daxis, int saxis,
            char *A_data, np.npy_intp *A_shape, np.npy_intp *A_strides,
            char *B_data, np.npy_intp *B_shape, np.npy_intp *B_strides,
            char *r_data, np.npy_intp *r_shape, np.npy_intp *r_strides,
            char *theta_data, np.npy_intp *theta_shape, np.npy_intp *theta_strides,
            char *phi_data, np.npy_intp *phi_shape, np.npy_intp *phi_strides,
            char *r_ge_d_data, np.npy_intp *phi_shape, np.npy_intp *phi_strides) nogil



