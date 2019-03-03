cimport numpy as np

ctypedef double complex cdouble

from libc.stdint cimport *

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
    struct qpms_quat_t:
        cdouble a
        cdouble b
    struct qpms_quat4d_t:
        double c1
        double ci
        double cj
        double ck
    struct qpms_irot3_t:
        qpms_quat_t rot
        short det
    ctypedef np.ulonglong_t qpms_uvswfi_t
    struct qpms_vswf_set_spec_t:
        size_t n
        qpms_uvswfi_t *ilist
        qpms_l_t lMax
        qpms_l_t lMax_M
        qpms_l_t lMax_N
        qpms_l_t lMax_L
        size_t capacity
        qpms_normalisation_t norm
    ctypedef enum qpms_errno_t:
        QPMS_SUCCESS
        QPMS_ERROR
        # more if needed
    ctypedef enum qpms_vswf_type_t:
        QPMS_VSWF_ELECTRIC
        QPMS_VSWF_MAGNETIC
        QPMS_VSWF_LONGITUDINAL
    ctypedef int32_t qpms_ss_tmi_t
    ctypedef int32_t qpms_ss_pi_t
    ctypedef int qpms_gmi_t
    ctypedef int qpms_iri_t
    ctypedef const char * qpms_permutation_t
    # maybe more if needed

cdef extern from "indexing.h":
    qpms_uvswfi_t qpms_tmn2uvswfi(qpms_vswf_type_t t, qpms_m_t m, qpms_l_t n)
    qpms_errno_t qpms_uvswfi2tmn(qpms_uvswfi_t u, qpms_vswf_type_t* t, qpms_m_t* m, qpms_l_t* n)
    qpms_m_t qpms_uvswfi2m(qpms_uvswfi_t u)
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


cdef extern from "wigner.h":
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
    qpms_irot3_t qpms_irot3_mult(qpms_irot3_t p, qpms_irot3_t q)
    qpms_irot3_t qpms_irot3_pow(qpms_irot3_t p, int n)

cdef extern from "groups.h":
    struct qpms_finite_group_irrep_t:
        int dim
        char *name
        cdouble *m
    struct qpms_finite_group_t:
        char *name
        qpms_gmi_t order
        qpms_gmi_t idi
        qpms_gmi_t *mt
        qpms_gmi_t *invi
        qpms_gmi_t *gens
        int ngens
        qpms_permutation_t *permrep
        char **elemlabels
        int permrep_nelem
        qpms_irot3_t *rep3d
        qpms_iri_t nirreps
        qpms_finite_group_irrep_t *irreps
    qpms_finite_group_t QPMS_FINITE_GROUP_TRIVIAL
    qpms_finite_group_t QPMS_FINITE_GROUP_TRIVIAL_G

cdef extern from "symmetries.h":
    cdouble *qpms_zflip_uvswi_dense(cdouble *target, const qpms_vswf_set_spec_t *bspec)
    cdouble *qpms_yflip_uvswi_dense(cdouble *target, const qpms_vswf_set_spec_t *bspec)
    cdouble *qpms_xflip_uvswi_dense(cdouble *target, const qpms_vswf_set_spec_t *bspec)
    cdouble *qpms_zrot_uvswi_dense(cdouble *target, const qpms_vswf_set_spec_t *bspec, double phi)
    cdouble *qpms_zrot_rational_uvswi_dense(cdouble *target, const qpms_vswf_set_spec_t *bspec, int N, int w)
    cdouble *qpms_irot3_uvswfi_dense(cdouble *target, const qpms_vswf_set_spec_t *bspec, qpms_irot3_t transf)

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

cdef extern from "scatsystem.h":
    struct qpms_tmatrix_t:
        qpms_vswf_set_spec_t *spec
        cdouble *m
        bint owns_m # FIXME in fact bool
    struct qpms_particle_t:
        cart3_t pos
        const qpms_tmatrix_t *tmatrix
    struct qpms_particle_tid_t:
        cart3_t pos
        qpms_ss_tmi_t tmatrix_id
    struct qpms_tmatrix_interpolator_t:
        const qpms_vswf_set_spec_t *bspec
    struct qpms_scatsys_t:
        qpms_tmatrix_t **tm
        qpms_ss_tmi_t tm_count
        qpms_particle_tid_t *p
        qpms_ss_pi_t p_count
        # We shouldn't need more to construct a symmetric scatsystem
    qpms_scatsys_t *qpms_scatsys_apply_symmetry(const qpms_scatsys_t *orig, const qpms_finite_group_t *sym)
    void qpms_scatsys_free(qpms_scatsys_t *s)
    qpms_errno_t qpms_scatsys_dump(qpms_scatsys_t *ss, char *path) #NI
    qpms_scatsys_t *qpms_scatsys_load(char *path) #NI
    qpms_tmatrix_isclose(const qpms_tmatrix_t *A, const qpms_tmatrix_t *B,
                const double rtol, const double atol)
        
        



