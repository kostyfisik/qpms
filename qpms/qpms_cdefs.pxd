cimport numpy as np

ctypedef double complex cdouble

from libc.stdint cimport *

cdef extern from "gsl/gsl_errno.h":
    ctypedef void gsl_error_handler_t (const char *reason, const char *file,
            int line, int gsl_errno)
    gsl_error_handler_t *gsl_set_error_handler(gsl_error_handler_t *new_handler)
    gsl_error_handler_t *gsl_set_error_handler_off();

cdef extern from "gsl/gsl_const_mksa.h":
    const double GSL_CONST_MKSA_SPEED_OF_LIGHT

cdef extern from "qpms_types.h":
    cdef struct cart3_t:
        double x
        double y
        double z
    cdef struct ccart3_t:
        cdouble x
        cdouble y
        cdouble z
    cdef struct cart2_t:
        double x
        double y
    cdef struct sph_t:
        double r
        double theta
        double phi
    cdef struct csph_t:
        cdouble r
        double theta
        double phi
    cdef struct csphvec_t:
        cdouble rc
        cdouble thetac
        cdouble phic
    cdef struct pol_t:
        double r
        double phi
    cdef union anycoord_point_t:
        double z
        cart3_t cart3
        cart2_t cart2
        pol_t pol
    ctypedef enum qpms_normalisation_t:
        QPMS_NORMALISATION_UNDEF
        QPMS_NORMALISATION_INVERSE
        QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE
        QPMS_NORMALISATION_SPHARM_REAL
        QPMS_NORMALISATION_CSPHASE
        QPMS_NORMALISATION_M_I
        QPMS_NORMALISATION_M_MINUS
        QPMS_NORMALISATION_N_I
        QPMS_NORMALISATION_N_MINUS
        QPMS_NORMALISATION_L_I
        QPMS_NORMALISATION_L_MINUS
        QPMS_NORMALISATION_NORM_BITSTART
        QPMS_NORMALISATION_NORM_POWER
        QPMS_NORMALISATION_NORM_SPHARM
        QPMS_NORMALISATION_NORM_NONE
        QPMS_NORMALISATION_NORM_BITS
        QPMS_NORMALISATION_CONVENTION_KRISTENSSON_REAL
        QPMS_NORMALISATION_CONVENTION_KRISTENSSON
        QPMS_NORMALISATION_CONVENTION_SCUFF
        QPMS_NORMALISATION_DEFAULT
    ctypedef enum qpms_bessel_t:
        QPMS_BESSEL_REGULAR
        QPMS_BESSEL_SINGULAR
        QPMS_HANKEL_PLUS
        QPMS_HANKEL_MINUS
        QPMS_BESSEL_UNDEF
    ctypedef int qpms_lm_t
    ctypedef int qpms_l_t
    ctypedef int qpms_m_t
    ctypedef size_t qpms_y_t
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
    ctypedef int32_t qpms_ss_tmgi_t
    ctypedef int32_t qpms_ss_tmi_t
    ctypedef int32_t qpms_ss_pi_t
    ctypedef int qpms_gmi_t
    ctypedef int qpms_iri_t
    qpms_iri_t QPMS_NO_IRREP
    ctypedef const char * qpms_permutation_t
    struct qpms_tmatrix_t:
        qpms_vswf_set_spec_t *spec
        cdouble *m
        bint owns_m # FIXME in fact bool
    ctypedef enum qpms_pointgroup_class:
        QPMS_PGS_CN
        QPMS_PGS_S2N
        QPMS_PGS_CNH
        QPMS_PGS_CNV
        QPMS_PGS_DN
        QPMS_PGS_DND
        QPMS_PGS_DNH
        QPMS_PGS_T
        QPMS_PGS_TD
        QPMS_PGS_TH
        QPMS_PGS_O
        QPMS_PGS_OH
        QPMS_PGS_I
        QPMS_PGS_IH
        QPMS_PGS_CINF
        QPMS_PGS_CINFH
        QPMS_PGS_CINFV
        QPMS_PGS_DINF
        QPMS_PGS_DINFH
        QPMS_PGS_SO3
        QPMS_PGS_O3
    struct qpms_pointgroup_t:
        qpms_pointgroup_class c
        qpms_gmi_t n
        qpms_irot3_t orientation
    struct qpms_epsmu_t:
        cdouble eps
        cdouble mu
    # maybe more if needed

cdef extern from "qpms_error.h":
    ctypedef enum qpms_dbgmsg_flags:
        QPMS_DBGMSG_MISC
        QPMS_DBGMSG_THREADS
        QPMS_DBGMSG_INTEGRATION
    qpms_dbgmsg_flags qpms_dbgmsg_enable(qpms_dbgmsg_flags types)
    qpms_dbgmsg_flags qpms_dbgmsg_disable(qpms_dbgmsg_flags types)

cdef extern from "tolerances.h":
    struct qpms_tolerance_spec_t:
        pass # TODO
    const qpms_tolerance_spec_t QPMS_TOLERANCE_DEFAULT


# This is due to the fact that cython apparently cannot nest the unnamed struct/unions in an obvious way
ctypedef union qpms_incfield_planewave_params_k:
    ccart3_t cart
    csph_t sph
ctypedef union qpms_incfield_planewave_params_E:
    ccart3_t cart
    csphvec_t sph

cdef extern from "vswf.h":
    ctypedef qpms_errno_t (*qpms_incfield_t)(cdouble *target, const qpms_vswf_set_spec_t *bspec,
            const cart3_t evalpoint, const void *args, bint add)
    ctypedef struct qpms_incfield_planewave_params_t:
        bint use_cartesian
        qpms_incfield_planewave_params_k k
        qpms_incfield_planewave_params_E E
    qpms_errno_t qpms_incfield_planewave(cdouble *target, const qpms_vswf_set_spec_t *bspec,
            const cart3_t evalpoint, const void *args, bint add)
    csphvec_t qpms_vswf_single_el_csph(qpms_m_t m, qpms_l_t n, csph_t kdlj, qpms_bessel_t btyp, qpms_normalisation_t norm)
    csphvec_t qpms_vswf_single_mg_csph(qpms_m_t m, qpms_l_t n, csph_t kdlj, qpms_bessel_t btyp, qpms_normalisation_t norm)

cdef extern from "indexing.h":
    qpms_y_t qpms_lMax2nelem(qpms_l_t lMax)
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
        void *stateData
    void PGen_destroy(PGen *g)
    
    int l2d_reciprocalBasis2pi(cart2_t b1, cart2_t b2, cart2_t *rb1, cart2_t *rb2);
    double l2d_unitcell_area(cart2_t b1, cart2_t b2)
    void l2d_reduceBasis(cart2_t in1, cart2_t in2, cart2_t *out1, cart2_t *out2)

    const double BASIS_RTOL
    size_t qpms_emptylattice2_modes_maxfreq(double **target_freqs, cart2_t b1_rec, cart2_t b2_rec,
            double rtol, cart2_t k, double wave_speed, double maxomega)
    void qpms_emptylattice2_modes_nearest(double *target_freqs, cart2_t b1_rec, cart2_t b2_rec,
            double rtol, cart2_t k, double wave_speed, double omega)

    # now the individual PGen implementations:
    # FIXME Is bint always guaranteed to be equivalent to _Bool? (I dont't think so.)
    PGen PGen_xyWeb_new(cart2_t b1, cart2_t b2, double rtol, cart2_t offset,
            double minR, bint inc_minR, double maxR, bint inc_maxR)
    ctypedef enum PGen_1D_incrementDirection:
        PGEN_1D_INC_FROM_ORIGIN
        PGEN_1D_INC_TOWARDS_ORIGIN
    PGen PGen_1D_new_minMaxR(double period, double offset, double minR, bint inc_minR,
            double maxR, bint inc_maxR, PGen_1D_incrementDirection incdir)
    int qpms_reduce_lattice_basis(double *b, size_t bsize, size_t ndim, double delta)
    ctypedef enum LatticeDimensionality:
        LAT1D
        LAT2D
        LAT3D
        SPACE1D
        SPACE2D
        SPACE3D
        LAT_1D_IN_3D
        LAT_2D_IN_3D
        LAT_3D_IN_3D
        LAT_ZONLY
        LAT_XYONLY
        LAT_1D_IN_3D_ZONLY
        LAT_2D_IN_3D_XYONLY


cdef extern from "vectors.h":
    cart2_t cart2_substract(cart2_t a, cart2_t b)
    cart2_t cart2_scale(const double c, cart2_t b)
    double cart2norm(cart2_t a)
    const cart2_t CART2_ZERO


cdef extern from "quaternions.h":
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
    qpms_quat_t qpms_quat_from_rotvector(cart3_t v)

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
    size_t qpms_zero_roundoff_clean(double *arr, size_t nmemb, double atol)
    size_t qpms_czero_roundoff_clean(cdouble *arr, size_t nmemb, double atol)

#cdef extern from "numpy/arrayobject.h":
#    cdef enum NPY_TYPES:
#        NPY_DOUBLE
#        NPY_CDOUBLE # complex double
#        NPY_LONG # int
#    ctypedef int npy_intp


cdef extern from "translations.h":
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
            int m, int n, int mu, int nu, cdouble kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J) nogil
    cdouble qpms_trans_calculator_get_B_ext(const qpms_trans_calculator* c,
            int m, int n, int mu, int nu, cdouble kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J) nogil
    int qpms_trans_calculator_get_AB_p_ext(const qpms_trans_calculator* c,
            cdouble *Adest, cdouble *Bdest,
            int m, int n, int mu, int nu, cdouble kdlj_r, double kdlj_th, double kdlj_phi,
            int r_ge_d, int J) nogil
    int qpms_trans_calculator_get_AB_arrays_ext(const qpms_trans_calculator *c,
            cdouble *Adest, cdouble *Bdest,
            size_t deststride, size_t srcstride,
            cdouble kdlj_r, double kdlj_theta, double kdlj_phi,
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

    int qpms_trans_calculator_get_trans_array(const qpms_trans_calculator *c,
                cdouble *target,
                const qpms_vswf_set_spec_t *destspec, size_t deststride,
                const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
                csph_t kdlj, bint r_ge_d, qpms_bessel_t J);

    int qpms_trans_calculator_get_trans_array_lc3p(
                const qpms_trans_calculator *c,
                cdouble *target,
                const qpms_vswf_set_spec_t *destspec, size_t deststride,
                const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
                cdouble k, cart3_t destpos, cart3_t srcpos,
                qpms_bessel_t J
                );

cdef extern from "qpms_specfunc.h":
    struct qpms_pitau_t:
        qpms_l_t lMax
        double *leg
        double *pi
        double *tau
    qpms_pitau_t qpms_pitau_get(double theta, qpms_l_t lMax, double csphase)
    qpms_errno_t qpms_pitau_fill(double *leg, double *pi, double *tau, double theta, qpms_l_t lMax, double csphase)
    void qpms_pitau_free(qpms_pitau_t pitau)


cdef extern from "gsl/gsl_interp.h":
    struct gsl_interp_type:
        pass
    const gsl_interp_type *gsl_interp_linear
    const gsl_interp_type *gsl_interp_cspline
    # ^^^ These are probably the only relevant ones.

cdef extern from "materials.h":
    struct qpms_epsmu_generator_t:
        qpms_epsmu_t (*function)(cdouble omega, const void *params)
        const void *params
    qpms_epsmu_t qpms_epsmu_const_g(cdouble omega, const void *params)
    qpms_epsmu_t qpms_permittivity_interpolator_epsmu_g(cdouble omega, const void *epsmu)
    qpms_epsmu_t qpms_lorentzdrude_epsmu_g(cdouble omega, const void *ldparams)

    struct qpms_permittivity_interpolator_t:
        pass
    qpms_permittivity_interpolator_t *qpms_permittivity_interpolator_create(const size_t incount,
            cdouble *wavelength_m, cdouble *n, cdouble *k, const gsl_interp_type *iptype)
    qpms_permittivity_interpolator_t *qpms_permittivity_interpolator_from_yml(const char *path,
            const gsl_interp_type *iptype)
    cdouble qpms_permittivity_interpolator_eps_at_omega(const qpms_permittivity_interpolator_t *interp, double omega_SI)
    double qpms_permittivity_interpolator_omega_max(const qpms_permittivity_interpolator_t *interp)
    double qpms_permittivity_interpolator_omega_min(const qpms_permittivity_interpolator_t *interp)
    void qpms_permittivity_interpolator_free(qpms_permittivity_interpolator_t *interp)
    struct qpms_ldparams_triple_t:
        double f
        double omega
        double gamma
    struct qpms_ldparams_t:
        cdouble eps_inf
        double omega_p
        size_t n
        qpms_ldparams_triple_t data[0]
    cdouble qpms_lorentzdrude_eps(cdouble, const qpms_ldparams_t *)

cdef extern from "tmatrices.h":
    bint qpms_load_scuff_tmatrix_crash_on_failure
    struct qpms_tmatrix_generator_t:
        qpms_errno_t (*function)(qpms_tmatrix_t *t, cdouble omega, const void *params)
        const void *params
    qpms_errno_t qpms_tmatrix_generator_axialsym(qpms_tmatrix_t *t, cdouble omega, const void *params)
    qpms_errno_t qpms_tmatrix_generator_interpolator(qpms_tmatrix_t *t, cdouble omega, const void *params)
    qpms_errno_t qpms_tmatrix_generator_sphere(qpms_tmatrix_t *t, cdouble omega, const void *params)
    qpms_errno_t qpms_tmatrix_generator_constant(qpms_tmatrix_t *t, cdouble omega, const void *params)
    struct qpms_tmatrix_generator_sphere_param_t:
        qpms_epsmu_generator_t outside
        qpms_epsmu_generator_t inside
        double radius
    struct qpms_arc_function_retval_t:
        double r
        double beta
    struct qpms_arc_function_t:
        qpms_arc_function_retval_t (*function)(double theta, const void *params)
        const void *params
    struct qpms_tmatrix_generator_axialsym_param_t:
        qpms_epsmu_generator_t outside
        qpms_epsmu_generator_t inside
        qpms_arc_function_t shape
        qpms_l_t lMax_extend
    struct qpms_arc_cylinder_params_t:
        double R
        double h
    qpms_arc_function_retval_t qpms_arc_cylinder(double theta, const void *params)
    qpms_arc_function_retval_t qpms_arc_sphere(double theta, const void *params)
    struct qpms_tmatrix_interpolator_t:
        const qpms_vswf_set_spec_t *bspec
    void qpms_tmatrix_interpolator_free(qpms_tmatrix_interpolator_t *interp)
    qpms_tmatrix_t *qpms_tmatrix_interpolator_eval(const qpms_tmatrix_interpolator_t *interp, double freq)
    qpms_tmatrix_interpolator_t *qpms_tmatrix_interpolator_create(size_t n, double *freqs, 
            const qpms_tmatrix_t *tmatrices_array, const gsl_interp_type *iptype)
    void qpms_tmatrix_free(qpms_tmatrix_t *tmatrix)
    qpms_tmatrix_isclose(const qpms_tmatrix_t *A, const qpms_tmatrix_t *B,
                const double rtol, const double atol)
    qpms_errno_t qpms_symmetrise_tmdata_irot3arr(
            cdouble *tmdata, const size_t tmcount,
            const qpms_vswf_set_spec_t *bspec,
            size_t n_symops,
            const qpms_irot3_t *symops
            )
    qpms_errno_t qpms_symmetrise_tmdata_finite_group(
            cdouble *tmdata, const size_t tmcount,
            const qpms_vswf_set_spec_t *bspec,
            const qpms_finite_group_t *pointgroup
            )
    qpms_tmatrix_t *qpms_tmatrix_symmetrise_irot3arr_inplace(
            qpms_tmatrix_t *T,
            size_t n_symops,
            const qpms_irot3_t *symops
            )
    qpms_tmatrix_t *qpms_tmatrix_symmetrise_finite_group_inplace(
            qpms_tmatrix_t *T,
            const qpms_finite_group_t *pointgroup
            )
    qpms_errno_t qpms_load_scuff_tmatrix(const char *path, const qpms_vswf_set_spec_t *bspec,
            size_t *n, double **freqs, double **freqs_su, qpms_tmatrix_t **tmatrices_array,
            cdouble **tmdata)
    cdouble *qpms_mie_coefficients_reflection(cdouble *target, const qpms_vswf_set_spec_t *bspec,
            double a, cdouble k_i, cdouble k_e, cdouble mu_i, cdouble mu_e, qpms_bessel_t J_ext, qpms_bessel_t J_scat)
    qpms_tmatrix_t *qpms_tmatrix_spherical(const qpms_vswf_set_spec_t *bspec, double a, 
            cdouble k_i, cdouble k_e, cdouble mu_i, cdouble mu_e)
    qpms_errno_t qpms_tmatrix_spherical_fill(qpms_tmatrix_t *t, double a, 
            cdouble k_i, cdouble k_e, cdouble mu_i, cdouble mu_e)
    qpms_tmatrix_t *qpms_tmatrix_spherical(const qpms_vswf_set_spec_t *bspec,
            double a, cdouble k_i, cdouble k_e, cdouble mu_i, cdouble mu_e)
    cdouble qpms_drude_epsilon(cdouble eps_inf, cdouble omega_p, cdouble gamma_p, cdouble omega)
    qpms_errno_t qpms_tmatrix_spherical_mu0_fill(qpms_tmatrix_t *t, double a, double omega,
            cdouble epsilon_fg, cdouble epsilon_bg)
    qpms_tmatrix_t *qpms_tmatrix_spherical_mu0(const qpms_vswf_set_spec_t *bspec, double a,
            double omega, cdouble epsilon_fg, cdouble epsilon_bg)
    qpms_errno_t qpms_tmatrix_generator_axialsym_RQ_transposed_fill(cdouble *target, cdouble omega,
            const qpms_tmatrix_generator_axialsym_param_t *param, qpms_normalisation_t norm, qpms_bessel_t J)
    struct qpms_tmatrix_function_t:
        const qpms_vswf_set_spec_t *spec
        const qpms_tmatrix_generator_t *gen
    ctypedef enum qpms_tmatrix_operation_kind_t:
        QPMS_TMATRIX_OPERATION_NOOP
        QPMS_TMATRIX_OPERATION_LRMATRIX
        QPMS_TMATRIX_OPERATION_IROT3
        QPMS_TMATRIX_OPERATION_IROT3ARR
        QPMS_TMATRIX_OPERATION_COMPOSE_SUM
        QPMS_TMATRIX_OPERATION_COMPOSE_CHAIN
        QPMS_TMATRIX_OPERATION_SCMULZ
        QPMS_TMATRIX_OPERATION_FINITE_GROUP_SYMMETRISE
    
    struct qpms_tmatrix_operation_t:
        qpms_tmatrix_operation_kind_t typ
        pass # TODO add the op union later if needed
    const qpms_tmatrix_operation_t qpms_tmatrix_operation_noop
    void qpms_tmatrix_operation_clear(qpms_tmatrix_operation_t *)

cdef extern from "pointgroups.h":
    bint qpms_pg_is_finite_axial(qpms_pointgroup_class cls)
    double qpms_pg_quat_cmp_atol
    int qpms_pg_irot3_cmp(const qpms_irot3_t *, const qpms_irot3_t *);
    int qpms_pg_irot3_cmp_v(const void *, const void *);
    int qpms_pg_irot3_approx_cmp(const qpms_irot3_t *a, const qpms_irot3_t *b, double atol)
    int qpms_pg_irot3_approx_cmp_v(const void *a, const void *b)

    qpms_gmi_t qpms_pg_order(qpms_pointgroup_class cls, qpms_gmi_t n)
    qpms_irot3_t *qpms_pg_canonical_elems( qpms_irot3_t *target, qpms_pointgroup_class cls, qpms_gmi_t n)
    qpms_gmi_t qpms_pg_genset_size(qpms_pointgroup_class cls, qpms_gmi_t n)
    qpms_gmi_t qpms_pg_genset(qpms_pointgroup_class cls, qpms_gmi_t n, qpms_irot3_t *gen)
    qpms_irot3_t *qpms_pg_elems(qpms_irot3_t *target, qpms_pointgroup_t g)
    bint qpms_pg_is_subgroup(qpms_pointgroup_t a, qpms_pointgroup_t b);

cdef extern from "scatsystem.h":
    void qpms_scatsystem_set_nthreads(long n)
    struct qpms_particle_t:
        cart3_t pos
        const qpms_tmatrix_function_t *tmg
        qpms_tmatrix_operation_t op
    struct qpms_particle_tid_t:
        cart3_t pos
        qpms_ss_tmi_t tmatrix_id
    struct qpms_ss_derived_tmatrix_t:
        qpms_ss_tmgi_t tmgi
        qpms_tmatrix_operation_t op
    struct qpms_scatsys_periodic_info_t:
        cart3_t lattice_basis[3]
        double unitcell_volume
        #etc.
    struct qpms_scatsys_t:
        int lattice_dimension
        qpms_epsmu_generator_t medium
        qpms_tmatrix_function_t *tmg
        qpms_ss_tmgi_t tmg_count
        qpms_ss_derived_tmatrix_t *tm
        qpms_ss_tmi_t tm_count
        qpms_particle_tid_t *p
        qpms_ss_pi_t p_count
        # We shouldn't need more to construct a symmetric scatsystem ^^^
        size_t fecv_size
        size_t *saecv_sizes
        const qpms_finite_group_t *sym
        qpms_scatsys_periodic_info_t per

        # per[] and other stuff not currently needed in cython
    void qpms_scatsys_free(qpms_scatsys_t *s)
    qpms_errno_t qpms_scatsys_dump(qpms_scatsys_t *ss, char *path) #NI
    qpms_scatsys_t *qpms_scatsys_load(char *path) #NI
    struct qpms_scatsys_at_omega_t:
        const qpms_scatsys_t *ss
        qpms_tmatrix_t **tm,
        cdouble omega
        qpms_epsmu_t medium
        cdouble wavenumber
    qpms_scatsys_at_omega_t *qpms_scatsys_apply_symmetry(const qpms_scatsys_t *orig, const qpms_finite_group_t *sym,
            cdouble omega, const qpms_tolerance_spec_t *tol)
    qpms_scatsys_at_omega_t *qpms_scatsys_at_omega(const qpms_scatsys_t *ss, cdouble omega)
    void qpms_scatsys_at_omega_free(qpms_scatsys_at_omega_t *ssw)
    cdouble *qpms_scatsys_irrep_pack_matrix(cdouble *target_packed,
            const cdouble *orig_full, const qpms_scatsys_t *ss, qpms_iri_t iri)
    cdouble *qpms_scatsys_irrep_unpack_matrix(cdouble *target_full, 
            const cdouble *orig_packed, const qpms_scatsys_t *ss, qpms_iri_t iri, bint add)
    cdouble *qpms_scatsys_irrep_pack_vector(cdouble *target_packed,
            const cdouble *orig_full, const qpms_scatsys_t *ss, qpms_iri_t iri)
    cdouble *qpms_scatsys_irrep_unpack_vector(cdouble *target_full,
            const cdouble *orig_packed, const qpms_scatsys_t *ss, qpms_iri_t iri, bint add)
    cdouble *qpms_scatsysw_build_modeproblem_matrix_full(cdouble *target,
            const qpms_scatsys_at_omega_t *ssw)
    cdouble *qpms_scatsys_build_translation_matrix_full(cdouble *target,
            const qpms_scatsys_t *ss, cdouble k)
    cdouble *qpms_scatsys_build_translation_matrix_e_full(cdouble *target,
            const qpms_scatsys_t *ss, cdouble k, qpms_bessel_t J)
    cdouble *qpms_scatsysw_build_modeproblem_matrix_irrep_packed(cdouble *target,
            const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri) nogil
    cdouble *qpms_scatsys_build_translation_matrix_e_irrep_packed(cdouble *target,
            const qpms_scatsys_t *ss, qpms_iri_t iri, cdouble k, qpms_bessel_t J) nogil
    cdouble *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_orbitorderR(
            cdouble *target, const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri) nogil
    cdouble *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_serial(
            cdouble *target, const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri) nogil
    cdouble *qpms_scatsys_incident_field_vector_full(cdouble *target_full,
            const qpms_scatsys_t *ss, qpms_incfield_t field_at_point, 
            const void *args, bint add)
    cdouble *qpms_scatsysw_apply_Tmatrices_full(cdouble *target_full, const cdouble *inc_full,
            const qpms_scatsys_at_omega_t *ssw)
    struct qpms_ss_LU:
        const qpms_scatsys_at_omega_t *ssw
        const qpms_scatsys_at_omega_k_t *sswk
        bint full
        qpms_iri_t iri
        cdouble *a
        int *ipiv
    void qpms_ss_LU_free(qpms_ss_LU lu)
    qpms_ss_LU qpms_scatsysw_build_modeproblem_matrix_full_LU(cdouble *target,
            int *target_piv, const qpms_scatsys_at_omega_t *ssw)
    qpms_ss_LU qpms_scatsysw_build_modeproblem_matrix_irrep_packed_LU(cdouble *target,
            int *target_piv, const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri)
    qpms_ss_LU qpms_scatsysw_modeproblem_matrix_full_factorise(cdouble *modeproblem_matrix_full,
            int *target_piv, const qpms_scatsys_at_omega_t *ssw)
    qpms_ss_LU qpms_scatsysw_modeproblem_matrix_irrep_packed_factorise(cdouble *modeproblem_matrix_full,
            int *target_piv, const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri)
    cdouble *qpms_scatsys_scatter_solve(cdouble *target_f, const cdouble *a_inc, qpms_ss_LU ludata)
    const qpms_vswf_set_spec_t *qpms_ss_bspec_tmi(const qpms_scatsys_t *ss, qpms_ss_tmi_t tmi)
    const qpms_vswf_set_spec_t *qpms_ss_bspec_pi(const qpms_scatsys_t *ss, qpms_ss_pi_t pi)
    beyn_result_t *qpms_scatsys_finite_find_eigenmodes(const qpms_scatsys_t *ss, qpms_iri_t iri,
            cdouble omega_centre, double omega_rr, double omega_ri, size_t contour_npoints, 
            double rank_tol, size_t rank_sel_min, double res_tol)
    # periodic-related funs
    struct qpms_scatsys_at_omega_k_t:
        const qpms_scatsys_at_omega_t *ssw
        double k[3]
    cdouble *qpms_scatsyswk_build_modeproblem_motrix_full(cdouble *target, const qpms_scatsys_at_omega_k_t *sswk)
    cdouble *qpms_scatsys_periodic_build_translation_matrix_full(cdouble *target, const qpms_scatsys_t *ss, cdouble wavenumber, const cart3_t *wavevector)
    qpms_ss_LU qpms_scatsyswk_build_modeproblem_matrix_full_LU(cdouble *target, int *target_piv, const qpms_scatsys_at_omega_k_t *sswk)
    beyn_result_t *qpms_scatsys_periodic_find_eigenmodes(const qpms_scatsys_t *ss, const double *k,
            cdouble omega_centre, double omega_rr, double omega_ri, size_t contour_npoints, 
            double rank_tol, size_t rank_sel_min, double res_tol)
    const qpms_vswf_set_spec_t *qpms_ss_bspec_pi(const qpms_scatsys_t *ss, qpms_ss_pi_t pi) 
    ccart3_t qpms_scatsys_scattered_E(const qpms_scatsys_t *ss, cdouble wavenumber,
            const cdouble *f_excitation_vector_full, cart3_t where)
    ccart3_t qpms_scatsysw_scattered_E(const qpms_scatsys_at_omega_t *ssw,
            const cdouble *f_excitation_vector_full, cart3_t where)
    ccart3_t qpms_scatsys_scattered_E__alt(const qpms_scatsys_t *ss, cdouble wavenumber,
            const cdouble *f_excitation_vector_full, cart3_t where)
    ccart3_t qpms_scatsysw_scattered_E__alt(const qpms_scatsys_at_omega_t *ssw,
            const cdouble *f_excitation_vector_full, cart3_t where)

cdef extern from "ewald.h":
    struct qpms_csf_result:
        cdouble val
        double err

    ctypedef enum qpms_ewald_part:
        QPMS_EWALD_LONG_RANGE
        QPMS_EWALD_SHORT_RANGE
        QPMS_EWALD_FULL
        QPMS_EWALD_0TERM

    struct qpms_ewald3_constants_t:
        qpms_l_t lMax
        qpms_y_t nelem_sc
    qpms_ewald3_constants_t *qpms_ewald3_constants_init(qpms_l_t lMax, int csphase)
    void qpms_ewald3_constants_free(qpms_ewald3_constants_t *c)
    
    cdouble lilgamma(double t)
    cdouble clilgamma(cdouble z)
    int cx_gamma_inc_series_e(double a, cdouble x, qpms_csf_result *result)
    int cx_gamma_inc_CF_e(double a, cdouble x, qpms_csf_result *result)
    int complex_gamma_inc_e(double a, cdouble x, int m, qpms_csf_result *result)

    int ewald3_sigma0(cdouble *target, double *err, const qpms_ewald3_constants_t *c, double eta, cdouble wavenumber)
    int ewald3_sigma_long(cdouble *target_sigmalr_y, double *target_sigmalr_y_err, const qpms_ewald3_constants_t *c, 
            double eta, cdouble wavenumber, double unitcell_volume, LatticeDimensionality latdim, PGen *pgen_K,
            bint pgen_generates_shifted_points, cart3_t k, cart3_t particle_shift)
    int ewald3_sigma_short(cdouble *target_sigmasr_y, double *target_sigmasr_y_err, const qpms_ewald3_constants_t *c, 
            double eta, cdouble wavenumber, LatticeDimensionality latdim, PGen *pgen_R, 
            bint pgen_generates_shifted_points, cart3_t k, cart3_t particle_shift)


cdef extern from "gsl/gsl_complex.h":
    ctypedef struct gsl_complex:
        double dat[2]

cdef extern from "gsl/gsl_matrix.h":
    ctypedef struct gsl_matrix_complex:
        pass
    ctypedef struct gsl_vector:
        pass
    ctypedef struct gsl_vector_complex:
        pass

cdef extern from "beyn.h":
    ctypedef struct beyn_contour_t:
        bint (*inside_test)(beyn_contour_t *, cdouble z)
        pass
    ctypedef struct beyn_result_gsl_t:
        pass
    ctypedef struct beyn_result_t:
        size_t neig
        size_t vlen
        cdouble *eigval
        cdouble *eigval_err
        double *residuals
        cdouble *eigvec
        double *ranktest_SV
        beyn_result_gsl_t *gsl
    ctypedef enum beyn_contour_halfellipse_orientation:
        BEYN_CONTOUR_HALFELLIPSE_RE_PLUS
        BEYN_CONTOUR_HALFELLIPSE_IM_PLUS
        BEYN_CONTOUR_HALFELLIPSE_RE_MINUS
        BEYN_CONTOUR_HALFELLIPSE_IM_MINUS

    ctypedef int (*beyn_function_M_gsl_t)(gsl_matrix_complex *target_M, cdouble z, void *params)
    ctypedef int (*beyn_function_M_inv_Vhat_gsl_t)(gsl_matrix_complex *target, const gsl_matrix_complex *Vhat, cdouble z, void *params)
    ctypedef int (*beyn_function_M_t)(cdouble *target_M, size_t m, cdouble z, void *params)
    ctypedef int (*beyn_function_M_inv_Vhat_t)(cdouble *target, size_t m, size_t l, const cdouble *Vhat, cdouble z, void *params)

    void beyn_result_gsl_free(beyn_result_gsl_t *result)
    void beyn_result_free(beyn_result_t *result)

    beyn_result_gsl_t *beyn_solve_gsl(size_t m, size_t l, beyn_function_M_gsl_t M,
            beyn_function_M_inv_Vhat_gsl_t M_inv_Vhat, void *params, const beyn_contour_t *contour,
            double rank_tol, size_t rank_min_sel, double res_tol)

    beyn_result_t *beyn_solve(size_t m, size_t l, beyn_function_M_t M,
            beyn_function_M_inv_Vhat_t M_inv_Vhat, void *params, const beyn_contour_t *contour,
            double rank_tol, size_t rank_min_sel, double res_tol)

    beyn_contour_t *beyn_contour_ellipse(cdouble centre, double halfax_re, double halfax_im, size_t npoints)
    beyn_contour_t *beyn_contour_halfellipse(cdouble centre, double halfax_re, double halfax_im, size_t npoints,
            beyn_contour_halfellipse_orientation ori)
    beyn_contour_t *beyn_contour_kidney(cdouble centre, double halfax_re, double halfax_im, size_t npoints,
            double rounding, beyn_contour_halfellipse_orientation ori)


    cdouble gsl_comlpex_tostd(gsl_complex z)
    gsl_complex gsl_complex_fromstd(cdouble z)

