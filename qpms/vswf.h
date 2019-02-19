/*! \file vswf.h
 * \brief Vector spherical wavefunctions.
 *
 * N.B. for the Legendre polynomial norm definitions, see
 * <a href="https://www.gnu.org/software/gsl/doc/html/specfunc.html#associated-legendre-polynomials-and-spherical-harmonics">the corresponding section of GSL docs</a>
 * or <a href="http://git.savannah.gnu.org/cgit/gsl.git/tree/specfunc/legendre_source.c">gsl/specfunc/legendre_source.c</a>.
 */
#ifndef QPMS_VSWF_H
#define QPMS_VSWF_H
#include "qpms_types.h"
#include <gsl/gsl_sf_legendre.h>

/// Specifies a finite set of VSWFs.
/**
 * When for example not all the M and N -type waves up to a degree lMax
 * need to be computed, this will specify the subset. 
 *
 * A typical use case would be when only even/odd fields wrt. z-plane
 * mirror symmetry are considered.
 */
typedef struct qpms_vswf_set_spec_t {
	size_t n; ///< Actual number of VSWF indices included in ilist.
	qpms_uvswfi_t *ilist; ///< List of wave indices.
	qpms_l_t lMax; ///< Maximum degree of the waves specified in ilist.
	qpms_l_t lMax_M, ///< Maximum degree of the magnetic (M-type) waves.
		 lMax_N, ///< Maximum degree of the electric (N-type) waves.
		 lMax_L; ///< Maximum degree of the longitudinal (L-type) waves.
	size_t capacity; ///< Allocated capacity of ilist.
	qpms_normalisation_t norm; ///< Normalisation convention. To be set manually if needed.
} qpms_vswf_set_spec_t;

/// Creates a qpms_vswf_set_spec_t structure with an empty list of wave indices.
qpms_vswf_set_spec_t *qpms_vswf_set_spec_init();

/// Appends a VSWF index to a \ref qpms_vswf_set_spec_t, also updating metadata.
qpms_errno_t qpms_vswf_set_spec_append(qpms_vswf_set_spec_t *self, qpms_uvswfi_t u);

/// Destroys a \ref qpms_vswf_set_spec_t.
void qpms_vswf_set_spec_free(qpms_vswf_set_spec_t *);

/// Electric wave N.
csphvec_t qpms_vswf_single_el(int m, int n, sph_t kdlj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);
/// Magnetic wave M.
csphvec_t qpms_vswf_single_mg(int m, int n, sph_t kdlj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);

/// Set of electric and magnetic VSWF values in spherical coordinate basis.
/* This is supposed to contain all the waves up to $l = lMax$.
 * for a custom set of waves, use \ref qpms_uvswfset_sph_t instead.
 */

typedef struct qpms_vswfset_sph_t {
	//qpms_normalisation_t norm;
	qpms_l_t lMax;
	//qpms_y_t nelem;
	//sph_t kdlj
	csphvec_t *el, *mg;
} qpms_vswfset_sph_t;



qpms_errno_t qpms_legendre_deriv_y_get(double **result, double **result_deriv, double x, qpms_l_t lMax, 
		gsl_sf_legendre_t lnorm, double csphase); // free() result and result_deriv yourself!
qpms_errno_t qpms_legendre_deriv_y_fill(double *where, double *where_deriv, double x, 
		qpms_l_t lMax, gsl_sf_legendre_t lnorm, double csphase); 

/* some of the result targets may be NULL */
qpms_errno_t qpms_vswf_fill(csphvec_t *resultL, csphvec_t *resultM, csphvec_t *resultN, qpms_l_t lMax, sph_t kdrj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);
// Should give the same results: for consistency checks
qpms_errno_t qpms_vswf_fill_alternative(csphvec_t *resultL, csphvec_t *resultM, csphvec_t *resultN, qpms_l_t lMax, sph_t kdrj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);

qpms_errno_t qpms_vecspharm_fill(csphvec_t *const a1target, csphvec_t *const a2target, csphvec_t *const a3target,
		qpms_l_t lMax, sph_t dir, qpms_normalisation_t norm);
qpms_errno_t qpms_vecspharm_dual_fill(csphvec_t *const a1target, csphvec_t *const a2target, csphvec_t *const a3target,
		qpms_l_t lMax, sph_t dir, qpms_normalisation_t norm);

qpms_errno_t qpms_planewave2vswf_fill_cart(cart3_t wavedir, ccart3_t amplitude,
		complex double *targt_longcoeff, complex double *target_mgcoeff, complex double *target_elcoeff,
		qpms_l_t lMax, qpms_normalisation_t norm);
qpms_errno_t qpms_planewave2vswf_fill_sph(sph_t wavedir, csphvec_t amplitude,
		complex double *targt_longcoeff, complex double *target_mgcoeff, complex double *target_elcoeff,
		qpms_l_t lMax, qpms_normalisation_t norm);


csphvec_t qpms_eval_vswf(sph_t where,
		complex double *longcoeffs, complex double *mgcoeffs, complex double *elcoeffs,
		qpms_l_t lMax, qpms_bessel_t btyp, qpms_normalisation_t norm);


qpms_vswfset_sph_t *qpms_vswfset_make(qpms_l_t lMax, sph_t kdlj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);//NI
void qpms_vswfset_sph_pfree(qpms_vswfset_sph_t *);//NI

#endif // QPMS_VSWF_H
