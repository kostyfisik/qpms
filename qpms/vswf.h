/*! \file vswf.h
 * \brief Vector spherical wavefunctions.
 *
 * N.B. for the Legendre polynomial norm definitions, see
 * <a href="https://www.gnu.org/software/gsl/doc/html/specfunc.html#associated-legendre-polynomials-and-spherical-harmonics">the corresponding section of GSL docs</a>
 * or <a href="http://git.savannah.gnu.org/cgit/gsl.git/tree/specfunc/legendre_source.c">gsl/specfunc/legendre_source.c</a>.
 */
#ifndef QPMS_VSWF_H
#define QPMS_VSWF_H
#include <unistd.h> // ssize_t 
#include "qpms_types.h"
#include <gsl/gsl_sf_legendre.h>

// -------------- Typedefs (function prototypes) for qpms_vswf_spec_t ----------------------

/// Calculates the (regular VSWF) expansion coefficients of an external incident field.
typedef qpms_errno_t (*qpms_incfield_t)(
        /// Target non-NULL array of the regular VSWF expansion coefficients of length bspec->n.
        complex double *target,
        const qpms_vswf_set_spec_t *bspec,
        const cart3_t evalpoint, ///< Point at which the VSWF expansion is made.
        const void *args, ///< Pointer to additional function-specific arguments.
        bool add ///< If true, add to target; rewrite target if false.
);

// ---------------Methods for qpms_vswf_spec_t-----------------------
//
/// Creates a qpms_vswf_set_spec_t structure with an empty list of wave indices.
qpms_vswf_set_spec_t *qpms_vswf_set_spec_init(void);
/// Appends a VSWF index to a \ref qpms_vswf_set_spec_t, also updating metadata.
qpms_errno_t qpms_vswf_set_spec_append(qpms_vswf_set_spec_t *self, qpms_uvswfi_t u);
/// Destroys a \ref qpms_vswf_set_spec_t.
void qpms_vswf_set_spec_free(qpms_vswf_set_spec_t *);
/// Compares two vswf basis specs.
/**
 * Checks whether ilist is the same and of the same length.
 * If yes, returns true, else returns false.
 */
bool qpms_vswf_set_spec_isidentical(const qpms_vswf_set_spec_t *a,
		const qpms_vswf_set_spec_t *b);
/// Copies an instance of qpms_vswf_set_spec_t
qpms_vswf_set_spec_t *qpms_vswf_set_spec_copy(const qpms_vswf_set_spec_t *orig);
/// Creates an instance of qpms_vswf_set_spec_t in the 'traditional' layout.
qpms_vswf_set_spec_t *qpms_vswf_set_spec_from_lMax(qpms_l_t lMax,
		qpms_normalisation_t norm);

/// Finds the position of a given index in the bspec's ilist.
/** If not found, returns -1. */
// TODO more consistency in types (here size_t vs. ptrdiff_t).
static inline ssize_t qpms_vswf_set_spec_find_uvswfi(const qpms_vswf_set_spec_t *bspec,
		const qpms_uvswfi_t index) {
	for(size_t i = 0; i < bspec->n; ++i) 
		if (bspec->ilist[i] == index)
			return i;
	return -1;
}

/// Evaluates a set of VSWF basis functions at a given point.
/** The list of basis wave indices is specified in \a setspec; 
 *  \a setspec->norm must be set as well.
 */
qpms_errno_t qpms_uvswf_fill(
		csphvec_t *const target, ///< Target array of size at least setspec->n.
		const qpms_vswf_set_spec_t *setspec,
		csph_t kr, ///< Evaluation point.
	       	qpms_bessel_t btyp);

/// Evaluates field specified by SVWF coefficients at a given point.
/** SVWF coefficients in \a coeffs must be ordered according to \a setspec->ilist.
 */
csphvec_t qpms_eval_uvswf(const qpms_vswf_set_spec_t *setspec,
		const complex double *coeffs, ///< SVWF coefficient vector of size setspec->n.
		csph_t kr, ///< Evaluation point.
		qpms_bessel_t btyp);


// --- qpms_incfield_t instances and their arguments

/// Parameter structure for qpms_incfield_planewave()
typedef struct qpms_incfield_planewane_params_t {
	bool use_cartesian; ///< If true, wave direction k and amplitude E are specified in cartesian coordinates (via k.cart, E.cart). If false, k is specified in spherical coordinates and E are specified in the corresponding geographical coordinates (via k.sph, E.sph).
	union {
		ccart3_t cart;
		csph_t sph;
	} k; ///< Wave vector.
	union {
		ccart3_t cart;
		csphvec_t sph;
	} E; ///< Electric field amplitude at origin.
} qpms_incfield_planewave_params_t;
	

/// Calculates the (regular VSWF) expansion coefficients of a plane wave.
/** The wave amplitude and wave vector is defined by struct qpms_incfield_planewave_params_t.
 *
 *  If the wave vector and amplitude are not orthogonal (i.e. the plane wave is not
 *  fully transversal) and bspec->lMaxL is non-negative, 
 *  the corresponding longitudinal components are calculated as well.
 *  
 *  For complex k vectors, the implementation is not completely correct right now.
 *  Locally, it corresponds to decomposition of a plane wave with a real \a k
 *  (using the real part of the \a k supplied), just the whole decomposition
 *  is modulated by the origin-dependent factor \f$ \vect E e^{i \vect k \cdot \vect r} \f$
 *  with \f$ \vect k \f$ complex.
 */
qpms_errno_t qpms_incfield_planewave(
        /// Target non-NULL array of the regular VSWF expansion coefficients of length bspec->n.
        complex double *target,
        const qpms_vswf_set_spec_t *bspec,
        const cart3_t evalpoint, ///< Point at which the VSWF expansion is made.
        const void *args, ///< Pointer to additional function-specific arguments (converted to (const qpms_incfield_planewave_params_t *)).
        bool add ///< If true, add to target; rewrite target if false.
);

// -----------------------------------------------------------------------

/// Electric wave N.
csphvec_t qpms_vswf_single_el(int m, int n, sph_t kdlj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);
/// Magnetic wave M.
csphvec_t qpms_vswf_single_mg(int m, int n, sph_t kdlj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);

/// Set of electric and magnetic VSWF values in spherical coordinate basis.
/** This is supposed to contain all the waves up to $l = lMax$.
 * 
 *  For a completely custom set of waves, use \ref qpms_uvswfset_sph_t instead.
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


/// Evaluate the zeroth-degree longitudinal VSWF \f$ \mathbf{L}_0^0 \f$.
/**
 * Any `norm` is being ignored right now.
 */
csphvec_t qpms_vswf_L00(
		csph_t kdrj, //< VSWF evaluation point.
		qpms_bessel_t btyp,
		qpms_normalisation_t norm //< Ignored!
		);

/// Evaluate VSWFs at a given point from \a l = 1 up to a given degree \a lMax.
/**
 * The target arrays \a resultL, \a resultM, \a resultN have to be large enough to contain
 * \a lMax * (\a lMax + 2) elements. If NULL is passed instead, the corresponding SVWF type
 * is not evaluated.
 * 
 * Does not evaluate the zeroth-order wave \f$ \mathbf{L}_0^0 \f$. 
 * If you need that, use qpms_vswf_L00().
 */
qpms_errno_t qpms_vswf_fill(
	csphvec_t *resultL, //< Target array for longitudinal VSWFs.
       	csphvec_t *resultM, //< Target array for magnetic VSWFs.
       	csphvec_t *resultN, //< Target array for electric VSWFs.
	qpms_l_t lMax, //< Maximum multipole degree to be calculated.
	sph_t kdrj, //< VSWF evaluation point.
	qpms_bessel_t btyp, qpms_normalisation_t norm);

// Should give the same results: for consistency checks
qpms_errno_t qpms_vswf_fill_alternative(csphvec_t *resultL, csphvec_t *resultM, csphvec_t *resultN, qpms_l_t lMax, sph_t kdrj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);

/// Evaluate VSWFs at a given point from \a l = 1 up to a given degree \a lMax (complex \a kr version).
/**
 * The target arrays \a resultL, \a resultM, \a resultN have to be large enough to contain
 * \a lMax * (\a lMax + 2) elements. If NULL is passed instead, the corresponding SVWF type
 * is not evaluated.
 * 
 * Does not evaluate the zeroth-order wave \f$ \mathbf{L}_0^0 \f$. 
 * If you need that, use qpms_vswf_L00().
 */
qpms_errno_t qpms_vswf_fill_csph(
	csphvec_t *resultL, //< Target array for longitudinal VSWFs.
       	csphvec_t *resultM, //< Target array for magnetic VSWFs.
       	csphvec_t *resultN, //< Target array for electric VSWFs.
	qpms_l_t lMax, //< Maximum multipole degree to be calculated.
	csph_t kdrj, //< VSWF evaluation point.
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

csphvec_t qpms_eval_vswf_csph(csph_t where,
		complex double *longcoeffs, complex double *mgcoeffs, complex double *elcoeffs,
		qpms_l_t lMax, qpms_bessel_t btyp, qpms_normalisation_t norm);

qpms_vswfset_sph_t *qpms_vswfset_make(qpms_l_t lMax, sph_t kdlj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);//NI
void qpms_vswfset_sph_pfree(qpms_vswfset_sph_t *);//NI

#endif // QPMS_VSWF_H
