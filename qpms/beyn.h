/** \file beyn.h
 * \brief Beyn's algorithm for nonlinear eigenvalue problems.
 */
#ifndef BEYN_H
#define BEYN_H

#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>

/// User-supplied function that provides the operator M(z) whose "roots" are to be found.
/** GSL matrix version */
typedef int (*beyn_function_M_gsl_t)(gsl_matrix_complex *target_M, complex double z, void *params);

/// (optional) User-supplied function that, given \f$ \hat V \f$, calculates \f$ M(z)^{-1} \hat V \f$.
/** GSL matrix version */
typedef int (*beyn_function_M_inv_Vhat_gsl_t)(gsl_matrix_complex *target_M_inv_Vhat,
	       	const gsl_matrix_complex *Vhat, complex double z, void *params); 

/// User-supplied function that provides the operator M(z) whose "roots" are to be found.
/** Pure C array version */
typedef int (*beyn_function_M_t)(complex double *target_M, size_t m, complex double z, void *params);

/// (optional) User-supplied function that, given \f$ \hat V \f$, calculates \f$ M(z)^{-1} \hat V \f$.
/** Pure C array version */
typedef int (*beyn_function_M_inv_Vhat_t)(complex double *target_M_inv_Vhat, size_t m, size_t l,
	       	const complex double *Vhat, complex double z, void *params); 

/// Complex plane integration contour structure.
typedef struct beyn_contour_t {
	size_t n; ///< Number of discretisation points.
	/// "Centre" of the contour.
	/** 
	 * This point is used in the rescaling of the \f$ A_1 \f$ matrix as in
	 * Beyn's Remark 3.2 (b) in order to improve the numerical stability.
	 * It does not have to be a centre in some strictly defined sense,
	 * but it should be "somewhere around" where the contour is.
	 */
	complex double centre;
	/// Function testing that a point \a z lies inside the contour (optional).
	_Bool (*inside_test)(struct beyn_contour_t *, complex double z);
	complex double z_dz[][2]; ///< Pairs of contour points and derivatives in that points.
} beyn_contour_t;

/// Complex plane elliptic integration contour with axes parallel to the real, imaginary axes.
/** Free using free(). */
beyn_contour_t *beyn_contour_ellipse(complex double centre, double halfax_re, double halfax_im, size_t npoints);


typedef enum {
	BEYN_CONTOUR_HALFELLIPSE_RE_PLUS = 3,
	BEYN_CONTOUR_HALFELLIPSE_RE_MINUS = 1,
	BEYN_CONTOUR_HALFELLIPSE_IM_PLUS = 0,
	BEYN_CONTOUR_HALFELLIPSE_IM_MINUS = 2,
} beyn_contour_halfellipse_orientation;


/// Complex plane "half-elliptic" integration contour with axes parallel to the real, imaginary axes.
/** Free using free(). */
beyn_contour_t *beyn_contour_halfellipse(complex double centre, double halfax_re, double halfax_im, size_t npoints,
		beyn_contour_halfellipse_orientation or);

/// Beyn algorithm result structure (GSL matrix/vector version).
typedef struct beyn_result_gsl_t {
	size_t neig; ///< Number of eigenvalues found (a bit redundant?).
	gsl_vector_complex *eigval;
	gsl_vector_complex *eigval_err;
	gsl_vector *residuals;
	gsl_matrix_complex *eigvec;
	gsl_vector *ranktest_SV;
} beyn_result_gsl_t;

void beyn_result_gsl_free(beyn_result_gsl_t *result);

/// Beyn algorithm result structure (pure C array version).
typedef struct beyn_result_t {
	size_t neig; ///< Number of eigenvalues found.
	size_t vlen; ///< Vector space dimension (also the leading dimension of eigvec).
	complex double *eigval;
	complex double *eigval_err;
	double *residuals;
	complex double *eigvec;
	double *ranktest_SV;

	/// Private, we wrap it around beyn_result_gsl_t for now.
	beyn_result_gsl_t *gsl;
	
} beyn_result_t;

/// Converts beyn_result_gsl_t from beyn_result_t.
/** After calling this, use beyn_result_free() on the returned pointer;
 *  do NOT run beyn_result_gsl_free() anymore.
 */
beyn_result_t *beyn_result_from_beyn_result_gsl(beyn_result_gsl_t *g);

void beyn_result_free(beyn_result_t *result);

beyn_result_gsl_t *beyn_solve_gsl(
		size_t m, ///< Dimension of the matrix \a M.
		size_t l, ///< Number of columns of the random matrix  \f$ \hat V \f$ (larger than the expected number of solutions).
		beyn_function_M_gsl_t M, ///< Function providing the matrix \f$ M(z) \f$.
		beyn_function_M_inv_Vhat_gsl_t M_inv_Vhat, ///< Fuction providing the matrix \f$ M^{-1}(z) \hat V \f$ (optional).
		void *params, ///< Parameter pointer passed to M() and M_inv_Vhat().
		const beyn_contour_t *contour, ///< Integration contour.
		double rank_tol, ///< (default: `1e-4`) TODO DOC.
		double res_tol ///< (default: `0.0`) TODO DOC.
	      );

beyn_result_t *beyn_solve(
		size_t m, ///< Dimension of the matrix \a M.
		size_t l, ///< Number of columns of the random matrix  \f$ \hat V \f$ (larger than the expected number of solutions).
		beyn_function_M_t M, ///< Function providing the matrix \f$ M(z) \f$.
		beyn_function_M_inv_Vhat_t M_inv_Vhat, ///< Fuction providing the matrix \f$ M^{-1}(z) \hat V \f$ (optional).
		void *params, ///< Parameter pointer passed to M() and M_inv_Vhat().
		const beyn_contour_t *contour, ///< Integration contour.
		double rank_tol, ///< (default: `1e-4`) TODO DOC.
		double res_tol ///< (default: `0.0`) TODO DOC.
	      );

static inline complex double gsl_complex_tostd(gsl_complex z) { return GSL_REAL(z) + I*GSL_IMAG(z); }
static inline gsl_complex gsl_complex_fromstd(complex double z) { return gsl_complex_rect(creal(z), cimag(z)); }

#endif // BEYN_H
