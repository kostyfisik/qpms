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
	       	const gsl_matrix_complex *Vhat, complex double z, void *params); 

/// Complex plane integration contour structure.
typedef struct beyn_contour_t {
	size_t n; ///< Number of discretisation points.
	complex double centre; ///< TODO what is the exact purpose of this?
	complex double z_dz[][2]; ///< Pairs of contour points and derivatives in that points.
} beyn_contour_t;

/// Complex plane elliptic integration contour with axes parallel to the real, imaginary axes.
/** Free using free(). */
beyn_contour_t *beyn_contour_ellipse(complex double centre, double halfax_re, double halfax_im, size_t npoints);

/// Beyn algorithm result structure (GSL matrix/vector version).
typedef struct beyn_result_gsl_t {
	size_t neig; ///< Number of eigenvalues found (a bit redundant?).
	gsl_vector_complex *eigval;
	gsl_vector_complex *eigval_err;
	gsl_vector *residuals;
	gsl_matrix_complex *eigvec;
} beyn_result_gsl_t;

void beyn_result_gsl_free(beyn_result_gsl_t *result);

/// Beyn algorithm result structure (pure C array version).
typedef struct beyn_result_t {
	size_t neig; ///< Number of eigenvalues found.
	complex double *eigval;
	complex double *eigval_err;
	double *residuals;
	complex double *eigvec;
} beyn_result_t;

void beyn_result_free(beyn_result_t *result);

int beyn_solve_gsl(beyn_result_gsl_t **result,
		size_t m, ///< Dimension of the matrix \a M.
		size_t l, ///< Number of columns of the random matrix  \f$ \hat V \f$ (larger than the expected number of solutions).
		beyn_function_M_gsl_t M, ///< Function providing the matrix \f$ M(z) \f$.
		beyn_function_M_inv_Vhat_gsl_t M_inv_Vhat, ///< Fuction providing the matrix \f$ M^{-1}(z) \hat V \f$ (optional).
		void *params, ///< Parameter pointer passed to M() and M_inv_Vhat().
		const beyn_contour_t *contour, ///< Integration contour.
		double rank_tol, ///< (default: `1e-4`) TODO DOC.
		double res_tol ///< (default: `0.0`) TODO DOC.
	      );

int beyn_solve(beyn_result_t **result,
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
