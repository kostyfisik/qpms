/*! \file ewald.h
 * \brief Lattice sums of spherical waves.
 *
 * Implementation of two-dimensional lattice sum in three dimensions
 * according to:
 * - [1] C.M. Linton, I. Thompson
 *     Journal of Computational Physics 228 (2009) 1815–1829
 * - [2] C.M.Linton
 *     SIAM Review Vol 52, No. 4, pp. 630–674
 *
 * N.B.!!! currently, the long-range parts are calculated 
 * not according to [1,(4.5)], but rather
 * according to the spherical-harmonic-normalisation-independent 
 * formulation in my notes notes/ewald.lyx.
 * Both parts of lattice sums are then calculated with 
 * the \f$ P_n^{|m|} e^{im\phi} \f$
 * (N.B. or \f$ P_n^{|m|} e^{imf} (-1)^m \f$ for negative m) 
 * substituted in place  of  \f$ Y_n^m \f$ 
 * (this is quite a weird normalisation especially 
 * for negative \f$ |m| \f$, but it is consistent
 * with the current implementation of the translation coefficients in
 * @ref translations.c;
 * in the long run, it might make more sense to replace it everywhere with normalised
 * Legendre polynomials).
 */

#ifndef EWALD_H
#define EWALD_H
#include <gsl/gsl_sf_result.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <math.h> // for inlined lilgamma
#include <complex.h>
#include "qpms_types.h"
#include "lattices.h"


typedef enum {
	QPMS_EWALD_LONG_RANGE = 1,
	QPMS_EWALD_SHORT_RANGE = 2,
	QPMS_EWALD_0TERM = 4,
	QPMS_EWALD_FULL = QPMS_EWALD_LONG_RANGE | QPMS_EWALD_SHORT_RANGE | QPMS_EWALD_0TERM,
} qpms_ewald_part;


/// Use this handler to ignore underflows of incomplete gamma.
gsl_error_handler_t IgnoreUnderflowsGSLErrorHandler;


/// Object holding the Ewald sum constant factors.
/**
 * Used internally by qpms_translation_calculator_t. 
 * Initialised by qpms_ewald3_constants_init() and freed by qpms_ewald3_constants_free().
 */
typedef struct qpms_ewald3_constants_t { 
	qpms_l_t lMax;
	qpms_y_t nelem_sc;
	/// The values of maximum \a j's in the long-range part summation, `[(l-|m|/2)]`.
	qpms_l_t *s1_jMaxes;
	/// The constant factors for the long range part of a 2D Ewald sum.
	complex double **s1_constfacs; // indices [y][j] where j is same as in [1, (4.5)]
	/* These are the actual numbers now: (in the EWALD32_CONSTANTS_AGNOSTIC version)
	 * for m + n EVEN:
	 *
	 * s1_constfacs[y(m,n)][j] = 
	 *
	 *   -2 * I**(n+1) * sqrt(π) * ((n-m)/2)! * ((n+m)/2)! * (-1)**j 
	 *   -----------------------------------------------------------
	 *              j! * ((n-m)/2 - j)! * ((n+m)/2 + j)!
	 *
	 * for m + n ODD:
	 *
	 * s1_constfacs[y(m,n)][j] = 0
	 */
	complex double *s1_constfacs_base; ///< Internal pointer holding memory for the 2D Ewald sum constant factors.
	// similarly for the 1D z-axis aligned case; now the indices are [n][j] (as m == 0)
	/// The constant factors for the long range part of a 1D Ewald sum along the \a z axis.
	/** If the summation points lie along a different direction, use the formula for
	 * 2D sum with additional factor of 
	 * \f$ \sqrt{pi} \kappa \gamma(\abs{\vect{k}+\vect{K}}/\kappa) \f$.
	 */
	complex double **s1_constfacs_1Dz; 
	/* These are the actual numbers now:
	 * s1_constfacs_1Dz[n][j] =
	 *   
	 *     -I**(n+1) (-1)**j * n!
	 *   --------------------------
	 *   j! * 2**(2*j) * (n - 2*j)!
	 */
	complex double *s1_constfacs_1Dz_base; ///<Internal pointer holding memory for the 1D Ewald sum constant factors.

	double *legendre0; /* now with GSL_SF_LEGENDRE_NONE normalisation, because this is what is
			    * what the multipliers from translations.c count with.
			    */
	double *legendre_plus1; // needed? TODO; in any case, nonzero only for m=0
	double *legendre_minus1; // needed? TODO; in any case, nonzero only for m=0
	gsl_sf_legendre_t legendre_normconv;
	int legendre_csphase;       /* 1 or -1; csphase of the Legendre polynomials saved in legendre0 etc.
					This is because I dont't actually consider this fixed in
					translations.c */

} qpms_ewald3_constants_t;

/// Constructor for qpms_ewald3_constants_t.
qpms_ewald3_constants_t *qpms_ewald3_constants_init(qpms_l_t lMax, int csphase);
/// Destructor for qpms_ewald3_constants_t.
void qpms_ewald3_constants_free(qpms_ewald3_constants_t *);


/// Structure for holding complex-valued result of computation and an error estimate.
/** Similar to gsl_sf_result, but with complex val. */
typedef struct qpms_csf_result { 
  complex double val; ///< Calculation result.
  double err; ///< Error estimate.
} qpms_csf_result;


// [1, (A.9)]
static inline complex double lilgamma(double t) {
  t = fabs(t);
  if (t >= 1) 
    return sqrt(t*t - 1);
  else
    return -I * sqrt(1 - t*t);
}

// [1, (A.8)], complex version of lilgamma()
static inline complex double clilgamma(complex double z) {
	complex double a1 = z - 1, a2 = z + 1;
	// ensure  -pi/2 < arg(z + 1) < 3*pi/2
	if (creal(a2) < 0 && cimag(a2) <= 0) 
		a2 = -csqrt(a2);
	else 
		a2 = csqrt(a2);
	// ensure -3*pi/2 < arg(z - 1) < pi/2
	if (creal(a1) < 0 && cimag(a1) >= 0)
		a1 = -csqrt(a1);
	else
		a1 = csqrt(a1);
	return a1 * a2;
}

/// Incomplete Gamma function as a series.
/** DLMF 8.7.3 (latter expression) for complex second argument.
 *
 * The principal value is calculated. On the negative real axis
 * (where the function has branch cut), the sign of the imaginary
 * part is what matters (even if it is zero). Therefore one
 * can have
 * `cx_gamma_inc_series_e(a, z1) != cx_gamma_inc_series_e(a, z2)`
 * even if `z1 == z2`, because `-0 == 0` according to IEEE 754.
 * The side of the branch cut can be determined using `signbit(creal(z))`.
 */
int cx_gamma_inc_series_e(double a, complex z, qpms_csf_result * result);

/// Incomplete Gamma function as continued fractions.
/** 
 * The principal value is calculated. On the negative real axis
 * (where the function has branch cut), the sign of the imaginary
 * part is what matters (even if it is zero). Therefore one
 * can have
 * `cx_gamma_inc_CF_e(a, z1) != cx_gamma_inc_CF_e(a, z2)`
 * even if `z1 == z2`, because `-0 == 0` according to IEEE 754.
 * The side of the branch cut can be determined using `signbit(creal(z))`.
 */
int cx_gamma_inc_CF_e(double a, complex z, qpms_csf_result * result);

/// Incomplete gamma for complex second argument.
/** 
 * If x is (almost) real, it just uses gsl_sf_gamma_inc_e(). 
 * 
 * On the negative real axis
 * (where the function has branch cut), the sign of the imaginary
 * part is what matters (even if it is zero). Therefore one
 * can have
 * `complex_gamma_inc_e(a, z1, m) != complex_gamma_inc_e(a, z2, m)`
 * even if `z1 == z2`, because `-0 == 0` according to IEEE 754.
 * The side of the branch cut can be determined using `signbit(creal(z))`.
 *
 * Another than principal branch can be selected using non-zero \a m 
 * argument.
 */
int complex_gamma_inc_e(double a, complex double x,
	/// Branch index.
	/** If zero, the principal value is calculated.	 
	 * Other branches might be chosen using non-zero \a m.
	 * In such case, the returned value corresponds to \f[
	 * \Gamma(a,ze^{2\pi mi})=e^{2\pi mia} \Gamma(a,z) 
	 *   + (1-e^{2\pi mia}) \Gamma(a).
	 * \f]
	 *
	 * If \a a is non-positive integer, the limiting value should
	 * be used, but this is not yet implemented!
	 */
	int m,
	qpms_csf_result *result);

/// Exponential integral for complex second argument.
/** If x is (almost) positive real, it just uses gsl_sf_expint_En_e(). */
int complex_expint_n_e(int n, complex double x, qpms_csf_result *result);


/// Hypergeometric 2F2, used to calculate some errors.
int hyperg_2F2_series(double a, double b, double c, double d, double x, 
		gsl_sf_result *result);

#if 0
// The integral from (4.6); maybe should be static and not here.
int ewald32_sr_integral(double r, double k, double n, double eta, double *result, double *err, gsl_integration_workspace *workspace);
#endif


// General functions acc. to [2], sec. 4.6 – currently valid for 2D and 1D lattices in 3D space

/// The Ewald sum "self-interaction" term that appears in the lattice sums with zero (direct-space) Bravais lattice displacement.
int ewald3_sigma0(complex double *result, ///< Pointer to save the result (single complex double).
	       	double *err, ///< Pointer to save the result error estimate (single double).
		const qpms_ewald3_constants_t *c, ///< Constant factors structure initialised by qpms_ewald3_constants_init().
		double eta, ///< Ewald parameter.
		complex double wavenumber ///< Wavenumber of the background medium.
);

/// Short-range part of outgoing scalar spherical wavefunctions' lattice sum \f$ \sigma_{l,m}^\mathrm{S}(\vect k,\vect s)\f$.
int ewald3_sigma_short(
		complex double *target_sigmasr_y, ///< Target array for \f$ \sigma_{l,m}^\mathrm{S} \f$, must be `c->nelem_sc` long.
		double *target_sigmasr_y_err, ///< Target array for error estimates, must be `c->nelem_sc` long or `NULL`.
		const qpms_ewald3_constants_t *c, ///< Constant factors structure initialised by qpms_ewald3_constants_init().
		double eta, ///< Ewald parameter.
	       	complex double wavenumber, ///< Wavenumber of the background medium. 
		/// Lattice dimensionality. 
		/** Ignored apart from asserts and possible optimisations, as the SR formula stays the same. */
		LatticeDimensionality latdim,
		/// Lattice point generator for the direct Bravais lattice.
		/** There is a possibility that the whole PGen is not consumed
		 *  (this might happen if the summand start to be consistently smaller
		 *  than the (partial) sums * DBL_EPSILON.
		 *  In such case, it is the responsibility of the caller to deallocate
		 *  the generator.
		 */
		PGen *pgen_R,
		/// Indicates whether pgen_R already generates shifted points.
		/** If false, the behaviour corresponds to the old ewald32_sigma_short_points_and_shift(),
		 * so the function assumes that the generated points correspond to the unshifted Bravais lattice,
		 * and adds particle_shift to the generated points before calculations.
		 * If true, it assumes that they are already shifted (if calculating interaction between
		 * different particles in the unit cell).
		 */
		bool pgen_generates_shifted_points,
		/// Wave vector \f$\vect k\f$.
		cart3_t k,
		/// Lattice offset \f$\vect s\f$ wrt. the Bravais lattice.
		cart3_t particle_shift
		);

/// Long-range part of outgoing scalar spherical wavefunctions' lattice sum \f$ \sigma_{l,m}^\mathrm{L}(\vect k,\vect s)\f$.
int ewald3_sigma_long( // calls ewald3_21_sigma_long or ewald3_3_sigma_long, depending on latdim
		complex double *target_sigmalr_y, ///< Target array for \f$ \sigma_{l,m}^\mathrm{L} \f$, must be `c->nelem_sc` long.
		double *target_sigmalr_y_err, ///< Target array for error estimates, must be `c->nelem_sc` long or `NULL`.
		const qpms_ewald3_constants_t *c, ///< Constant factors structure initialised by qpms_ewald3_constants_init().
		double eta, ///< Ewald parameter.
	       	complex double wavenumber, ///< Wavenumber of the background medium. 
	        double unitcell_volume, ///< Volume of the (direct lattice) unit cell (with dimension corresponding to the lattice dimensionality).
		/// Lattice dimensionality. 
		LatticeDimensionality latdim,
		/// Lattice point generator for the reciprocal lattice.
		/** There is a possibility that the whole PGen is not consumed
		 *  (this might happen if the summand start to be consistently smaller
		 *  than the (partial) sums * DBL_EPSILON.
		 *  In such case, it is the responsibility of the caller to deallocate
		 *  the generator.
		 */
		PGen *pgen_K, 
		/// Indicates whether pgen_K already generates shifted points.
		/** If false, the behaviour corresponds to the old ewald32_sigma_long_points_and_shift(),
		 * so the function assumes that the generated points correspond to the unshifted reciprocal Bravais lattice,
		 * and adds beta to the generated points before calculations.
		 * If true, it assumes that they are already shifted.
		 */
		bool pgen_generates_shifted_points,
		/// Wave vector \f$\vect k\f$.
		cart3_t k,
		/// Lattice offset \f$\vect s\f$ wrt. the Bravais lattice.
		cart3_t particle_shift
		);
		
#endif //EWALD_H
