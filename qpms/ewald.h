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

/* 
 * Implementation of two-dimensional lattice sum in three dimensions
 * according to:
 * [1] C.M. Linton, I. Thompson
 *     Journal of Computational Physics 228 (2009) 1815–1829
 * [2] C.M.Linton
 *     SIAM Review Vol 52, No. 4, pp. 630–674
 */

/*
 * N.B.!!! currently, the long-range parts are calculated not according to [1,(4.5)], but rather
 * according to the spherical-harmonic-normalisation-independent formulation in my notes
 * notes/ewald.lyx.
 * Both parts of lattice sums are then calculated with the P_n^|m| e^imf substituted in place
 *                                             // (N.B. or P_n^|m| e^imf (-1)^m for negative m) 
 * of Y_n^m (this is quite a weird normalisation especially for negative |m|, but it is consistent
 * with the current implementation of the translation coefficients in translations.c;
 * in the long run, it might make more sense to replace it everywhere with normalised
 * Legendre polynomials).
 */


// Use this handler to ignore underflows of incomplete gamma.
gsl_error_handler_t IgnoreUnderflowsGSLErrorHandler;


/* Object holding the constant factors from [1, (4.5)] */ 
typedef struct { 
	qpms_l_t lMax;
	qpms_y_t nelem_sc;
	qpms_l_t *s1_jMaxes;
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
	complex double *s1_constfacs_base; // internal pointer holding the memory for the constants
	// similarly for the 1D z-axis aligned case; now the indices are [n][j] (as m == 0)
	complex double **s1_constfacs_1Dz; 
	/* These are the actual numbers now:
	 * s1_consstfacs_1Dz[n][j] =
	 *   
	 *     -I**(n+1) (-1)**j * n!
	 *   --------------------------
	 *   j! * 2**(2*j) * (n - 2*j)!
	 */
	complex double *s1_constfacs_1Dz_base;

	double *legendre0; /* now with GSL_SF_LEGENDRE_NONE normalisation, because this is what is
			    * what the multipliers from translations.c count with.
			    */
	double *legendre_plus1; // needed? TODO; in any case, nonzero only for m=0
	double *legendre_minus1; // needed? TODO; in any case, nonzero only for m=0
	gsl_sf_legendre_t legendre_normconv;
	int legendre_csphase;       /* 1 or -1; csphase of the Legendre polynomials saved in legendre0 etc.
					This is because I dont't actually consider this fixed in
					translations.c */

} qpms_ewald32_constants_t;

qpms_ewald32_constants_t *qpms_ewald32_constants_init(qpms_l_t lMax, int csphase);
void qpms_ewald32_constants_free(qpms_ewald32_constants_t *);


typedef struct { // as gsl_sf_result, but with complex val
  complex double val;
  double err;
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


 


// Incomplete Gamma function as series
// DLMF 8.7.3 (latter expression) for complex second argument
int cx_gamma_inc_series_e(double a, complex z, qpms_csf_result * result);

// Incomplete gamma for complex second argument
// if x is (almost) real, it just uses gsl_sf_gamma_inc_e
int complex_gamma_inc_e(double a, complex double x, qpms_csf_result *result);

// Exponential integral for complex second argument
// If x is (almost) positive real, it just uses gsl_sf_expint_En_e
int complex_expint_n_e(int n, complex double x, qpms_csf_result *result);


// hypergeometric 2F2, used to calculate some errors
int hyperg_2F2_series(const double a, const double b, const double c, const double d,
    const double x, gsl_sf_result *result);

#if 0
// The integral from (4.6); maybe should be static and not here.
int ewald32_sr_integral(double r, double k, double n, double eta, double *result, double *err, gsl_integration_workspace *workspace);
#endif

#include "lattices.h"

// General functions acc. to [2], sec. 4.6 – currently valid for 2D and 1D lattices in 3D space

int ewald3_sigma0(complex double *result, double *err,
		const qpms_ewald32_constants_t *c,
		double eta, complex double k
);

int ewald3_sigma_short(
		complex double *target_sigmasr_y, // must be c->nelem_sc long
		double *target_sigmasr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c,
		const double eta, const complex double k,
		const LatticeDimensionality latdim, // apart from asserts and possible optimisations ignored, as the SR formula stays the same
		PGen *pgen_R, const bool pgen_generates_shifted_points 
		/* If false, the behaviour corresponds to the old ewald32_sigma_short_points_and_shift,
		 * so the function assumes that the generated points correspond to the unshifted Bravais lattice,
		 * and adds particle_shift to the generated points before calculations.
		 * If true, it assumes that they are already shifted (if calculating interaction between
		 * different particles in the unit cell).
		 */,
		const cart3_t beta,
		const cart3_t particle_shift
		);

int ewald3_sigma_long( // calls ewald3_21_sigma_long or ewald3_3_sigma_long, depending on latdim
		complex double *target_sigmalr_y, // must be c->nelem_sc long
		double *target_sigmalr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c,
		const double eta, const complex double k,
	        const double unitcell_volume /* with the corresponding lattice dimensionality */,
		const LatticeDimensionality latdim,
		PGen *pgen_K, const bool pgen_generates_shifted_points 
		/* If false, the behaviour corresponds to the old ewald32_sigma_long_points_and_shift,
		 * so the function assumes that the generated points correspond to the unshifted reciprocal Bravais lattice,
		 * and adds beta to the generated points before calculations.
		 * If true, it assumes that they are already shifted.
		 */,
		const cart3_t beta,
		const cart3_t particle_shift
		);
		
/// !!!!!!!!!!!!!!! ZDE JSEM SKONČIL !!!!!!!!!!!!!!!!!!!!!!.


int ewald32_sigma0(complex double *result, double *err, // actually, this should be only alias for ewald3_sigma0
		const qpms_ewald32_constants_t *c,
		double eta, double k
);

// TODO make "compressed versions" where the (m+n)-odd terms (which are zero)
// are not included.

int ewald32_sigma_long_shiftedpoints ( 
		complex double *target_sigmalr_y, // must be c->nelem_sc long
		double *target_sigmalr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c,
		double eta, double k, double unitcell_area,
		size_t npoints, const point2d *Kpoints_plus_beta,
		point2d beta,
		point2d particle_shift
);
int ewald32_sigma_long_points_and_shift (
		complex double *target_sigmalr_y, // must be c->nelem_sc long
		double *target_sigmalr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c,
		double eta, double k, double unitcell_area,
		size_t npoints, const point2d *Kpoints,
		point2d beta,
		point2d particle_shift
		);
int ewald32_sigma_long_shiftedpoints_rordered(//NI
		complex double *target_sigmalr_y, // must be c->nelem_sc long
		double *target_sigmalr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c,
		double eta, double k, double unitcell_area,
		const points2d_rordered_t *Kpoints_plus_beta_rordered,
		point2d particle_shift
		);

int ewald32_sigma_short_shiftedpoints(
		complex double *target_sigmasr_y, // must be c->nelem_sc long
		double *target_sigmasr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c, // N.B. not too useful here
		double eta, double k,
		size_t npoints, const point2d *Rpoints_plus_particle_shift,
		point2d beta,
		point2d particle_shift           // used only in the very end to multiply it by the phase
		);
int ewald32_sigma_short_points_and_shift(
		complex double *target_sigmasr_y, // must be c->nelem_sc long
		double *target_sigmasr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c, // N.B. not too useful here
		double eta, double k,
		size_t npoints, const point2d *Rpoints, 
		point2d beta,
		point2d particle_shift
		);
int ewald32_sigma_short_points_rordered(//NI
		complex double *target_sigmasr_y, // must be c->nelem_sc long
		double *target_sigmasr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c, // N.B. not too useful here
		double eta, double k,
		const points2d_rordered_t *Rpoints_plus_particle_shift_rordered,
		point2d particle_shift    // used only in the very end to multiply it by the phase
		);


// 1D sums aligned along z-axis
int ewald31z_sigma_long_points_and_shift (
		complex double *target_sigmalr_y, // must be c->nelem_sc long
		double *target_sigmalr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c,
		double eta, double k, double unitcell_area,
		size_t npoints, const double *Kpoints,
		double beta,
		double particle_shift
		);
int ewald31z_sigma_short_points_and_shift(
		complex double *target_sigmasr_y, // must be c->nelem_sc long
		double *target_sigmasr_y_err, // must be c->nelem_sc long or NULL
		const qpms_ewald32_constants_t *c, // N.B. not too useful here
		double eta, double k,
		size_t npoints, const double *Rpoints, 
		double beta,
		double particle_shift
		);
int ewald31z_sigma0(complex double *result, double *err, 
		const qpms_ewald32_constants_t *c,
		double eta, double k
		); // exactly the same as ewald32_sigma0


#endif //EWALD_H
