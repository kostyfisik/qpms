#ifndef EWALD_H
#define EWALD_H
#include <gsl/gsl_sf_result.h>
#include <math.h> // for inlined lilgamma
#include <complex.h>
#include "qpms_types.h"

/* 
 * Implementation of two-dimensional lattice sum in three dimensions
 * according to:
 * [1] C.M. Linton, I. Thompson
 *     Journal of Computational Physics 228 (2009) 1815–1829
 */

/* Object holding the constant factors from [1, (4.5)] */ 
typedef struct { 
	qpms_l_t lMax;
	qpms_y_t nelem;
	qpms_l_t *s1_jMaxes;
	complex double **s1_constfacs; // indices [y][j] where j is same as in [1, (4.5)]
	complex double *s1_constfacs_base; // internal pointer holding the memory for the constants
} qpms_ewald32_constants_t;

qpms_ewald32_constants_t *qpms_ewald32_constants_init(qpms_l_t lMax);
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


// Incomplete Gamma function as series
// DLMF 8.7.3 (latter expression) for complex second argument
int cx_gamma_inc_series_e(double a, complex z, qpms_csf_result * result);

// Incomplete gamma for complex second argument
// if x is (almost) real, it just uses gsl_sf_gamma_inc_e
int complex_gamma_inc_e(double a, complex double x, qpms_csf_result *result);


#include "lattices.h"

// TODO make "compressed versions" where the (m+n)-odd terms (which are zero)
// are not included.

int ewald32_sigma_long_shiftedpoints_e ( 
		complex double *target_sigmalr_y, // must be c->nelem long
		double *target_sigmalr_y_err, // must be c->nelem long or NULL
		const qpms_ewald32_constants_t *c,
		double eta, double k, double unitcell_area,
		size_t npoints, const point2d *Kpoints_plus_beta,
		point2d particle_shift
);
int ewald32_sigma_long_points_and_shift (
		complex double *target_sigmalr_y, // must be c->nelem long
		double *target_sigmalr_y_err, // must be c->nelem long or NULL
		const qpms_ewald32_constants_t *c,
		double eta, double k, double unitcell_area,
		size_t npoints, const point2d *Kpoints,
		point2d beta,
		point2d particle_shift
		);
int ewald32_sigma_long_shiftedpoints_rordered(
		complex double *target_sigmalr_y, // must be c->nelem long
		double *target_sigmalr_y_err, // must be c->nelem long or NULL
		const qpms_ewald32_constants_t *c,
		double eta, double k, double unitcell_area,
		const points2d_rordered_t *Kpoints_plus_beta_rordered,
		point2d particle_shift
		);

int ewald32_sigma_short_shiftedpoints(
		complex double *target_sigmasr_y, // must be c->nelem long
		double *target_sigmasr_y_err, // must be c->nelem long or NULL
		const qpms_ewald32_constants_t *c, // N.B. not too useful here
		double eta, double k,
		size_t npoints, const point2d *Rpoints_plus_particle_shift,
		point2d particle_shift           // used only in the very end to multiply it by the phase
		);
int ewald32_sigma_short_points_and_shift(
		complex double *target_sigmasr_y, // must be c->nelem long
		double *target_sigmasr_y_err, // must be c->nelem long or NULL
		const qpms_ewald32_constants_t *c, // N.B. not too useful here
		double eta, double k,
		size_t npoints, const point2d *Rpoints, 
		point2d particle_shift
		);
int ewald32_sigma_short_points_rordered(
		complex double *target_sigmasr_y, // must be c->nelem long
		double *target_sigmasr_y_err, // must be c->nelem long or NULL
		const qpms_ewald32_constants_t *c, // N.B. not too useful here
		double eta, double k,
		const points2d_rordered_t *Rpoints_plus_particle_shift_rordered,
		point2d particle_shift    // used only in the very end to multiply it by the phase
		);



#endif //EWALD_H
