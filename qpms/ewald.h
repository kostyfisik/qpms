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

#endif //EWALD_H
