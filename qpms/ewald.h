#ifndef EWALD_H
#define EWALD_H
#include <gsl/gsl_sf_result.h>
#include <math.h>
#include <complex.h>

typedef struct { // as gsl_sf_result, but with complex val
  complex double val;
  double err;
} qpms_csf_result;


// Linton&Thompson (A.9)
// TODO put into a header file as inline
static inline complex double lilgamma(double t) {
  t = fabs(t);
  if (t >= 1) 
    return sqrt(t*t - 1);
  else
    return -I * sqrt(1 - t*t);
}


// DLMF 8.7.3 (latter expression) for complex second argument
int cx_gamma_inc_series_e(double a, complex z, qpms_csf_result * result);

// incomplete gamma for complex second argument
// if x is (almost) real, it just uses gsl_sf_gamma_inc_e
int complex_gamma_inc_e(double a, complex double x, qpms_csf_result *result);

#endif //EWALD_H
