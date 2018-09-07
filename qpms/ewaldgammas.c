#include "ewald.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include "kahansum.h"
#include <math.h>
#include <complex.h>
//#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <float.h>
#include <stdbool.h>

#ifndef COMPLEXPART_REL_ZERO_LIMIT
#define COMPLEXPART_REL_ZERO_LIMIT 1e-14
#endif

gsl_error_handler_t IgnoreUnderflowsGSLErrorHandler;

void IgnoreUnderflowsGSLErrorHandler (const char * reason, 
                    const char * file,
                    const int line, 
                    const int gsl_errno) {
  if (gsl_errno == GSL_EUNDRFLW)
    return;

  gsl_stream_printf ("ERROR", file, line, reason);

  fflush(stdout);
  fprintf (stderr, "Underflow-ignoring error handler invoked.\n");
  fflush(stderr);

  abort();
}


// DLMF 8.7.3 (latter expression) for complex second argument
// BTW if a is large negative, it might take a while to evaluate.
int cx_gamma_inc_series_e(double a, complex z, qpms_csf_result * result) {
  if (a <= 0 && a == (int) a) {
    result->val = NAN + NAN*I;
    result->err = NAN;
    GSL_ERROR("Undefined for non-positive integer values", GSL_EDOM);
  }
  gsl_sf_result fullgamma;
  int retval = gsl_sf_gamma_e(a, &fullgamma);
  if (GSL_EUNDRFLW == retval)
    result->err += DBL_MIN;
  else if (GSL_SUCCESS != retval){
    result->val = NAN + NAN*I; result->err = NAN;
    return retval;
  }

  complex double sumprefac = cpow(z, a) * cexp(-z);
  double sumprefac_abs = cabs(sumprefac);
  complex double sum, sumc; ckahaninit(&sum, &sumc);
  double err, errc; kahaninit(&err, &errc);

  bool breakswitch = false;
  for (int k = 0; (!breakswitch) && (a + k + 1 <= GSL_SF_GAMMA_XMAX); ++k) {
    gsl_sf_result fullgamma_ak;
    if (GSL_SUCCESS != (retval = gsl_sf_gamma_e(a+k+1, &fullgamma_ak))) {
      result->val = NAN + NAN*I; result->err = NAN;
      return retval;
    }
    complex double summand = - cpow(z, k) / fullgamma_ak.val;
    ckahanadd(&sum, &sumc, summand);
    double summanderr = fabs(fullgamma_ak.err * cabs(summand / fullgamma_ak.val));
    // TODO add also the rounding error
    kahanadd(&err, &errc, summanderr);
    // TODO put some smarter cutoff break here?
    if (a + k >= 18 && (cabs(summand) < err || cabs(summand) < DBL_EPSILON))
      breakswitch = true;
  }
  sum *= sumprefac; // Not sure if not breaking the Kahan summation here
  sumc *= sumprefac;
  err *= sumprefac_abs;
  errc *= sumprefac_abs;
  ckahanadd(&sum, &sumc, 1.);
  kahanadd(&err, &errc, DBL_EPSILON);
  result->err = cabs(sum) * fullgamma.err + err * fabs(fullgamma.val);
  result->val = sum * fullgamma.val; // + sumc*fullgamma.val???
  if (breakswitch)
    return GSL_SUCCESS;
  else GSL_ERROR("Overflow; the absolute value of the z argument is probably too large.", GSL_ETOL); // maybe different error code...
}


// incomplete gamma for complex second argument
int complex_gamma_inc_e(double a, complex double x, qpms_csf_result *result) {
  if (creal(x) >= 0 &&
      (0 == fabs(cimag(x)) || // x is real positive; just use the real fun
      fabs(cimag(x)) < fabs(creal(x)) * COMPLEXPART_REL_ZERO_LIMIT)) {
    gsl_sf_result real_gamma_inc_result;
    int retval = gsl_sf_gamma_inc_e(a, creal(x), &real_gamma_inc_result);
    result->val = real_gamma_inc_result.val;
    result->err = real_gamma_inc_result.err;
    return retval;
  } else 
    return cx_gamma_inc_series_e(a, x, result);
}

