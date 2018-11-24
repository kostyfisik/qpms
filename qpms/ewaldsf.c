#include "ewald.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_machine.h> // Maybe I should rather use DBL_EPSILON instead of GSL_DBL_EPSILON.
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
int cx_gamma_inc_series_e(const double a, const complex double z, qpms_csf_result * result) {
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
    complex double summand = - cpow(z, k) / fullgamma_ak.val; // TODO test branch selection here with cimag(z) = -0.0
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

// Exponential integral for complex argument; !UNTESTED! and probably not needed, as I expressed everything in terms of inc. gammas anyways.
int complex_expint_n_e(int n, complex double x, qpms_csf_result *result) {
  if (creal(x) >= 0 &&
    (0 == fabs(cimag(x)) || // x is real positive; just use the real fun
    fabs(cimag(x)) < fabs(creal(x)) * COMPLEXPART_REL_ZERO_LIMIT)) {
    gsl_sf_result real_expint_result;
    int retval = gsl_sf_expint_En_e(n, creal(x), &real_expint_result);
    result->val = real_expint_result.val;
    result->err = real_expint_result.err;
    return retval;
  } else {
    int retval = complex_gamma_inc_e(-n+1, x, result);
    complex double f = cpow(x, 2*n-2);
    retval.val *= f;
    retval.err *= cabs(f);
    return retval;
  }
}

// inspired by GSL's hyperg_2F1_series
int hyperg_2F2_series(const double a, const double b, const double c, const double d,
    const double x, gsl_sf_result *result
    )
{
  double sum_pos = 1.0;
  double sum_neg = 0.0;
  double del_pos = 1.0;
  double del_neg = 0.0;
  double del = 1.0;
  double del_prev;
  double k = 0.0;
  int i = 0;

  if(fabs(c) < GSL_DBL_EPSILON || fabs(d) < GSL_DBL_EPSILON) {
    result->val = NAN; 
    result->err = INFINITY;
    GSL_ERROR ("error", GSL_EDOM);
  }

  do {
    if(++i > 30000) {
      result->val  = sum_pos - sum_neg;
      result->err  = del_pos + del_neg;
      result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
      result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k)+1.0) * fabs(result->val);
      GSL_ERROR ("error", GSL_EMAXITER);
    }
    del_prev = del;
    del *= (a+k)*(b+k) * x / ((c+k) * (d+k) * (k+1.0));  /* Gauss series */

    if(del > 0.0) {
      del_pos  =  del;
      sum_pos +=  del;
    }
    else if(del == 0.0) {
      /* Exact termination (a or b was a negative integer).
       */
      del_pos = 0.0;
      del_neg = 0.0;
      break;
    }
    else {
      del_neg  = -del;
      sum_neg -=  del;
    }

    /*
     * This stopping criteria is taken from the thesis
     * "Computation of Hypergeometic Functions" by J. Pearson, pg. 31
     * (http://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf)
     * and fixes bug #45926
     */
    if (fabs(del_prev / (sum_pos - sum_neg)) < GSL_DBL_EPSILON &&
        fabs(del / (sum_pos - sum_neg)) < GSL_DBL_EPSILON)
      break;

    k += 1.0;
  } while(fabs((del_pos + del_neg)/(sum_pos-sum_neg)) > GSL_DBL_EPSILON);

  result->val  = sum_pos - sum_neg;
  result->err  = del_pos + del_neg;
  result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
  result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k) + 1.0) * fabs(result->val);

  return GSL_SUCCESS;
}



