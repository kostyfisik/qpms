#include "ewald.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>
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
// This can't be used for non-positive integer a due to
// Г(a) in the formula.
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

/* Continued fraction which occurs in evaluation
 * of Q(a,z) or Gamma(a,z).
 * Borrowed from GSL and adapted for complex z.
 *
 *              1   (1-a)/z  1/z  (2-a)/z   2/z  (3-a)/z
 *   F(a,z) =  ---- ------- ----- -------- ----- -------- ...
 *             1 +   1 +     1 +   1 +      1 +   1 +
 *
 */
static int
cx_gamma_inc_F_CF(const double a, const complex double z, qpms_csf_result * result)
{
  const int    nmax  =  5000;
  const double small =  DBL_EPSILON * DBL_EPSILON * DBL_EPSILON;

  complex double hn = 1.0;           /* convergent */
  complex double Cn = 1.0 / small;
  complex double Dn = 1.0;
  int n;

  /* n == 1 has a_1, b_1, b_0 independent of a,z,
     so that has been done by hand                */
  for ( n = 2 ; n < nmax ; n++ )
  {
    complex double an;
    complex double delta;

    if(n % 2)
      an = 0.5*(n-1)/z;
    else
      an = (0.5*n-a)/z;

    Dn = 1.0 + an * Dn;
    if ( cabs(Dn) < small )
      Dn = small;
    Cn = 1.0 + an/Cn;
    if ( cabs(Cn) < small )
      Cn = small;
    Dn = 1.0 / Dn;
    delta = Cn * Dn;
    hn *= delta;
    if(cabs(delta-1.0) < DBL_EPSILON) break;
  }

  result->val = hn;
  result->err = 2.0*GSL_DBL_EPSILON * cabs(hn);
  result->err += GSL_DBL_EPSILON * (2.0 + 0.5*n) * cabs(result->val);

  if(n == nmax)
    GSL_ERROR ("error in CF for F(a,x)", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}

// Incomplete gamma fuction with complex second argument as continued fraction.
int cx_gamma_inc_CF_e(const double a, const complex double z, qpms_csf_result *result) 
{
  qpms_csf_result F;
  gsl_sf_result pre;
  const complex double am1lgz = (a-1.0)*clog(z); // TODO check branches
  const int stat_F = cx_gamma_inc_F_CF(a, z, &F);
  const int stat_E = gsl_sf_exp_err_e(creal(am1lgz - z), GSL_DBL_EPSILON*cabs(am1lgz), &pre);
  complex double cpre = pre.val * cexp(I*cimag(am1lgz - z));// TODO add the error estimate for this
  //complex double cpre = cpow(z, a-1) * cexp(-z);


  result->val = F.val * cpre;
  result->err = fabs(F.err * pre.val) + fabs(F.val * pre.err);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_ERROR_SELECT_2(stat_F, stat_E);
}


// Incomplete gamma function for complex second argument.
int complex_gamma_inc_e(double a, complex double x, int m, qpms_csf_result *result) {
  int retval;
  if (creal(x) >= 0 &&
      (0 == fabs(cimag(x)) || // x is real positive; just use the real fun
      fabs(cimag(x)) < fabs(creal(x)) * COMPLEXPART_REL_ZERO_LIMIT)) {
    gsl_sf_result real_gamma_inc_result;
    retval = gsl_sf_gamma_inc_e(a, creal(x), &real_gamma_inc_result);
    result->val = real_gamma_inc_result.val;
    result->err = real_gamma_inc_result.err;
  } else if (creal(x) >= 0 && cabs(x) > 0.5)
    retval = cx_gamma_inc_CF_e(a, x, result);
  else if (QPMS_LIKELY(a > 0 || fmod(a, 1.0)))
    retval = cx_gamma_inc_series_e(a, x, result);
  else
  /* FIXME cx_gamma_inc_series_e() probably fails for non-positive integer a.
   * This does not matter for 2D lattices in 3D space, 
   * but it might cause problems in the other cases.
   */
    QPMS_NOT_IMPLEMENTED("Incomplete Gamma function with non-positive integer a.");
  if (m) { // Non-principal branch.
    /* This might be sub-optimal, as Γ(a) has probably been already evaluated
     * somewhere in the functions called above. */
    gsl_sf_result fullgamma;
    int retval_fg = gsl_sf_gamma_e(a, &fullgamma);
    if (GSL_EUNDRFLW == retval_fg)
      fullgamma.err += DBL_MIN;
    else if (GSL_SUCCESS != retval_fg){
      result->val = NAN + NAN*I; result->err = NAN;
      return GSL_ERROR_SELECT_2(retval_fg, retval);
    }
    complex double f = cexp(2*m*M_PI*a*I);
    result->val *= f;
    f = 1 - f;
    result->err += cabs(f) * fullgamma.err;
    result->val += f * fullgamma.val;
  }
  return retval;
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
    int retval = complex_gamma_inc_e(-n+1, x, 0, result);
    complex double f = cpow(x, 2*n-2);
    result->val *= f;
    result->err *= cabs(f);
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



