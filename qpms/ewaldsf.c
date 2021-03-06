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
#include <Faddeeva.h>
#include "tiny_inlines.h"

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
    f = -f + 1;
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

// Complex square root with branch selection
static inline complex double csqrt_branch(complex double x, int xbranch) {
  return csqrt(x) * min1pow(xbranch);
}

// The Delta_n factor from [Kambe II], Appendix 3
// \f[ \Delta_n = \int_n^\infty t^{-1/2 - n} \exp(-t + z^2/(4t))\ud t \f]
void ewald3_2_sigma_long_Delta_recurrent(complex double *target, double *err,
    int maxn, complex double x, int xbranch, complex double z) {
  complex double expfac = cexp(-x + 0.25 * z*z / x);
  complex double sqrtx = csqrt_branch(x, xbranch); // TODO check carefully, which branch is needed
  // These are used in the first two recurrences
  complex double w_plus  = Faddeeva_w(+z/(2*sqrtx) + I*sqrtx, 0);
  complex double w_minus = Faddeeva_w(-z/(2*sqrtx) + I*sqrtx, 0);
  QPMS_ASSERT(maxn >= 0); 
  if (maxn >= 0) 
    target[0] = 0.5 * M_SQRTPI * expfac * (w_minus + w_plus);
  if (maxn >= 1) 
    target[1] = I / z * M_SQRTPI * expfac * (w_minus - w_plus);
  for(int n = 1; n < maxn; ++n) { // The rest via recurrence
    // TODO The cpow(x, 0.5 - n) might perhaps better be replaced with a recurrently computed variant
    target[n+1] = -(4 / (z*z)) * (-(0.5 - n) * target[n] + target[n-1] - sqrtx * cpow(x, -n) * expfac);
  }
  if (err) {
    // The error estimates for library math functions are based on 
    // https://www.gnu.org/software/libc/manual/html_node/Errors-in-Math-Functions.html
    // and are not guaranteed to be extremely precise.
    // The error estimate for Faddeeva's functions is based on the package's authors at
    // http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
    // "we find that the accuracy is typically at at least 13 significant digits in both the real and imaginary parts"
    // FIXME the error estimate seems might be off by several orders of magnitude (try parameters x = -3, z = 0.5, maxn=20)
    double w_plus_abs = cabs(w_plus), w_minus_abs = cabs(w_minus), expfac_abs = cabs(expfac);
    double w_plus_err = w_plus_abs * 1e-13, w_minus_err = w_minus_abs * 1e-13; // LPTODO argument error contrib.
    double expfac_err = expfac_abs * (4 * DBL_EPSILON); // LPTODO add argument error contrib.
    double z_abs = cabs(z);
    double z_err = z_abs * DBL_EPSILON;
    double x_abs = cabs(x);
    if (maxn >= 0) 
      err[0] = 0.5 * M_SQRTPI * (expfac_abs * (w_minus_err + w_plus_err) + (w_minus_abs + w_plus_abs) * expfac_err);
    if (maxn >= 1) 
      err[1] = 2 * err[0] / z_abs + cabs(target[1]) * z_err / (z_abs*z_abs);
    for(int n = 1; n < maxn; ++n) {
      err[n+1] = (2 * cabs(target[n+1]) / z_abs + 4 * ((0.5+n) * err[n] + err[n-1] +
            pow(x_abs, 0.5 - n) * (2*DBL_EPSILON * expfac_abs + expfac_err))  // LPTODO not ideal, pow() call is an overkill
            ) * z_err / (z_abs*z_abs);
    }
  }
}


void ewald3_2_sigma_long_Delta_series(complex double *target, double *err,
    int maxn, complex double x, int xbranch, complex double z) {
  complex double w = 0.25*z*z;
  double w_abs = cabs(w);
  int maxk;
  if (w_abs == 0)
    maxk = 0; // Delta is equal to the respective incomplete Gamma functions
  else {
    // Estimate a suitable maximum k, using Stirling's formula, so that w**maxk / maxk! is less than DBL_EPSILON
    // This implementation is quite stupid, but it is still cheap compared to the actual computation, so LPTODO better one
    maxk = 1;
    double log_w_abs = log(w_abs);
    while (maxk * (log_w_abs - log(maxk) + 1) >= -DBL_MANT_DIG)
      ++maxk;
  }
  // TODO asserts on maxn, maxk

  complex double *Gammas;
  double *Gammas_err = NULL, *Gammas_abs = NULL;
  QPMS_CRASHING_CALLOC(Gammas, maxk+maxn+1, sizeof(*Gammas));
  if(err) {
    QPMS_CRASHING_CALLOC(Gammas_err, maxk+maxn+1, sizeof(*Gammas_err));
    QPMS_CRASHING_MALLOC(Gammas_abs, (maxk+maxn+1) * sizeof(*Gammas_abs));
  }

  for(int j = 0; j <= maxn+maxk; ++j) {
    qpms_csf_result g;
    QPMS_ENSURE_SUCCESS(complex_gamma_inc_e(0.5-j, x, xbranch, &g));
    Gammas[j] = g.val;
    if(err) {
      Gammas_abs[j] = cabs(g.val);
      Gammas_err[j] = g.err;
    }
  }

  for(int n = 0; n <= maxn; ++n) target[n] = 0;
  if(err) for(int n = 0; n <= maxn; ++n) err[n] = 0;

  complex double wpowk_over_fack = 1.;
  double wpowk_over_fack_abs = 1.;
  for(int k = 0; k <= maxk; ++k, wpowk_over_fack *= w/k) { // TODO? Kahan sum, Horner's method?
    // Also TODO? for small n, continue for higher k if possible/needed
    for(int n = 0; n <= maxn; ++n) {
      target[n] += Gammas[n+k] * wpowk_over_fack;
      if(err) {
        // DBL_EPSILON might not be the best estimate here, but...
        err[n] += wpowk_over_fack_abs * Gammas_err[n+k] + DBL_EPSILON * Gammas_abs[n+k];
        wpowk_over_fack_abs *= w_abs / (k+1);
      }
    }
  }

  // TODO add an error estimate for the k-cutoff!!!

  free(Gammas);
  free(Gammas_err);
  free(Gammas_abs);
}


void ewald3_2_sigma_long_Delta(complex double *target, double *err,
    int maxn, complex double x, int xbranch, complex double z) {
  double absz = cabs(z);
  if (absz < 2.) // TODO take into account also the other parameters
    ewald3_2_sigma_long_Delta_series(target, err, maxn, x, xbranch, z);
  else
    ewald3_2_sigma_long_Delta_recurrent(target, err, maxn, x, xbranch, z);
}

