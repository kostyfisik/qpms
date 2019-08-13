#include "qpms_specfunc.h"
#include "qpms_types.h"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include "indexing.h"
#include <string.h>
#include "qpms_error.h"

// Legendre functions also for negative m, see DLMF 14.9.3
qpms_errno_t qpms_legendre_deriv_y_fill(double *target, double *target_deriv, double x, qpms_l_t lMax,
    gsl_sf_legendre_t lnorm, double csphase)
{
  const size_t n = gsl_sf_legendre_array_n(lMax);
  double *legendre_tmp, *legendre_deriv_tmp;
  QPMS_CRASHING_MALLOC(legendre_tmp, n * sizeof(double));
  QPMS_CRASHING_MALLOC(legendre_deriv_tmp, n * sizeof(double));
  int gsl_errno = gsl_sf_legendre_deriv_array_e(
      lnorm, (size_t)lMax, x, csphase, legendre_tmp,legendre_deriv_tmp);
  for (qpms_l_t l = 1; l <= lMax; ++l)
    for (qpms_m_t m = 0; m <= l; ++m) {
      qpms_y_t y = qpms_mn2y(m,l);
      size_t i = gsl_sf_legendre_array_index(l,m);
      target[y] = legendre_tmp[i];
      target_deriv[y] = legendre_deriv_tmp[i];
    }

  // Fill negative m's.
  for (qpms_l_t l = 1; l <= lMax; ++l)
    for (qpms_m_t m = 1; m <= l; ++m) {
      qpms_y_t y = qpms_mn2y(-m,l);
      size_t i = gsl_sf_legendre_array_index(l,m);
      // cf. DLMF 14.9.3, but we're normalised.
      double factor = ((m%2)?-1:1); 
      target[y] = factor * legendre_tmp[i];
      target_deriv[y] = factor * legendre_deriv_tmp[i];
    }

  free(legendre_tmp);
  free(legendre_deriv_tmp);
  return gsl_errno;
}

qpms_errno_t qpms_legendre_deriv_y_get(double **target, double **dtarget, double x, qpms_l_t lMax, gsl_sf_legendre_t lnorm,
    double csphase)
{

  const qpms_y_t nelem = qpms_lMax2nelem(lMax);
  QPMS_CRASHING_MALLOC(target, nelem * sizeof(double));
  QPMS_CRASHING_MALLOC(dtarget, nelem * sizeof(double));
  return qpms_legendre_deriv_y_fill(*target, *dtarget, x, lMax, lnorm, csphase);
}

qpms_pitau_t qpms_pitau_get(double theta, qpms_l_t lMax, const double csphase)
{
  qpms_pitau_t res;
  const qpms_y_t nelem = qpms_lMax2nelem(lMax);
  QPMS_CRASHING_MALLOC(res.leg, nelem * sizeof(double));
  QPMS_CRASHING_MALLOC(res.pi, nelem * sizeof(double));
  QPMS_CRASHING_MALLOC(res.tau, nelem * sizeof(double));
  qpms_pitau_fill(res.leg, res.pi, res.tau, theta, lMax, csphase);
  return res;
}

qpms_errno_t qpms_pitau_fill(double *target_leg, double *pi, double *tau, double theta, qpms_l_t lMax, double csphase)
{
  QPMS_ENSURE(fabs(csphase) == 1, "The csphase argument must be either 1 or -1, not %g.", csphase);
  const qpms_y_t nelem = qpms_lMax2nelem(lMax);
 
  double ct = cos(theta), st = sin(theta);
  if (1 == fabs(ct)) { // singular case, use DLMF 14.8.2
    if(pi) memset(pi, 0, nelem * sizeof(double));
    if(tau) memset(tau, 0, nelem * sizeof(double));
    if(target_leg) memset(target_leg, 0, nelem * sizeof(double));
    for (qpms_l_t l = 1; l <= lMax; ++l) {
      if(target_leg) target_leg[qpms_mn2y(0, l)] = ((l%2)?ct:1.)*sqrt((2*l+1)/(4*M_PI *l*(l+1)));
      double fl = 0.25 * sqrt((2*l+1)*M_1_PI);
      int lpar = (l%2)?-1:1;
      if(pi) {
        pi [qpms_mn2y(+1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
        pi [qpms_mn2y(-1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
      }
      if(tau) {
        tau[qpms_mn2y(+1, l)] = ((ct>0) ? +1 : lpar) * fl * csphase;
        tau[qpms_mn2y(-1, l)] = -((ct>0) ? +1 : lpar) * fl * csphase;
      }
    }
  }
  else { // cos(theta) in (-1,1), use normal calculation
    double *leg, *legder;
    if (target_leg)
      leg = target_leg;
    else
      QPMS_CRASHING_MALLOC(leg, nelem*sizeof(double));
    QPMS_CRASHING_MALLOC(legder, nelem * sizeof(double));
    QPMS_ENSURE_SUCCESS(qpms_legendre_deriv_y_fill(leg, legder, ct, lMax,
          GSL_SF_LEGENDRE_SPHARM, csphase));
    // Multiply by the "power normalisation" factor
    for (qpms_l_t l = 1; l <= lMax; ++l) {
      double prefac = 1./sqrt(l*(l+1));
      for (qpms_m_t m = -l; m <= l; ++m) {
        leg[qpms_mn2y(m,l)] *= prefac;
        legder[qpms_mn2y(m,l)] *= prefac;
      }
    }
    for (qpms_l_t l = 1; l <= lMax; ++l) {
      for (qpms_m_t m = -l; m <= l; ++m) {
        if(pi)  pi [qpms_mn2y(m,l)] = m / st * leg[qpms_mn2y(m,l)];
        if(tau) tau[qpms_mn2y(m,l)] = - st * legder[qpms_mn2y(m,l)];
      }
    }
    free(legder);
    if(!target_leg)
      free(leg);
  }
  return QPMS_SUCCESS;
}

void qpms_pitau_free(qpms_pitau_t x) {
  free(x.leg);
  free(x.pi);
  free(x.tau);
}

