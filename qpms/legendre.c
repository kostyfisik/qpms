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
  QPMS_ENSURE(fabs(csphase) == 1, "The csphase argument must be either 1 or -1, not %g.", csphase);
  qpms_pitau_t res;
  const qpms_y_t nelem = qpms_lMax2nelem(lMax);
  QPMS_CRASHING_MALLOC(res.leg, nelem * sizeof(double));
  QPMS_CRASHING_MALLOC(res.pi, nelem * sizeof(double));
  QPMS_CRASHING_MALLOC(res.tau, nelem * sizeof(double));
  double ct = cos(theta), st = sin(theta);
  if (1 == fabs(ct)) { // singular case, use DLMF 14.8.2
    memset(res.pi, 0, nelem * sizeof(double));
    memset(res.tau, 0, nelem * sizeof(double));
    memset(res.leg, 0, nelem * sizeof(double));
    for (qpms_l_t l = 1; l <= lMax; ++l) {
      res.leg[qpms_mn2y(0, l)] = ((l%2)?ct:1.)*sqrt((2*l+1)/(4*M_PI *l*(l+1)));
      double fl = 0.25 * sqrt((2*l+1)*M_1_PI);
      int lpar = (l%2)?-1:1;
      res.pi [qpms_mn2y(+1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
      res.pi [qpms_mn2y(-1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
      res.tau[qpms_mn2y(+1, l)] = ((ct>0) ? +1 : lpar) * fl * csphase;
      res.tau[qpms_mn2y(-1, l)] = -((ct>0) ? +1 : lpar) * fl * csphase;
    }
  }
  else { // cos(theta) in (-1,1), use normal calculation
    double *legder;
    QPMS_CRASHING_MALLOC(legder, nelem * sizeof(double));
    QPMS_ENSURE_SUCCESS(qpms_legendre_deriv_y_fill(res.leg, legder, ct, lMax,
          GSL_SF_LEGENDRE_SPHARM, csphase));
    // Multiply by the "power normalisation" factor
    for (qpms_l_t l = 1; l <= lMax; ++l) {
      double prefac = 1./sqrt(l*(l+1));
      for (qpms_m_t m = -l; m <= l; ++m) {
        res.leg[qpms_mn2y(m,l)] *= prefac;
        legder[qpms_mn2y(m,l)] *= prefac;
      }
    }
    for (qpms_l_t l = 1; l <= lMax; ++l) {
      for (qpms_m_t m = -l; m <= l; ++m) {
        res.pi [qpms_mn2y(m,l)] = m / st * res.leg[qpms_mn2y(m,l)];
        res.tau[qpms_mn2y(m,l)] = - st * legder[qpms_mn2y(m,l)];
      }
    }
    free(legder);
  }
  res.lMax = lMax;
  return res;
}

void qpms_pitau_free(qpms_pitau_t x) {
  free(x.leg);
  free(x.pi);
  free(x.tau);
}

