#include <math.h>
#include "qpms_types.h"
#include "gaunt.h"
#include "translations.h"
#include "indexing.h" // TODO replace size_t and int with own index types here
#include <stdbool.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include "assert_cython_workaround.h"
#include "kahansum.h"
#include <stdlib.h> //abort()
#include <gsl/gsl_sf_coupling.h>


/*
 * Define macros with additional factors that "should not be there" according
 * to the "original" formulae but are needed to work with my vswfs.
 * (actually, I don't know whether the error is in using "wrong" implementation
 * of vswfs, "wrong" implementation of Xu's translation coefficient formulae,
 * error/inconsintency in Xu's paper or something else)
 * Anyway, the zeroes give the correct _numerical_ values according to Xu's
 * paper tables (without Xu's typos, of course), while 
 * the predefined macros give the correct translations of the VSWFs for the 
 * QPMS_NORMALIZATION_TAYLOR_CS norm.
 */
#if !(defined AN0 || defined AN1 || defined AN2 || defined AN3)
#pragma message "using AN1 macro as default"
#define AN1
#endif
//#if !(defined AM0 || defined AM2)
//#define AM1
//#endif
#if !(defined BN0 || defined BN1 || defined BN2 || defined BN3)
#pragma message "using BN1 macro as default"
#define BN1
#endif
//#if !(defined BM0 || defined BM2)
//#define BM1
//#endif
//#if !(defined BF0 || defined BF1 || defined BF2 || defined BF3)
//#define BF1
//#endif

// if defined, the pointer B_multipliers[y] corresponds to the q = 1 element;
// otherwise, it corresponds to the q = 0 element, which should be identically zero
#ifdef QPMS_PACKED_B_MULTIPLIERS
#define BQ_OFFSET 1
#else
#define BQ_OFFSET 0
#endif


/*
 * References:
 * [Xu_old] Yu-Lin Xu, Journal of Computational Physics 127, 285–298 (1996)
 * [Xu] Yu-Lin Xu, Journal of Computational Physics 139, 137–165 (1998)
 */

/*
 * GENERAL TODO: use normalised Legendre functions for Kristensson and Taylor conventions directly
 * instead of normalising them here (the same applies for csphase).
 */

static const double sqrtpi = 1.7724538509055160272981674833411451827975494561223871;
//static const double ln2 = 0.693147180559945309417232121458176568075500134360255254120;

// Associated Legendre polynomial at zero argument (DLMF 14.5.1)
double qpms_legendre0(int m, int n) {
  return pow(2,m) * sqrtpi / tgamma(.5*n - .5*m + .5) / tgamma(.5*n-.5*m);
}

static inline int min1pow(int x) {
  return (x % 2) ? -1 : 1;
}

static inline complex double ipow(int x) {
  return cpow(I, x);
}

// Derivative of associated Legendre polynomial at zero argument (DLMF 14.5.2)
double qpms_legendreD0(int m, int n) {
  return -2 * qpms_legendre0(m, n);
}


static inline int imin(int x, int y) {
  return x > y ? y : x;
}

// The uppermost value of q index for the B coefficient terms from [Xu](60).
// N.B. this is different from [Xu_old](79) due to the n vs. n+1 difference.
// However, the trailing terms in [Xu_old] are analytically zero (although
// the numerical values will carry some non-zero rounding error).
static inline int gauntB_Q_max(int M, int n, int mu, int nu) {
  return imin(n, imin(nu, (n+nu+1-abs(M+mu))/2));
}

int qpms_sph_bessel_fill(qpms_bessel_t typ, int lmax, double x, complex double *result_array) {
  int retval;
  double tmparr[lmax+1];
  switch(typ) {
    case QPMS_BESSEL_REGULAR:
      retval = gsl_sf_bessel_jl_steed_array(lmax, x, tmparr);
      for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
      return retval;
      break;
    case QPMS_BESSEL_SINGULAR: //FIXME: is this precise enough? Would it be better to do it one-by-one?
      retval = gsl_sf_bessel_yl_array(lmax,x,tmparr);
      for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
      return retval;
      break;
    case QPMS_HANKEL_PLUS:
    case QPMS_HANKEL_MINUS:
      retval = gsl_sf_bessel_jl_steed_array(lmax, x, tmparr);
      for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
      if(retval) return retval;
      retval = gsl_sf_bessel_yl_array(lmax, x, tmparr);
      if (typ==QPMS_HANKEL_PLUS)
        for (int l = 0; l <= lmax; ++l) result_array[l] += I * tmparr[l];
      else 
        for (int l = 0; l <= lmax; ++l) result_array[l] +=-I * tmparr[l];
      return retval;
      break;
    default:
      abort();
      //return GSL_EDOM;
  }
  assert(0);
}

static inline double qpms_trans_normlogfac(qpms_normalisation_t norm,
    int m, int n, int mu, int nu) {
  //int csphase = qpms_normalisation_t csphase(norm); // probably not needed here
  norm = qpms_normalisation_t_normonly(norm);
  switch(norm) {
    case QPMS_NORMALISATION_KRISTENSSON:
    case QPMS_NORMALISATION_TAYLOR:
      return -0.5*(lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1));
      break;
    case QPMS_NORMALISATION_NONE:
      return -(lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1));
      break;
#ifdef USE_XU_ANTINORMALISATION
    case QPMS_NORMALISATION_XU:
      return 0;
      break;
#endif
    default:
      abort();
  }
}

static inline double qpms_trans_normfac(qpms_normalisation_t norm,
    int m, int n, int mu, int nu) {
  int csphase = qpms_normalisation_t_csphase(norm);
  norm = qpms_normalisation_t_normonly(norm);
  /* Account for csphase here. Alternatively, this could be done by
   * using appropriate csphase in the legendre polynomials when calculating
   * the translation operator.
   */
  double normfac = (1 == csphase) ? min1pow(m-mu) : 1.;
  switch(norm) {
    case QPMS_NORMALISATION_KRISTENSSON:
      normfac *= sqrt((n*(n+1.))/(nu*(nu+1.)));
      normfac *= sqrt((2.*n+1)/(2.*nu+1));
      break;
    case QPMS_NORMALISATION_TAYLOR:
      normfac *= sqrt((2.*n+1)/(2.*nu+1));
      break;
    case QPMS_NORMALISATION_NONE:
      normfac *= (2.*n+1)/(2.*nu+1);
      break;
#ifdef USE_XU_ANTINORMALISATION
    case QPMS_NORMALISATION_XU:
      break;
#endif
    default:
      abort();
  }

  return normfac;
}

complex double qpms_trans_single_A(qpms_normalisation_t norm,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  if(r_ge_d) J = QPMS_BESSEL_REGULAR;

  double costheta = cos(kdlj.theta);

  int qmax = gaunt_q_max(-m,n,mu,nu); // nemá tu být +m?
  // N.B. -m !!!!!!
  double a1q[qmax+1];
  int err;
  gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
  double a1q0 = a1q[0];
  if (err) abort();

  int csphase = qpms_normalisation_t_csphase(norm);

  double leg[gsl_sf_legendre_array_n(n+nu)];
  if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu,costheta,csphase,leg)) abort();
  complex double bes[n+nu+1];
  if (qpms_sph_bessel_fill(J, n+nu, kdlj.r, bes)) abort();
  complex double sum = 0;
  for(int q = 0; q <= qmax; ++q) {
    int p = n+nu-2*q;
    int Pp_order = mu-m;
    //if(p < abs(Pp_order)) continue; // FIXME raději nastav lépe meze
    assert(p >= abs(Pp_order));
    double a1q_n = a1q[q] / a1q0;
    double Pp = leg[gsl_sf_legendre_array_index(p, abs(Pp_order))];
    if (Pp_order < 0) Pp *= min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
    complex double zp = bes[p];
    complex double summandq = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n * zp * Pp;
    sum += summandq; // TODO KAHAN
  }

  double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+1) - lgamma(2*(n+nu)+1));
  complex double presum = exp(exponent);
  presum *= cexp(I*(mu-m)*kdlj.phi) * min1pow(m) * ipow(nu+n) / (4*n);

  double normlogfac = qpms_trans_normlogfac(norm,m,n,mu,nu);
  double normfac = qpms_trans_normfac(norm,m,n,mu,nu);

  // ipow(n-nu) is the difference from the Taylor formula!
  presum *= /*ipow(n-nu) * */
    (normfac * exp(normlogfac))
#ifdef AN1
      * ipow(n-nu)
#elif defined AN2
      * min1pow(-n+nu)
#elif defined AN3
      * ipow (nu - n)
#endif
#ifdef AM2
      * min1pow(-m+mu)
#endif //NNU
  ;
  return presum * sum;
}

complex double qpms_trans_single_A_Taylor(int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  if(r_ge_d) J = QPMS_BESSEL_REGULAR;

  double costheta = cos(kdlj.theta);

  int qmax = gaunt_q_max(-m,n,mu,nu); // nemá tu být +m?
  // N.B. -m !!!!!!
  double a1q[qmax+1];
  int err;
  gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
  double a1q0 = a1q[0];
  if (err) abort();

  double leg[gsl_sf_legendre_array_n(n+nu)];
  if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu,costheta,-1,leg)) abort();
  complex double bes[n+nu+1];
  if (qpms_sph_bessel_fill(J, n+nu, kdlj.r, bes)) abort();
  complex double sum = 0;
  for(int q = 0; q <= qmax; ++q) {
    int p = n+nu-2*q;
    int Pp_order = mu-m;
    //if(p < abs(Pp_order)) continue; // FIXME raději nastav lépe meze
    assert(p >= abs(Pp_order));
    double a1q_n = a1q[q] / a1q0;
    double Pp = leg[gsl_sf_legendre_array_index(p, abs(Pp_order))];
    if (Pp_order < 0) Pp *= min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
    complex double zp = bes[p];
    complex double summandq = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n * zp * Pp;
    sum += summandq; // TODO KAHAN
  }

  double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+1) - lgamma(2*(n+nu)+1));
  complex double presum = exp(exponent);
  presum *= cexp(I*(mu-m)*kdlj.phi) * min1pow(m) * ipow(nu+n) / (4*n);

  // N.B. ipow(nu-n) is different from the general formula!
  complex double prenormratio = ipow(nu-n) *  sqrt(((2.*nu+1)/(2.*n+1))* exp(
        lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1)));
  return (presum / prenormratio) * sum;
}

// [Xu_old], eq. (83)
complex double qpms_trans_single_B_Xu(int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) { 
  assert(0); // FIXME probably gives wrong values, do not use.
  if(r_ge_d) J = QPMS_BESSEL_REGULAR;
  double costheta = cos(kdlj.theta);

  // TODO Qmax cleanup: can I replace Qmax with realQmax???
  int q2max = gaunt_q_max(-m-1,n+1,mu+1,nu);
  int Qmax = gaunt_q_max(-m,n+1,mu,nu);
  int realQmax = gauntB_Q_max(-m, n, mu, nu);
  double a2q[q2max+1], a3q[Qmax+1], a2q0, a3q0;
  int err;
  if (mu == nu) {
    for (int q = 0; q <= q2max; ++q) 
      a2q[q] = 0;
    a2q0 = 1;
  }
  else { 
    gaunt_xu(-m-1,n+1,mu+1,nu,q2max,a2q,&err); if (err) abort();
    a2q0 = a2q[0];
  }
  gaunt_xu(-m,n+1,mu,nu,Qmax,a3q,&err); if (err) abort();
  a3q0 = a3q[0];

  double leg[gsl_sf_legendre_array_n(n+nu+1)];
  if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,costheta,-1,leg)) abort();
  complex double bes[n+nu+2];
  if (qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bes)) abort();

  complex double sum = 0;
  for (int q = 0; q <= realQmax; ++q) {
    int p = n+nu-2*q;
    double a2q_n = a2q[q]/a2q0;
    double a3q_n = a3q[q]/a3q0;
    complex double zp_ = bes[p+1];
    int Pp_order_ = mu-m;
    //if(p+1 < abs(Pp_order_)) continue; // FIXME raději nastav lépe meze
    assert(p+1 >= abs(Pp_order_));
    double Pp_ = leg[gsl_sf_legendre_array_index(p+1, abs(Pp_order_))];
    if (Pp_order_ < 0) Pp_ *= min1pow(mu-m) * exp(lgamma(1+1+p+Pp_order_)-lgamma(1+1+p-Pp_order_));
    complex double summandq = ((2*(n+1)*(nu-mu)*a2q_n
          -(-nu*(nu+1) - n*(n+3) - 2*mu*(n+1)+p*(p+3))* a3q_n)
        *min1pow(q) * zp_ * Pp_);
    sum += summandq; // TODO KAHAN
  }

  double exponent=(lgamma(2*n+3)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+2)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+2) - lgamma(2*(n+nu)+3));
  complex double presum = exp(exponent);
  presum *= cexp(I*(mu-m)*kdlj.phi) * min1pow(m) * ipow(nu+n+1) / (
      (4*n)*(n+1)*(n+m+1));

  // Taylor normalisation v2, proven to be equivalent
  complex double prenormratio = ipow(nu-n);

  return (presum / prenormratio) * sum;
}

complex double qpms_trans_single_B(qpms_normalisation_t norm,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) { 
#ifndef USE_BROKEN_SINGLETC
  assert(0); // FIXME probably gives wrong values, do not use.
#endif
  if(r_ge_d) J = QPMS_BESSEL_REGULAR;
  double costheta = cos(kdlj.theta);

  int q2max = gaunt_q_max(-m-1,n+1,mu+1,nu);
  int Qmax = gaunt_q_max(-m,n+1,mu,nu);
  int realQmax = gauntB_Q_max(-m,n,mu,nu);
  double a2q[q2max+1], a3q[Qmax+1], a2q0, a3q0;
  int err;
  if (mu == nu) {
    for (int q = 0; q <= q2max; ++q) 
      a2q[q] = 0;
    a2q0 = 1;
  }
  else { 
    gaunt_xu(-m-1,n+1,mu+1,nu,q2max,a2q,&err); if (err) abort();
    a2q0 = a2q[0];
  }
  gaunt_xu(-m,n+1,mu,nu,Qmax,a3q,&err); if (err) abort();
  a3q0 = a3q[0];
  
  int csphase = qpms_normalisation_t_csphase(norm);

  double leg[gsl_sf_legendre_array_n(n+nu+1)];
  if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,costheta,csphase,leg)) abort();
  complex double bes[n+nu+2];
  if (qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bes)) abort();

  complex double sum = 0;
  for (int q = 0; q <= realQmax; ++q) {
    int p = n+nu-2*q;
    double a2q_n = a2q[q]/a2q0;
    double a3q_n = a3q[q]/a3q0;
    complex double zp_ = bes[p+1];
    int Pp_order_ = mu-m;
    //if(p+1 < abs(Pp_order_)) continue; // FIXME raději nastav lépe meze
    assert(p+1 >= abs(Pp_order_));
    double Pp_ = leg[gsl_sf_legendre_array_index(p+1, abs(Pp_order_))];
    if (Pp_order_ < 0) Pp_ *= min1pow(mu-m) * exp(lgamma(1+1+p+Pp_order_)-lgamma(1+1+p-Pp_order_));
    complex double summandq = ((2*(n+1)*(nu-mu)*a2q_n
          -(-nu*(nu+1) - n*(n+3) - 2*mu*(n+1)+p*(p+3))* a3q_n)
        *min1pow(q) * zp_ * Pp_);
    sum += summandq; //TODO KAHAN
  }

  double exponent=(lgamma(2*n+3)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+2)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+2) - lgamma(2*(n+nu)+3));
  complex double presum = exp(exponent);
  presum *= cexp(I*(mu-m)*kdlj.phi) * min1pow(m) * ipow(nu+n+1) / (
      (4*n)*(n+1)*(n+m+1));

  double normlogfac = qpms_trans_normlogfac(norm,m,n,mu,nu);
  double normfac = qpms_trans_normfac(norm,m,n,mu,nu);

  // ipow(n-nu) is the difference from the "old Taylor" formula	
  presum *= /*ipow(n-nu) * */(exp(normlogfac) * normfac)
#ifdef AN1
      * ipow(n-nu)
#elif defined AN2
      * min1pow(-n+nu)
#elif defined AN3
      * ipow (nu - n)
#endif
#ifdef AM2
      * min1pow(-m+mu)
#endif //NNU
  ;

  return presum * sum;
}

complex double qpms_trans_single_B_Taylor(int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) { 
  assert(0); // FIXME probably gives wrong values, do not use.
  if(r_ge_d) J = QPMS_BESSEL_REGULAR;
  double costheta = cos(kdlj.theta);

  int q2max = gaunt_q_max(-m-1,n+1,mu+1,nu);
  int Qmax = gaunt_q_max(-m,n+1,mu,nu);
  int realQmax = gauntB_Q_max(-m,n,mu,nu);
  double a2q[q2max+1], a3q[Qmax+1], a2q0, a3q0;
  int err;
  if (mu == nu) {
    for (int q = 0; q <= q2max; ++q) 
      a2q[q] = 0;
    a2q0 = 1;
  }
  else { 
    gaunt_xu(-m-1,n+1,mu+1,nu,q2max,a2q,&err); if (err) abort();
    a2q0 = a2q[0];
  }
  gaunt_xu(-m,n+1,mu,nu,Qmax,a3q,&err); if (err) abort();
  a3q0 = a3q[0];

  double leg[gsl_sf_legendre_array_n(n+nu+1)];
  if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,costheta,-1,leg)) abort();
  complex double bes[n+nu+2];
  if (qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bes)) abort();

  complex double sum = 0;
  for (int q = 0; q <= realQmax; ++q) {
    int p = n+nu-2*q;
    double a2q_n = a2q[q]/a2q0;
    double a3q_n = a3q[q]/a3q0;
    complex double zp_ = bes[p+1];
    int Pp_order_ = mu-m;
    //if(p+1 < abs(Pp_order_)) continue; // FIXME raději nastav lépe meze
    assert(p+1 >= abs(Pp_order_));
    double Pp_ = leg[gsl_sf_legendre_array_index(p+1, abs(Pp_order_))];
    if (Pp_order_ < 0) Pp_ *= min1pow(mu-m) * exp(lgamma(1+1+p+Pp_order_)-lgamma(1+1+p-Pp_order_));
    complex double summandq = ((2*(n+1)*(nu-mu)*a2q_n
          -(-nu*(nu+1) - n*(n+3) - 2*mu*(n+1)+p*(p+3))* a3q_n)
        *min1pow(q) * zp_ * Pp_);
    sum += summandq; //TODO KAHAN
  }

  double exponent=(lgamma(2*n+3)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+2)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+2) - lgamma(2*(n+nu)+3));
  complex double presum = exp(exponent);
  presum *= cexp(I*(mu-m)*kdlj.phi) * min1pow(m) * ipow(nu+n+1) / (
      (4*n)*(n+1)*(n+m+1));

  // Taylor normalisation v2, proven to be equivalent
  // ipow(nu-n) is different from the new general formula!!!
  complex double prenormratio = ipow(nu-n) * sqrt(((2.*nu+1)/(2.*n+1))* exp(
        lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1)));

  return (presum / prenormratio) * sum;
}

complex double qpms_trans_single_A_Taylor_ext(int m, int n, int mu, int nu, 
    double kdlj_r, double kdlj_theta, double kdlj_phi, int r_ge_d, int J) {
  sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_single_A_Taylor(m,n,mu,nu,kdlj,r_ge_d,J);
}

complex double qpms_trans_single_B_Taylor_ext(int m, int n, int mu, int nu, 
    double kdlj_r, double kdlj_theta, double kdlj_phi, int r_ge_d, int J) {
  sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_single_B_Taylor(m,n,mu,nu,kdlj,r_ge_d,J);
}

void qpms_trans_calculator_free(qpms_trans_calculator *c) {
  free(c->A_multipliers[0]);
  free(c->A_multipliers);
  free(c->B_multipliers[0]);
  free(c->B_multipliers);
#ifdef LATTICESUMS
  free(c->hct);
  free(c->legendre0);
#endif
  free(c);
}

static inline size_t qpms_trans_calculator_index_mnmunu(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu){
  return c->nelem * qpms_mn2y(m,n) + qpms_mn2y(mu,nu);
}

static inline size_t qpms_trans_calculator_index_yyu(const qpms_trans_calculator *c,
    size_t y, size_t yu) {
  return c->nelem * y + yu;
}


#define SQ(x) ((x)*(x))

static inline int isq(int x) { return x * x; }
static inline double fsq(double x) {return x * x; }

static void qpms_trans_calculator_multipliers_A_general(
    qpms_normalisation_t norm,
    complex double *dest, int m, int n, int mu, int nu, int qmax) {
  assert(qmax == gaunt_q_max(-m,n,mu,nu));
  double a1q[qmax+1];
  int err;
  gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
  if (err) abort();
  double a1q0 = a1q[0];

  double normlogfac = qpms_trans_normlogfac(norm,m,n,mu,nu);
  double normfac = qpms_trans_normfac(norm,m,n,mu,nu);

  normfac *= min1pow(m); //different from old Taylor

  double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+1) - lgamma(2*(n+nu)+1))
    + normlogfac;
  complex double presum = exp(exponent);
  presum *= normfac / (4.*n);
  presum *= ipow(n+nu); // different from old Taylor

  for(int q = 0; q <= qmax; q++) {
    int p = n+nu-2*q;
    int Pp_order = mu - m;
    assert(p >= abs(Pp_order));
    double a1q_n = a1q[q] / a1q0;
    // Assuming non_normalized legendre polynomials (normalisation done here by hand)!
    double Ppfac = (Pp_order >= 0) ? 1 :
      min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
    double summandfac = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n;
    dest[q] = presum * summandfac * Ppfac
#ifdef AN1
      * ipow(n-nu)
#elif defined AN2
      * min1pow(-n+nu)
#elif defined AN3
      * ipow (nu - n)
#endif
#ifdef AM2
      * min1pow(-m+mu)
#endif //NNU
      ;
    // FIXME I might not need complex here
  }
}


// as in [Xu](61)
double cruzan_bfactor(int M, int n, int mu, int nu, int p) {
  double logprefac = lgamma(n+M+1) - lgamma(n-M+1) + lgamma(nu+mu+1) - lgamma(nu-mu+1) 
    + lgamma(p-M-mu+2) - lgamma(p+M+mu+2);
  logprefac *= 0.5;
  return min1pow(mu+M) * (2*p+3) * exp(logprefac) 
    * gsl_sf_coupling_3j(2*n, 2*nu, 2*(p+1), 2*M, 2*mu, 2*(-M-mu))
    * gsl_sf_coupling_3j(2*n, 2*nu, 2*p, 0, 0, 0);
}


void qpms_trans_calculator_multipliers_B_general(
    qpms_normalisation_t norm,
    complex double *dest, int m, int n, int mu, int nu, int Qmax){
  // This is according to the Cruzan-type formula [Xu](59)
  assert(Qmax == gauntB_Q_max(-m,n,mu,nu));



  double normlogfac= qpms_trans_normlogfac(norm,m,n,mu,nu);
  double normfac = qpms_trans_normfac(norm,m,n,mu,nu);

  double presum = min1pow(1-m) * (2*nu+1)/(2.*(n*(n+1))) 
    * exp(lgamma(n+m+1) - lgamma(n-m+1) + lgamma(nu-mu+1) - lgamma(nu+mu+1)
        + normlogfac)
    * normfac;

  for(int q = BQ_OFFSET; q <= Qmax; ++q) {
    int p = n+nu-2*q;
    int Pp_order = mu - m;
    // Assuming non-normalised Legendre polynomials, normalise here by hand.
    // Ppfac_ differs from Ppfac in the A-case by the substitution p->p+1
    double Ppfac_ = (Pp_order >= 0)? 1 :
      min1pow(mu-m) * exp(lgamma(1+1+p+Pp_order)-lgamma(1+1+p-Pp_order));
    double t = sqrt(
        (isq(p+1)-isq(n-nu))
        * (isq(n+nu+1)-isq(p+1))
        );
    dest[q-BQ_OFFSET] = presum * t * Ppfac_ 
      * cruzan_bfactor(-m,n,mu,nu,p) * ipow(p+1) 
#ifdef BN1
      * ipow(n-nu)
#elif defined BN2
      * min1pow(-n+nu)
#elif defined BN3
      * ipow (nu - n)
#endif
#ifdef BM2
      * min1pow(-m+mu)
#endif  
#ifdef BF1
      * I
#elif defined BF2
      * (-1)
#elif defined BF3
      * (-I)
#endif
      ;// NNU
  }
}

/*static*/ void qpms_trans_calculator_multipliers_B_general_oldXu(
    qpms_normalisation_t norm,
    complex double *dest, int m, int n, int mu, int nu, int Qmax) {
  assert(0); // FIXME probably gives wrong values, do not use.
  assert(Qmax == gauntB_Q_max(-m,n,mu,nu));
  int q2max = gaunt_q_max(-m-1,n+1,mu+1,nu);
  // assert(Qmax == q2max);
  // FIXME is it safe to replace q2max with Qmax in gaunt_xu??
  double a2q[q2max+1], a3q[Qmax+1], a2q0, a3q0;
  int err;
  if (mu == nu) {
    for (int q = 0; q <= q2max; ++q)
      a2q[q] = 0;
    a2q0 = 1;
  }
  else {
    gaunt_xu(-m-1,n+1,mu+1,nu,q2max,a2q,&err); if (err) abort();
    a2q0 = a2q[0];
  }
  gaunt_xu(-m,n+1,mu,nu,q2max,a3q,&err); if (err) abort(); // FIXME this should probably go away
  a3q0 = a3q[0];


  int csphase = qpms_normalisation_t_csphase(norm); //TODO FIXME use this
  norm = qpms_normalisation_t_normonly(norm);
  double normlogfac= qpms_trans_normlogfac(norm,m,n,mu,nu);
  double normfac = qpms_trans_normfac(norm,m,n,mu,nu);
  // TODO use csphase to modify normfac here!!!!
  // normfac = xxx ? -normfac : normfac;
  normfac *= min1pow(m);//different from old taylor



  double exponent=(lgamma(2*n+3)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+2)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+2) - lgamma(2*(n+nu)+3))
    +normlogfac;
  complex double presum = exp(exponent);
  presum *= I * ipow(nu+n) /*different from old Taylor */ * normfac / (
      (4*n)*(n+1)*(n+m+1));

  for (int q = BQ_OFFSET; q <= Qmax; ++q) {
    int p = n+nu-2*q;
    double a2q_n = a2q[q]/a2q0;
    double a3q_n = a3q[q]/a3q0;
    int Pp_order_ = mu-m;
    //if(p+1 < abs(Pp_order_)) continue; // FIXME raději nastav lépe meze
    assert(p+1 >= abs(Pp_order_));
    double Ppfac = (Pp_order_ >= 0) ? 1 :		

      min1pow(mu-m) * exp(lgamma(1+1+p+Pp_order_)-lgamma(1+1+p-Pp_order_));
    double summandq = ((2*(n+1)*(nu-mu)*a2q_n
          -(-nu*(nu+1) - n*(n+3) - 2*mu*(n+1)+p*(p+3))* a3q_n)
        *min1pow(q));
    dest[q-BQ_OFFSET] = Ppfac * summandq * presum;
  }
}

//#if 0
static void qpms_trans_calculator_multipliers_A_Taylor(
    complex double *dest, int m, int n, int mu, int nu, int qmax) {
  assert(qmax == gaunt_q_max(-m,n,mu,nu));
  double a1q[qmax+1];
  int err;
  gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
  if (err) abort();
  double a1q0 = a1q[0];

  double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+1) - lgamma(2*(n+nu)+1)) - 0.5*( // ex-prenormratio
      lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1));
  double presum = exp(exponent);
  presum *=  min1pow(m+n) * sqrt((2.*n+1)/(2.*nu+1)) / (4*n);

  for(int q = 0; q <= qmax; q++) {
    int p = n+nu-2*q;
    int Pp_order = mu - m;
    assert(p >= abs(Pp_order));
    double a1q_n = a1q[q] / a1q0;
    // Assuming non_normalized legendre polynomials!
    double Ppfac = (Pp_order >= 0) ? 1 :
      min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
    double summandfac = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n;
    dest[q] = presum * summandfac * Ppfac;
    // FIXME I might not need complex here
  }
}
//#endif
#if 0
static void qpms_trans_calculator_multipliers_A_Taylor(
    complex double *dest, int m, int n, int mu, int nu, int qmax) {
  assert(qmax == gaunt_q_max(-m,n,mu,nu));
  double a1q[qmax+1];
  int err;
  gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
  if (err) abort();
  double a1q0 = a1q[0];
  for(int q = 0; q <= qmax; ++q) {
    int p = n+nu-2*q;
    int Pp_order = mu-m;
    //if(p < abs(Pp_order)) continue; // FIXME raději nastav lépe meze
    assert(p >= abs(Pp_order));
    double a1q_n = a1q[q] / a1q0;
    //double Pp = leg[gsl_sf_legendre_array_index(p, abs(Pp_order))];
    //complex double zp = bes[p];
    dest[q] = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n /* * zp * Pp*/;
    if (Pp_order < 0) dest[q] *= min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
    //sum += summandq;
  }

  double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+1) - lgamma(2*(n+nu)+1));
  complex double presum = exp(exponent);
  presum *=/* cexp(I*(mu-m)*kdlj.phi) * */  min1pow(m) * ipow(nu+n) / (4*n);

  complex double prenormratio = ipow(nu-n) *  sqrt(((2.*nu+1)/(2.*n+1))* exp(
        lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1)));
  //return (presum / prenormratio) * sum;
  for(int q=0;q<=qmax;++q) dest[q] *= presum / prenormratio;
}
#endif	



static void qpms_trans_calculator_multipliers_B_Taylor(
    complex double *dest, int m, int n, int mu, int nu, int Qmax) {
  assert(0); // FIXME probably gives wrong values, do not use.
  assert(Qmax == gauntB_Q_max(-m,n,mu,nu));
  int q2max = gaunt_q_max(-m-1,n+1,mu+1,nu);
  //assert(Qmax == q2max);
  // FIXME remove the q2max variable altogether, as it is probably equal
  // to Qmax
  double a2q[q2max+1], a3q[q2max+1], a2q0, a3q0;
  int err;
  if (mu == nu) {
    for (int q = 0; q <= q2max; ++q)
      a2q[q] = 0;
    a2q0 = 1;
  }
  else {
    gaunt_xu(-m-1,n+1,mu+1,nu,q2max,a2q,&err); if (err) abort();
    a2q0 = a2q[0];
  }
  gaunt_xu(-m,n+1,mu,nu,q2max,a3q,&err); if (err) abort();
  a3q0 = a3q[0];

  double exponent=(lgamma(2*n+3)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
      +lgamma(n+nu+m-mu+2)-lgamma(n-m+1)-lgamma(nu+mu+1)
      +lgamma(n+nu+2) - lgamma(2*(n+nu)+3)) - 0.5 * (
      lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)
      -lgamma(nu+mu+1));
  complex double presum = exp(exponent);
  presum *= I * (min1pow(m+n) *sqrt((2.*n+1)/(2.*nu+1)) / (
        (4*n)*(n+1)*(n+m+1)));

  for (int q = BQ_OFFSET; q <= Qmax; ++q) {
    int p = n+nu-2*q;
    double a2q_n = a2q[q]/a2q0;
    double a3q_n = a3q[q]/a3q0;
    int Pp_order_ = mu-m;
    //if(p+1 < abs(Pp_order_)) continue; // FIXME raději nastav lépe meze
    assert(p+1 >= abs(Pp_order_));
    double Ppfac = (Pp_order_ >= 0) ? 1 :		

      min1pow(mu-m) * exp(lgamma(1+1+p+Pp_order_)-lgamma(1+1+p-Pp_order_));
    double summandq = ((2*(n+1)*(nu-mu)*a2q_n
          -(-nu*(nu+1) - n*(n+3) - 2*mu*(n+1)+p*(p+3))* a3q_n)
        *min1pow(q));
    dest[q-BQ_OFFSET] = Ppfac * summandq * presum;
  }
}

int qpms_trans_calculator_multipliers_A(qpms_normalisation_t norm, complex double *dest, int m, int n, int mu, int nu, int qmax) {
  switch (qpms_normalisation_t_normonly(norm)) {
    case QPMS_NORMALISATION_TAYLOR:
#ifdef USE_SEPARATE_TAYLOR
      qpms_trans_calculator_multipliers_A_Taylor(dest,m,n,mu,nu,qmax);
      return 0;
      break;
#endif
    case QPMS_NORMALISATION_NONE:
#ifdef USE_XU_ANTINORMALISATION
    case QPMS_NORMALISATION_XU:
#endif   
    case QPMS_NORMALISATION_KRISTENSSON:
      qpms_trans_calculator_multipliers_A_general(norm, dest, m, n, mu, nu, qmax);
      return 0;
      break;
    default:
      abort();
  }
  assert(0);
}

int qpms_trans_calculator_multipliers_B(qpms_normalisation_t norm, complex double *dest, int m, int n, int mu, int nu, int Qmax) {
  switch (qpms_normalisation_t_normonly(norm)) {
    case QPMS_NORMALISATION_TAYLOR:
#ifdef USE_SEPARATE_TAYLOR
      qpms_trans_calculator_multipliers_B_Taylor(dest,m,n,mu,nu,Qmax);
      return 0;
      break;
#endif
    case QPMS_NORMALISATION_NONE:
#ifdef USE_XU_ANTINORMALISATION
    case QPMS_NORMALISATION_XU:
#endif   
    case QPMS_NORMALISATION_KRISTENSSON:
      qpms_trans_calculator_multipliers_B_general(norm, dest, m, n, mu, nu, Qmax);
      return 0;
      break;
    default:
      abort();
  }
  assert(0);
}

qpms_trans_calculator
*qpms_trans_calculator_init (int lMax, qpms_normalisation_t normalisation) {
  assert(lMax > 0);
  qpms_trans_calculator *c = malloc(sizeof(qpms_trans_calculator));
  c->lMax = lMax;
  c->nelem = lMax * (lMax+2);
  c->A_multipliers = malloc((1+SQ(c->nelem)) * sizeof(complex double *));
  c->B_multipliers = malloc((1+SQ(c->nelem)) * sizeof(complex double *));
  c->normalisation = normalisation;
  size_t *qmaxes = malloc(SQ(c->nelem) * sizeof(size_t));
  size_t qmaxsum = 0;
  for(size_t y = 0; y < c->nelem; y++)
    for(size_t yu = 0; yu < c->nelem; yu++) {
      int m,n, mu, nu;
      qpms_y2mn_p(y,&m,&n);
      qpms_y2mn_p(yu,&mu,&nu);
      qmaxsum += 1 + (
          qmaxes[qpms_trans_calculator_index_yyu(c,y,yu)] 
          = gaunt_q_max(-m,n,mu,nu));
    }
  c->A_multipliers[0] = malloc(qmaxsum * sizeof(complex double));
  // calculate multiplier beginnings
  for(size_t i = 0; i < SQ(c->nelem); ++i) 
    c->A_multipliers[i+1] = c->A_multipliers[i] + qmaxes[i] + 1;
  // calculate the multipliers
  for(size_t y = 0; y < c->nelem; ++y)
    for(size_t yu = 0; yu < c->nelem; ++yu) {
      size_t i = y * c->nelem + yu;
      int m, n, mu, nu;
      qpms_y2mn_p(y, &m, &n);
      qpms_y2mn_p(yu, &mu, &nu);
      qpms_trans_calculator_multipliers_A(normalisation,
          c->A_multipliers[i], m, n, mu, nu, qmaxes[i]);
    }

  qmaxsum = 0;
  for(size_t y=0; y < c->nelem; y++)
    for(size_t yu = 0; yu < c->nelem; yu++) {
      int m, n, mu, nu;
      qpms_y2mn_p(y,&m,&n);
      qpms_y2mn_p(yu,&mu,&nu);
      qmaxsum += (1 - BQ_OFFSET) + (
          qmaxes[qpms_trans_calculator_index_yyu(c,y,yu)] 
          = gauntB_Q_max(-m,n,mu,nu));
    }
  c->B_multipliers[0] = malloc(qmaxsum * sizeof(complex double));
  // calculate multiplier beginnings
  for(size_t i = 0; i < SQ(c->nelem); ++i) 
    c->B_multipliers[i+1] = c->B_multipliers[i] + qmaxes[i] + (1 - BQ_OFFSET);
  // calculate the multipliers
  for(size_t y = 0; y < c->nelem; ++y)
    for(size_t yu = 0; yu < c->nelem; ++yu) {
      size_t i = y * c->nelem + yu;
      int m, n, mu, nu;
      qpms_y2mn_p(y, &m, &n);
      qpms_y2mn_p(yu, &mu, &nu);
      qpms_trans_calculator_multipliers_B(normalisation,
          c->B_multipliers[i], m, n, mu, nu, qmaxes[i]);
    }

  free(qmaxes);
#ifdef LATTICESUMS
  c->hct = hankelcoefftable_init(2*lMax+1);
  c->legendre0 = malloc(gsl_sf_legendre_array_n(2*lMax+1) * sizeof(double));
  if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,2*lMax+1,
              0,-1,c->legendre0)) abort(); // TODO maybe use some "precise" analytical formula instead?
#endif
  return c;
}

static inline complex double qpms_trans_calculator_get_A_precalcbuf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J,
    const complex double *bessel_buf, const double *legendre_buf) {
  size_t i = qpms_trans_calculator_index_mnmunu(c, m, n, mu, nu);
  size_t qmax = c->A_multipliers[i+1] - c->A_multipliers[i] - 1;
  assert(qmax == gaunt_q_max(-m,n,mu,nu));
  complex double sum, kahanc;
  ckahaninit(&sum, &kahanc);
  for(size_t q = 0; q <= qmax; ++q) {
    int p = n+nu-2*q;
    double Pp = legendre_buf[gsl_sf_legendre_array_index(p, abs(mu-m))];
    complex double zp = bessel_buf[p];
    complex double multiplier = c->A_multipliers[i][q];
    ckahanadd(&sum, &kahanc, Pp * zp *  multiplier);
  }
  complex double eimf =  cexp(I*(mu-m)*kdlj.phi);
  return sum * eimf;
}

complex double qpms_trans_calculator_get_A_buf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J,
    complex double *bessel_buf, double *legendre_buf) {
  // This functions gets preallocated memory for bessel and legendre functions, but computes them itself
  if (r_ge_d) J = QPMS_BESSEL_REGULAR;
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) 
    // TODO warn? 
    return NAN+I*NAN;
  int csphase = qpms_normalisation_t_csphase(c->normalisation);
  switch(qpms_normalisation_t_normonly(c->normalisation)) {
    // TODO use normalised legendre functions for Taylor and Kristensson
    case QPMS_NORMALISATION_TAYLOR:
    case QPMS_NORMALISATION_KRISTENSSON:
    case QPMS_NORMALISATION_NONE:
#ifdef USE_XU_ANTINORMALISATION
    case QPMS_NORMALISATION_XU:   
#endif
      {
        double costheta = cos(kdlj.theta);
        if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu,
              costheta,csphase,legendre_buf)) abort();
        if (qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bessel_buf)) abort();
        return qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
            kdlj,r_ge_d,J,bessel_buf,legendre_buf);
      }
      break;
    default:
      abort();
  }
  assert(0);
}

static inline complex double qpms_trans_calculator_get_B_precalcbuf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J,
    const complex double *bessel_buf, const double *legendre_buf) {
  size_t i = qpms_trans_calculator_index_mnmunu(c, m, n, mu, nu);
  size_t qmax = c->B_multipliers[i+1] - c->B_multipliers[i] - (1 - BQ_OFFSET);
  assert(qmax == gauntB_Q_max(-m,n,mu,nu));
  complex double sum, kahanc;
  ckahaninit(&sum, &kahanc);
  for(int q = BQ_OFFSET; q <= qmax; ++q) {
    int p = n+nu-2*q;
    double Pp_ = legendre_buf[gsl_sf_legendre_array_index(p+1, abs(mu-m))];
    complex double zp_ = bessel_buf[p+1];
    complex double multiplier = c->B_multipliers[i][q-BQ_OFFSET];
    ckahanadd(&sum, &kahanc, Pp_ * zp_ * multiplier);
  }
  complex double eimf =  cexp(I*(mu-m)*kdlj.phi);
  return sum * eimf;
}

complex double qpms_trans_calculator_get_B_buf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J,
    complex double *bessel_buf, double *legendre_buf) {
  // This functions gets preallocated memory for bessel and legendre functions, but computes them itself
  if (r_ge_d) J = QPMS_BESSEL_REGULAR;
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) 
    // TODO warn? 
    return NAN+I*NAN;
  int csphase = qpms_normalisation_t_csphase(c->normalisation);
  switch(qpms_normalisation_t_normonly(c->normalisation)) {
    case QPMS_NORMALISATION_TAYLOR:
    case QPMS_NORMALISATION_KRISTENSSON:
    case QPMS_NORMALISATION_NONE:
#ifdef USE_XU_ANTINORMALISATION
    case QPMS_NORMALISATION_XU:   
#endif
      {
        double costheta = cos(kdlj.theta);
        if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,
              costheta,csphase,legendre_buf)) abort();
        if (qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bessel_buf)) abort();
        return qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
            kdlj,r_ge_d,J,bessel_buf,legendre_buf);
      }
      break;
    default:
      abort();
  }
  assert(0);
}

int qpms_trans_calculator_get_AB_buf_p(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J,
    complex double *bessel_buf, double *legendre_buf) {
  if (r_ge_d) J = QPMS_BESSEL_REGULAR;
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) {
    *Adest = NAN+I*NAN;
    *Bdest = NAN+I*NAN;
    // TODO warn? different return value?
    return 0;
  }
  switch(qpms_normalisation_t_normonly(c->normalisation)) {
    case QPMS_NORMALISATION_TAYLOR:
    case QPMS_NORMALISATION_KRISTENSSON:
    case QPMS_NORMALISATION_NONE:
#ifdef USE_XU_ANTINORMALISATION
    case QPMS_NORMALISATION_XU:   
#endif
      {
        double costheta = cos(kdlj.theta);
        if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,
              costheta,-1,legendre_buf)) abort();
        if (qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bessel_buf)) abort();
        *Adest = qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
            kdlj,r_ge_d,J,bessel_buf,legendre_buf);
        *Bdest = qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
            kdlj,r_ge_d,J,bessel_buf,legendre_buf);
        return 0;
      }
      break;
    default:
      abort();
  }
  assert(0);
}


int qpms_trans_calculator_get_AB_arrays_buf(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride,
    sph_t kdlj, bool r_ge_d, qpms_bessel_t J,
    complex double *bessel_buf, double *legendre_buf) {
  if (r_ge_d) J = QPMS_BESSEL_REGULAR;
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) {
    for (size_t i = 0; i < c->nelem; ++i)
      for (size_t j = 0; j < c->nelem; ++j) {
        *(Adest + i*srcstride + j*deststride) = NAN+I*NAN;
        *(Bdest + i*srcstride + j*deststride) = NAN+I*NAN;
      }
    // TODO warn? different return value?
    return 0;
  }
  switch(qpms_normalisation_t_normonly(c->normalisation)) {
    case QPMS_NORMALISATION_TAYLOR:
    case QPMS_NORMALISATION_POWER:
    case QPMS_NORMALISATION_NONE:
      {
        double costheta = cos(kdlj.theta);
        if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,2*c->lMax+1,
              costheta,-1,legendre_buf)) abort();
        if (qpms_sph_bessel_fill(J, 2*c->lMax+1, kdlj.r, bessel_buf)) abort();
        size_t desti = 0, srci = 0;
        for (int n = 1; n <= c->lMax; ++n) for (int m = -n; m <= n; ++m) {
          for (int nu = 1; nu <= c->lMax; ++nu) for (int mu = -nu; mu <= nu; ++mu) {
            size_t assertindex = qpms_trans_calculator_index_mnmunu(c,m,n,mu,nu);
            assert(assertindex == desti*c->nelem + srci);
            *(Adest + deststride * desti + srcstride * srci) = 
              qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
                  kdlj,r_ge_d,J,bessel_buf,legendre_buf);
            *(Bdest + deststride * desti + srcstride * srci) = 
              qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
                  kdlj,r_ge_d,J,bessel_buf,legendre_buf);
            ++srci;
          }
          ++desti;
          srci = 0;
        }
        return 0;
      }
      break;
    default:
      abort();
  }
  assert(0);
}

complex double qpms_trans_calculator_get_A(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  double leg[gsl_sf_legendre_array_n(n+nu)];
  complex double bes[n+nu+1]; // maximum order is 2n for A coeffs, plus the zeroth.
  return qpms_trans_calculator_get_A_buf(c,m,n,mu,nu,kdlj,r_ge_d,J,
      bes,leg);
}

complex double qpms_trans_calculator_get_B(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  double leg[gsl_sf_legendre_array_n(n+nu+1)];
  complex double bes[n+nu+2]; // maximum order is 2n+1 for B coeffs, plus the zeroth.
  return qpms_trans_calculator_get_B_buf(c,m,n,mu,nu,kdlj,r_ge_d,J,
      bes,leg);
}

int qpms_trans_calculator_get_AB_p(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    int m, int n, int mu, int nu, sph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  double leg[gsl_sf_legendre_array_n(2*c->lMax+1)];
  complex double bes[2*c->lMax+2]; // maximum order is 2n+1 for B coeffs, plus the zeroth.
  return qpms_trans_calculator_get_AB_buf_p(c,Adest, Bdest,m,n,mu,nu,kdlj,r_ge_d,J,
      bes,leg);
}

int qpms_trans_calculator_get_AB_arrays(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride,
    sph_t kdlj, bool r_ge_d, qpms_bessel_t J) {
  double leg[gsl_sf_legendre_array_n(c->lMax+c->lMax+1)];
  complex double bes[2*c->lMax+2]; // maximum order is 2n+1 for B coeffs, plus the zeroth.
  return qpms_trans_calculator_get_AB_arrays_buf(c, 
      Adest, Bdest, deststride, srcstride,
      kdlj, r_ge_d, J, 
      bes, leg);
}


#ifdef LATTICESUMS

int qpms_trans_calculator_get_shortrange_AB_arrays_buf(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride,
    sph_t kdlj, qpms_bessel_t J,
    qpms_l_t lrcutoff, unsigned kappa, double cc, // regularisation params
    complex double *bessel_buf, double *legendre_buf
    ) {
  assert(J == QPMS_HANKEL_PLUS); // support only J == 3 for now
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) {
    for (size_t i = 0; i < c->nelem; ++i)
      for (size_t j = 0; j < c->nelem; ++j) {
        *(Adest + i*srcstride + j*deststride) = NAN+I*NAN;
        *(Bdest + i*srcstride + j*deststride) = NAN+I*NAN;
      }
    // TODO warn? different return value?
    return 0;
  }
  switch(qpms_normalisation_t_normonly(c->normalisation)) {
    case QPMS_NORMALISATION_TAYLOR:
    case QPMS_NORMALISATION_POWER:
    case QPMS_NORMALISATION_NONE: 
      {
        double costheta = cos(kdlj.theta);
        if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,2*c->lMax+1,
              costheta,-1,legendre_buf)) abort();
        // if (qpms_sph_bessel_fill(J, 2*c->lMax+1, kdlj.r, bessel_buf)) abort(); // original
        hankelparts_fill(NULL, bessel_buf, 2*c->lMax+1, lrcutoff, c->hct, kappa, cc, kdlj.r);
        size_t desti = 0, srci = 0;
        for (int n = 1; n <= c->lMax; ++n) for (int m = -n; m <= n; ++m) {
          for (int nu = 1; nu <= c->lMax; ++nu) for (int mu = -nu; mu <= nu; ++mu) {
            size_t assertindex = qpms_trans_calculator_index_mnmunu(c,m,n,mu,nu);
            assert(assertindex == desti*c->nelem + srci);
            *(Adest + deststride * desti + srcstride * srci) = 
              qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
                  kdlj,false,J,bessel_buf,legendre_buf);
            *(Bdest + deststride * desti + srcstride * srci) = 
              qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
                  kdlj,false,J,bessel_buf,legendre_buf);
            ++srci;
          }
          ++desti;
          srci = 0;
        }
        return 0;
      }
      break;
    default:
      abort();
  }
  assert(0);
}

int qpms_trans_calculator_get_shortrange_AB_buf_p(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    int m, int n, int mu, int nu, sph_t kdlj,
    qpms_bessel_t J,
    qpms_l_t lrcutoff, unsigned kappa, double cc, // regularisation params
    complex double *bessel_buf, double *legendre_buf) {
  assert(J == QPMS_HANKEL_PLUS); // support only J == 3 for now
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) {
    *Adest = NAN+I*NAN;
    *Bdest = NAN+I*NAN;
    // TODO warn? different return value?
    return 0;
  }
  switch(qpms_normalisation_t_normonly(c->normalisation)) {
    case QPMS_NORMALISATION_TAYLOR:
    case QPMS_NORMALISATION_KRISTENSSON:
    case QPMS_NORMALISATION_NONE:
      {
        double costheta = cos(kdlj.theta);
        if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,
              costheta,-1,legendre_buf)) abort();
        //if (qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bessel_buf)) abort(); // original
        hankelparts_fill(NULL, bessel_buf, 2*c->lMax+1, lrcutoff, c->hct, kappa, cc, kdlj.r);
        
        *Adest = qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
            kdlj,false,J,bessel_buf,legendre_buf);
        *Bdest = qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
            kdlj,false,J,bessel_buf,legendre_buf);
        return 0;
      }
      break;
    default:
      abort();
  }
  assert(0);
}

// Short-range parts of the translation coefficients
int qpms_trans_calculator_get_shortrange_AB_p(const qpms_trans_calculator *c,
                complex double *Adest, complex double *Bdest,
                qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
                qpms_bessel_t J /* Only J=3 valid for now */,
                qpms_l_t lrcutoff, unsigned kappa, double cc) {
  double leg[gsl_sf_legendre_array_n(2*c->lMax+1)];
  complex double bes[2*c->lMax+2]; // maximum order is 2n+1 for B coeffs, plus the zeroth.
  return qpms_trans_calculator_get_shortrange_AB_buf_p(c,Adest, Bdest,m,n,mu,nu,kdlj,J,
      lrcutoff, kappa, cc,
      bes, leg);
}

int qpms_trans_calculator_get_shortrange_AB_arrays(const qpms_trans_calculator *c,
                complex double *Adest, complex double *Bdest,
                size_t deststride, size_t srcstride,
                sph_t kdlj, qpms_bessel_t J /* Only J=3 valid for now */,
                qpms_l_t lrcutoff, unsigned kappa, double cc) {
  double leg[gsl_sf_legendre_array_n(c->lMax+c->lMax+1)];
  complex double bes[2*c->lMax+2]; // maximum order is 2n+1 for B coeffs, plus the zeroth.
  return qpms_trans_calculator_get_shortrange_AB_arrays_buf(c, 
      Adest, Bdest, deststride, srcstride,
      kdlj, J,
      lrcutoff, kappa, cc,
      bes, leg);
}


// Long-range parts
static inline complex double qpms_trans_calculator_get_2DFT_longrange_A_precalcbuf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, sph_t k_sph /* theta must be M_PI_2 */,
    qpms_bessel_t J /* must be 3 for now */,
    const complex double *lrhankel_recparts_buf) {
  assert(J == QPMS_HANKEL_PLUS);
  //assert(k_sph.theta == M_PI_2); CHECK IN ADVANCE INSTEAD
  //assert(k_sph.r > 0);
  size_t i = qpms_trans_calculator_index_mnmunu(c, m, n, mu, nu);
  size_t qmax = c->A_multipliers[i+1] - c->A_multipliers[i] - 1;
  assert(qmax == gaunt_q_max(-m,n,mu,nu));
  complex double sum, kahanc;
  ckahaninit(&sum, &kahanc);
  for(size_t q = 0; q <= qmax; ++q) {
    int p = n+nu-2*q;
    double Pp = c->legendre0[gsl_sf_legendre_array_index(p, abs(mu-m))];
    complex double zp = trindex_cd(lrhankel_recparts_buf, p)[abs(mu-m)]; // orig: bessel_buf[p];
    if (mu - m < 0) zp *= min1pow(mu-m); // DLMF 10.4.1
    complex double multiplier = c->A_multipliers[i][q];
    ckahanadd(&sum, &kahanc, Pp * zp *  multiplier);
  }
  complex double eimf =  cexp(I*(mu-m)*k_sph.phi);
  return sum * eimf * ipow(mu-m);
}

static inline complex double qpms_trans_calculator_get_2DFT_longrange_B_precalcbuf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, sph_t k_sph /* theta must be M_PI_2 */,
    qpms_bessel_t J /* must be 3 for now */,
    const complex double *lrhankel_recparts_buf) {
  assert(J == QPMS_HANKEL_PLUS);
  size_t i = qpms_trans_calculator_index_mnmunu(c, m, n, mu, nu);
  size_t qmax = c->B_multipliers[i+1] - c->B_multipliers[i] - (1 - BQ_OFFSET);
  assert(qmax == gauntB_Q_max(-m,n,mu,nu));
  complex double sum, kahanc;
  ckahaninit(&sum, &kahanc);
  for(int q = BQ_OFFSET; q <= qmax; ++q) {
    int p = n+nu-2*q;
    double Pp_ = c->legendre0[gsl_sf_legendre_array_index(p+1, abs(mu-m))];
    complex double zp_ = trindex_cd(lrhankel_recparts_buf, p+1)[abs(mu-m)]; // orig: bessel_buf[p+1];
    if (mu - m < 0) zp_ *= min1pow(mu-m); // DLMF 10.4.1 
    complex double multiplier = c->B_multipliers[i][q-BQ_OFFSET];
    ckahanadd(&sum, &kahanc, Pp_ * zp_ * multiplier);
  }
  complex double eimf =  cexp(I*(mu-m)*k_sph.phi);
  return sum * eimf * ipow(mu-m);
}

int qpms_trans_calculator_get_2DFT_longrange_AB_buf_p(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    int m, int n, int mu, int nu, sph_t k_sph,
    qpms_bessel_t J,
    qpms_l_t lrk_cutoff, unsigned kappa, double cv, double k0,
    complex double *lrhankel_recparts_buf) {
  
  assert (J == QPMS_HANKEL_PLUS);
  assert(k_sph.theta == M_PI_2);
  
  switch(qpms_normalisation_t_normonly(c->normalisation)) {
    case QPMS_NORMALISATION_TAYLOR:
    case QPMS_NORMALISATION_KRISTENSSON:
    case QPMS_NORMALISATION_NONE:
#ifdef USE_XU_ANTINORMALISATION
    case QPMS_NORMALISATION_XU:   
#endif
      {
        //double costheta = cos(kdlj.theta);
        //if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,
        //      costheta,-1,legendre_buf)) abort();
        //if (qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bessel_buf)) abort();
        lrhankel_recpart_fill(lrhankel_recparts_buf, 2*c->lMax+1 /* TODO n+nu+1 might be enough */,
            lrk_cutoff, c->hct, kappa, cv, k0, k_sph.r);
        *Adest = qpms_trans_calculator_get_2DFT_longrange_A_precalcbuf(c,m,n,mu,nu,
            k_sph,J,lrhankel_recparts_buf);
        *Bdest = qpms_trans_calculator_get_2DFT_longrange_B_precalcbuf(c,m,n,mu,nu,
            k_sph,J,lrhankel_recparts_buf);
        return 0;
      }
      break;
    default:
      abort();
  }
  assert(0);
}

// Fourier transforms of the long-range parts of the translation coefficients
int qpms_trans_calculator_get_2DFT_longrange_AB_p(const qpms_trans_calculator *c,
                complex double *Adest, complex double *Bdest,
                qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t k_sph,
                qpms_bessel_t J /* Only J=3 valid for now */,
                qpms_l_t lrcutoff, unsigned kappa, double cv, double k0) {
  int maxp = 2*c->lMax+1; // TODO this may not be needed here, n+nu+1 could be enough instead
  complex double lrhankel_recpart[maxp * (maxp+1) / 2];
  return qpms_trans_calculator_get_2DFT_longrange_AB_buf_p(c, Adest, Bdest,m,n,mu,nu,k_sph,
      J, lrcutoff, kappa, cv, k0, lrhankel_recpart);
}

int qpms_trans_calculator_get_2DFT_longrange_AB_arrays_buf(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride,
    sph_t k_sph, qpms_bessel_t J /* must be 3 for now */,
    qpms_l_t lrk_cutoff, unsigned kappa, double cv, double k0,
    complex double *lrhankel_recparts_buf) {
  assert(J == QPMS_HANKEL_PLUS);
  assert(k_sph.theta == M_PI_2);
#if 0
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) {
    for (size_t i = 0; i < c->nelem; ++i)
      for (size_t j = 0; j < c->nelem; ++j) {
        *(Adest + i*srcstride + j*deststride) = NAN+I*NAN;
        *(Bdest + i*srcstride + j*deststride) = NAN+I*NAN;
      }
    // TODO warn? different return value?
    return 0;
  }
#endif 
  switch(qpms_normalisation_t_normonly(c->normalisation)) {
    case QPMS_NORMALISATION_TAYLOR:
    case QPMS_NORMALISATION_POWER:
    case QPMS_NORMALISATION_NONE:
      {
        lrhankel_recpart_fill(lrhankel_recparts_buf, 2*c->lMax+1,
            lrk_cutoff, c->hct, kappa, cv, k0, k_sph.r);
        // if (qpms_sph_bessel_fill(J, 2*c->lMax+1, kdlj.r, bessel_buf)) abort();
        size_t desti = 0, srci = 0;
        for (int n = 1; n <= c->lMax; ++n) for (int m = -n; m <= n; ++m) {
          for (int nu = 1; nu <= c->lMax; ++nu) for (int mu = -nu; mu <= nu; ++mu) {
            size_t assertindex = qpms_trans_calculator_index_mnmunu(c,m,n,mu,nu);
            assert(assertindex == desti*c->nelem + srci);
            *(Adest + deststride * desti + srcstride * srci) = 
              qpms_trans_calculator_get_2DFT_longrange_A_precalcbuf(c,m,n,mu,nu,
                  k_sph,J,lrhankel_recparts_buf);
            *(Bdest + deststride * desti + srcstride * srci) = 
              qpms_trans_calculator_get_2DFT_longrange_B_precalcbuf(c,m,n,mu,nu,
                  k_sph,J,lrhankel_recparts_buf);
            ++srci;
          }
          ++desti;
          srci = 0;
        }
        return 0;
      }
      break;
    default:
      abort();
  }
  assert(0);
}


int qpms_trans_calculator_get_2DFT_longrange_AB_arrays(const qpms_trans_calculator *c,
                complex double *Adest, complex double *Bdest,
                size_t deststride, size_t srcstride,
                sph_t k_sph, qpms_bessel_t J /* Only J=3 valid for now */,
                qpms_l_t lrcutoff, unsigned kappa, double cv, double k0) {
  int maxp = 2*c->lMax+1;
  complex double lrhankel_recpart[maxp * (maxp+1) / 2];
  return qpms_trans_calculator_get_2DFT_longrange_AB_arrays_buf(c,
     Adest, Bdest, deststride, srcstride, k_sph, J, 
     lrcutoff, kappa, cv, k0,
     lrhankel_recpart);
}

#endif // LATTICESUMS


complex double qpms_trans_calculator_get_A_ext(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu,
    double kdlj_r, double kdlj_theta, double kdlj_phi,
    int r_ge_d, int J) {
  sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_calculator_get_A(c,m,n,mu,nu,kdlj,r_ge_d,J);
}

complex double qpms_trans_calculator_get_B_ext(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu,
    double kdlj_r, double kdlj_theta, double kdlj_phi,
    int r_ge_d, int J) {
  sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_calculator_get_B(c,m,n,mu,nu,kdlj,r_ge_d,J);
}

int qpms_trans_calculator_get_AB_p_ext(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    int m, int n, int mu, int nu,
    double kdlj_r, double kdlj_theta, double kdlj_phi,
    int r_ge_d, int J) {
  sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_calculator_get_AB_p(c,Adest,Bdest,m,n,mu,nu,kdlj,r_ge_d,J);
}

int qpms_trans_calculator_get_AB_arrays_ext(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride,
    double kdlj_r, double kdlj_theta, double kdlj_phi,
    int r_ge_d, int J) {
  sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_calculator_get_AB_arrays(c,Adest,Bdest,deststride,srcstride,
      kdlj, r_ge_d, J);
}
#ifdef QPMS_COMPILE_PYTHON_EXTENSIONS
#include <string.h>

#ifdef QPMS_USE_OMP
#include <omp.h>
#endif

int qpms_cython_trans_calculator_get_AB_arrays_loop(
    const qpms_trans_calculator *c, const qpms_bessel_t J, const int resnd,
    const int daxis, const int saxis,
    char *A_data, const npy_intp *A_shape, const npy_intp *A_strides,
    char *B_data, const npy_intp *B_shape, const npy_intp *B_strides,
    const char *r_data, const npy_intp *r_shape, const npy_intp *r_strides,
    const char *theta_data, const npy_intp *theta_shape, const npy_intp *theta_strides,
    const char *phi_data, const npy_intp *phi_shape, const npy_intp *phi_strides,
    const char *r_ge_d_data, const npy_intp *r_ge_d_shape, const npy_intp *r_ge_d_strides){
  assert(daxis != saxis);
  assert(resnd >= 2);
  int longest_axis = 0;
  int longestshape = 1;
  const npy_intp *resultshape = A_shape, *resultstrides = A_strides;
  // TODO put some restrict's everywhere?
  for (int ax = 0; ax < resnd; ++ax){
    assert(A_shape[ax] == B_shape[ax]);
    assert(A_strides[ax] == B_strides[ax]);
    if (daxis == ax || saxis == ax) continue;
    if (A_shape[ax] > longestshape) {
      longest_axis = ax;
      longestshape = 1;
    }
  }
  const npy_intp longlen = resultshape[longest_axis];

  npy_intp innerloop_shape[resnd];
  for (int ax = 0; ax < resnd; ++ax) {
    innerloop_shape[ax] = resultshape[ax];
  }
  /* longest axis will be iterated in the outer (parallelized) loop. 
   * Therefore, longest axis, together with saxis and daxis, 
   * will not be iterated in the inner loop:
   */  
  innerloop_shape[longest_axis] = 1;
  innerloop_shape[daxis] = 1;
  innerloop_shape[saxis] = 1;

  // these are the 'strides' passed to the qpms_trans_calculator_get_AB_arrays_ext
  // function, which expects 'const double *' strides, not 'char *' ones.
  const npy_intp dstride = resultstrides[daxis] / sizeof(complex double);
  const npy_intp sstride = resultstrides[saxis] / sizeof(complex double);

  int errval = 0;
  // TODO here start parallelisation
  //#pragma omp parallel 
  {
    npy_intp local_indices[resnd];
    memset(local_indices, 0, sizeof(local_indices));
    int errval_local = 0;
    size_t longi;
    //#pragma omp for
    for(longi = 0; longi < longlen; ++longi) {
      // this might be done also in the inverse order, but this is more 
      // 'c-contiguous' way of incrementing the indices
      int ax = resnd - 1;
      while(ax >= 0) {
        /* calculate the correct index/pointer for each array used. 
         * This can be further optimized from O(resnd * total size of 
         * the result array) to O(total size of the result array), but 
         * fick that now
         */
        const char *r_p = r_data + r_strides[longest_axis] * longi;
        const char *theta_p = theta_data + theta_strides[longest_axis] * longi;
        const char *phi_p = phi_data + phi_strides[longest_axis] * longi;
        const char *r_ge_d_p = r_ge_d_data + r_ge_d_strides[longest_axis] * longi;
        char *A_p = A_data + A_strides[longest_axis] * longi;
        char *B_p = B_data + B_strides[longest_axis] * longi;
        for(int i = 0; i < resnd; ++i) {
          // following two lines are probably not needed, as innerloop_shape is there 1 anyway
          // so if i == daxis, saxis, or longest_axis, local_indices[i] is zero.
          if (i == longest_axis) continue;
          if (daxis == i || saxis == i) continue;
          r_p += r_strides[i] * local_indices[i];
          theta_p += theta_strides[i] * local_indices[i];
          phi_p += phi_strides[i] * local_indices[i];
          A_p += A_strides[i] * local_indices[i];
          B_p += B_strides[i] * local_indices[i];
        }

        // perform the actual task here
        errval_local |= qpms_trans_calculator_get_AB_arrays_ext(c, (complex double *)A_p, 
            (complex double *)B_p,
            dstride, sstride,
            // FIXME change all the _ext function types to npy_... so that
            // these casts are not needed
            *((double *) r_p), *((double *) theta_p), *((double *)phi_p),
            (int)(*((npy_bool *) r_ge_d_p)), J);
        if (errval_local) abort();

        // increment the last index 'digit' (ax is now resnd-1; we don't have do-while loop in python)
        ++local_indices[ax];
        while(local_indices[ax] == innerloop_shape[ax] && ax >= 0) {
          // overflow to the next digit but stop when reached below the last one
          local_indices[ax] = 0;
          local_indices[--ax]++;
        }
        if (ax >= 0) // did not overflow, get back to the lowest index
          ax = resnd - 1;
      }
    }
    errval |= errval_local;
  }
  // FIXME when parallelizing
  // TODO Here end parallelisation
  return errval;
}


#endif // QPMS_COMPILE_PYTHON_EXTENSIONS

