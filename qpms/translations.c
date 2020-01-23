#include <math.h>
#include "qpms_types.h"
#include "qpms_specfunc.h"
#include "gaunt.h"
#include "translations.h"
#include "indexing.h" // TODO replace size_t and int with own index types here
#include <stdbool.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include "tiny_inlines.h"
#include "assert_cython_workaround.h"
#include "kahansum.h"
#include <gsl/gsl_sf_coupling.h>
#include "qpms_error.h"
#include "normalisation.h"

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


// Translation operators for real sph. harm. based waves are not yet implemented...
static inline void TROPS_ONLY_EIMF_IMPLEMENTED(qpms_normalisation_t norm) {
  if (norm & (QPMS_NORMALISATION_SPHARM_REAL | QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE)) 
    QPMS_NOT_IMPLEMENTED("Translation operators for real or inverse complex spherical harmonics based waves are not implemented.");
}

// Use if only the symmetric form [[A, B], [B, A]] (without additional factors) of translation operator is allowed.
static inline void TROPS_ONLY_AB_SYMMETRIC_NORMS_IMPLEMENTED(qpms_normalisation_t norm) {
  switch (norm & QPMS_NORMALISATION_NORM_BITS) {
    case QPMS_NORMALISATION_NORM_SPHARM:
    case QPMS_NORMALISATION_NORM_POWER:
      break; // OK
    default:
      QPMS_NOT_IMPLEMENTED("Only spherical harmonic and power normalisation supported.");
  }
  if (
      ( !(norm & QPMS_NORMALISATION_N_I) != !(norm & QPMS_NORMALISATION_M_I)  )
      ||
      ( !(norm & QPMS_NORMALISATION_N_MINUS) != !(norm & QPMS_NORMALISATION_M_MINUS) )
     )
    QPMS_NOT_IMPLEMENTED("Only normalisations without a phase factors between M and N waves are supported.");
}


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

static inline double qpms_trans_normlogfac(qpms_normalisation_t norm,
    int m, int n, int mu, int nu) {
      return -0.5*(lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1));
}

static inline double qpms_trans_normfac(qpms_normalisation_t norm,
    int m, int n, int mu, int nu) {
  int csphase = qpms_normalisation_t_csphase(norm);
  /* Account for csphase here. Alternatively, this could be done by
   * using appropriate csphase in the legendre polynomials when calculating
   * the translation operator.
   */
  double normfac = (1 == csphase) ? min1pow(m-mu) : 1.;
  normfac *= sqrt((n*(n+1.))/(nu*(nu+1.)));
  normfac *= sqrt((2.*n+1)/(2.*nu+1));
  return normfac;
}

complex double qpms_trans_single_A(qpms_normalisation_t norm,
    int m, int n, int mu, int nu, csph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  TROPS_ONLY_EIMF_IMPLEMENTED(norm);
  if(r_ge_d) J = QPMS_BESSEL_REGULAR;

  double costheta = cos(kdlj.theta);

  int qmax = gaunt_q_max(-m,n,mu,nu); // nemá tu být +m?
  // N.B. -m !!!!!!
  double a1q[qmax+1];
  int err;
  gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
  QPMS_ENSURE_SUCCESS(err);
  double a1q0 = a1q[0];

  double leg[gsl_sf_legendre_array_n(n+nu)];
  QPMS_ENSURE_SUCCESS(gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu,costheta,-1,leg));
  complex double bes[n+nu+1];
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(J, n+nu, kdlj.r, bes));
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
  /// N<-N type coefficients w.r.t. Kristensson's convention. Csphase has been already taken into acct ^^^.
  normfac *=  qpms_normalisation_factor_N_noCS(norm, nu, mu)
              / qpms_normalisation_factor_N_noCS(norm, n, m);
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
#ifdef AM1
      * ipow(-m+mu)
#endif //NNU
#ifdef AM2
      * min1pow(-m+mu)
#endif //NNU
#ifdef AM3
      * ipow(m-mu)
#endif //NNU
  ;
  return presum * sum;
}


complex double qpms_trans_single_B(qpms_normalisation_t norm,
    int m, int n, int mu, int nu, csph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) { 
  TROPS_ONLY_EIMF_IMPLEMENTED(norm);
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
    gaunt_xu(-m-1,n+1,mu+1,nu,q2max,a2q,&err);
    QPMS_ENSURE_SUCCESS(err);
    a2q0 = a2q[0];
  }
  gaunt_xu(-m,n+1,mu,nu,Qmax,a3q,&err);
  QPMS_ENSURE_SUCCESS(err);
  a3q0 = a3q[0];
  
  double leg[gsl_sf_legendre_array_n(n+nu+1)];
  QPMS_ENSURE_SUCCESS(gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,costheta,-1,leg));
  complex double bes[n+nu+2];
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bes));

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
  /// N<-M type coefficients w.r.t. Kristensson's convention. Csphase has been already taken into acct ^^^.
  normfac *=  qpms_normalisation_factor_M_noCS(norm, nu, mu)
              / qpms_normalisation_factor_N_noCS(norm, n, m);

  // ipow(n-nu) is the difference from the "old Taylor" formula	
  presum *= /*ipow(n-nu) * */(exp(normlogfac) * normfac)
#ifdef AN1
      * ipow(n-nu)
#elif defined AN2
      * min1pow(-n+nu)
#elif defined AN3
      * ipow (nu - n)
#endif
#ifdef AM1
      * ipow(-m+mu)
#endif //NNU
#ifdef AM2
      * min1pow(-m+mu)
#endif //NNU
#ifdef AM3
      * ipow(m-mu)
#endif //NNU
  ;

  return presum * sum;
}

void qpms_trans_calculator_free(qpms_trans_calculator *c) {
  free(c->A_multipliers[0]);
  free(c->A_multipliers);
  free(c->B_multipliers[0]);
  free(c->B_multipliers);
#ifdef LATTICESUMS
  qpms_ewald3_constants_free(e3c);
#endif
  free(c->legendre0);
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

static inline double fsq(double x) {return x * x; }

static void qpms_trans_calculator_multipliers_A(
    qpms_normalisation_t norm,
    complex double *dest, int m, int n, int mu, int nu, int qmax) {
  assert(qmax == gaunt_q_max(-m,n,mu,nu));
  double a1q[qmax+1];
  int err;
  gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
  QPMS_ENSURE_SUCCESS(err);
  double a1q0 = a1q[0];

  double normlogfac = qpms_trans_normlogfac(norm,m,n,mu,nu);
  double normfac = qpms_trans_normfac(norm,m,n,mu,nu);
  /// N<-N type coefficients w.r.t. Kristensson's convention. Csphase has been already taken into acct ^^^.
  normfac *=  qpms_normalisation_factor_N_noCS(norm, nu, mu)
              / qpms_normalisation_factor_N_noCS(norm, n, m);

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
#ifdef AM1
      * ipow(-m+mu)
#endif //NNU
#ifdef AM2
      * min1pow(-m+mu)
#endif //NNU
#ifdef AM3
      * ipow(m-mu)
#endif //NNU
      ;
    // FIXME I might not need complex here
  }
}


// as in [Xu](61)
static double cruzan_bfactor(int M, int n, int mu, int nu, int p) {
  double logprefac = lgamma(n+M+1) - lgamma(n-M+1) + lgamma(nu+mu+1) - lgamma(nu-mu+1) 
    + lgamma(p-M-mu+2) - lgamma(p+M+mu+2);
  logprefac *= 0.5;
  return min1pow(mu+M) * (2*p+3) * exp(logprefac) 
    * gsl_sf_coupling_3j(2*n, 2*nu, 2*(p+1), 2*M, 2*mu, 2*(-M-mu))
    * gsl_sf_coupling_3j(2*n, 2*nu, 2*p, 0, 0, 0);
}


void qpms_trans_calculator_multipliers_B(
    qpms_normalisation_t norm,
    complex double *dest, int m, int n, int mu, int nu, int Qmax){
  // This is according to the Cruzan-type formula [Xu](59)
  assert(Qmax == gauntB_Q_max(-m,n,mu,nu));

  double normlogfac= qpms_trans_normlogfac(norm,m,n,mu,nu);
  double normfac = qpms_trans_normfac(norm,m,n,mu,nu);
  /// N<-M type coefficients w.r.t. Kristensson's convention. Csphase has been already taken into acct ^^^.
  normfac *=  qpms_normalisation_factor_M_noCS(norm, nu, mu)
              / qpms_normalisation_factor_N_noCS(norm, n, m);

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
#ifdef BM1
      * ipow(-m+mu)
#endif  
#ifdef BM2
      * min1pow(-m+mu)
#endif  
#ifdef BM3
      * ipow(m-mu)
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

qpms_trans_calculator
*qpms_trans_calculator_init (const int lMax, const qpms_normalisation_t normalisation) {
  TROPS_ONLY_EIMF_IMPLEMENTED(normalisation);
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
#ifdef LATTICESUMS32
  c->e3c = qpms_ewald3_constants_init(2 * lMax + 1, /*csphase*/ qpms_normalisation_t_csphase(normalisation));
#endif
  c->legendre0 = malloc(gsl_sf_legendre_array_n(2*lMax+1) * sizeof(double));
  QPMS_ENSURE_SUCCESS(gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,2*lMax+1,
              0,-1,c->legendre0)); // TODO maybe use some "precise" analytical formula instead?
  return c;
}

static inline complex double qpms_trans_calculator_get_A_precalcbuf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, double kdlj_phi,
    const complex double *bessel_buf, const double *legendre_buf) {
  TROPS_ONLY_EIMF_IMPLEMENTED(c->normalisation);
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
    ckahanadd(&sum, &kahanc, Pp * zp * multiplier);
  }
  complex double eimf =  cexp(I*(mu-m)*kdlj_phi);
  return sum * eimf;
}

complex double qpms_trans_calculator_get_A_buf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, csph_t kdlj,
    bool r_ge_d, qpms_bessel_t J,
    complex double *bessel_buf, double *legendre_buf) {
  // This functions gets preallocated memory for bessel and legendre functions, but computes them itself
  if (r_ge_d) J = QPMS_BESSEL_REGULAR;
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) 
    // TODO warn? 
    return NAN+I*NAN;
  int csphase = qpms_normalisation_t_csphase(c->normalisation);
  
  double costheta = cos(kdlj.theta);
  QPMS_ENSURE_SUCCESS(gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu,
        costheta,csphase,legendre_buf));
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bessel_buf));
  return qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
      kdlj.phi,bessel_buf,legendre_buf);
}

static inline complex double qpms_trans_calculator_get_B_precalcbuf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, double kdlj_phi,
    const complex double *bessel_buf, const double *legendre_buf) {
  TROPS_ONLY_EIMF_IMPLEMENTED(c->normalisation);
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
  complex double eimf =  cexp(I*(mu-m)*kdlj_phi);
  return sum * eimf;
}

complex double qpms_trans_calculator_get_B_buf(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, csph_t kdlj,
    bool r_ge_d, qpms_bessel_t J,
    complex double *bessel_buf, double *legendre_buf) {
  // This functions gets preallocated memory for bessel and legendre functions, but computes them itself
  if (r_ge_d) J = QPMS_BESSEL_REGULAR;
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) 
    // TODO warn? 
    return NAN+I*NAN;
  int csphase = qpms_normalisation_t_csphase(c->normalisation);
  double costheta = cos(kdlj.theta);
  QPMS_ENSURE_SUCCESS(gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,
        costheta,csphase,legendre_buf));
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bessel_buf));
  return qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
      kdlj.phi,bessel_buf,legendre_buf);
}

int qpms_trans_calculator_get_AB_buf_p(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    int m, int n, int mu, int nu, csph_t kdlj,
    bool r_ge_d, qpms_bessel_t J,
    complex double *bessel_buf, double *legendre_buf) {
  if (r_ge_d) J = QPMS_BESSEL_REGULAR;
  if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) {
    *Adest = NAN+I*NAN;
    *Bdest = NAN+I*NAN;
    // TODO warn? different return value?
    return 0;
  }
  double costheta = cos(kdlj.theta);
  QPMS_ENSURE_SUCCESS(gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,
        costheta,-1,legendre_buf));
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(J, n+nu+1, kdlj.r, bessel_buf));
  *Adest = qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
      kdlj.phi,bessel_buf,legendre_buf);
  *Bdest = qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
      kdlj.phi,bessel_buf,legendre_buf);
  return 0;
}

int qpms_trans_calculator_get_AB_arrays_precalcbuf(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride, double kdlj_phi,
    const complex double *bessel_buf, const double *legendre_buf) {
  size_t desti = 0, srci = 0;
  for (int n = 1; n <= c->lMax; ++n) for (int m = -n; m <= n; ++m) {
    for (int nu = 1; nu <= c->lMax; ++nu) for (int mu = -nu; mu <= nu; ++mu) {
#ifndef NDEBUG
      size_t assertindex = qpms_trans_calculator_index_mnmunu(c,m,n,mu,nu);
#endif
      assert(assertindex == desti*c->nelem + srci);
      *(Adest + deststride * desti + srcstride * srci) = 
        qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
            kdlj_phi, bessel_buf, legendre_buf);
      *(Bdest + deststride * desti + srcstride * srci) = 
        qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
            kdlj_phi,bessel_buf,legendre_buf);
      ++srci;
    }
    ++desti;
    srci = 0;
  }
  return 0;
}

int qpms_trans_calculator_get_AB_arrays_buf(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride,
    csph_t kdlj, bool r_ge_d, qpms_bessel_t J,
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
  {
    double costheta = cos(kdlj.theta);
    QPMS_ENSURE_SUCCESS(gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,2*c->lMax+1,
          costheta,-1,legendre_buf));
    QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(J, 2*c->lMax+1, kdlj.r, bessel_buf));
  }
  return qpms_trans_calculator_get_AB_arrays_precalcbuf(c, Adest, Bdest,
      deststride, srcstride, kdlj.phi, bessel_buf, legendre_buf);
}

complex double qpms_trans_calculator_get_A(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, csph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  double leg[gsl_sf_legendre_array_n(n+nu)];
  complex double bes[n+nu+1]; // maximum order is 2n for A coeffs, plus the zeroth.
  return qpms_trans_calculator_get_A_buf(c,m,n,mu,nu,kdlj,r_ge_d,J,
      bes,leg);
}

complex double qpms_trans_calculator_get_B(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu, csph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  double leg[gsl_sf_legendre_array_n(n+nu+1)];
  complex double bes[n+nu+2]; // maximum order is 2n+1 for B coeffs, plus the zeroth.
  return qpms_trans_calculator_get_B_buf(c,m,n,mu,nu,kdlj,r_ge_d,J,
      bes,leg);
}

int qpms_trans_calculator_get_AB_p(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    int m, int n, int mu, int nu, csph_t kdlj,
    bool r_ge_d, qpms_bessel_t J) {
  double leg[gsl_sf_legendre_array_n(2*c->lMax+1)];
  complex double bes[2*c->lMax+2]; // maximum order is 2n+1 for B coeffs, plus the zeroth.
  return qpms_trans_calculator_get_AB_buf_p(c,Adest, Bdest,m,n,mu,nu,kdlj,r_ge_d,J,
      bes,leg);
}

int qpms_trans_calculator_get_AB_arrays(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride,
    csph_t kdlj, bool r_ge_d, qpms_bessel_t J) {
  double leg[gsl_sf_legendre_array_n(c->lMax+c->lMax+1)];
  complex double bes[2*c->lMax+2]; // maximum order is 2n+1 for B coeffs, plus the zeroth.
  return qpms_trans_calculator_get_AB_arrays_buf(c, 
      Adest, Bdest, deststride, srcstride,
      kdlj, r_ge_d, J, 
      bes, leg);
}

// Convenience functions using VSWF base specs
qpms_errno_t qpms_trans_calculator_get_trans_array(const qpms_trans_calculator *c,
    complex double *target,
    /// Must be destspec->lMax <= c-> lMax && destspec->norm == c->norm.
    const qpms_vswf_set_spec_t *destspec, size_t deststride,
    /// Must be srcspec->lMax <= c-> lMax && srcspec->norm == c->norm.
    const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
    csph_t kdlj, bool r_ge_d, qpms_bessel_t J) 
{
  TROPS_ONLY_AB_SYMMETRIC_NORMS_IMPLEMENTED(c->normalisation);
  assert(c->normalisation == destspec->norm && c->normalisation == srcspec->norm);
  assert(c->lMax >= destspec->lMax && c->lMax >= srcspec->lMax);
  assert(destspec->lMax_L < 0 && srcspec->lMax_L < 0);
  // TODO don't use c->lMax etc. if both destspec->lMax and srcspec->lMax are smaller
  complex double A[c->nelem][c->nelem];
  complex double B[c->nelem][c->nelem];
  qpms_errno_t retval = qpms_trans_calculator_get_AB_arrays(c,
      A[0], B[0], c->nelem, 1,
      kdlj, r_ge_d, J);
  for (size_t desti = 0; desti < destspec->n; ++desti) {
    qpms_y_t desty; qpms_vswf_type_t destt;
    if(QPMS_SUCCESS != qpms_uvswfi2ty(destspec->ilist[desti], &destt, &desty))
        qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,
          "Invalid u. vswf index %llx.", destspec->ilist[desti]);
    for (size_t srci = 0; srci < srcspec->n; ++srci){
      qpms_y_t srcy; qpms_vswf_type_t srct;
      if(QPMS_SUCCESS != qpms_uvswfi2ty(srcspec->ilist[srci], &srct, &srcy))
          qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,
            "Invalid u. vswf index %llx.", srcspec->ilist[srci]);
      target[srci * srcstride + desti * deststride]
        = (srct == destt) ? A[desty][srcy] : B[desty][srcy];
    }
  }
  return retval;
}

qpms_errno_t qpms_trans_calculator_get_trans_array_e32_e(const qpms_trans_calculator *c,
    complex double *target, double *err,
    /// Must be destspec->lMax <= c-> lMax && destspec->norm == c->norm.
    const qpms_vswf_set_spec_t *destspec, size_t deststride,
    /// Must be srcspec->lMax <= c-> lMax && srcspec->norm == c->norm.
    const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
    const double eta, const complex double k,
    cart2_t b1, cart2_t b2,
    const cart2_t beta,
    const cart2_t particle_shift,
    double maxR, double maxK,
    const qpms_ewald_part parts
    )
{
  TROPS_ONLY_AB_SYMMETRIC_NORMS_IMPLEMENTED(c->normalisation);
  QPMS_ENSURE(c->normalisation == destspec->norm && c->normalisation == srcspec->norm,
      "The normalisation conventions must be the same");
  assert(c->lMax >= destspec->lMax && c->lMax >= srcspec->lMax);
  assert(destspec->lMax_L < 0 && srcspec->lMax_L < 0);
  // TODO don't use c->lMax etc. if both destspec->lMax and srcspec->lMax are smaller
  const ptrdiff_t ldAB = c->nelem;
  complex double *A, *B;
  double *Aerr = NULL, *Berr = NULL;
  QPMS_CRASHING_MALLOC(A, c->nelem*c->nelem*sizeof(complex double));
  QPMS_CRASHING_MALLOC(B, c->nelem*c->nelem*sizeof(complex double));
  if(err) {
    QPMS_CRASHING_MALLOC(Aerr, c->nelem*c->nelem*sizeof(double));
    QPMS_CRASHING_MALLOC(Berr, c->nelem*c->nelem*sizeof(double));
  }
  qpms_errno_t retval = qpms_trans_calculator_get_AB_arrays_e32_e(c,
      A, Aerr, B, Berr, ldAB, 1, 
      eta, k, b1, b2, beta, particle_shift, maxR, maxK, parts);
  for (size_t desti = 0; desti < destspec->n; ++desti) {
    qpms_y_t desty; qpms_vswf_type_t destt;
    if(QPMS_SUCCESS != qpms_uvswfi2ty(destspec->ilist[desti], &destt, &desty))
        qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,
          "Invalid u. vswf index %llx.", destspec->ilist[desti]);
    for (size_t srci = 0; srci < srcspec->n; ++srci){
      qpms_y_t srcy; qpms_vswf_type_t srct;
      if(QPMS_SUCCESS != qpms_uvswfi2ty(srcspec->ilist[srci], &srct, &srcy))
          qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,
            "Invalid u. vswf index %llx.", srcspec->ilist[srci]);
      target[srci * srcstride + desti * deststride]
        = (srct == destt) ? A[ldAB*desty + srcy] : B[ldAB*desty + srcy];
      if(err) err[srci * srcstride + desti * deststride]
        = (srct == destt) ? Aerr[ldAB*desty + srcy] : Berr[ldAB*desty + srcy];
    }
  }
  free(A); free(B);
  if (err) { free(Aerr); free(Berr); }
  return retval;
}

qpms_errno_t qpms_trans_calculator_get_trans_array_e32(const qpms_trans_calculator *c,
    complex double *target, double *err,
    /// Must be destspec->lMax <= c-> lMax && destspec->norm == c->norm.
    const qpms_vswf_set_spec_t *destspec, size_t deststride,
    /// Must be srcspec->lMax <= c-> lMax && srcspec->norm == c->norm.
    const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
    const double eta, const complex double k,
    cart2_t b1, cart2_t b2,
    const cart2_t beta,
    const cart2_t particle_shift,
    double maxR, double maxK
    )
{
  return qpms_trans_calculator_get_trans_array_e32_e(c, target, err, destspec, deststride,
      srcspec, srcstride, eta, k, b1, b2, beta, particle_shift, maxR, maxK, QPMS_EWALD_FULL);
}


qpms_errno_t qpms_trans_calculator_get_trans_array_lc3p( 
    const qpms_trans_calculator *c,
    complex double *target,
    /// Must be destspec->lMax <= c-> lMax && destspec->norm == c->norm.
    const qpms_vswf_set_spec_t *destspec, size_t deststride,
    /// Must be srcspec->lMax <= c-> lMax && srcspec->norm == c->norm.
    const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
    complex double k, cart3_t destpos, cart3_t srcpos, qpms_bessel_t J
    /// Workspace has to be at least 2 * c->neleme**2 long
    )
{
  csph_t kdlj = cart2csph(cart3_substract(destpos, srcpos));
  kdlj.r *= k;
  return qpms_trans_calculator_get_trans_array(c, target,
      destspec, deststride, srcspec, srcstride, kdlj, 
      false, J);
}

#ifdef LATTICESUMS31
int qpms_trans_calculator_get_AB_arrays_e31z_both_points_and_shift(const qpms_trans_calculator *c,
    complex double * const Adest, double * const Aerr,
    complex double * const Bdest, double * const Berr,
    const ptrdiff_t deststride, const ptrdiff_t srcstride,
    /* qpms_bessel_t J*/ // assume QPMS_HANKEL_PLUS
    const double eta, const double k, const double unitcell_area,
    const size_t nRpoints, const cart2_t *Rpoints, // n.b. can't contain 0; TODO automatic recognition and skip
    const size_t nKpoints, const cart2_t *Kpoints,
    const double beta,//DIFF21
    const double particle_shift//DIFF21
    )
{

  const qpms_y_t nelem2_sc = qpms_lMax2nelem_sc(c->e3c->lMax);
  //const qpms_y_t nelem = qpms_lMax2nelem(c->lMax);
  const bool doerr = Aerr || Berr;
  const bool do_sigma0 = (particle_shift == 0)//DIFF21((particle_shift.x == 0) && (particle_shift.y == 0)); // FIXME ignoring the case where particle_shift equals to lattice vector

  complex double *sigmas_short = malloc(sizeof(complex double)*nelem2_sc);
  complex double *sigmas_long = malloc(sizeof(complex double)*nelem2_sc);
  complex double *sigmas_total = malloc(sizeof(complex double)*nelem2_sc);
  double *serr_short, *serr_long, *serr_total;
  if(doerr) {
    serr_short = malloc(sizeof(double)*nelem2_sc);
    serr_long = malloc(sizeof(double)*nelem2_sc);
    serr_total = malloc(sizeof(double)*nelem2_sc);
  } else serr_short = serr_long = serr_total = NULL;

  QPMS_ENSURE_SUCCESS(ewald31z_sigma_long_points_and_shift(sigmas_long, serr_long, //DIFF21
      c->e3c, eta, k, unitcell_area, nKpoints, Kpoints, beta, particle_shift));

  QPMS_ENSURE_SUCCESS(ewald31z_sigma_short_points_and_shift(sigmas_short, serr_short, //DIFF21
      c->e3c, eta, k, nRpoints, Rpoints, beta, particle_shift));

  for(qpms_y_t y = 0; y < nelem2_sc; ++y)
    sigmas_total[y] = sigmas_short[y] + sigmas_long[y];
  if (doerr) for(qpms_y_t y = 0; y < nelem2_sc; ++y)
    serr_total[y]  = serr_short[y] + serr_long[y];

  complex double sigma0 = 0; double sigma0_err = 0;
  if (do_sigma0) {
    QPMS_ENSURE_SUCCESS(ewald31z_sigma0(&sigma0, &sigma0_err, c->e3c, eta, k));
    const qpms_l_t y = qpms_mn2y_sc(0,0);
    sigmas_total[y] += sigma0;
    if(doerr) serr_total[y] += sigma0_err;
  }
  
      {
        ptrdiff_t desti = 0, srci = 0;
        for (qpms_l_t n = 1; n <= c->lMax; ++n) for (qpms_m_t m = -n; m <= n; ++m) {
          for (qpms_l_t nu = 1; nu <= c->lMax; ++nu) for (qpms_m_t mu = -nu; mu <= nu; ++mu){
            const size_t i = qpms_trans_calculator_index_mnmunu(c, m, n, mu, nu);
            const size_t qmax = c->A_multipliers[i+1] - c->A_multipliers[i] - 1;
            complex double Asum, Asumc; ckahaninit(&Asum, &Asumc);
            double Asumerr, Asumerrc; if(Aerr) kahaninit(&Asumerr, &Asumerrc);
            
            const qpms_m_t mu_m  = mu - m;
            // TODO skip if ... (N.B. skip will be different for 31z and 32)
            for(qpms_l_t q = 0; q <= qmax; ++q) {
              const qpms_l_t p = n + nu - 2*q;
              const qpms_y_t y_sc = qpms_mn2y_sc(mu_m, p); 
              const complex double multiplier = c->A_multipliers[i][q];
              complex double sigma = sigmas_total[y_sc];
              ckahanadd(&Asum, &Asumc, multiplier * sigma);
              if (Aerr) kahanadd(&Asumerr, &Asumerrc, multiplier * serr_total[y_sc]);
            }

            *(Adest + deststride * desti + srcstride * srci) = Asum;
            if (Aerr) *(Aerr + deststride * desti + srcstride * srci) = Asumerr;
            
            // TODO skip if ...
            complex double Bsum, Bsumc; ckahaninit(&Bsum, &Bsumc);
            double Bsumerr, Bsumerrc; if(Berr) kahaninit(&Bsumerr, &Bsumerrc);
             for(qpms_l_t q = 0; q <= qmax; ++q) {
              const qpms_l_t p_ = n + nu - 2*q + 1;
              const qpms_y_t y_sc = qpms_mn2y_sc(mu_m, p_); 
              const complex double multiplier = c->B_multipliers[i][q-BQ_OFFSET];
              complex double sigma = sigmas_total[y_sc];
              ckahanadd(&Bsum, &Bsumc, multiplier * sigma);
              if (Berr) kahanadd(&Bsumerr, &Bsumerrc, multiplier * serr_total[y_sc]);
            }

            *(Bdest + deststride * desti + srcstride * srci) = Bsum;
            if (Berr) *(Berr + deststride * desti + srcstride * srci) = Bsumerr;

            ++srci;
          }
          ++desti;
          srci = 0;
        }
      }

  free(sigmas_short);
  free(sigmas_long);
  free(sigmas_total);
  if(doerr) {
    free(serr_short);
    free(serr_long);
    free(serr_total);
  }
  return 0;
}
#endif // LATTICESUMS_31


#ifdef LATTICESUMS32

// N.B. alternative point generation strategy toggled by macro GEN_RSHIFTEDPOINTS
// and GEN_KSHIFTEDPOINTS.
// The results should be the same. The performance can slightly differ (especially
// if some optimizations in the point generators are implemented.)
int qpms_trans_calculator_get_AB_arrays_e32_e(const qpms_trans_calculator *c,
    complex double * const Adest, double * const Aerr,
    complex double * const Bdest, double * const Berr,
    const ptrdiff_t deststride, const ptrdiff_t srcstride,
    /* qpms_bessel_t J*/ // assume QPMS_HANKEL_PLUS
    const double eta, const complex double k,
    const cart2_t b1, const cart2_t b2,
    const cart2_t beta,
    const cart2_t particle_shift,
    double maxR, double maxK,
    const qpms_ewald_part parts
    )
{

  const qpms_y_t nelem2_sc = qpms_lMax2nelem_sc(c->e3c->lMax);
  //const qpms_y_t nelem = qpms_lMax2nelem(c->lMax);
  const bool doerr = Aerr || Berr;
  const bool do_sigma0 = ((particle_shift.x == 0) && (particle_shift.y == 0)); // FIXME ignoring the case where particle_shift equals to lattice vector

  complex double *sigmas_short = malloc(sizeof(complex double)*nelem2_sc);
  complex double *sigmas_long = malloc(sizeof(complex double)*nelem2_sc);
  complex double *sigmas_total = malloc(sizeof(complex double)*nelem2_sc);
  double *serr_short, *serr_long, *serr_total;
  if(doerr) {
    serr_short = malloc(sizeof(double)*nelem2_sc);
    serr_long = malloc(sizeof(double)*nelem2_sc);
    serr_total = malloc(sizeof(double)*nelem2_sc);
  } else serr_short = serr_long = serr_total = NULL;

  const double unitcell_area = l2d_unitcell_area(b1, b2);
  cart2_t rb1, rb2; // reciprocal basis
  QPMS_ENSURE_SUCCESS(l2d_reciprocalBasis2pi(b1, b2, &rb1, &rb2));

  if (parts & QPMS_EWALD_LONG_RANGE) {
    PGen Kgen = PGen_xyWeb_new(rb1, rb2, BASIS_RTOL,
#ifdef GEN_KSHIFTEDPOINTS
        beta,
#else
        CART2_ZERO,
#endif
        0, true, maxK, false);

    QPMS_ENSURE_SUCCESS(ewald3_sigma_long(sigmas_long, serr_long, c->e3c, eta, k, 
          unitcell_area, LAT_2D_IN_3D_XYONLY, &Kgen,
#ifdef GEN_KSHIFTEDPOINTS
          true,
#else
          false,
#endif
          cart22cart3xy(beta), cart22cart3xy(particle_shift)));
    if(Kgen.stateData) // PGen not consumed entirely (converged earlier)
      PGen_destroy(&Kgen);
  }

  if (parts & QPMS_EWALD_SHORT_RANGE) {
    PGen Rgen = PGen_xyWeb_new(b1, b2, BASIS_RTOL, 
#ifdef GEN_RSHIFTEDPOINTS
        cart2_scale(-1 /*CHECKSIGN*/, particle_shift),
#else
        CART2_ZERO,
#endif
        0, !do_sigma0, maxR, false);

    QPMS_ENSURE_SUCCESS(ewald3_sigma_short(sigmas_short, serr_short, c->e3c, eta, k,
          LAT_2D_IN_3D_XYONLY, &Rgen, 
#ifdef GEN_RSHIFTEDPOINTS
          true,
#else
          false,
#endif
          cart22cart3xy(beta), cart22cart3xy(particle_shift)));

    if(Rgen.stateData) // PGen not consumed entirely (converged earlier)
      PGen_destroy(&Rgen);
  }

  for(qpms_y_t y = 0; y < nelem2_sc; ++y)
    sigmas_total[y] = ((parts & QPMS_EWALD_SHORT_RANGE) ? sigmas_short[y] : 0)
                    + ((parts & QPMS_EWALD_LONG_RANGE)  ? sigmas_long[y] : 0);
  if (doerr) for(qpms_y_t y = 0; y < nelem2_sc; ++y)
    serr_total[y]  = ((parts & QPMS_EWALD_SHORT_RANGE) ? serr_short[y] : 0)
                    + ((parts & QPMS_EWALD_LONG_RANGE)  ? serr_long[y] : 0);

  complex double sigma0 = 0; double sigma0_err = 0;
  if (do_sigma0 && (parts & QPMS_EWALD_0TERM)) {
    QPMS_ENSURE_SUCCESS(ewald3_sigma0(&sigma0, &sigma0_err, c->e3c, eta, k));
    const qpms_l_t y = qpms_mn2y_sc(0,0);
    sigmas_total[y] += sigma0;
    if(doerr) serr_total[y] += sigma0_err;
  }
  
      {
        ptrdiff_t desti = 0, srci = 0;
        for (qpms_l_t n = 1; n <= c->lMax; ++n) for (qpms_m_t m = -n; m <= n; ++m) {
          for (qpms_l_t nu = 1; nu <= c->lMax; ++nu) for (qpms_m_t mu = -nu; mu <= nu; ++mu){
            const size_t i = qpms_trans_calculator_index_mnmunu(c, m, n, mu, nu);
            const size_t qmax = c->A_multipliers[i+1] - c->A_multipliers[i] - 1;
            complex double Asum, Asumc; ckahaninit(&Asum, &Asumc);
            double Asumerr, Asumerrc; if(Aerr) kahaninit(&Asumerr, &Asumerrc);
            
            const qpms_m_t mu_m  = mu - m;
            // TODO skip if ...
            for(qpms_l_t q = 0; q <= qmax; ++q) {
              const qpms_l_t p = n + nu - 2*q;
              const qpms_y_t y_sc = qpms_mn2y_sc(mu_m, p); 
              const complex double multiplier = c->A_multipliers[i][q];
              complex double sigma = sigmas_total[y_sc];
              ckahanadd(&Asum, &Asumc, multiplier * sigma);
              if (Aerr) kahanadd(&Asumerr, &Asumerrc, multiplier * serr_total[y_sc]);
            }

            *(Adest + deststride * desti + srcstride * srci) = Asum;
            if (Aerr) *(Aerr + deststride * desti + srcstride * srci) = Asumerr;
            
            // TODO skip if ...
            complex double Bsum, Bsumc; ckahaninit(&Bsum, &Bsumc);
            double Bsumerr, Bsumerrc; if(Berr) kahaninit(&Bsumerr, &Bsumerrc);
             for(qpms_l_t q = 0; q <= qmax; ++q) {
              const qpms_l_t p_ = n + nu - 2*q + 1;
              const qpms_y_t y_sc = qpms_mn2y_sc(mu_m, p_); 
              const complex double multiplier = c->B_multipliers[i][q-BQ_OFFSET];
              complex double sigma = sigmas_total[y_sc];
              ckahanadd(&Bsum, &Bsumc, multiplier * sigma);
              if (Berr) kahanadd(&Bsumerr, &Bsumerrc, multiplier * serr_total[y_sc]);
            }

            *(Bdest + deststride * desti + srcstride * srci) = Bsum;
            if (Berr) *(Berr + deststride * desti + srcstride * srci) = Bsumerr;

            ++srci;
          }
          ++desti;
          srci = 0;
        }
      }

  free(sigmas_short);
  free(sigmas_long);
  free(sigmas_total);
  if(doerr) {
    free(serr_short);
    free(serr_long);
    free(serr_total);
  }
  return 0;
}

int qpms_trans_calculator_get_AB_arrays_e32(const qpms_trans_calculator *c,
    complex double * const Adest, double * const Aerr,
    complex double * const Bdest, double * const Berr,
    const ptrdiff_t deststride, const ptrdiff_t srcstride,
    /* qpms_bessel_t J*/ // assume QPMS_HANKEL_PLUS
    const double eta, const complex double k,
    const cart2_t b1, const cart2_t b2,
    const cart2_t beta,
    const cart2_t particle_shift,
    double maxR, double maxK) 
{
  return qpms_trans_calculator_get_AB_arrays_e32_e(
      c, Adest, Aerr, Bdest, Berr, deststride, srcstride, 
      eta, k, b1, b2, beta, particle_shift, maxR, maxK, QPMS_EWALD_FULL);
}

#endif // LATTICESUMS32


complex double qpms_trans_calculator_get_A_ext(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu,
    complex double kdlj_r, double kdlj_theta, double kdlj_phi,
    int r_ge_d, int J) {
  csph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_calculator_get_A(c,m,n,mu,nu,kdlj,r_ge_d,J);
}

complex double qpms_trans_calculator_get_B_ext(const qpms_trans_calculator *c,
    int m, int n, int mu, int nu,
    complex double kdlj_r, double kdlj_theta, double kdlj_phi,
    int r_ge_d, int J) {
  csph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_calculator_get_B(c,m,n,mu,nu,kdlj,r_ge_d,J);
}

int qpms_trans_calculator_get_AB_p_ext(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    int m, int n, int mu, int nu,
    complex double kdlj_r, double kdlj_theta, double kdlj_phi,
    int r_ge_d, int J) {
  csph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_calculator_get_AB_p(c,Adest,Bdest,m,n,mu,nu,kdlj,r_ge_d,J);
}

int qpms_trans_calculator_get_AB_arrays_ext(const qpms_trans_calculator *c,
    complex double *Adest, complex double *Bdest,
    size_t deststride, size_t srcstride,
    complex double kdlj_r, double kdlj_theta, double kdlj_phi,
    int r_ge_d, int J) {
  csph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
  return qpms_trans_calculator_get_AB_arrays(c,Adest,Bdest,deststride,srcstride,
      kdlj, r_ge_d, J);
}

