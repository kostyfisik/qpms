#include "ewald.h"
#include <stdlib.h>
#include "indexing.h"
#include "kahansum.h"
#include <assert.h>
#include <string.h>
#include <complex.h>
#include "tiny_inlines.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

// parameters for the quadrature of integral in (4.6)
#ifndef INTEGRATION_WORKSPACE_LIMIT
#define INTEGRATION_WORKSPACE_LIMIT 30000
#endif

#ifndef INTEGRATION_EPSABS
#define INTEGRATION_EPSABS 0
#endif

#ifndef INTEGRATION_EPSREL
#define INTEGRATION_EPSREL 1e-13
#endif

#ifndef M_SQRTPI
#define M_SQRTPI 1.7724538509055160272981674833411452
#endif



// sloppy implementation of factorial
static inline double factorial(const int n) {
  assert(n >= 0);
  if (n < 0)
    return 0; // should not happen in the functions below. (Therefore the assert above)
  else if (n <= 20) {
    double fac = 1;
    for (int i = 1; i <= n; ++i)
      fac *= i;
    return fac;
  }
  else 
    return tgamma(n + 1); // hope it's precise and that overflow does not happen
}

static inline complex double csq(complex double x) { return x * x; }
static inline double sq(double x) { return x * x; }


typedef enum {
  EWALD32_CONSTANTS_ORIG, // As in [1, (4,5)], NOT USED right now.
  EWALD32_CONSTANTS_AGNOSTIC /* Not depending on the spherical harmonic sign/normalisation
                              * convention – the $e^{im\alpha_pq}$ term in [1,(4.5)] being
                              * replaced by the respective $Y_n^m(\pi/2,\alpha)$ 
                              * spherical harmonic. See notes/ewald.lyx.
                              */

} ewald32_constants_option;

static const ewald32_constants_option type = EWALD32_CONSTANTS_AGNOSTIC;

qpms_ewald32_constants_t *qpms_ewald32_constants_init(const qpms_l_t lMax /*, const ewald32_constants_option type */,
    const int csphase)
{
  qpms_ewald32_constants_t *c = malloc(sizeof(qpms_ewald32_constants_t));
  //if (c == NULL) return NULL; // Do I really want to do this?
  c->lMax = lMax;
  c->nelem = qpms_lMax2nelem(lMax);
  c->s1_jMaxes = malloc(c->nelem * sizeof(qpms_l_t));
  c->s1_constfacs = malloc(c->nelem * sizeof(complex double *));
  //if (c->s1_jMaxes == NULL) return NULL;

  //determine sizes
  size_t s1_constfacs_sz = 0;
  for (qpms_y_t y = 0; y < c->nelem; ++y) {
    qpms_l_t n; qpms_m_t m; qpms_y2mn_p(y, &m, &n);
    if ((m + n) % 2 == 0) 
      s1_constfacs_sz += 1 + (c->s1_jMaxes[y] = (n-abs(m))/2);
    else
      c->s1_jMaxes[y] = -1;
  }

  c->s1_constfacs[0]; //WTF???
  c->s1_constfacs_base = malloc(c->nelem * sizeof(complex double));
  size_t s1_constfacs_sz_cumsum = 0;
  for (qpms_y_t y = 0; y < c->nelem; ++y) {
    qpms_l_t n; qpms_m_t m; qpms_y2mn_p(y, &m, &n);
    if ((m + n) % 2 == 0) {
      c->s1_constfacs[y] = c->s1_constfacs_base + s1_constfacs_sz_cumsum;
      // and here comes the actual calculation
      for (qpms_l_t j = 0; j <= c->s1_jMaxes[y]; ++j){
        switch(type) {
          case EWALD32_CONSTANTS_ORIG: // NOT USED
            c->s1_constfacs[y][j] = -0.5 * ipow(n+1) * min1pow((n+m)/2) 
              * sqrt((2*n + 1) * factorial(n-m) * factorial(n+m))
              * min1pow(j) * pow(0.5, n-2*j)
              / (factorial(j) * factorial((n-m)/2-j) * factorial((n+m)/2-j))
              * pow(0.5, 2*j-1);
            break;
          case EWALD32_CONSTANTS_AGNOSTIC:
            c->s1_constfacs[y][j] = -2 * ipow(n+1) * M_SQRTPI   
              * factorial((n-m)/2) * factorial((n+m)/2)
              * min1pow(j) 
              / (factorial(j) * factorial((n-m)/2-j) * factorial((n+m)/2-j));
            break;
          default:
            abort();
        }
      }
      s1_constfacs_sz_cumsum += 1 + c->s1_jMaxes[y];
    }
    else
      c->s1_constfacs[y] = NULL;
  }

  c->legendre0_csphase = csphase;
  c->legendre0 = malloc(gsl_sf_legendre_array_n(lMax) * sizeof(double));
  // N.B. here I use the GSL_SF_LEGENRE_NONE, in order to be consistent with translations.c
  // Moreover, using this approach (i.e. gsl) takes about 64kB extra memory
  if(GSL_SUCCESS != gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE, lMax, 0, csphase, c->legendre0))
    abort();

  return c;
}

void qpms_ewald32_constants_free(qpms_ewald32_constants_t *c) {
  free(c->legendre0);
  free(c->s1_constfacs);
  free(c->s1_constfacs_base);
  free(c->s1_jMaxes);
  free(c);
}



int ewald32_sigma0(complex double *result, double *err,
    const qpms_ewald32_constants_t *c,
    const double eta, const double k)
{
  qpms_csf_result gam;
  int retval = complex_gamma_inc_e(-0.5, -k/sq(2*eta), &gam);
  if (0 != retval)
    abort();
  *result = gam.val * c->legendre0[gsl_sf_legendre_array_index(0,0)] / 2 / M_SQRTPI;
  if(err) 
    *err = gam.err * c->legendre0[gsl_sf_legendre_array_index(0,0)] / 2 / M_SQRTPI;
  return 0;
}



int ewald32_sigma_long_shiftedpoints (
    complex double *target, // must be c->nelem long
    double *err,
    const qpms_ewald32_constants_t *c,
    const double eta, const double k, const double unitcell_area,
    const size_t npoints, const point2d *Kpoints_plus_beta,
    const point2d particle_shift // target - src
    ) 
{
  const qpms_y_t nelem = c->nelem;
  const qpms_l_t lMax = c->lMax;
  
  // Manual init of the ewald summation targets
  complex double *target_c = calloc(nelem, sizeof(complex double));
  memset(target, 0, nelem * sizeof(complex double));
  double *err_c = NULL;
  if (err) {
    err_c = calloc(nelem, sizeof(double));
    memset(err, 0, nelem * sizeof(double));
  }

  const double commonfac = 1/(k*k*unitcell_area); // used in the very end (CFC)
  assert(commonfac > 0);

  // space for Gamma_pq[j]'s
  qpms_csf_result Gamma_pq[lMax/2];

  // CHOOSE POINT BEGIN
  for (size_t i = 0; i < npoints; ++i) { // BEGIN POINT LOOP
    point2d beta_pq = Kpoints_plus_beta[i];
    double rbeta_pq = cart2norm(beta_pq);
  // CHOOSE POINT END

    complex double phasefac = cexp(I*cart2_dot(beta_pq,particle_shift)); // POINT-DEPENDENT (PFC) // !!!CHECKSIGN!!!
    double arg_pq = atan2(beta_pq.y, beta_pq.x); // POINT-DEPENDENT

    // R-DEPENDENT BEGIN
    complex double gamma_pq = lilgamma(rbeta_pq/k);
    complex double z = csq(gamma_pq*k/(2*eta)); // Když o tom tak přemýšlím, tak tohle je vlastně vždy reálné
    for(qpms_l_t j = 0; j < lMax/2; ++j) {
      int retval = complex_gamma_inc_e(0.5-j, z, Gamma_pq+j);
      if(retval) abort();
    }
    // R-DEPENDENT END
    
    // TODO optimisations: all the j-dependent powers can be done for each j only once, stored in array
    // and just fetched for each n, m pair
    for(qpms_l_t n = 1; n <= lMax; ++n)
      for(qpms_m_t m = -n; m <= n; ++m) {
        if((m+n) % 2 != 0) // odd coefficients are zero.
          continue;
        qpms_y_t y = qpms_mn2y(m, n);
        complex double e_imalpha_pq = cexp(I*m*arg_pq);
        complex double jsum, jsum_c; ckahaninit(&jsum, &jsum_c);
        double jsum_err, jsum_err_c; kahaninit(&jsum_err, &jsum_err_c); // TODO do I really need to kahan sum errors?
        for(qpms_l_t j = 0; j < (n-abs(m))/2; ++j) {
          complex double summand = pow(rbeta_pq/k, n-2*j) 
            * e_imalpha_pq  * c->legendre0[gsl_sf_legendre_array_index(n,abs(m))] // This line can actually go outside j-loop
            * cpow(gamma_pq, 2*j-1) // * Gamma_pq[j] bellow (GGG) after error computation
            * c->s1_constfacs[y][j];
          if(err) {
            // FIXME include also other errors than Gamma_pq's relative error
             kahanadd(&jsum_err, &jsum_err_c, Gamma_pq[j].err * cabs(summand));
          }
          summand *= Gamma_pq[j].val; // GGG
          ckahanadd(&jsum, &jsum_c, summand);
        }
        jsum *= phasefac; // PFC
        ckahanadd(target + y, target_c + y, jsum);
        if(err) kahanadd(err + y, err_c + y, jsum_err);
      }
  } // END POINT LOOP
 
  free(err_c);
  free(target_c);
  for(qpms_y_t y = 0; y < nelem; ++y) // CFC common factor from above
    target[y] *= commonfac;
  if(err)
    for(qpms_y_t y = 0; y < nelem; ++y)
      err[y] *= commonfac;
  return 0;
}


struct sigma2_integrand_params {
  int n;
  double k, R;
};

static double sigma2_integrand(double ksi, void *params) {
  struct sigma2_integrand_params *p = (struct sigma2_integrand_params *) params;
  return exp(-sq(p->R * ksi) + sq(p->k / ksi / 2)) * pow(ksi, 2*p->n);
}

static int ewald32_sr_integral(double r, double k, int n, double eta,
    double *result, double *err, gsl_integration_workspace *workspace)
{
  struct sigma2_integrand_params p;


  const double R0 = r; // Maybe could be chosen otherwise, but fuck it for now.
  p.n = n;
  eta *= R0;
  p.k = k * R0;
  p.R = r / R0;  // i.e. p.R = 1

  gsl_function F;
  F.function = sigma2_integrand;
  F.params = &p;

  int retval = gsl_integration_qagiu(&F, eta, INTEGRATION_EPSABS, 
      INTEGRATION_EPSREL, INTEGRATION_WORKSPACE_LIMIT,
      workspace, result, err);
  double normfac = pow(R0, -2*p.n - 1);
  *result *= normfac;
  *err *= normfac;
  return retval;
}

int ewald32_sigma_short_shiftedpoints(
    complex double *target, // must be c->nelem long
    double *err,
    const qpms_ewald32_constants_t *c, // N.B. not too useful here
    const double eta, const double k,
    const size_t npoints, const point2d *Rpoints_plus_particle_shift,
    const point2d beta,
    const point2d particle_shift           // used only in the very end to multiply it by the phase
    )
{
  const qpms_y_t nelem = c->nelem;
  const qpms_l_t lMax = c->lMax;
  gsl_integration_workspace *workspace = 
    gsl_integration_workspace_alloc(INTEGRATION_WORKSPACE_LIMIT);
  
  // Manual init of the ewald summation targets
  complex double *target_c = calloc(nelem, sizeof(complex double));
  memset(target, 0, nelem * sizeof(complex double));
  double *err_c = NULL;
  if (err) {
    err_c = calloc(nelem, sizeof(double));
    memset(err, 0, nelem * sizeof(double));
  }
  

// CHOOSE POINT BEGIN
  for (size_t i = 0; i < npoints; ++i) { // BEGIN POINT LOOP
    point2d Rpq_shifted = Rpoints_plus_particle_shift[i];
    double r_pq_shifted = cart2norm(Rpq_shifted);
// CHOOSE POINT END

    // This should be imaginary, but gcc does not support:
    double Rpq_shifted_arg = atan2(Rpq_shifted.x, Rpq_shifted.y); // POINT-DEPENDENT
    complex double e_beta_Rpq = cexp(I*cart2_dot(beta, Rpq_shifted)); // POINT-DEPENDENT
    
    for(qpms_l_t n = 1; n <= lMax; ++n) {
      double complex prefacn = - I * pow(2./k, n+1) * M_2_SQRTPI / 2; // TODO put outside the R-loop and multiply in the end
      double R_pq_pown = pow(r_pq_shifted, n);
      // TODO the integral here
      double intres, interr;
      int retval = ewald32_sr_integral(r_pq_shifted, k, n, eta,
          &intres, &interr, workspace);
      if (retval) abort();
      for (qpms_m_t m = -n; m <= n; ++m){
        if((m+n) % 2 != 0) // odd coefficients are zero.
          continue; // nothing needed, already done by memset
        complex double e_imf = cexp(I*m*Rpq_shifted_arg);
        double leg = c->legendre0[gsl_sf_legendre_array_index(n, m)];
        qpms_y_t y = qpms_mn2y(m,n);
        if(err)
          kahanadd(err + y, err_c + y, cabs(leg * (prefacn / I) * R_pq_pown
              * interr)); // TODO include also other errors
        ckahanadd(target + y, target_c + y,
            prefacn * R_pq_pown * leg * intres * e_beta_Rpq * e_imf);
      }

    }
  }

  gsl_integration_workspace_free(workspace);
  if(err) free(err_c);
  free(target_c);
  return 0;
}

 
#if 0


int ewald32_sigma_long_points_and_shift (
    complex double *target_sigmalr_y, // must be c->nelem long
    const qpms_ewald32_constants_t *c,
    double eta, double k, double unitcell_area,
    size_t npoints, const point2d *Kpoints,
    point2d beta,
    point2d particle_shift
    );
int ewald32_sigma_long_shiftedpoints_rordered(
    complex double *target_sigmalr_y, // must be c->nelem long
    const qpms_ewald32_constants_t *c,
    double eta, double k, double unitcell_area,
    const points2d_rordered_t *Kpoints_plus_beta_rordered,
    point2d particle_shift
    );
int ewald32_sigma_short_points_and_shift(
    complex double *target_sigmasr_y, // must be c->nelem long
    const qpms_ewald32_constants_t *c, // N.B. not too useful here
    double eta, double k,
    size_t npoints, const point2d *Rpoints,
    point2d particle_shift
    );
int ewald32_sigma_short_points_rordered(
    complex double *target_sigmasr_y, // must be c->nelem long
    const qpms_ewald32_constants_t *c, // N.B. not too useful here
    double eta, double k,
    const points2d_rordered_t *Rpoints_plus_particle_shift_rordered,
    point2d particle_shift    // used only in the very end to multiply it by the phase
    );

#endif
