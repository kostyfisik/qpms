#include "ewald.h"
#include <stdlib.h>
#include "indexing.h"
#include "kahansum.h"
#include <assert.h>
#include <string.h>
#include "tiny_inlines.h"

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


qpms_ewald32_constants_t *qpms_ewald32_constants_init(const qpms_l_t lMax)
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
        c->s1_constfacs[y][j] = -0.5 * ipow(n+1) * min1pow((n+m)/2) 
          * sqrt((2*n + 1) * factorial(n-m) * factorial(n+m))
          * min1pow(j) * pow(0.5, n-2*j)
          / (factorial(j) * factorial((n-m)/2-j) * factorial((n+m)/2-j))
          * pow(0.5, 2*j-1);
      }
      s1_constfacs_sz_cumsum += 1 + c->s1_jMaxes[y];
    }
    else
      c->s1_constfacs[y] = NULL;
  }
  return c;
}

void qpms_ewald32_constants_free(qpms_ewald32_constants_t *c) {
  free(c->s1_constfacs);
  free(c->s1_constfacs_base);
  free(c->s1_jMaxes);
  free(c);
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
    for(qpms_l_t n = 0; n <= lMax; ++n)
      for(qpms_m_t m = -n; m <= n; ++m) {
        if((m+n) % 2 != 0) // odd coefficients are zero.
          continue;
        qpms_y_t y = qpms_mn2y(m, n);
        complex double e_imalpha_pq = cexp(I*m*arg_pq);
        complex double jsum, jsum_c; ckahaninit(&jsum, &jsum_c);
        double jsum_err, jsum_err_c; kahaninit(&jsum_err, &jsum_err_c); // TODO do I really need to kahan sum errors?
        for(qpms_l_t j = 0; j < (n-abs(m))/2; ++j) {
          complex double summand = pow(rbeta_pq/k, n-2*j) 
            * e_imalpha_pq * cpow(gamma_pq, 2*j-1) // * Gamma_pq[j] bellow (GGG) after error computation
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
 
int ewald32_sigma_short_shiftedpoints(
    complex double *target, // must be c->nelem long
    double *err,
    const qpms_ewald32_constants_t *c, // N.B. not too useful here
    const double eta, const double k,
    const size_t npoints, const point2d *Rpoints_plus_particle_shift,
    const point2d particle_shift           // used only in the very end to multiply it by the phase
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
  

  //// Zde jsem skončil







  free(err_c);
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
