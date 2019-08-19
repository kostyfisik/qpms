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
#include <gsl/gsl_sf_expint.h>

// parameters for the quadrature of integral in (4.6)
#ifndef INTEGRATION_WORKSPACE_LIMIT
#define INTEGRATION_WORKSPACE_LIMIT 30000
#endif

#ifndef INTEGRATION_EPSABS
#define INTEGRATION_EPSABS 1e-13
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

/// Metadata describing the normalisation conventions used in ewald32_constants_t.
typedef enum {
  EWALD32_CONSTANTS_ORIG, // As in [1, (4,5)], NOT USED right now.
  EWALD32_CONSTANTS_AGNOSTIC /* Not depending on the spherical harmonic sign/normalisation
                              * convention – the $e^{im\alpha_pq}$ term in [1,(4.5)] being
                              * replaced by the respective $Y_n^m(\pi/2,\alpha)$ 
                              * spherical harmonic. See notes/ewald.lyx.
                              */

} ewald3_constants_option;

static const ewald3_constants_option type = EWALD32_CONSTANTS_AGNOSTIC;

qpms_ewald3_constants_t *qpms_ewald3_constants_init(const qpms_l_t lMax /*, const ewald3_constants_option type */,
    const int csphase)
{
  qpms_ewald3_constants_t *c = malloc(sizeof(qpms_ewald3_constants_t));
  //if (c == NULL) return NULL; // Do I really want to do this?
  c->lMax = lMax;
  c->nelem_sc = qpms_lMax2nelem_sc(lMax);
  c->s1_jMaxes = malloc(c->nelem_sc * sizeof(qpms_l_t));
  c->s1_constfacs = malloc(c->nelem_sc * sizeof(complex double *));
  //if (c->s1_jMaxes == NULL) return NULL;

  // ----- the "xy-plane constants" ------
  //determine sizes
  size_t s1_constfacs_sz = 0;
  for (qpms_y_t y = 0; y < c->nelem_sc; ++y) {
    qpms_l_t n; qpms_m_t m; qpms_y2mn_sc_p(y, &m, &n);
    if ((m + n) % 2 == 0) 
      s1_constfacs_sz += 1 + (c->s1_jMaxes[y] = (n-abs(m))/2);
    else
      c->s1_jMaxes[y] = -1;
  }

  c->s1_constfacs_base = malloc(s1_constfacs_sz * sizeof(complex double));
  size_t s1_constfacs_sz_cumsum = 0;
  for (qpms_y_t y = 0; y < c->nelem_sc; ++y) {
    qpms_l_t n; qpms_m_t m; qpms_y2mn_sc_p(y, &m, &n);
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

  // ------ the "z-axis constants" -----
  // determine sizes
  size_t s1_constfacs_1Dz_sz;
  {
    const qpms_l_t sz_n_lMax = lMax/2 + 1;  // number of different j's for n = lMax
    s1_constfacs_1Dz_sz = (lMax % 2) ? isq(sz_n_lMax) + sz_n_lMax
                                     : isq(sz_n_lMax);
  }
  c->s1_constfacs_1Dz_base = malloc(s1_constfacs_1Dz_sz * sizeof(complex double));
  c->s1_constfacs_1Dz = malloc((lMax+1)*sizeof(complex double *));

  size_t s1_constfacs_1Dz_sz_cumsum = 0;
  for (qpms_l_t n = 0; n <= lMax; ++n) {
    c->s1_constfacs_1Dz[n] = c->s1_constfacs_1Dz_base + s1_constfacs_1Dz_sz_cumsum;
    for (qpms_l_t j = 0; j <= n/2; ++j) {
      switch(type) {
        case EWALD32_CONSTANTS_AGNOSTIC:
          c->s1_constfacs_1Dz[n][j] = -ipow(n+1) * min1pow(j) * factorial(n)
            / (factorial(j) * pow(2, 2*j) * factorial(n - 2*j));
          break;
        default:
          abort(); // wrong type argument or not implemented
      }
    }
    s1_constfacs_1Dz_sz_cumsum += 1 + n / 2;
  }
  assert(s1_constfacs_1Dz_sz == s1_constfacs_1Dz_sz_cumsum);


  c->legendre_csphase = csphase;
  c->legendre0 = malloc(gsl_sf_legendre_array_n(lMax) * sizeof(double));
  c->legendre_plus1 = malloc(gsl_sf_legendre_array_n(lMax) * sizeof(double));
  c->legendre_minus1 = malloc(gsl_sf_legendre_array_n(lMax) * sizeof(double));
  // N.B. here I use the GSL_SF_LEGENRE_NONE, in order to be consistent with translations.c
  c->legendre_normconv = GSL_SF_LEGENDRE_NONE;
  // Moreover, using this approach (i.e. gsl) takes about 64kB extra memory
  if(GSL_SUCCESS != gsl_sf_legendre_array_e(c->legendre_normconv, lMax, 0, csphase, c->legendre0))
    abort();
  if(GSL_SUCCESS != gsl_sf_legendre_array_e(c->legendre_normconv, lMax, +1, csphase, c->legendre_plus1))
    abort();
  if(GSL_SUCCESS != gsl_sf_legendre_array_e(c->legendre_normconv, lMax, -1, csphase, c->legendre_minus1))
    abort();
  return c;
}

void qpms_ewald3_constants_free(qpms_ewald3_constants_t *c) {
  free(c->legendre0);
  free(c->legendre_plus1);
  free(c->legendre_minus1);
  free(c->s1_constfacs);
  free(c->s1_constfacs_base);
  free(c->s1_constfacs_1Dz_base);
  free(c->s1_constfacs_1Dz);
  free(c->s1_jMaxes);
  free(c);
}



int ewald3_sigma0(complex double *result, double *err,
    const qpms_ewald3_constants_t *c,
    const double eta, const complex double k)
{
  qpms_csf_result gam;
  int retval = complex_gamma_inc_e(-0.5, -csq(k/(2*eta)), &gam);
  // FIXME DO THIS CORRECTLY gam.val = conj(gam.val); // We take the other branch, cf. [Linton, p. 642 in the middle]
  if (0 != retval)
    abort();
  *result = gam.val * c->legendre0[gsl_sf_legendre_array_index(0,0)] / 2 / M_SQRTPI;
  if(err) 
    *err = gam.err * fabs(c->legendre0[gsl_sf_legendre_array_index(0,0)] / 2 / M_SQRTPI);
  return 0;
}

int ewald3_21_xy_sigma_long (
    complex double *target, // must be c->nelem_sc long
    double *err,
    const qpms_ewald3_constants_t *c,
    const double eta, const complex double k,
    const double unitcell_volume /* with the corresponding lattice dimensionality */,
    const LatticeDimensionality latdim,
    PGen *pgen_K, const bool pgen_generates_shifted_points
    /* If false, the behaviour corresponds to the old ewald32_sigma_long_points_and_shift,
     * so the function assumes that the generated points correspond to the unshifted reciprocal Bravais lattice,
     * and adds beta to the generated points before calculations.
     * If true, it assumes that they are already shifted.
     */,
    const cart3_t beta,
    const cart3_t particle_shift
    )
{
  const bool k_is_real = (cimag(k) == 0);
  assert((latdim & LAT_XYONLY) && (latdim & SPACE3D));
  assert((latdim & LAT1D) || (latdim & LAT2D));
  const qpms_y_t nelem_sc = c->nelem_sc;
  const qpms_l_t lMax = c->lMax;
  
  // Manual init of the ewald summation targets
  complex double *target_c = calloc(nelem_sc, sizeof(complex double));
  memset(target, 0, nelem_sc * sizeof(complex double));
  double *err_c = NULL;
  if (err) {
    err_c = calloc(nelem_sc, sizeof(double));
    memset(err, 0, nelem_sc * sizeof(double));
  }

  const complex double commonfac = 1/(k*k*unitcell_volume); // used in the very end (CFC)
  if (k_is_real)
    assert(creal(commonfac) > 0);

  PGenSphReturnData pgen_retdata;
#ifndef NDEBUG
  double rbeta_pq_prev;
#endif
  // recycleable values if rbeta_pq stays the same:
  complex double gamma_pq;
  complex double z;
  complex double factor1d = 1; // the "additional" factor for the 1D case (then it is not 1)
  // space for Gamma_pq[j]'s
  qpms_csf_result Gamma_pq[lMax/2+1];

  // CHOOSE POINT BEGIN
  // TODO mayby PGen_next_sph is not the best coordinate system choice here
  while ((pgen_retdata = PGen_next_sph(pgen_K)).flags & PGEN_NOTDONE) { // BEGIN POINT LOOP
    cart3_t K_pq_cart;
    sph_t beta_pq_sph;
    if (pgen_generates_shifted_points) {
      beta_pq_sph = pgen_retdata.point_sph;
      const cart3_t beta_pq_cart = sph2cart(beta_pq_sph);
      K_pq_cart = cart3_add(cart3_scale(-1, beta), beta_pq_cart);
    } else { // as in the old _points_and_shift functions
      const sph_t K_pq_sph = pgen_retdata.point_sph;
      K_pq_cart = sph2cart(K_pq_sph);
      const cart3_t beta_pq_cart = cart3_add(K_pq_cart, beta);
      beta_pq_sph = cart2sph(beta_pq_cart);
    }

    const double rbeta_pq = beta_pq_sph.r;
    const double arg_pq = beta_pq_sph.phi;
    //const double beta_pq_theta = beta_pq_sph.theta; //unused

  // CHOOSE POINT END

    const complex double phasefac = cexp(I*cart3_dot(K_pq_cart,particle_shift)); // POINT-DEPENDENT (PFC) // !!!CHECKSIGN!!!

    const bool new_rbeta_pq = (!pgen_generates_shifted_points) || (pgen_retdata.flags & !PGEN_OLD_R);
    if (!new_rbeta_pq) assert(rbeta_pq == rbeta_pq_prev);


    // R-DEPENDENT BEGIN
    if (new_rbeta_pq) {
      gamma_pq = clilgamma(rbeta_pq/k);
      z = csq(gamma_pq*k/(2*eta)); 
      for(qpms_l_t j = 0; j <= lMax/2; ++j) {
        // TODO COMPLEX FIXME check the branches in the old lilgamma case
        int retval = complex_gamma_inc_e(0.5-j, z, Gamma_pq+j);
        // we take the other branch, cf. [Linton, p. 642 in the middle]: FIXME instead use the C11 CMPLX macros and fill in -O*I part to z in the line above
        //if(creal(z) < 0) 
        //  Gamma_pq[j].val = conj(Gamma_pq[j].val); //FIXME as noted above
        if(!(retval==0 || retval==GSL_EUNDRFLW)) abort();
      }
      if (latdim & LAT1D)
        factor1d =  M_SQRT1_2 * .5 * k * gamma_pq;
    }
    // R-DEPENDENT END
    
    // TODO optimisations: all the j-dependent powers can be done for each j only once, stored in array
    // and just fetched for each n, m pair
    for(qpms_l_t n = 0; n <= lMax; ++n)
      for(qpms_m_t m = -n; m <= n; ++m) {
        if((m+n) % 2 != 0) // odd coefficients are zero.
          continue;
        const qpms_y_t y = qpms_mn2y_sc(m, n);
        const complex double e_imalpha_pq = cexp(I*m*arg_pq);
        complex double jsum, jsum_c; ckahaninit(&jsum, &jsum_c);
        double jsum_err, jsum_err_c; kahaninit(&jsum_err, &jsum_err_c); // TODO do I really need to kahan sum errors?
        assert((n-abs(m))/2 == c->s1_jMaxes[y]);
        for(qpms_l_t j = 0; j <= c->s1_jMaxes[y]/*(n-abs(m))/2*/; ++j) { // FIXME </<= ?
          complex double summand = cpow(rbeta_pq/k, n-2*j) 
            * e_imalpha_pq  * c->legendre0[gsl_sf_legendre_array_index(n,abs(m))] * min1pow_m_neg(m) // This line can actually go outside j-loop
            * cpow(gamma_pq, 2*j-1) // * Gamma_pq[j] bellow (GGG) after error computation
            * c->s1_constfacs[y][j];
          if(err) {
            // FIXME include also other errors than Gamma_pq's relative error
             kahanadd(&jsum_err, &jsum_err_c, Gamma_pq[j].err * cabs(summand));
          }
          summand *= Gamma_pq[j].val; // GGG
          ckahanadd(&jsum, &jsum_c, summand);
        }
        jsum *= phasefac * factor1d; // PFC
        ckahanadd(target + y, target_c + y, jsum);
        if(err) kahanadd(err + y, err_c + y, jsum_err);
      }
#ifndef NDEBUG
    rbeta_pq_prev = rbeta_pq;
#endif
  } // END POINT LOOP
 
  free(err_c);
  free(target_c);
  for(qpms_y_t y = 0; y < nelem_sc; ++y) // CFC common factor from above
    target[y] *= commonfac;
  if(err)
    for(qpms_y_t y = 0; y < nelem_sc; ++y)
      err[y] *= commonfac;
  return 0;
}


// 1D chain along the z-axis; not many optimisations here as the same
// shifted beta radius could be recycled only once anyways
int ewald3_1_z_sigma_long (
    complex double *target, // must be c->nelem_sc long
    double *err,
    const qpms_ewald3_constants_t *c,
    const double eta, const complex double k,
    const double unitcell_volume /* length (periodicity) in this case */,
    const LatticeDimensionality latdim,
    PGen *pgen_K, const bool pgen_generates_shifted_points
    /* If false, the behaviour corresponds to the old ewald32_sigma_long_points_and_shift,
     * so the function assumes that the generated points correspond to the unshifted reciprocal Bravais lattice,
     * and adds beta to the generated points before calculations.
     * If true, it assumes that they are already shifted.
     */,
    const cart3_t beta,
    const cart3_t particle_shift
    )
{
  assert(LatticeDimensionality_checkflags(latdim, LAT_1D_IN_3D_ZONLY));
  assert(beta.x == 0 && beta.y == 0);
  assert(particle_shift.x == 0 && particle_shift.y == 0);
  const double beta_z = beta.z;
  const double particle_shift_z = particle_shift_z;
  const qpms_y_t nelem_sc = c->nelem_sc;
  const qpms_l_t lMax = c->lMax;
  
  // Manual init of the ewald summation targets
  complex double *target_c = calloc(nelem_sc, sizeof(complex double));
  memset(target, 0, nelem_sc * sizeof(complex double));
  double *err_c = NULL;
  if (err) {
    err_c = calloc(nelem_sc, sizeof(double));
    memset(err, 0, nelem_sc * sizeof(double));
  }

  const double commonfac = 1/(k*unitcell_volume); // multiplied in the very end (CFC)
  assert(commonfac > 0);

  // space for Gamma_pq[j]'s (I rewrite the exp. ints. E_n in terms of Gamma fns., cf. my ewald.lyx notes, (eq:1D_z_LRsum).
  qpms_csf_result Gamma_pq[lMax/2+1];

  PGenSphReturnData pgen_retdata;
  // CHOOSE POINT BEGIN
  // TODO FIXME USE PGen_next_z
  while ((pgen_retdata = PGen_next_sph(pgen_K)).flags & PGEN_NOTDONE) { // BEGIN POINT LOOP
    assert(pgen_retdata.flags & PGEN_AT_Z);
    double K_z, beta_mu_z;
    if (pgen_generates_shifted_points) {
      beta_mu_z = ((pgen_retdata.point_sph.theta == 0) ?
        pgen_retdata.point_sph.r : -pgen_retdata.point_sph.r); //!!!CHECKSIGN!!!
      K_z = beta_mu_z - beta_z;
    } else { // as in the old _points_and_shift functions
      K_z = ((pgen_retdata.point_sph.theta == 0) ?
          pgen_retdata.point_sph.r : -pgen_retdata.point_sph.r); // !!!CHECKSIGN!!!
      beta_mu_z = K_z + beta_z;
    }
    double rbeta_mu = fabs(beta_mu_z);
  // CHOOSE POINT END

    const complex double phasefac = cexp(I * K_z * particle_shift_z); // POINT-DEPENDENT (PFC) // !!!CHECKSIGN!!!

    // R-DEPENDENT BEGIN
    complex double gamma_pq = clilgamma(rbeta_mu/k);  // For real beta and k this is real or pure imaginary ...
    const complex double z = csq(gamma_pq*k/(2*eta));// ... so the square (this) is in fact real.
    for(qpms_l_t j = 0; j <= lMax/2; ++j) {
      int retval = complex_gamma_inc_e(0.5-j, z, Gamma_pq+j);
      // we take the other branch, cf. [Linton, p. 642 in the middle]: FIXME instead use the C11 CMPLX macros and fill in -O*I part to z in the line above
      //if(creal(z) < 0) 
      //  Gamma_pq[j].val = conj(Gamma_pq[j].val); //FIXME as noted above
      if(!(retval==0 || retval==GSL_EUNDRFLW)) abort();
    }
    // R-DEPENDENT END
    // TODO optimisations: all the j-dependent powers can be done for each j only once, stored in array
    // and just fetched for each n
    for(qpms_l_t n = 0; n <= lMax; ++n) {
        const qpms_y_t y = qpms_mn2y_sc(0, n);
        complex double jsum, jsum_c; ckahaninit(&jsum, &jsum_c);
        double jsum_err, jsum_err_c; kahaninit(&jsum_err, &jsum_err_c); // TODO do I really need to kahan sum errors?
        for(qpms_l_t j = 0; j <= n/2; ++j) {
          complex double summand = pow(rbeta_mu/k, n-2*j) 
            * ((beta_mu_z > 0) ? // TODO this can go outsize the j-loop
                c->legendre_plus1[gsl_sf_legendre_array_index(n,0)] :
                (c->legendre_minus1[gsl_sf_legendre_array_index(n,0)] * min1pow(n))
              )
            // * min1pow_m_neg(m) // not needed as m == 0
            * cpow(gamma_pq, 2*j) // * Gamma_pq[j] bellow (GGG) after error computation
            * c->s1_constfacs_1Dz[n][j];
              /* s1_consstfacs_1Dz[n][j] =
               *
               *     -I**(n+1) (-1)**j * n!
               *   --------------------------
               *   j! * 2**(2*j) * (n - 2*j)!
               */   

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
  for(qpms_y_t y = 0; y < nelem_sc; ++y) // CFC common factor from above
    target[y] *= commonfac;
  if(err)
    for(qpms_y_t y = 0; y < nelem_sc; ++y)
      err[y] *= commonfac;
  return 0;
}

int ewald3_sigma_long (
    complex double *target, // must be c->nelem_sc long
    double *err,
    const qpms_ewald3_constants_t *c,
    const double eta, const complex double k,
    const double unitcell_volume /* with the corresponding lattice dimensionality */,
    const LatticeDimensionality latdim,
    PGen *pgen_K, const bool pgen_generates_shifted_points
    /* If false, the behaviour corresponds to the old ewald32_sigma_long_points_and_shift,
     * so the function assumes that the generated points correspond to the unshifted reciprocal Bravais lattice,
     * and adds beta to the generated points before calculations.
     * If true, it assumes that they are already shifted.
     */,
    const cart3_t beta,
    const cart3_t particle_shift
    )
{
  assert(latdim & SPACE3D);
  if (latdim & LAT_XYONLY)
    return ewald3_21_xy_sigma_long(target, err, c, eta, k, unitcell_volume,
        latdim, pgen_K, pgen_generates_shifted_points, beta, particle_shift);
  else if (latdim & LAT_ZONLY)
    return ewald3_1_z_sigma_long(target, err, c, eta, k, unitcell_volume,
        latdim, pgen_K, pgen_generates_shifted_points, beta, particle_shift);
  // TODO 3D case and general 2D case
  else abort(); // NOT IMPLEMENTED
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

// a version allowing complex k

struct sigma2_integrand_params_ck {
  int n;
  complex double k;
  double R;
};

// TODO ther might be some space for optimisation if needed, as now we calculate everything twice
// including the whole complex exponential (throwing the imaginary/real part away)
static double sigma2_integrand_ck_real(double ksi, void *params) {
  struct sigma2_integrand_params_ck *p = (struct sigma2_integrand_params_ck *) params;
  return creal(cexp(-csq(p->R * ksi) + csq(p->k / ksi / 2))) * pow(ksi, 2*p->n);
}
static double sigma2_integrand_ck_imag(double ksi, void *params) {
  struct sigma2_integrand_params_ck *p = (struct sigma2_integrand_params_ck *) params;
  return cimag(cexp(-csq(p->R * ksi) + csq(p->k / ksi / 2))) * pow(ksi, 2*p->n);
}

static int ewald32_sr_integral_ck(double r, complex double k, int n, double eta,
    complex double *result, double *err, gsl_integration_workspace *workspace)
{
  struct sigma2_integrand_params_ck p;


  const double R0 = r; // Maybe could be chosen otherwise, but fuck it for now.
  p.n = n;
  eta *= R0;
  p.k = k * R0;
  p.R = r / R0;  // i.e. p.R = 1

  gsl_function F;
  F.params = &p;
  double result_real, result_imag, err_real, err_imag;
  
  F.function = sigma2_integrand_ck_real;
  // TODO check return values
  int retval = gsl_integration_qagiu(&F, eta, INTEGRATION_EPSABS, 
      INTEGRATION_EPSREL, INTEGRATION_WORKSPACE_LIMIT,
      workspace, &result_real, &err_real);
  
  F.function = sigma2_integrand_ck_imag;
  // TODO check return values
  retval = gsl_integration_qagiu(&F, eta, INTEGRATION_EPSABS, 
      INTEGRATION_EPSREL, INTEGRATION_WORKSPACE_LIMIT,
      workspace, &result_imag, &err_imag);

  *result = result_real + I*result_imag;
  *err = sqrt(sq(err_real) + sq(err_imag));

  double normfac = pow(R0, -2*p.n - 1);
  *result *= normfac;
  *err *= normfac;
  return retval;
}

int ewald3_sigma_short( 
                complex double *target, // must be c->nelem_sc long
                double *err, // must be c->nelem_sc long or NULL
                const qpms_ewald3_constants_t *c,
                const double eta, const complex double k /* TODO COMPLEX */,
                const LatticeDimensionality latdim, // apart from asserts and possible optimisations ignored, as the SR formula stays the same
                PGen *pgen_R, const bool pgen_generates_shifted_points
                /* If false, the behaviour corresponds to the old ewald32_sigma_short_points_and_shift,
                 * so the function assumes that the generated points correspond to the unshifted Bravais lattice,
                 * and adds particle_shift to the generated points before calculations.
                 * If true, it assumes that they are already shifted (if calculating interaction between
                 * different particles in the unit cell).
                 */,
                const cart3_t beta,
                const cart3_t particle_shift
                )
{
  const bool k_is_real = (cimag(k) == 0); // TODO check how the compiler optimises the loops
  const double kreal = creal(k);
  const qpms_y_t nelem_sc = c->nelem_sc;
  const qpms_l_t lMax = c->lMax;
  gsl_integration_workspace *workspace = 
    gsl_integration_workspace_alloc(INTEGRATION_WORKSPACE_LIMIT);

  double legendre_buf[gsl_sf_legendre_array_n(c->lMax)]; // work space for the legendre polynomials (used only in the general case)
  
  // Manual init of the ewald summation targets
  complex double * const target_c = calloc(nelem_sc, sizeof(complex double));
  memset(target, 0, nelem_sc * sizeof(complex double));
  double *err_c = NULL;
  if (err) {
    err_c = calloc(nelem_sc, sizeof(double));
    memset(err, 0, nelem_sc * sizeof(double));
  }
  

  PGenSphReturnData pgen_retdata;
#ifndef NDEBUG
  double r_pq_shifted_prev;
#endif
  // recyclable variables if r_pq_shifted stays the same:
  double intres[lMax+1], interr[lMax+1];
  complex double cintres[lMax+1];

// CHOOSE POINT BEGIN
// TODO check whether _next_sph is the optimal coordinate system choice here
  while ((pgen_retdata = PGen_next_sph(pgen_R)).flags & PGEN_NOTDONE) { // BEGIN POINT LOOP
// CHOOSE POINT END
    cart3_t Rpq_shifted_cart; // I will need both sph and cart representations anyway...
    sph_t Rpq_shifted_sph;
    if (pgen_generates_shifted_points) {
      Rpq_shifted_sph = pgen_retdata.point_sph;
      Rpq_shifted_cart = sph2cart(Rpq_shifted_sph);
    } else { // as in the old _points_and_shift functions
      //const point2d Rpq_shifted = Rpoints_plus_particle_shift[i];
      const sph_t bravais_point_sph = pgen_retdata.point_sph;
      const cart3_t bravais_point_cart = sph2cart(bravais_point_sph);
      Rpq_shifted_cart = cart3_add(bravais_point_cart, cart3_scale(-1, particle_shift)); // CHECKSIGN!!!
      Rpq_shifted_sph = cart2sph(Rpq_shifted_cart);
    }

    // TODO eliminate and use the Rpq_shifted_sph members directly (but in compiler optimizations we trust)
    const double Rpq_shifted_arg = Rpq_shifted_sph.phi; //atan2(Rpq_shifted.y, Rpq_shifted.x); // POINT-DEPENDENT
    const double Rpq_shifted_theta = Rpq_shifted_sph.theta; // POINT-DEPENDENT
    const double r_pq_shifted = Rpq_shifted_sph.r;

    // if the radius is the same as in previous cycle, most of the calculations can be recycled
    const bool new_r_pq_shifted = (!pgen_generates_shifted_points) || (pgen_retdata.flags & !PGEN_OLD_R);
    if (!new_r_pq_shifted) assert(r_pq_shifted_prev == r_pq_shifted);

    const complex double e_beta_Rpq = cexp(I*cart3_dot(beta, Rpq_shifted_cart)); // POINT-DEPENDENT
    LatticeDimensionality speedup_regime = 0;
    if ((latdim & LAT_2D_IN_3D_XYONLY) == LAT_2D_IN_3D_XYONLY) speedup_regime = LAT_2D_IN_3D_XYONLY;
    if ((latdim & LAT_1D_IN_3D_ZONLY)  == LAT_1D_IN_3D_ZONLY)  speedup_regime = LAT_1D_IN_3D_ZONLY;

    const double * legendre_array;

    switch(speedup_regime) {
      // speedup checks for special geometries and Legendre polynomials
      case LAT_1D_IN_3D_ZONLY:
        assert((pgen_retdata.flags & PGEN_AT_Z) == PGEN_AT_Z);
        assert(Rpq_shifted_theta == M_PI || Rpq_shifted_theta == 0);
        legendre_array = (Rpq_shifted_theta == 0) ? c->legendre_plus1 : c->legendre_minus1; // CHECKSIGN
        break;
      case LAT_2D_IN_3D_XYONLY:
        assert((pgen_retdata.flags &PGEN_AT_XY) == PGEN_AT_XY);
        assert(fabs(Rpq_shifted_theta - M_PI_2) < DBL_EPSILON * 1024);
        // assert(Rpq_shifted_theta == M_PI_2); // FIXME this should work as well
        legendre_array = c->legendre0;
        break;
      default:
        if(GSL_SUCCESS != gsl_sf_legendre_array_e(c->legendre_normconv, lMax, cos(Rpq_shifted_theta), c->legendre_csphase, legendre_buf))
        abort();
        legendre_array = legendre_buf;
        break;
    }
    
    for(qpms_l_t n = 0; n <= lMax; ++n) {
      const double complex prefacn = - I * (k_is_real ? pow(2./creal(k),n+1) : cpow(2./k, n+1)) * M_2_SQRTPI / 2; // profiling TODO put outside the R-loop and multiply in the end?
      const double R_pq_pown = pow(r_pq_shifted, n); // profiling TODO: maybe put this into the new_r_pq_shifted condition as well?
      if (new_r_pq_shifted) {
        int retval;
        if (k_is_real) {
          double intres_real;
          retval = ewald32_sr_integral(r_pq_shifted, kreal, n, eta,
               &intres_real, interr + n, workspace);
          cintres[n] = intres_real;
        } else
          retval = ewald32_sr_integral_ck(r_pq_shifted, k, n, eta,
               cintres+n, interr + n, workspace);
        if (retval) abort();
      } // otherwise recycle the integrals
      for (qpms_m_t m = -n; m <= n; ++m){
        complex double e_imf;
        // SPEEDUPS for some special geometries
        if(speedup_regime == LAT_2D_IN_3D_XYONLY) { //2D lattice inside the xy plane
          if((m+n) % 2 != 0) // odd coefficients are zero.
            continue; // nothing needed, already done by memset
          e_imf = cexp(I*m*Rpq_shifted_arg); // profiling TODO: calculate outside the n-loop?
        } else if (speedup_regime == LAT_1D_IN_3D_ZONLY) { // 1D lattice along the z axis
          if (m != 0) // m-non-zero coefficients are zero
            continue; // nothing needed, already done by memset
          e_imf = 1;
        } else { // general 1D/2D/3D lattice in 3D space
          e_imf = cexp(I*m*Rpq_shifted_arg);
        }

        const double leg = legendre_array[gsl_sf_legendre_array_index(n, abs(m))];

        const qpms_y_t y = qpms_mn2y_sc(m,n);
        if(err)
          kahanadd(err + y, err_c + y, cabs(leg * (prefacn / I) * R_pq_pown
              * interr[n])); // TODO include also other errors
        ckahanadd(target + y, target_c + y,
            prefacn * R_pq_pown * leg * cintres[n] * e_beta_Rpq * e_imf * min1pow_m_neg(m));
      }

    }
#ifndef NDEBUG
    r_pq_shifted_prev = r_pq_shifted;
#endif
  }

  gsl_integration_workspace_free(workspace);
  if(err) free(err_c);
  free(target_c);
  return 0;
}

