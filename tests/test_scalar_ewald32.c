// Perform Ewald summation (2D xy-lattice in 3D space ) of SSWFs with different Ewald parameters and check whether the difference is inside the tolerance range.
// run as
// test_scalar_ewald32 lMax a1.x a1.y a2.x a2.y wavenum.real wavenum.imag k.x k.y particle_shift.x particle_shift.y csphase rtol atol maxR maxK eta1 [eta2 [eta 3 ...]]
// c99 -o ewaldshift3g_vargeom -ggdb -Wall -I ../ ewaldshift3g_vargeom.c ../qpms/ewald.c ../qpms/ewaldsf.c  ../qpms/lattices2d.c ../qpms/latticegens.c -lgsl -lm -lblas



// implementation of the [LT(4.16)] test
#include <math.h>
#define M_SQRTPI 1.7724538509055160272981674833411452 
#define M_SQRT3 1.7320508075688772935274463415058724
#include <qpms/ewald.h>
#include <qpms/tiny_inlines.h>
#include <qpms/indexing.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_sf_legendre.h>
#include <fenv.h>

typedef struct ewaldtest2d_params {
  qpms_l_t lMax;
  point2d b1, b2;
  point2d beta; 
  point2d particle_shift;
  complex double k; 
  //double a;
  double eta;
  double maxR;
  double maxK;
  double csphase;
} ewaldtest2d_params;

typedef struct ewaldtest2d_results {
  ewaldtest2d_params p;
  complex double *sigmas_short,
          *sigmas_long,
          sigma0,
          *sigmas_total;
  double *err_sigmas_short,
    *err_sigmas_long,
    err_sigma0,
    *err_sigmas_total;
  complex double *regsigmas_416;
} ewaldtest2d_results;


void ewaldtest2d_results_free(ewaldtest2d_results *r) {
  free(r->sigmas_short);
  free(r->sigmas_long);
  free(r->sigmas_total);
  free(r->err_sigmas_long);
  free(r->err_sigmas_total);
  free(r->err_sigmas_short);
  free(r->regsigmas_416);
  free(r);
}

static inline double san(double x) {
  return fabs(x) < 1e-13 ? 0 : x;
}

int isclose_cmplx(complex double a, complex double b, double rtol, double atol) {
  return cabs(a-b) <= atol + rtol * .5 * (cabs(a) + cabs(b));
}

ewaldtest2d_results *ewaldtest2d(const ewaldtest2d_params p);

int main(int argc, char **argv) {
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  bool verbose = !!getenv("QPMS_VERBOSE_TESTS");
  gsl_set_error_handler(IgnoreUnderflowsGSLErrorHandler);
  QPMS_ENSURE(argc >= 18, "At least 16 arguments expected, I see only %d.", argc-1);
  int netas = argc - 17;
  ewaldtest2d_params plist[netas];
  double atol, rtol;
  plist[0].lMax = atoi(argv[1]);
  plist[0].b1.x = strtod(argv[2], NULL);
  plist[0].b1.y = strtod(argv[3], NULL);
  plist[0].b2.x = strtod(argv[4], NULL);
  plist[0].b2.y = strtod(argv[5], NULL);
  plist[0].k = strtod(argv[6], NULL) + I*strtod(argv[7], NULL);
  plist[0].beta.x = strtod(argv[8], NULL);
  plist[0].beta.y = strtod(argv[9], NULL);
  plist[0].particle_shift.x = strtod(argv[10], NULL);
  plist[0].particle_shift.y = strtod(argv[11], NULL);
  plist[0].csphase = strtod(argv[12], NULL);
  atol = strtod(argv[13], NULL);
  rtol = strtod(argv[14], NULL);
  plist[0].maxR = strtod(argv[15], NULL);
  plist[0].maxK = strtod(argv[16], NULL);
  plist[0].eta  = strtod(argv[17], NULL);
  for(int i = 1; i < netas; ++i) {
    plist[i] = plist[0];
    plist[i].eta = strtod(argv[17+i], NULL);
  }
  
  ewaldtest2d_results *r[netas];

  int fails = 0;

  for (size_t i = 0; i < netas; ++i) {
    ewaldtest2d_params p = plist[i];
    r[i] = ewaldtest2d(p);
    // TODO print per-test header here
    printf("===============================\n");
    printf("b1 = (%g, %g), b2 = (%g, %g)," /* "K1 = (%g, %g), K2 = (%g, %g),"*/ " Kmax = %g, Rmax = %g, lMax = %d, eta = %g, k = %g%+gj, beta = (%g,%g), ps = (%g,%g), csphase = %g\n",
        p.b1.x, p.b1.y, p.b2.x, p.b2.y,/*TODO K1, K2*/ p.maxK, p.maxR, p.lMax, p.eta, creal(p.k), cimag(p.k), p.beta.x, p.beta.y, p.particle_shift.x, p.particle_shift.y, p.csphase);
    printf("sigma0: %.16g%+.16gj\n", creal(r[i]->sigma0), cimag(r[i]->sigma0));    
    for (qpms_l_t n = 0; n <= p.lMax; ++n) {
      for (qpms_m_t m = -n; m <= n; ++m){
        if ((m+n)%2) continue;
        qpms_y_t y = qpms_mn2y_sc(m,n);
        qpms_y_t y_conj = qpms_mn2y_sc(-m,n);
        // y n m sigma_total (err), regsigmas_416 regsigmas_415_recon
        if (verbose) printf("%zd %d %d: T:%.16g%+.16gj(%.3g) L:%.16g%+.16gj(%.3g) S:%.16g%+.16gj(%.3g) \n"
            //"| predict %.16g%+.16gj \n| actual  %.16g%+.16gj\n"
            ,
            y, n, m, creal(san(r[i]->sigmas_total[y])), san(cimag(r[i]->sigmas_total[y])),
            r[i]->err_sigmas_total[y],
            san(creal(r[i]->sigmas_long[y])), san(cimag(r[i]->sigmas_long[y])),
            r[i]->err_sigmas_long[y],
            san(creal(r[i]->sigmas_short[y])), san(cimag(r[i]->sigmas_short[y])),
            r[i]->err_sigmas_short[y]
            // TODO and count big differences as failures.
            //san(creal(r[i]->regsigmas_416[y])), san(cimag(r[i]->regsigmas_416[y])),
            //san(creal(r[i]->sigmas_total[y]) + creal(r[i]->sigmas_total[y_conj])),
            //san(cimag(r[i]->sigmas_total[y]) - cimag(r[i]->sigmas_total[y_conj]))
        );
      }
    }
  }

  bool toprint[netas];

  for (qpms_l_t n = 0; n <= plist[0].lMax; ++n) {
    for (qpms_m_t m = -n; m <= n; ++m){
      memset(toprint, 0, netas*sizeof(bool));
      if ((m+n)%2) continue;
      qpms_y_t y = qpms_mn2y_sc(m,n);
      qpms_y_t y_conj = qpms_mn2y_sc(-m,n);
      for (size_t i = 0; i < netas; ++i) {
        for (size_t j = i+1; j < netas; ++j){ 
          if (!isclose_cmplx(r[i]->sigmas_total[y], r[j]->sigmas_total[y], rtol, atol)) 
            toprint[i] = toprint[j] = true;
        }
        if (toprint[i]) {
          ++fails;
          printf("with eta = %.16g:\n", plist[i].eta);
          printf("%zd %d %d: T:%.16g%+.16gj(%.3g) L:%.16g%+.16gj(%.3g) S:%.16g%+.16gj(%.3g) \n"
              //"| predict %.16g%+.16gj \n| actual  %.16g%+.16gj\n"
              ,
              y, n, m, creal(san(r[i]->sigmas_total[y])), san(cimag(r[i]->sigmas_total[y])),
              r[i]->err_sigmas_total[y],
              san(creal(r[i]->sigmas_long[y])), san(cimag(r[i]->sigmas_long[y])),
              r[i]->err_sigmas_long[y],
              san(creal(r[i]->sigmas_short[y])), san(cimag(r[i]->sigmas_short[y])),
              r[i]->err_sigmas_short[y]
              // TODO and count big differences as failures.
              //san(creal(r[i]->regsigmas_416[y])), san(cimag(r[i]->regsigmas_416[y])),
              //san(creal(r[i]->sigmas_total[y]) + creal(r[i]->sigmas_total[y_conj])),
              //san(cimag(r[i]->sigmas_total[y]) - cimag(r[i]->sigmas_total[y_conj]))
              );
          //if(!y) printf("0:%.16g%+.16g\n", san(creal(r[i]->sigma0)), san(cimag(r[i]->sigma0)));
        }
      }
    }
  }

  for (size_t i = 0; i < netas; ++i) 
    ewaldtest2d_results_free(r[i]);

  return fails;
}


int ewaldtest_counter = 0;


ewaldtest2d_results *ewaldtest2d(const ewaldtest2d_params p) {
  cart3_t beta3 = cart22cart3xy(p.beta);
  cart3_t particle_shift3 = cart22cart3xy(p.particle_shift);
 
  cart2_t b1 = p.b1, b2 = p.b2, rb1, rb2;
  if (QPMS_SUCCESS != l2d_reciprocalBasis2pi(b1, b2, &rb1, &rb2))
    abort();

  const double A = l2d_unitcell_area(b1, b2); // sqrt(3) * a * a / 2.; // unit cell size
  const double K_len = cart2norm(rb1)+cart2norm(rb2); //4*M_PI/a/sqrt(3); // reciprocal vector length

  ewaldtest2d_results *results = malloc(sizeof(ewaldtest2d_results));
  results->p = p;
 
  // skip zeroth point if it coincides with origin
  bool include_origin = !(fabs(p.particle_shift.x) == 0 
      && fabs(p.particle_shift.y) == 0);

  PGen Rlgen = PGen_xyWeb_new(b1, b2, BASIS_RTOL, CART2_ZERO, 0, include_origin, p.maxR, false);
  //PGen Rlgen_plus_shift = PGen_xyWeb_new(b1, b2, BASIS_RTOL, cart2_scale(-1 /* CHECKSIGN */, particle_shift2), 0, include_origin, p.maxR + a, false);
  PGen Klgen = PGen_xyWeb_new(rb1, rb2, BASIS_RTOL, CART2_ZERO, 0, true, p.maxK + K_len, false);
  //PGen Klgen_plus_beta = PGen_xyWeb_new(rb1, rb2, BASIS_RTOL, beta2, 0, true, p.maxK + K_len, false); 

  qpms_y_t nelem_sc = qpms_lMax2nelem_sc(p.lMax);

  results->sigmas_short = malloc(sizeof(complex double)*nelem_sc);
  results->sigmas_long = malloc(sizeof(complex double)*nelem_sc);
  results->sigmas_total = malloc(sizeof(complex double)*nelem_sc);
  results->err_sigmas_short = malloc(sizeof(double)*nelem_sc);
  results->err_sigmas_long = malloc(sizeof(double)*nelem_sc);
  results->err_sigmas_total = malloc(sizeof(double)*nelem_sc);

  qpms_ewald3_constants_t *c = qpms_ewald3_constants_init(p.lMax, p.csphase);

  if (0!=ewald3_sigma_long(results->sigmas_long,
        results->err_sigmas_long, c, p.eta, p.k, A,
        LAT_2D_IN_3D_XYONLY, &Klgen, false, beta3, particle_shift3))
      abort();
  if (0!=ewald3_sigma_short(
    results->sigmas_short, results->err_sigmas_short, c,
    p.eta, p.k, LAT_2D_IN_3D_XYONLY, &Rlgen, false,  beta3, particle_shift3))
    abort();
  if (0!=ewald3_sigma0(&(results->sigma0), &(results->err_sigma0), c, p.eta, p.k))
    abort();
  for(qpms_y_t y = 0; y < nelem_sc; ++y) {
    results->sigmas_total[y] = results->sigmas_short[y] + results->sigmas_long[y];
    results->err_sigmas_total[y] = results->err_sigmas_short[y] + results->err_sigmas_long[y];
  }
  if(!include_origin) { // "Renormalised" contribution of origin point
    results->sigmas_total[0] += results->sigma0;
    results->err_sigmas_total[0] += results->err_sigma0;
  }

  // Now calculate the reference values [LT(4.16)]
  results->regsigmas_416 = calloc(nelem_sc, sizeof(complex double));
  results->regsigmas_416[0] = -2 * c->legendre0[gsl_sf_legendre_array_index(0,0)];
  
#if 0 // not yet implemented for the new API
  {
    double legendres[gsl_sf_legendre_array_n(p.lMax)];
    points2d_rordered_t sel = 
      points2d_rordered_annulus(Kpoints_plus_beta, 0, true, p.k, false);
    if (0 != sel.nrs) 
    { 
      point2d *beta_pq_lessthan_k = sel.base + sel.r_offsets[0];
      size_t beta_pq_lessthan_k_count = sel.r_offsets[sel.nrs] - sel.r_offsets[0];
      for(size_t i = 0; i < beta_pq_lessthan_k_count; ++i) {
        point2d beta_pq = beta_pq_lessthan_k[i];
        double rbeta_pq = cart2norm(beta_pq);
        double arg_pq = atan2(beta_pq.y, beta_pq.x);
        double denom = sqrt(p.k*p.k - rbeta_pq*rbeta_pq);
        if( gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,
              p.lMax, denom/p.k, p.csphase, legendres) != 0)
          abort();
        for (qpms_y_t y = 0; y < nelem_sc; ++y) {
          qpms_l_t n; qpms_m_t m; 
          qpms_y2mn_sc_p(y, &m, &n);
          if ((m+n)%2 != 0)
            continue;
          complex double eimf = cexp(I*m*arg_pq);
          results->regsigmas_416[y] +=
            4*M_PI*ipow(n)/p.k/A 
            * eimf * legendres[gsl_sf_legendre_array_index(n,abs(m))] * min1pow_m_neg(m)
            / denom;
        }
      }
    }
  }
#else
  for(qpms_y_t y = 0; y < nelem_sc; ++y) {
          qpms_l_t n; qpms_m_t m; 
          qpms_y2mn_sc_p(y, &m, &n);
          if ((m+n)%2 != 0)
            continue;
          results->regsigmas_416[y] = NAN;
  }
#endif


  qpms_ewald3_constants_free(c);
  ++ewaldtest_counter;
  return results;
}
