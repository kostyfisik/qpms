// c99 -o ewaldshift_3g -ggdb -Wall -I ../ ewaldshift_3g.c ../qpms/ewald.c ../qpms/ewaldsf.c  ../qpms/lattices2d.c ../qpms/latticegens.c -lgsl -lm -lblas

// implementation of the [LT(4.16)] test
#include <math.h>
#define M_SQRTPI 1.7724538509055160272981674833411452 
#define M_SQRT3 1.7320508075688772935274463415058724
#include <qpms/ewald.h>
#include <qpms/tiny_inlines.h>
#include <qpms/indexing.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <gsl/gsl_sf_legendre.h>
typedef struct ewaldtest_triang_params {
  qpms_l_t lMax;
  point2d beta; 
  point2d particle_shift;
  double k; 
  double a;
  double eta;
  double maxR;
  double maxK;
  double csphase;
  TriangularLatticeOrientation orientation;
} ewaldtest_triang_params;

typedef struct ewaldtest_triang_results {
  ewaldtest_triang_params p;
  complex double *sigmas_short,
          *sigmas_long,
          sigma0,
          *sigmas_total;
  double *err_sigmas_short,
    *err_sigmas_long,
    err_sigma0,
    *err_sigmas_total;
  complex double *regsigmas_416;
} ewaldtest_triang_results;


/*
const double a = 582e-9;
const double inv_a = 4*M_PI/a/M_SQRT3;
const double Klen = 4*M_PI/a/3;
*/

#define AA (582.e-9)
#define INV_A (4*M_PI/AA/M_SQRT3)
#define KLEN (4*M_PI/AA/3)

ewaldtest_triang_params paramslist[] = {
// lMax, beta, shift, k, a, eta, maxR, maxK, csphase, orientation
/*
 { 2, {2.7, 1}, {0.5,0.1325}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {2.7, 1}, {0.5,0.1325}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {2.7, 1}, {0.5,0.1325}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {2.7, 1}, {0.5,0.1325}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 2, {1.1, 1}, {0.5,0.1325}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 1}, {0.5,0.1325}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 1}, {0.5,0.1325}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 1}, {0.5,0.1325}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 2, {1.1, 1}, {0.5,0.}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 1}, {0.5,0.}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 1}, {0.5,0.}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 1}, {0.5,0.}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},
*/

//
  { 3, {0,4198609.6394310603}, {0,0}, 11255786.828366444, 9.9766126515967311e-07, 29088820.866572164, 20*9.9766126515967311e-07, 20*7272205.21664304, 1., TRIANGULAR_VERTICAL},
  { 3, {0,4198609.6394310603}, {0,0}, 11255786.828366444, 9.9766126515967311e-07, 0.18*29088820.866572164, 20*9.9766126515967311e-07, 20*7272205.21664304, 1., TRIANGULAR_VERTICAL},
  { 3, {0,4198609.6394310603}, {0,0}, 11255786.828366444, 9.9766126515967311e-07, 0.13*29088820.866572164, 20*9.9766126515967311e-07, 20*7272205.21664304, 1., TRIANGULAR_VERTICAL},
  { 3, {0,4198609.6394310603}, {0,0}, 11255786.828366444, 9.9766126515967311e-07, 0.07*29088820.866572164, 20*9.9766126515967311e-07, 20*7272205.21664304, 1., TRIANGULAR_VERTICAL},
  { 3, {0,4198609.6394310603}, {0,0}, 11255786.828366444, 9.9766126515967311e-07, 0.03*29088820.866572164, 20*9.9766126515967311e-07, 20*7272205.21664304, 1., TRIANGULAR_VERTICAL},

//  { 3, {0,KLEN}, {0,0}, 2.62 * 4 * M_PI/3/AA, AA, 0.5 / AA, 20*AA, 160/AA, 1., TRIANGULAR_VERTICAL},
  { 3, {0,KLEN}, {0,0}, 2.62 * 4 * M_PI/3/AA, AA, 2.5 / AA, 20*AA, 160/AA, 1., TRIANGULAR_VERTICAL},
  { 3, {0,KLEN}, {0,0}, 2.62 * 4 * M_PI/3/AA, AA, 4.5 / AA, 20*AA, 160/AA, 1., TRIANGULAR_VERTICAL},
  { 3, {0,KLEN}, {0,0}, 2.62 * 4 * M_PI/3/AA, AA, 6.5 / AA, 20*AA, 160/AA, 1., TRIANGULAR_VERTICAL},
  { 3, {0,KLEN}, {0,0}, 2.62 * 4 * M_PI/3/AA, AA, 8.5 / AA, 20*AA, 160/AA, 1., TRIANGULAR_VERTICAL},
/*
  { 3, {0,Klen}, {0,0}, 2.62 * 4 * M_PI/3/a, a, 0.5 / a, 20*a, 2*M_PI*160/a, 1., TRIANGULAR_VERTICAL},
  { 3, {0,Klen}, {0,0}, 2.62 * 4 * M_PI/3/a, a, 2.5 / a, 20*a, 2*M_PI*160/a, 1., TRIANGULAR_VERTICAL},
  { 3, {0,Klen}, {0,0}, 2.62 * 4 * M_PI/3/a, a, 4.5 / a, 20*a, 2*M_PI*160/a, 1., TRIANGULAR_VERTICAL},
  { 3, {0,Klen}, {0,0}, 2.62 * 4 * M_PI/3/a, a, 6.5 / a, 20*a, 2*M_PI*160/a, 1., TRIANGULAR_VERTICAL},
  { 3, {0,Klen}, {0,0}, 2.62 * 4 * M_PI/3/a, a, 8.5 / a, 20*a, 2*M_PI*160/a, 1., TRIANGULAR_VERTICAL},
*/


  { 2, {1.1, 2.1}, {0.5,0.1325}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 2.1}, {0.5,0.1325}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 2.1}, {0.5,0.1325}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1, 2.1}, {0.5,0.1325}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},
/*
  { 2, {0, 3.1}, {0.5,0}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 3.1}, {0.5,0}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 3.1}, {0.5,0}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 3.1}, {0.5,0}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 2, {0, 1.1}, {0.5,0}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 1.1}, {0.5,0}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 1.1}, {0.5,0}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 1.1}, {0.5,0}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 2, {3.1,0}, {0,0.5}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1,0}, {0,0.5}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1,0}, {0,0.5}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1,0}, {0,0.5}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 2, {1.1,0}, {0,0.5}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1,0}, {0,0.5}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1,0}, {0,0.5}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1,0}, {0,0.5}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 2, {3.1,0}, {0,0.5}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1,0}, {0,0.5}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1,0}, {0,0.5}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1,0}, {0,0.5}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},
*/
  { 2, {3.1*0.5,-3.1*0.8}, {0.8,0.5}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1*0.5,-3.1*0.8}, {0.8,0.5}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1*0.5,-3.1*0.8}, {0.8,0.5}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1*0.5,-3.1*0.8}, {0.8,0.5}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 2, {1.1*0.5,-1.1*0.8}, {0.8,0.5}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1*0.5,-1.1*0.8}, {0.8,0.5}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1*0.5,-1.1*0.8}, {0.8,0.5}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1*0.5,-1.1*0.8}, {0.8,0.5}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  // Poloviční posun oproti přodchozímu
  { 2, {3.1*0.5,-3.1*0.8}, {0.4,0.25}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1*0.5,-3.1*0.8}, {0.4,0.25}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1*0.5,-3.1*0.8}, {0.4,0.25}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {3.1*0.5,-3.1*0.8}, {0.4,0.25}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 2, {1.1*0.5,-1.1*0.8}, {0.4,0.25}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1*0.5,-1.1*0.8}, {0.4,0.25}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1*0.5,-1.1*0.8}, {0.4,0.25}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {1.1*0.5,-1.1*0.8}, {0.4,0.25}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  // miniposun
  { 3, {3.1*0.5,-3.1*0.8}, {0.004,0.0025}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {3.1*0.5,-3.1*0.8}, {0.004,0.0025}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {3.1*0.5,-3.1*0.8}, {0.004,0.0025}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {3.1*0.5,-3.1*0.8}, {0.004,0.0025}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 3, {1.1*0.5,-1.1*0.8}, {0.004,0.0025}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {1.1*0.5,-1.1*0.8}, {0.004,0.0025}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {1.1*0.5,-1.1*0.8}, {0.004,0.0025}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {1.1*0.5,-1.1*0.8}, {0.004,0.0025}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 3, {3.1*0.5,-3.1*0.8}, {-0.004,-0.0025}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {3.1*0.5,-3.1*0.8}, {-0.004,-0.0025}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {3.1*0.5,-3.1*0.8}, {-0.004,-0.0025}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {3.1*0.5,-3.1*0.8}, {-0.004,-0.0025}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},

  { 3, {1.1*0.5,-1.1*0.8}, {-0.004,-0.0025}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {1.1*0.5,-1.1*0.8}, {-0.004,-0.0025}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {1.1*0.5,-1.1*0.8}, {-0.004,-0.0025}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 3, {1.1*0.5,-1.1*0.8}, {-0.004,-0.0025}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},




/*
  { 2, {0, 3.1}, {0,0.5}, 2.3, 0.97, 0.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 3.1}, {0,0.5}, 2.3, 0.97, 1.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 3.1}, {0,0.5}, 2.3, 0.97, 2.5, 20, 160, 1., TRIANGULAR_VERTICAL},
  { 2, {0, 3.1}, {0,0.5}, 2.3, 0.97, 3.5, 20, 160, 1., TRIANGULAR_VERTICAL},
*/




// end:
//  { 0,  {0, 0}, 0, 0, 0, 0, 0, 0, 0}
};

void ewaldtest_triang_results_free(ewaldtest_triang_results *r) {
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

ewaldtest_triang_results *ewaldtest_triang_3g(const ewaldtest_triang_params p);

int main() {
  gsl_set_error_handler(IgnoreUnderflowsGSLErrorHandler);
  for (size_t i = 0; i < sizeof(paramslist)/sizeof(ewaldtest_triang_params); ++i) {
    ewaldtest_triang_params p = paramslist[i];
    ewaldtest_triang_results *r = ewaldtest_triang_3g(p);
    // TODO print per-test header here
    printf("===============================\n");
    printf("a = %g, K = %g, Kmax = %g, Rmax = %g, lMax = %d, eta = %g, k = %g, beta = (%g,%g), ps = (%g,%g), csphase = %g\n",
        p.a, 4*M_PI/sqrt(3)/p.a, p.maxK, p.maxR, p.lMax, p.eta, p.k, p.beta.x, p.beta.y, p.particle_shift.x, p.particle_shift.y, p.csphase);
    printf("sigma0: %.16g%+.16gj\n", creal(r->sigma0), cimag(r->sigma0));    
    for (qpms_l_t n = 0; n <= p.lMax; ++n) {
      for (qpms_m_t m = -n; m <= n; ++m){
        if ((m+n)%2) continue;
        qpms_y_t y = qpms_mn2y_sc(m,n);
        qpms_y_t y_conj = qpms_mn2y_sc(-m,n);
        // y n m sigma_total (err), regsigmas_416 regsigmas_415_recon
        printf("%zd %d %d: T:%.16g%+.16gj(%.3g) L:%.16g%+.16gj(%.3g) S:%.16g%+.16gj(%.3g) \n"
            //"| predict %.16g%+.16gj \n| actual  %.16g%+.16gj\n"
            ,
            y, n, m, creal(san(r->sigmas_total[y])), san(cimag(r->sigmas_total[y])),
            r->err_sigmas_total[y],
            san(creal(r->sigmas_long[y])), san(cimag(r->sigmas_long[y])),
            r->err_sigmas_long[y],
            san(creal(r->sigmas_short[y])), san(cimag(r->sigmas_short[y])),
            r->err_sigmas_short[y]
            //san(creal(r->regsigmas_416[y])), san(cimag(r->regsigmas_416[y])),
            //san(creal(r->sigmas_total[y]) + creal(r->sigmas_total[y_conj])),
            //san(cimag(r->sigmas_total[y]) - cimag(r->sigmas_total[y_conj]))
        );
      }
    }
    ewaldtest_triang_results_free(r);
  }
  return 0;
}


int ewaldtest_counter = 0;


ewaldtest_triang_results *ewaldtest_triang_3g(const ewaldtest_triang_params p) {
  const double a = p.a; //const double a = p.h * sqrt(3);
  cart3_t beta3 = cart22cart3xy(p.beta);
  cart3_t particle_shift3 = cart22cart3xy(p.particle_shift);
 
  cart2_t b1, b2, rb1, rb2;
  if (p.orientation == TRIANGULAR_VERTICAL) {
    b1.x = 0;
    b1.y = a;
    b2.x = a * M_SQRT3 * .5;
    b2.y = a * .5;
  } else {
    b1.x = a;
    b1.y = 0;
    b2.x = a * .5;
    b2.y = a * M_SQRT3 * .5;
  }
  if (QPMS_SUCCESS != l2d_reciprocalBasis2pi(b1, b2, &rb1, &rb2))
    abort();

  const double A = l2d_unitcell_area(b1, b2); // sqrt(3) * a * a / 2.; // unit cell size
  const double K_len = cart2norm(rb1); //4*M_PI/a/sqrt(3); // reciprocal vector length

  
  ewaldtest_triang_results *results = malloc(sizeof(ewaldtest_triang_results));
  results->p = p;

  //triangular_lattice_gen_t *Rlg = triangular_lattice_gen_init(a, p.orientation, true, 0); // N.B. orig is included 
  //triangular_lattice_gen_extend_to_r(Rlg, p.maxR + a);
  //triangular_lattice_gen_t *Klg = triangular_lattice_gen_init(K_len, reverseTriangularLatticeOrientation(p.orientation), true, 0);
  //triangular_lattice_gen_extend_to_r(Klg, p.maxK + K_len);

  //point2d *Rpoints = Rlg->ps.base; 
  //size_t nR = Rlg->ps.r_offsets[Rlg->ps.nrs];
  
  /*if (fabs(p.particle_shift.x) ==0 && fabs(p.particle_shift.y) == 0) {
    points2d_rordered_t Rpos = points2d_rordered_annulus(&(Rlg->ps), 0, false, INFINITY, false);
    Rpoints = Rpos.base + Rpos.r_offsets[0];
    nR = Rpos.r_offsets[Rpos.nrs] - Rpos.r_offsets[0];
  }*/

  //point2d  *Kpoints = Klg->ps.base;
  //size_t nK = Klg->ps.r_offsets[Klg->ps.nrs];
  
  /*
  point2d particle_shift = p.particle_shift;
  point2d minus_ps = {-particle_shift.x, -particle_shift.y};
  point2d Rpoints_plus_shift[nR];
  for(size_t i = 0; i < nR; ++i){
    Rpoints_plus_shift[i].x = Rpoints[i].x - particle_shift.x;
    Rpoints_plus_shift[i].y = Rpoints[i].y - particle_shift.y;
  }
  */

 
  // skip zeroth point if it coincides with origin
  bool include_origin = !(fabs(p.particle_shift.x) == 0 
      && fabs(p.particle_shift.y) == 0);

  PGen Rlgen = PGen_xyWeb_new(b1, b2, BASIS_RTOL, CART2_ZERO, 0, include_origin, p.maxR + a, false);
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

  qpms_ewald32_constants_t *c = qpms_ewald32_constants_init(p.lMax, p.csphase);

  //points2d_rordered_t *Kpoints_plus_beta = points2d_rordered_shift(&(Klg->ps), p.beta,
  //    8*DBL_EPSILON, 8*DBL_EPSILON);

  //char filename[BUFSIZ];
  //sprintf(filename, "betalattice_%d.out", ewaldtest_counter);
  //dump_points2d_rordered(Kpoints_plus_beta, filename);


  if (0!=ewald3_sigma_long(results->sigmas_long,
        results->err_sigmas_long, c, p.eta, p.k, A,
        LAT_2D_IN_3D_XYONLY, &Klgen, false, beta3, particle_shift3))
      abort();
#if 0
  if (0!=ewald32_sigma_long_points_and_shift(results->sigmas_long,
        results->err_sigmas_long, c, p.eta, p.k, A,
        nK, Kpoints,
        p.beta, 
        particle_shift /*minus_ps*/ )) 
    abort();
#endif
  if (0!=ewald3_sigma_short(
    results->sigmas_short, results->err_sigmas_short, c,
    p.eta, p.k, LAT_2D_IN_3D_XYONLY, &Rlgen, false,  beta3, particle_shift3))
    abort();
  /*if (0!=ewald32_sigma_short_points_and_shift(
        results->sigmas_short, results->err_sigmas_short, c,
        p.eta, p.k,
        nR, Rpoints, p.beta, particle_shift)) 
    abort();*/
  //if (0!=ewald32_sigma0(&(results->sigma0), &(results->err_sigma0), c, p.eta, p.k))
  if (0!=ewald3_sigma0(&(results->sigma0), &(results->err_sigma0), c, p.eta, p.k))
    abort();
  for(qpms_y_t y = 0; y < nelem_sc; ++y) {
    results->sigmas_total[y] = results->sigmas_short[y] + results->sigmas_long[y];
    results->err_sigmas_total[y] = results->err_sigmas_short[y] + results->err_sigmas_long[y];
  }
  results->sigmas_total[0] += results->sigma0;
  results->err_sigmas_total[0] += results->err_sigma0;

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


  //points2d_rordered_free(Kpoints_plus_beta);
  qpms_ewald32_constants_free(c);
  //triangular_lattice_gen_free(Klg);
  //triangular_lattice_gen_free(Rlg);
  ++ewaldtest_counter;
  return results;
}
