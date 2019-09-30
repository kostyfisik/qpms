#define _GNU_SOURCE
#include <qpms/beyn.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

static double randU(double a, double b) {return a + (b-a) * random() * (1. / RAND_MAX); }
// Normal distribution via Box-Muller transform
static double randN(double sigma, double mu) {
  double u1 = randU(0,1);
  double u2 = randU(0,1);
  return mu + sigma*sqrt(-2*log(u1))*cos(2.*M_PI*u2);
}

struct param {
  double *T0;
  double *T1;
  double *T2;
};

int M_function(complex double *target, const size_t m, const complex double z, void *params) {
  struct param *p = params;

  for(size_t i = 0; i < m*m; ++i)
    target[i] = p->T0[i] + z*p->T1[i] + cexp(
#ifdef VARIANTB // Also note that this case requires pretty large contour point number (>~ 3000)
        (1+3*I)
#else // VARIANTA or VARIANTC
        (1+1*I)
#endif
        *z*p->T2[i]) + 
#ifdef VARIANTC // Essential singularity at zero
      cexp(3/z);
#elif defined VARIANTD // Essential singularity outside the contour
      cexp(3/(z-1))
#elif defined VARIANTE // High-order pole at zero
      3/cpow(z,10);
#elif defined VARIANTF // High-order pole at zero, higher order than dim
      .0003/cpow(z,12);
#else // double pole at zero
      3/z/z
#endif
      ;
  return 0;
}

int main(int argc, char **argv) {
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  complex double z0 = 0+3e-1*I;
#ifdef RXSMALL
  double Rx = .1;
#else 
  double Rx = .3; // Variant B will fail in this case due to large number of eigenvalues (>30)
#endif
  double Ry = .25;
#ifdef VARIANTF
  int L = 10, N = 150, dim = 10;
#else
  int L = 30, N = 150, dim = 60;
#endif
  if (argc > 1) N = atoi(argv[1]);
  if (argc > 2) L = atoi(argv[2]);
#ifdef IMPLUS
  beyn_contour_t *contour = beyn_contour_halfellipse(z0, Rx, Ry, N, BEYN_CONTOUR_HALFELLIPSE_IM_PLUS);
#elif defined IMPLUS_KIDNEY
  beyn_contour_t *contour = beyn_contour_kidney(z0, Rx, Ry, 0.3, N, BEYN_CONTOUR_HALFELLIPSE_IM_PLUS);
#else
  beyn_contour_t *contour = beyn_contour_ellipse(z0, Rx, Ry, N);
#endif
  struct param p;
  p.T0 = malloc(dim*dim*sizeof(double));
  p.T1 = malloc(dim*dim*sizeof(double));
  p.T2 = malloc(dim*dim*sizeof(double));
  for(size_t i = 0; i < dim*dim; ++i) {
    p.T0[i] = randN(1,0);
    p.T1[i] = randN(1,0);
    p.T2[i] = randN(1,0);
  }

  beyn_result_t *result =
    beyn_solve(dim, L, M_function, NULL /*M_inv_Vhat_function*/, &p /*params*/,
      contour, 1e-4, 1e-4);
  printf("Found %zd eigenvalues:\n", result->neig);
  for (size_t i = 0; i < result->neig; ++i) {
    complex double eig = result->eigval[i];
    printf("%zd: %g%+gj\n", i, creal(eig), cimag(eig));
  }
  free(contour);
  beyn_result_free(result);
  free(p.T0);
  free(p.T1);
  free(p.T2);
  return 0;
}


