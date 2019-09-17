#include <qpms/beyn.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static double randU(double a, double b) {return a + (b-a) * random() * (1. / RAND_MAX); }
// Normal distribution via Box-Muller transform
static double randN(double sigma, double mu) {
  double u1 = randU(0,1);
  double u2 = randU(0,1);
  return mu + sigma*sqrt(-2*log(u1))*cos(2.*M_PI*u2);
}

struct param49 {
  double *T0;
  double *T1;
  double *T2;
};

// Matrix as in Beyn, example 4.9
int M_function(complex double *target, const size_t m, const complex double z, void *params) {
  struct param49 *p = params;

  for(size_t i = 0; i < m*m; ++i)
    target[i] = p->T0[i] + z*p->T1[i] + z*z*p->T2[i];
  return 0;
}

int main(int argc, char **argv) {
  complex double z0 = 0;
  double Rx = .3, Ry = .3;
  int L = 30, N = 150, dim = 60;
  if (argc > 1) N = atoi(argv[1]);
  beyn_contour_t *contour = beyn_contour_ellipse(z0, Rx, Ry, N);
  struct param49 p;
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


