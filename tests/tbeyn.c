#include <qpms/beyn.h>
#include <stdio.h>
#include <string.h>

// Matrix as in Beyn, section 4.11
int M_function(complex double *target, const size_t m, const complex double z, void *no_params) {
  complex double d =  2*m - 4*z / (6*m);
  complex double od =  -((double)m) - z / (6*m);

  memset(target, 0, m*m*sizeof(complex double));
  for (int i = 0; i < m; ++i) {
    target[i*m + i] = d;
    if(i > 0) target[i*m + i-1] = od;
    if(i < m - 1) target[i*m + i+1] = od;
  }
  target[m*(m-1) + m-1] = d/2 + z/(z-1);

  return 0;
}

int main() {
  complex double z0 = 150+2*I;
  double Rx = 148, Ry = 148;
  int L = 10, N = 50, dim = 400;
  beyn_contour_t *contour = beyn_contour_ellipse(z0, Rx, Ry, N);

  beyn_result_t *result =
    beyn_solve(dim, L, M_function, NULL /*M_inv_Vhat_function*/, NULL /*params*/,
      contour, 1e-4, 0);
  printf("Found %zd eigenvalues:\n", result->neig);
  for (size_t i = 0; i < result->neig; ++i) {
    complex double eig = result->eigval[i];
    printf("%zd: %g%+gj\n", i, creal(eig), cimag(eig));
  }
  free(contour);
  beyn_result_free(result);
  return 0;
}


