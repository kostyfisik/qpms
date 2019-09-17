#include <qpms/beyn.h>
#include <stdio.h>


// Matrix as in Beyn, section 4.11
int M_function(gsl_matrix_complex *target, complex double z, void *no_params) {
  int m = target->size1;

  gsl_complex d = gsl_complex_fromstd(  2*m - 4*z / (6*m) );
  gsl_complex od =  gsl_complex_fromstd( -(double)m - z / (6*m) );

  gsl_matrix_complex_set_zero(target);
  for (int i = 0; i < m; ++i) {
    gsl_matrix_complex_set(target, i, i, d);
    if(i > 0) gsl_matrix_complex_set(target, i, i-1, od);
    if(i < m - 1) gsl_matrix_complex_set(target, i, i+1, od);
  }
  gsl_matrix_complex_set(target, m-1, m-1, gsl_complex_fromstd(gsl_complex_tostd(d)/2 + z/(z-1)));

  return 0;
}

int main() {
  complex double z0 = 150+2*I;
  double Rx = 148, Ry = 148;
  int L = 10, N = 50, dim = 400;
  beyn_contour_t *contour = beyn_contour_ellipse(z0, Rx, Ry, N);

  beyn_result_gsl_t *result =
    beyn_solve_gsl(dim, L, M_function, NULL /*M_inv_Vhat_function*/, NULL /*params*/,
      contour, 1e-4, 1e-4);
  printf("Found %zd eigenvalues:\n", result->neig);
  for (size_t i = 0; i < result->neig; ++i) {
    gsl_complex eig = gsl_vector_complex_get(result->eigval, i);
    printf("%zd: %g%+gj\n", i, GSL_REAL(eig), GSL_IMAG(eig));
  }
  free(contour);
  beyn_result_gsl_free(result);
  return 0;
}


