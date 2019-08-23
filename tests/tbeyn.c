#include <qpms/beyn.h>
#include <stdio.h>


// Matrix as in Beyn, section 4.11
int M_function(gsl_matrix_complex *target, complex double z, void *no_params) {
  int m = target->size1;

  gsl_complex d = gsl_complex_fromstd(  2*m - 4*z / (6*m) );
  gsl_complex od =  gsl_complex_fromstd( -m - z / (6*m) );

  gsl_matrix_complex_set_zero(target);
  for (int i = 0; i < m; ++i) {
    gsl_matrix_complex_set(target, i, i, (gsl_complex) d);
    if(i > 0) gsl_matrix_complex_set(target, i, i-1, (gsl_complex) od);
    if(i < m - 1) gsl_matrix_complex_set(target, i, i+1, (gsl_complex) od);
  }

  return 0;
}

int main() {
  complex double z0 = 150+2*I;
  double Rx = 148, Ry = 148;
  int L = 10, N = 50, dim = 400;
  BeynSolver * solver = CreateBeynSolver(dim, L);

  int K = BeynSolve(solver, M_function, NULL /*M_inv_Vhat_function*/, NULL /*params*/,
      z0, Rx, Ry, N);
  printf("Found %d eigenvalues:\n", K);
  for (int i = 0; i < K; ++i) {
    gsl_complex eig = gsl_vector_complex_get(solver->Eigenvalues, i);
    printf("%d: %g%+gj\n", i, GSL_REAL(eig), GSL_IMAG(eig));
  }
  DestroyBeynSolver(solver);
  return 0;
}


