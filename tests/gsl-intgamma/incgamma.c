#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <qpms/ewald.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
//#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


int main(int argc, char **argv) {
  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
  int j;
  double x;
  while (scanf("%d %lf", &j, &x) == 2) {
    printf("%d %.16g", j, x);
    qpms_csf_result res;
    complex double argfac = 1;
    for (int i = 0; i < 4; ++i, argfac *= I) {
      int retval = complex_gamma_inc_e(0.5-j, argfac * x, &res);
      printf(" | %.16g+%.16gj %.4g %d", creal(res.val), cimag(res.val), res.err, retval);
    }
    putchar('\n');
  }
  return 0;
}


