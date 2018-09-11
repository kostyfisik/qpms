// c99 -I ../.. ewald2F2.c ../../qpms/ewaldsf.c -lm -lgsl -lblas
#include <gsl/gsl_sf_result.h>
#include <qpms/ewald.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>


int main(int argc, char **argv) {
  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
  double a, b, c, d, x;
  while (scanf("%lf %lf %lf %lf %lf", &a, &b, &c, &d, &x) == 5) {
    printf("%.16g %.16g %.16g %.16g %.16g", a, b, c, d, x);
    gsl_sf_result res;
    int retval = hyperg_2F2_series(a, b, c, d, x, &res);
    printf(" | %.16g (%.3g) %d\n", res.val,  res.err, retval);
  }
  return 0;
}


