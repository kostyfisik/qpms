#include <stddef.h>
#include <gsl/gsl_sf_legendre.h>
#include <stdio.h>

int main() {
  const size_t lmax = 2;
  const double x = .5;
  size_t arrsiz = gsl_sf_legendre_array_n(lmax);
  double csphase = 1;
  double target[arrsiz];
  printf("lmax = %zd, x = %g, csphase = %g:\n", lmax, x, csphase);
  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE, lmax, x, csphase, target);
  for(int l = 0; l <= lmax; ++l) for (int m = 0; m <= l; ++m)
    printf("P_%d^%d(%g)\t= %g\n", l, m, x, target[gsl_sf_legendre_array_index(l,m)]);
  csphase = -1;
  printf("lmax = %zd, x = %g, csphase = %g:\n", lmax, x, csphase);
  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE, lmax, x, csphase, target);
  for(int l = 0; l <= lmax; ++l) for (int m = 0; m <= l; ++m)
    printf("P_%d^%d(%g)\t= %g\n", l, m, x, target[gsl_sf_legendre_array_index(l,m)]);
  return 0;
}

