#include <qpms/beyn.h>
#include <stdio.h>

#define CPAIR(x) creal(x), cimag(x)

int main(int argc, char **argv)
{
  double rRe = 2e3, rIm = 1.5e3, rounding = 0.2;
  complex double centre = 1e3 * I;
  size_t n = 100;

  beyn_contour_t *c = beyn_contour_kidney(centre, rRe, rIm, rounding, n, BEYN_CONTOUR_HALFELLIPSE_IM_PLUS);

  for(size_t i = 0; i < n; ++i) 
    printf("%g\t%g\t%g\t%g\n", CPAIR(c->z_dz[i][0]), CPAIR(c->z_dz[i][1]));

  free(c);
  return 0;
}
