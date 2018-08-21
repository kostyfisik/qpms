#include <qpms/ewald.h>
#include <qpms/indexing.h>
#include <assert.h>
#include <stdio.h>
#define LMAX 10

int main()
{
  qpms_ewald32_constants_t *c = qpms_ewald32_constants_init(LMAX);
  for (qpms_l_t l = 1; l <= LMAX; ++l)
    for (qpms_m_t m = -l; m <= l; ++m) {
      printf("%d %d: (", l, m);
      qpms_y_t y = qpms_mn2y(m, l);
      if (0 == (l+m)%2) {
        for (int j = 0; j <= c->s1_jMaxes[y]; ++j)
          printf("%.16g %+.16gj, ", creal(c->s1_constfacs[y][j]), 
              cimag(c->s1_constfacs[y][j]));
      }
      else assert(c->s1_constfacs[y] == NULL);
      printf(")\n");
    }
  return 0;
}



