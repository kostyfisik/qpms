#include "ewald.h"
#include <stdlib.h>
#include "indexing.h"
#include <assert.h>
#include "tiny_inlines.h"

// sloppy implementation of factorial
static inline double factorial(const int n) {
  assert(n >= 0);
  if (n < 0)
    return 0; // should not happen in the functions below. (Therefore the assert above)
  else if (n <= 20) {
    double fac = 1;
    for (int i = 1; i <= n; ++i)
      fac *= i;
    return fac;
  }
  else 
    return tgamma(n + 1); // hope it's precise and that overflow does not happen
}



qpms_ewald32_constants_t *qpms_ewald32_constants_init(const qpms_l_t lMax)
{
  qpms_ewald32_constants_t *c = malloc(sizeof(qpms_ewald32_constants_t));
  //if (c == NULL) return NULL; // Do I really want to do this?
  c->lMax = lMax;
  c->nelem = qpms_lMax2nelem(lMax);
  c->s1_jMaxes = malloc(c->nelem * sizeof(qpms_l_t));
  c->s1_constfacs = malloc(c->nelem * sizeof(complex double *));
  //if (c->s1_jMaxes == NULL) return NULL;

  //determine sizes
  size_t s1_constfacs_sz = 0;
  for (qpms_y_t y = 0; y < c->nelem; ++y) {
    qpms_l_t n; qpms_m_t m; qpms_y2mn_p(y, &m, &n);
    if ((m + n) % 2 == 0) 
      s1_constfacs_sz += 1 + (c->s1_jMaxes[y] = (n-abs(m))/2);
    else
      c->s1_jMaxes[y] = -1;
  }

  c->s1_constfacs[0];
  c->s1_constfacs_base = malloc(c->nelem * sizeof(complex double));
  size_t s1_constfacs_sz_cumsum = 0;
  for (qpms_y_t y = 0; y < c->nelem; ++y) {
    qpms_l_t n; qpms_m_t m; qpms_y2mn_p(y, &m, &n);
    if ((m + n) % 2 == 0) {
      c->s1_constfacs[y] = c->s1_constfacs_base + s1_constfacs_sz_cumsum;
      // and here comes the actual calculation
      for (qpms_l_t j = 0; j <= c->s1_jMaxes[y]; ++j){
        c->s1_constfacs[y][j] = -0.5 * ipow(n+1) * min1pow((n+m)/2) 
          * sqrt((2*n + 1) * factorial(n-m) * factorial(n+m))
          * min1pow(j) * pow(0.5, n-2*j)
          / (factorial(j) * factorial((n-m)/2-j) * factorial((n+m)/2-j))
          * pow(0.5, 2*j-1);
      }
      s1_constfacs_sz_cumsum += 1 + c->s1_jMaxes[y];
    }
    else
      c->s1_constfacs[y] = NULL;
  }
  return c;
}

void qpms_ewald32_constants_free(qpms_ewald32_constants_t *c) {
  free(c->s1_constfacs);
  free(c->s1_constfacs_base);
  free(c->s1_jMaxes);
  free(c);
}

