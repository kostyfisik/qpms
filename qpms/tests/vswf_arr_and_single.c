// c99 -I .. vswf_single.c -lqpms -lgsl -lblas
#include "indexing.h"
#include "vswf.h"
#include <stdio.h>

int main(int argc, char **argv) {
  const qpms_l_t lMax = 2;
  csph_t point = {1.4554, 1.2424, 4.2545};

  qpms_vswf_set_spec_t *bspec = qpms_vswf_set_spec_from_lMax(lMax, QPMS_NORMALISATION_CONVENTION_KRISTENSSON);

  csphvec_t arr[bspec->n];
  qpms_uvswf_fill(arr, bspec, point, QPMS_BESSEL_REGULAR);

  for(size_t i = 0; i < bspec->n; i++)  {
    csphvec_t v = arr[i];
    qpms_vswf_type_t t; qpms_m_t m; qpms_l_t l;
    qpms_uvswfi2tmn(bspec->ilist[i], &t, &m, &l);
    printf("arr; l=%d,m=%+d,t=%d, @(%g,%g,%g): (%g%+gj, %g%+gj, %g%+gj)\n",
        l, m, t, creal(point.r), point.theta, point.phi, creal(v.rc), cimag(v.rc),
        creal(v.thetac), cimag(v.thetac), creal(v.phic), cimag(v.phic));
    if (t == QPMS_VSWF_ELECTRIC)
      v = qpms_vswf_single_el_csph(m, l, point,QPMS_BESSEL_REGULAR,bspec->norm);
    else
      v = qpms_vswf_single_mg_csph(m, l, point,QPMS_BESSEL_REGULAR,bspec->norm);
    printf("sgl; l=%d,m=%+d,t=%d, @(%g,%g,%g): (%g%+gj, %g%+gj, %g%+gj)\n",
        l, m, t, creal(point.r), point.theta, point.phi, creal(v.rc), cimag(v.rc),
        creal(v.thetac), cimag(v.thetac), creal(v.phic), cimag(v.phic));
  }
  qpms_vswf_set_spec_free(bspec);
  return 0;
}


