//  c99 -g -I.. bspectransl.c ../qpms/translations.c ../qpms/legendre.c ../qpms/vswf.c ../qpms/gaunt.c ../qpms/error.c -lm -lgsl -lblas
#include <qpms/translations.h>
#include <qpms/vswf.h>
#include <complex.h>
#include <stdio.h>

int main() {
  cart3_t pos1={0,1,2}, pos2={3,5,2};
  qpms_vswf_set_spec_t 
    *bspec1 = qpms_vswf_set_spec_from_lMax(1, QPMS_NORMALISATION_POWER),
    *bspec2 = qpms_vswf_set_spec_from_lMax(2, QPMS_NORMALISATION_POWER);

  qpms_trans_calculator *c = qpms_trans_calculator_init(3, QPMS_NORMALISATION_POWER);

  complex double s_2_1[bspec2->n][bspec1->n];
  complex double s_1_2[bspec1->n][bspec2->n];
  const double k = 1.8;

  qpms_trans_calculator_get_trans_array_lc3p(c, s_2_1[0], bspec2, bspec1->n,
      bspec1, 1, k, pos2, pos1);

  qpms_trans_calculator_get_trans_array_lc3p(c, s_1_2[0], bspec1, bspec2->n,
      bspec2, 1, k, pos1, pos2);

  for(size_t R = 0; R < bspec2->n; ++R) {
    for(size_t C = 0; C < bspec1->n; ++C)
      printf("%.3lg+%.3lgj\t", creal(s_2_1[R][C]), cimag(s_2_1[C][R]));
    putchar('\n');
  }
  putchar('\n');

  for(size_t R = 0; R < bspec1->n; ++R) {
    for(size_t C = 0; C < bspec2->n; ++C)
      printf("%.3lg+%.3lgj\t", creal(s_1_2[R][C]), cimag(s_1_2[R][C]));
    putchar('\n');
  }

  qpms_trans_calculator_free(c);
  qpms_vswf_set_spec_free(bspec1);
  qpms_vswf_set_spec_free(bspec2);
  return 0;
}
