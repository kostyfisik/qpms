#include <qpms/wigner.h>
#include <stdio.h>

int main() {
  cart3_t v = {1, 2, 3};
  qpms_quat4d_t q4[8] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
    {-1, 0, 0, 0},
    {0, -1, 0, 0},
    {0, 0, -1, 0},
    {0, 0, 0, -1},
  };
  printf("original: (%g, %g, %g)\n", v.x, v.y, v.z);
  for(int i = 0; i < 8; ++i) {
    for (int det = -1; det < 2; det += 2) {
      const qpms_quat4d_t qr = q4[i];
      const qpms_quat_t qc = qpms_quat_2c_from_4d(qr);
      const qpms_irot3_t r = {qc, det};
      const cart3_t w = qpms_irot3_apply_cart3(r, v);
      printf("[%+g %+g %+g %+g] -> [%+g%+gj %+g%+g] *%+d: ==> (%g, %g, %g)\n",
          qr.c1, qr.ci, qr.cj, qr.ck, creal(qc.a), cimag(qc.a), creal(qc.b), cimag(qc.b),
          det, w.x, w.y, w.z);
    }
  }
  return 0;
}
