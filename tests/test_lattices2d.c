// c99 -o l2dtest -ggdb -I ../ test_lattices2d.c ../qpms/lattices2d.c -lm
#include <qpms/lattices.h>
#include <stdio.h>

int main(int argc, char **argv) {
  cart2_t b1 = {1.2,0.1};
  cart2_t b2 = {0.5, 0.5};
  cart2_t b3;

  printf("original b1 = (%g, %g), b2 = (%g, %g)\n",
      b1.x, b1.y, b2.x, b2.y);
  printf("Inscribed circle radius: %g\n",
      l2d_hexWebInCircleRadius(b1, b2));
  l2d_shortestBase3(b1, b2, &b1, &b2, &b3);
  printf("shortestBase3: b1 = (%g, %g), b2 = (%g, %g), b3 = (%g, %g)\n",
      b1.x, b1.y, b2.x, b2.y, b3.x, b3.y);
  return 0;
}
