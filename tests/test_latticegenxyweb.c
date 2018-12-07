//c99 -o test_latticegenxyweb -ggdb -Wall -I ../ test_latticegenxyweb.c ../qpms/latticegens.c -lm
#define QPMS_VECTORS_NICE_TRANSFORMATIONS
#include <qpms/lattices.h>
#include <stdio.h>

void fprint_PGenCart2ReturnData(FILE *f, PGenCart2ReturnData d) {
  cart2_t c = d.point_cart2;
  fprintf(f, "%g\t%g\tflags: %#xd\n",
      c.x, c.y, d.flags);
}

void dump_PGenCart2(char *filename, PGen *g) {
  FILE *out = fopen(filename, "w");
  PGenCart2ReturnData d;
  while ((d = PGen_next_cart2(g)).flags & PGEN_NOTDONE) {
    fprint_PGenCart2ReturnData(out, d);
  } 
}

#if 0
void print_PGenSphReturnData(PGenSphReturnData d) {
  sph_t s = d.point_sph;
  cart3_t c = sph2cart(s);
  printf("sph: (%g, %g π, %g π), cart: (%g, %g, %g), flags: %#xd\n",
      s.r, s.theta / M_PI, s.phi / M_PI, c.x, c.y, c.z, d.flags);
}

void dump_PGenSph(PGen *g) {
  PGenSphReturnData d;
  do {
    d = PGen_next_sph(g);
    print_PGenSphReturnData(d);
  } while (d.flags & PGEN_NOTDONE);
}
#endif

#define DO_AND_PRINT(label, s) printf(#label ":\n" #s "\n"); s ;

int main(int argc, char **argv) {
  PGen g;
  cart2_t b1 = {1.1,0}, b2 = {0,-1.3}, offset = {0.1, -0.1};
  g = PGen_xyWeb_new(b1, b2, 1e-13, offset, 0, true, 5, false);
  dump_PGenCart2("rect1.xydump", &g);
  g = PGen_xyWeb_new(b1, b2, 1e-13, offset, 5, true, 8, false);
  dump_PGenCart2("rect2.xydump", &g);
  b1.x = 1.342, b1.y = 4.3121; b2.x = -1.83, b2.y = 1.4;
  g = PGen_xyWeb_new(b1, b2, 1e-13, offset, 0, true, 8, false);
  dump_PGenCart2("oblique1.xydump", &g);

  return 0;


}

