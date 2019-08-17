//c99 -o test_latticegens -ggdb -Wall -I ../ test_latticegens.c ../qpms/latticegens.c -lm
#define QPMS_VECTORS_NICE_TRANSFORMATIONS
#include <qpms/lattices.h>
#include <stdio.h>

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

#define DO_AND_PRINT(label, s) printf(#label ":\n" #s "\n"); s ;

int main(int argc, char **argv) {
  PGen g;
  DO_AND_PRINT(test1a, g = PGen_1D_new_minMaxR(0.2, 0.14, 5, true,  7, true, PGEN_1D_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test1b, g = PGen_1D_new_minMaxR(0.2, 0.14, 5, true,  7, true, PGEN_1D_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test2a, g = PGen_1D_new_minMaxR(0.2, 0.05, 5.05, true, 7.05, true, PGEN_1D_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test2b, g = PGen_1D_new_minMaxR(0.2, 0.05, 5.05, true, 7.05, true, PGEN_1D_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test3a, g = PGen_1D_new_minMaxR(0.2, 0.05, 5.05, false, 7.05, false, PGEN_1D_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test3b, g = PGen_1D_new_minMaxR(0.2, 0.05, 5.05, false, 7.05, false, PGEN_1D_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test4a, g = PGen_1D_new_minMaxR(0.2, 0.0, 0, false, 1, false, PGEN_1D_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test4b, g = PGen_1D_new_minMaxR(0.2, 0.0, 0, false, 1, false, PGEN_1D_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test5a, g = PGen_1D_new_minMaxR(0.2, 0.0, 0, true, 1, true, PGEN_1D_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test5b, g = PGen_1D_new_minMaxR(0.2, 0.0, 0, true, 1, true, PGEN_1D_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);


}

