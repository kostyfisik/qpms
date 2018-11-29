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

void dump_PGenSph(PGenSph *g) {
  PGenSphReturnData d;
  do {
    d = PGenSph_next(g);
    print_PGenSphReturnData(d);
  } while (d.flags & PGEN_NOTDONE);
}

#define DO_AND_PRINT(label, s) printf(#label ":\n" #s "\n"); s ;

int main(int argc, char **argv) {
  PGenSph g;
  DO_AND_PRINT(test1a, g = PGenSph_zAxis_new_minMaxR(0.2, 0.14, 5, true,  7, true, PGENSPH_ZAXIS_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test1b, g = PGenSph_zAxis_new_minMaxR(0.2, 0.14, 5, true,  7, true, PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test2a, g = PGenSph_zAxis_new_minMaxR(0.2, 0.05, 5.05, true, 7.05, true, PGENSPH_ZAXIS_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test2b, g = PGenSph_zAxis_new_minMaxR(0.2, 0.05, 5.05, true, 7.05, true, PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test3a, g = PGenSph_zAxis_new_minMaxR(0.2, 0.05, 5.05, false, 7.05, false, PGENSPH_ZAXIS_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test3b, g = PGenSph_zAxis_new_minMaxR(0.2, 0.05, 5.05, false, 7.05, false, PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test4a, g = PGenSph_zAxis_new_minMaxR(0.2, 0.0, 0, false, 1, false, PGENSPH_ZAXIS_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test4b, g = PGenSph_zAxis_new_minMaxR(0.2, 0.0, 0, false, 1, false, PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test5a, g = PGenSph_zAxis_new_minMaxR(0.2, 0.0, 0, true, 1, true, PGENSPH_ZAXIS_INC_FROM_ORIGIN))
  dump_PGenSph(&g);
  DO_AND_PRINT(test5b, g = PGenSph_zAxis_new_minMaxR(0.2, 0.0, 0, true, 1, true, PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN))
  dump_PGenSph(&g);


}

