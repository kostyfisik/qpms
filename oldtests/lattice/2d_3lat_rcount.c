#include <qpms/lattices.h>
#include <stdio.h>
#include <float.h>

void dump_points2d_rordered(const points2d_rordered_t *ps, char *filename) {
  FILE *f = fopen(filename, "w");
  for (size_t i = 0; i < ps->nrs; ++i) {
    fprintf(f, "# r = %.16g\n", ps->rs[i]);
    for (ptrdiff_t j = ps->r_offsets[i]; j < ps->r_offsets[i+1]; ++j)
      fprintf(f, "%.16g %.16g\n", ps->base[j].x, ps->base[j].y);
  }
  fclose(f);
}

int main() {
  triangular_lattice_gen_t *g = triangular_lattice_gen_init(1, TRIANGULAR_HORIZONTAL, false,0);
  triangular_lattice_gen_extend_to_steps(g, 1000);
  for(size_t i = 0; i < g->ps.nrs; ++i) {
    printf("%zd %.16g %td\n", i, g->ps.rs[i], g->ps.r_offsets[i+1]);
  }
  triangular_lattice_gen_free(g);

  return 0;
}

  

