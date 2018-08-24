#include <qpms/lattices.h>
#include <stdio.h>

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
  triangular_lattice_gen_t *g = triangular_lattice_gen_init(5, TRIANGULAR_HORIZONTAL, false,0);
  dump_points2d_rordered(&(g->ps), "triang_h_empty.out");

  points2d_rordered_t *p = points2d_rordered_scale(&(g->ps), 1.5464);
  dump_points2d_rordered(p, "triang_h_empty_scaled.out");

  points2d_rordered_free(p);
  triangular_lattice_gen_free(g);

  g = triangular_lattice_gen_init(5, TRIANGULAR_HORIZONTAL, false,0);
  triangular_lattice_gen_extend_to_steps(g, 5);
  dump_points2d_rordered(&(g->ps), "triang_h_s5.out");
  triangular_lattice_gen_extend_to_steps(g, 20);
  dump_points2d_rordered(&(g->ps), "triang_h_s20.out");
  triangular_lattice_gen_extend_to_steps(g, 160);
  dump_points2d_rordered(&(g->ps), "triang_h_s160.out");
  triangular_lattice_gen_extend_to_steps(g, 20);
  dump_points2d_rordered(&(g->ps), "triang_h_s160_20.out");
  p = points2d_rordered_scale(&(g->ps), -.2);
  dump_points2d_rordered(p, "triang_h_s160_20_scaled.out");
  points2d_rordered_free(p);
  triangular_lattice_gen_free(g);


  g = triangular_lattice_gen_init(7, TRIANGULAR_VERTICAL, true,1);
  triangular_lattice_gen_extend_to_steps(g, 7);
  dump_points2d_rordered(&(g->ps), "triang_v_plus_s7.out");
  p = points2d_rordered_scale(&(g->ps), 1/7.);
  triangular_lattice_gen_free(g);
  dump_points2d_rordered(p, "triang_v_plus_s7_scaled_out");
  points2d_rordered_free(p);
  g = triangular_lattice_gen_init(7, TRIANGULAR_VERTICAL, true,-1);
  triangular_lattice_gen_extend_to_steps(g, 7);
  dump_points2d_rordered(&(g->ps), "triang_v_minus_s7.out");
  p = points2d_rordered_scale(&(g->ps), 1/7.);
  triangular_lattice_gen_free(g);
  dump_points2d_rordered(p, "triang_v_minus_s7_scaled_out");
  points2d_rordered_free(p);

  return 0;
}

  

