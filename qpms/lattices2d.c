#include "lattices.h"
#include <assert.h>

typedef struct {
  int i, j;
} intcoord2_t;

static inline int sqi(int x) { return x*x; }


/*
 * EQUILATERAL TRIANGULAR LATTICE
 */


/*
 * N. B. the possible radii (distances from origin) of the lattice points can be described as
 *
 *     r**2 / a**2 == i**2 + j**2 + i*j ,
 *
 * where i, j are integer indices describing steps along two basis vectors (which have
 * 60 degree angle between them).
 *
 * The plane can be divided into six sextants, characterized as:
 *
 * 0) i >= 0 && j >= 0,
 *    [a] i > 0,
 *    [b] j > 0,
 * 1) i <= 0 && {j >= 0} && i + j >= 0,
 *    [a] i + j > 0,
 *    [b] i < 0,
 * 2) {i <= 0} && j >= 0 && i + j <= 0,
 *    [a] j > 0,
 *    [b] i + j < 0,
 * 3) i <= 0 && j <= 0,
 *    [a] i < 0,
 *    [b] j < 0,
 * 4) i >= 0 && {j <= 0} && i + j <= 0,
 *    [a] i + j < 0,
 *    [b] i > 0,
 * 5) {i >= 0} && j <= 0 && i + j >= 0,
 *    [a] j < 0,
 *    [b] i + j > 0.
 *
 * The [a], [b] are two variants that uniquely assign the points at the sextant boundaries.
 * The {conditions} in braces are actually redundant.
 *
 * In each sextant, the "minimum steps from the origin" value is calculated as:
 * 0) i + j,
 * 1) j
 * 2) -i
 * 3) -i - j,
 * 4) -j,
 * 5) i.
 *
 * The "spider web" generation for s steps from the origin (s-th layer) goes as following (variant [a]):
 * 0) for (i = s, j = 0; i > 0; --i, ++j)
 * 1) for (i = 0, j = s; i + j > 0; --i)
 * 2) for (i = -s, j = s; j > 0; --j)
 * 3) for (i = -s, j = 0; i < 0; ++i, --j)
 * 4) for (i = 0, j = -s; i + j < 0; ++i)
 * 5) for (i = s, j = -s; j < 0; ++j)
 *
 *
 * Length of the s-th layer is 6*s for s >= 1. Size (number of lattice points) of the whole s-layer "spider web"
 * is therefore 3*s*(s+1), excluding origin.
 * The real area inside the web is (a*s)**2 * 3 * sqrt(3) / 2.
 * Area of a unit cell is a**2 * sqrt(3)/2.
 * Inside the web, but excluding the circumscribed circle, there is no more
 * than 3/4.*s*(s+1) lattice cells (FIXME pretty stupid but safe estimate).
 *
 * s-th layer circumscribes a circle of radius a * s * sqrt(3)/2.
 *
 */

static inline int trilat_r2_ij(const int i, const int j) {
  return sqi(i) + sqi(j) + i*j;
}

static inline int trilat_r2_coord(const intcoord2_t c) {
  return trilat_r2_ij(c.i, c.j);
}

static int trilat_cmp_intcoord2_by_r2(const void *p1, const void *p2) {
  // CHECK the sign is right
  return trilat_r2_coord(*(const intcoord2_t *)p1) - trilat_r2_coord(*(const intcoord2_t *)p2);
}

// Classify points into sextants (variant [a])
static int trilat_sextant_ij_a(const int i, const int j) {
  const int w = i + j;
  if (i >  0 && j >= 0) return  0;
  if (i <= 0 && w >  0) return  1;
  if (w <= 0 && j >  0) return  2;
  if (i <  0 && j <= 0) return  3;
  if (i >= 0 && w <  0) return  4;
  if (w >= 0 && j <  0) return  5;
  if (i == 0 && j == 0) return -1; // origin
  assert(0); // other options should be impossible
}

typedef struct {
  intcoord2_t *pointlist_base; // allocated memory for the point "buffer"
  size_t pointlist_capacity; 
  // beginning and end of the point "buffer"
  // not 100% sure what type should I use here
  // (these are both relative to pointlist_base, due to possible realloc's)
  ptrdiff_t pointlist_beg, pointlist_end; 
  int maxs; // the highest layer of the spider web generated (-1 by init, 0 is only origin (if applicable))
  // capacities of the arrays in ps
  size_t ps_rs_capacity;
  size_t ps_points_capacity; // this is the "base" array
  // TODO anything else?
} triangular_lattice_gen_t_privstuff_t;

triangular_lattice_gen_t * triangular_lattice_gen_init(double a, TriangularLatticeOrientation ori, bool include_origin)
{
  triangular_lattice_gen_t *g = malloc(sizeof(triangular_latice_gen_t));
  g->orientation = ori;
  g->includes_origin = include_origin;
  g->ps.nrs = 0;
  g->ps.rs = NULL;
  g->ps.base = NULL;
  g->ps.r_offsets = NULL;
  g->priv = malloc(sizeof(triangular_lattice_gen_privstuff_t));
  g->priv->maxs = -1;
  g->priv->pointlist_capacity = 0;
  g->priv->pointlist_base = NULL;
  g->priv->pointlist_beg = 0;
  g->priv->pointlist_end = 0;
  g->priv->ps_rs_capacity = 0;
  g->priv->ps_points_capacity = 0;
  return g;
}

void triangular_lattice_gen_free(triangular_lattice_get_t *g) {
  free(g->ps.rs);
  free(g->ps.base);
  free(g->ps.r_offsets);
  free(g->priv->pointlist_base);
  free(g->priv);
  free(g);
}

const points2d_reordered_t * triangular_lattice_gen_getpoints(const triangular_lattice_generator_t *g) {
  return &(g->ps);
}

int triangular_lattice_gen_extend_to_steps(triangular_lattice_generator_t * g, int maxsteps) {
  if (maxsteps <= g->priv->maxs)  // nothing needed
    return 0;
  // TODO FIXME: check for maximum possible maxsteps (not sure what it is)
  int err;
  err = trilatgen_ensure_pointlist_capacity(g, maxsteps);
  if(err) return err;
  err = trilatgen_ensure_ps_rs_capacity(g, maxsteps);
  if(err) return err;
  err = trilatgen_ensure_ps_points_capacity(g, maxsteps);
  if(err) return err;
  
  if(g->includes_origin && g->priv->maxs < 0) // Add origin if not there yet
    trilatgen_append_ij(g, 0, 0);
  
  for (int s = g->priv->maxs + 1; s <= maxsteps; ++s) {
    int i, j; 
    // now go along the spider web layer as indicated in the lenghthy comment above
    for (i = s, j = 0; i > 0; --i; ++j) trilatgen_append_ij(g,i,j);
    for (i = 0, j = s; i + j > 0; --i) trilatgen_append_ij(g,i,j);
    for (i = -s, j = s; j > 0; --j) trilatgen_append_ij(g,i,j);
    for (i = -s, j = 0; i < 0; ++i, --j) trilatgen_append_ij(g,i,j);
    for (i = 0, j = -s; i + j < 0; ++i) trilatgen_append_ij(g,i,j);
    for (i = s, j = -s; j < 0; ++j) trilatgen_append_ij(g,i,j);
  }

  trilatgen_sort_pointlist(g);
  
  // TODO doma a ted je potřeba vytahat potřebný počet bodů z fronty a naflákat je do ps.




