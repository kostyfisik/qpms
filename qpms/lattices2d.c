#include "lattices.h"
#include <assert.h>

typedef struct {
  int i, j;
} intcoord2d;

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
  TODO;
} triangular_lattice_gen_t_privstuff_t;

triangular_lattice_gen_t * triangular_lattice_gen_init(double a, TriangularLatticeOrientation ori, bool include_origin)
{
  triangular_lattice_gen_t *g = malloc(sizeof(triangular_latice_gen_t));
  g->orientation = ori;
  g->includes_origin = include_origin;
  g->ps.nrs = 0;
  g->ps.rs = NULL;
  g->ps.points_at_r = NULL;
  g->priv = malloc(sizeof(triangular_lattice_gen_privstuff_t));
  TODO;

  return g;
}


