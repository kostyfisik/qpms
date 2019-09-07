#include "lattices.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int i, j;
} intcoord2_t;

static inline int sqi(int x) { return x*x; }

static inline double sqf(double x) { return x*x; }

void points2d_rordered_free(points2d_rordered_t *p) {
  free(p->rs);
  free(p->base);
  free(p->r_offsets);
  free(p);
}

points2d_rordered_t *points2d_rordered_scale(const points2d_rordered_t *orig, const double f)
{
  points2d_rordered_t *p = malloc(sizeof(points2d_rordered_t));
  if(0 == orig->nrs) { // orig is empty
    p->nrs = 0;
    p->rs = NULL;
    p->base = NULL;
    p->r_offsets = NULL;
    return p;
  }
  p->nrs = orig->nrs;
  p->rs = malloc(p->nrs*sizeof(double));
  p->r_offsets = malloc((p->nrs+1)*sizeof(ptrdiff_t));

  const double af = fabs(f);
  for(size_t i = 0; i < p->nrs; ++i) {
    p->rs[i] = orig->rs[i] * af;
    p->r_offsets[i] = orig->r_offsets[i];
  }
  p->r_offsets[p->nrs] = orig->r_offsets[p->nrs];
  p->base = malloc(sizeof(point2d) * p->r_offsets[p->nrs]);
  for(size_t i = 0; i < p->r_offsets[p->nrs]; ++i)
    p->base[i] = point2d_fromxy(orig->base[i].x * f, orig->base[i].y * f);
  return p;
}

ptrdiff_t points2d_rordered_locate_r(const points2d_rordered_t *p, const double r) {
  //if(p->r_rs[0] > r)
  //  return -1;
  //if(p->r_rs[p->nrs-1] < r)
  //  return p->nrs;
  ptrdiff_t lo = 0, hi = p->nrs-1, piv;
  while(lo < hi) {
    piv = (lo + hi + 1) / 2;
    if(p->rs[piv] > r) // the result will be less or equal
      hi = piv - 1;
    else
      lo = piv;
  }
  return lo;
}

points2d_rordered_t points2d_rordered_annulus(const points2d_rordered_t *orig,
    double minr, bool inc_minr, double maxr, bool inc_maxr) {
  points2d_rordered_t p;
  ptrdiff_t imin, imax;
  imin = points2d_rordered_locate_r(orig, minr);
  imax = points2d_rordered_locate_r(orig, maxr);
  // TODO check 
  if(imax >= orig->nrs) --imax;
  if(imax < 0) goto nothing;
  // END TODO
  if (!inc_minr && (orig->rs[imin] <= minr)) ++imin;
  if (!inc_maxr && (orig->rs[imax] >= maxr)) --imax;
  if (imax < imin) { // it's empty
nothing:
    p.nrs = 0;
    p.base = NULL;
    p.rs = NULL;
    p.r_offsets = NULL;
  } else {
    p.base = orig->base;
    p.nrs = imax - imin + 1;
    p.rs = orig->rs + imin;
    p.r_offsets = orig->r_offsets + imin;
  }
  return p;
}


static inline double pr2(const point2d p) {
  return sqf(p.x) + sqf(p.y);
}

static inline double prn(const point2d p) {
  return sqrt(pr2(p));
}

static int point2d_cmp_by_r2(const void *p1, const void *p2) {
  const point2d *z1 = (point2d *) p1, *z2 = (point2d *) p2; 
  double dif = pr2(*z1) - pr2(*z2);
  if(dif > 0) return 1;
  else if(dif < 0) return -1;
  else return 0;
}

static points2d_rordered_t *points2d_rordered_frompoints_c(point2d *orig_base, const size_t nmemb,
    const double rtol, const double atol, bool copybase)
{
  // TODO should the rtol and atol relate to |r| or r**2? (Currently: |r|)
  assert(rtol >= 0);
  assert(atol >= 0);
  points2d_rordered_t *p = malloc(sizeof(points2d_rordered_t));
  if(nmemb == 0) {
    p->nrs = 0;
    p->rs = NULL;
    p->base = NULL;
    p->r_offsets = NULL;
    return p;
  }

  if (copybase) {
    p->base = malloc(nmemb * sizeof(point2d));
    memcpy(p->base, orig_base, nmemb * sizeof(point2d));
  } else 
    p->base = orig_base;

  qsort(p->base, nmemb, sizeof(point2d), point2d_cmp_by_r2);

  // first pass: determine the number of "different" r's.
  size_t rcount = 0;
  double rcur = -INFINITY;
  double rcurmax = -INFINITY;
  for (size_t i = 0; i < nmemb; ++i)
    if ((rcur = prn(p->base[i])) > rcurmax) {
      ++rcount;
      rcurmax = rcur * (1 + rtol) + atol;
    }

  p->nrs = rcount;
  // TODO check malloc return values
  p->rs = malloc(rcount * sizeof(double));
  p->r_offsets = malloc((rcount+1) * sizeof(ptrdiff_t));


  // second pass: fill teh rs;
  size_t ri = 0;
  size_t rcurcount = 0;
  rcur = prn(p->base[0]);
  rcurmax = rcur * (1 + rtol) + atol;
  double rcursum = 0;
  p->r_offsets[0] = 0;
  for (size_t i = 0; i < nmemb; ++i) {
    rcur = prn(p->base[i]);
    if (rcur > rcurmax) {
      p->rs[ri] = rcursum / (double) rcurcount; // average of the accrued r's within tolerance
      ++ri;
      p->r_offsets[ri] = i; //r_offsets[ri-1] + rcurcount (is the same)
      rcurcount = 0;
      rcursum = 0;
      rcurmax = rcur * (1 + rtol) + atol;
    } 
    rcursum += rcur;
    ++rcurcount;
  }
  p->rs[ri] = rcursum / (double) rcurcount;
  p->r_offsets[rcount] = nmemb;

  return p;
}

points2d_rordered_t *points2d_rordered_frompoints(const point2d *orig_base, const size_t nmemb,
    const double rtol, const double atol) 
{
  return points2d_rordered_frompoints_c((point2d *)orig_base, nmemb, rtol, atol, true);
}

points2d_rordered_t *points2d_rordered_shift(const points2d_rordered_t *orig, const  point2d shift,
    double rtol, double atol)
{
  size_t n = (orig->nrs > 0) ?
      orig->r_offsets[orig->nrs] - orig->r_offsets[0] : 0;

  point2d * shifted = malloc(n * sizeof(point2d));

  for(size_t i = 0; i < n; ++i) 
    shifted[i] = cart2_add(orig->base[i+orig->r_offsets[0]], shift);

  return points2d_rordered_frompoints_c(shifted, n,rtol, atol, false);
}

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
 * than 3/4.*s*(s+1) + 6*s lattice cells (FIXME pretty stupid but safe estimate).
 *
 * s-th layer circumscribes a circle of radius a * s * sqrt(3)/2.
 *
 */


typedef struct triangular_lattice_gen_privstuff_t {
  intcoord2_t *pointlist_base; // allocated memory for the point "buffer"
  size_t pointlist_capacity; 
  // beginning and end of the point "buffer"
  // not 100% sure what type should I use here
  // (these are both relative to pointlist_base, due to possible realloc's)
  ptrdiff_t pointlist_beg, pointlist_n; // end of queue is at [(pointlist_beg+pointlist_n)%pointlist_capacity]
  int maxs; // the highest layer of the spider web generated (-1 by init, 0 is only origin (if applicable))
  // capacities of the arrays in ps
  size_t ps_rs_capacity;
  size_t ps_points_capacity; // this is the "base" array
  // TODO anything else?
} triangular_lattice_gen_privstuff_t;

static inline int trilat_r2_ij(const int i, const int j) {
  return sqi(i) + sqi(j) + i*j;
}

static inline int trilat_r2_coord(const intcoord2_t c) {
  return trilat_r2_ij(c.i, c.j);
}

// version with offset (n.b. this is includes a factor of 3)
static inline int trilat_3r2_ijs(const int i, const int j, const int s) {
  return 3*(sqi(i) + sqi(j) + i*j + j*s) + sqi(s);
}

static inline int trilat_3r2_coord_s(const intcoord2_t c, const int s) {
  return trilat_3r2_ijs(c.i, c.j, s);
}



// Classify points into sextants (variant [a] above)
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

static inline size_t tlgp_pl_end(const triangular_lattice_gen_privstuff_t *p) {
  return (p->pointlist_beg + p->pointlist_n) % p->pointlist_capacity;
}

#if 0
static inline void tlgpl_end_inc(triangular_lattice_gen_privstuff_t *p) {
  p->p_pointlist_n += 1;
}
#endif

// Puts a point to the end of the point queue
static inline void trilatgen_pointlist_append_ij(triangular_lattice_gen_t *g, int i, int j) {
  intcoord2_t thepoint = {i, j};
  triangular_lattice_gen_privstuff_t *p = g->priv;
  assert(p->pointlist_n < p->pointlist_capacity);
  
  // the actual addition
  p->pointlist_base[tlgp_pl_end(p)] = thepoint;
  p->pointlist_n += 1;
}

// Arange the pointlist queue into a continuous chunk of memory, so that we can qsort() it
static void trilatgen_pointlist_linearise(triangular_lattice_gen_t *g) {
  triangular_lattice_gen_privstuff_t *p = g->priv;
  assert(p->pointlist_n <= p->pointlist_capacity);
  if (p->pointlist_beg + p->pointlist_n <= p->pointlist_capacity)
    return;  // already linear, do nothing
  else if (p->pointlist_n == p->pointlist_capacity) { // full, therefore linear
    p->pointlist_beg = 0;
    return;
  } else { // non-linear; move "to the right"
    while (p->pointlist_beg < p->pointlist_capacity) {
      p->pointlist_base[tlgp_pl_end(p)] = p->pointlist_base[p->pointlist_beg];
      ++(p->pointlist_beg);
    }
    p->pointlist_beg = 0;
    return;
  }
}

static inline intcoord2_t trilatgen_pointlist_first(const triangular_lattice_gen_t *g) {
  return g->priv->pointlist_base[g->priv->pointlist_beg];
}

static inline void trilatgen_pointlist_deletefirst(triangular_lattice_gen_t *g) {
  triangular_lattice_gen_privstuff_t *p = g->priv;
  assert(p->pointlist_n > 0);
  ++p->pointlist_beg;
  if(p->pointlist_beg == p->pointlist_capacity) 
    p->pointlist_beg = 0;
  --(p->pointlist_n);
}

// TODO abort() and void or errorchecks and int?
static int trilatgen_pointlist_extend_capacity(triangular_lattice_gen_t *g, size_t newcapacity) {
  triangular_lattice_gen_privstuff_t *p = g->priv;
  if (newcapacity <= p->pointlist_capacity)
    return 0;

  trilatgen_pointlist_linearise(g);

  intcoord2_t *newmem = realloc(p->pointlist_base, newcapacity * sizeof(intcoord2_t));
  if (newmem != NULL) {
    p->pointlist_base = newmem;
    p->pointlist_capacity = newcapacity;
    return 0;
  } else 
    abort();
  
}

// lower estimate for the number of lattice points inside the circumscribed hexagon, but outside the circle
static inline size_t tlg_circumscribe_reserve(int maxs) {
  if (maxs <= 0) 
    return 0;
  return 3*maxs*(maxs+1)/4 + 6*maxs;
}

static inline size_t tlg_websize(int maxs) {
  if (maxs <= 0)
    return 0;
  else 
    return 3*maxs*(maxs+1); // does not include origin point!
}

static int trilatgen_ensure_pointlist_capacity(triangular_lattice_gen_t *g, int newmaxs) {
  return trilatgen_pointlist_extend_capacity(g,
        tlg_circumscribe_reserve(g->priv->maxs) // Space for those which are already in
      + tlg_websize(newmaxs) - tlg_websize(g->priv->maxs) // space for the new web layers
      + 1 // reserve for the origin
      );
}

static int trilatgen_ensure_ps_rs_capacity(triangular_lattice_gen_t *g, int maxs) {
  if (maxs < 0)
    return 0;

  size_t needed_capacity = 1 // reserve for origin
    + maxs*(maxs+1)/2; // stupid but safe estimate: number of points in a sextant of maxs-layered spider web

  if (needed_capacity <= g->priv->ps_rs_capacity)
    return 0; // probably does not happen, but fuck it

  double *newmem = realloc(g->ps.rs, needed_capacity * sizeof(double));
  if (newmem != NULL)  
    g->ps.rs = newmem; 
  else 
    abort();
  ptrdiff_t *newmem2 = realloc(g->ps.r_offsets, (needed_capacity + 1) * sizeof(ptrdiff_t));
  if (newmem2 != NULL)
    g->ps.r_offsets = newmem2;
  else 
    abort();
  
  g->priv->ps_rs_capacity = needed_capacity;
  return 0;
}

static int trilatgen_ensure_ps_points_capacity(triangular_lattice_gen_t *g, int maxs) {
  if (maxs < 0)
    return 0;
  size_t needed_capacity = 1 /*res. for origin */ + tlg_websize(maxs) /* stupid but safe */;
  if(needed_capacity <= g->priv->ps_points_capacity)
    return 0;

  point2d *newmem = realloc(g->ps.base, needed_capacity * sizeof(point2d));
  if (newmem != NULL) 
    g->ps.base = newmem;
  else
    abort();

  g->priv->ps_points_capacity = needed_capacity;
  return 0;
}

static int trilat_cmp_intcoord2_by_r2(const void *p1, const void *p2) {
  return trilat_r2_coord(*(const intcoord2_t *)p1) - trilat_r2_coord(*(const intcoord2_t *)p2);
}


static int trilat_cmp_intcoord2_by_3r2_plus1s(const void *p1, const void *p2) {
  return trilat_3r2_coord_s(*(const intcoord2_t *)p1, +1) - trilat_3r2_coord_s(*(const intcoord2_t *)p2, +1);
}

static int trilat_cmp_intcoord2_by_3r2_minus1s(const void *p1, const void *p2) {
  return trilat_3r2_coord_s(*(const intcoord2_t *)p1, -1) - trilat_3r2_coord_s(*(const intcoord2_t *)p2, -1);
}

static int trilat_cmp_intcoord2_by_3r2(const void *p1, const void *p2, void *sarg) {
  return trilat_3r2_coord_s(*(const intcoord2_t *)p1, *(int *)sarg) - trilat_3r2_coord_s(*(const intcoord2_t *)p2, *(int *)sarg);
}

static void trilatgen_sort_pointlist(triangular_lattice_gen_t *g) {
  trilatgen_pointlist_linearise(g);
  triangular_lattice_gen_privstuff_t *p = g->priv;
  int (*compar)(const void *, const void *);
  switch (g->hexshift) {
    case 0:
      compar = trilat_cmp_intcoord2_by_r2;
      break;
    case -1:
      compar = trilat_cmp_intcoord2_by_3r2_minus1s;
      break;
    case 1:
      compar = trilat_cmp_intcoord2_by_3r2_plus1s;
      break;
    default:
      abort();
  }
  qsort(p->pointlist_base + p->pointlist_beg, p->pointlist_n, sizeof(intcoord2_t), compar);
}

triangular_lattice_gen_t * triangular_lattice_gen_init(double a, TriangularLatticeOrientation ori, bool include_origin, 
    int hexshift)
{
  triangular_lattice_gen_t *g = malloc(sizeof(triangular_lattice_gen_t));
  g->a = a;
  g->hexshift = ((hexshift % 3)+3)%3; // reduce to the set {-1, 0, 1}
  if (2 == g->hexshift)
    g->hexshift = -1;
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
  g->priv->pointlist_n = 0;
  g->priv->ps_rs_capacity = 0;
  g->priv->ps_points_capacity = 0;
  return g;
}

void triangular_lattice_gen_free(triangular_lattice_gen_t *g) {
  free(g->ps.rs);
  free(g->ps.base);
  free(g->ps.r_offsets);
  free(g->priv->pointlist_base);
  free(g->priv);
  free(g);
}

const points2d_rordered_t * triangular_lattice_gen_getpoints(const triangular_lattice_gen_t *g) {
  return &(g->ps);
}

int triangular_lattice_gen_extend_to_r(triangular_lattice_gen_t * g, const double maxr) {
  return triangular_lattice_gen_extend_to_steps(g, maxr/g->a);
}

int triangular_lattice_gen_extend_to_steps(triangular_lattice_gen_t * g, int maxsteps) 
{
  if (maxsteps <= g->priv->maxs)  // nothing needed
    return 0;
  // TODO FIXME: check for maximum possible maxsteps (not sure what it is)
  int err;
  err = trilatgen_ensure_pointlist_capacity(g, maxsteps
      + abs(g->hexshift) /*FIXME this is quite brainless addition, probably not even needed.*/);
  if(err) return err;
  err = trilatgen_ensure_ps_rs_capacity(g, maxsteps
      + abs(g->hexshift) /*FIXME this is quite brainless addition, probably not even needed.*/);
  if(err) return err;
  err = trilatgen_ensure_ps_points_capacity(g, maxsteps
      + abs(g->hexshift) /*FIXME this is quite brainless addition, probably not even needed.*/);
  if(err) return err;
  
  if(g->includes_origin && g->priv->maxs < 0) // Add origin if not there yet
    trilatgen_pointlist_append_ij(g, 0, 0);
  
  for (int s = g->priv->maxs + 1; s <= maxsteps; ++s) {
    int i, j; 
    // now go along the spider web layer as indicated in the lenghthy comment above
    for (i = s, j = 0; i > 0; --i, ++j) trilatgen_pointlist_append_ij(g,i,j);
    for (i = 0, j = s; i + j > 0; --i) trilatgen_pointlist_append_ij(g,i,j);
    for (i = -s, j = s; j > 0; --j) trilatgen_pointlist_append_ij(g,i,j);
    for (i = -s, j = 0; i < 0; ++i, --j) trilatgen_pointlist_append_ij(g,i,j);
    for (i = 0, j = -s; i + j < 0; ++i) trilatgen_pointlist_append_ij(g,i,j);
    for (i = s, j = -s; j < 0; ++j) trilatgen_pointlist_append_ij(g,i,j);
  }

  trilatgen_sort_pointlist(g);
  
  // initialise first r_offset if needed
  if (0 == g->ps.nrs)
    g->ps.r_offsets[0] = 0;

  //ted je potřeba vytahat potřebný počet bodů z fronty a naflákat je do ps.
  // FIXME pohlídat si kapacitu datových typů
  //int maxr2i = sqi(maxsteps) * 3 / 4;
  int maxr2i3 = sqi(maxsteps) * 9 / 4 + sqi(g->hexshift) - abs(3*maxsteps*g->hexshift);
  while (g->priv->pointlist_n > 0) { // This condition should probably be always true anyways.
    intcoord2_t coord = trilatgen_pointlist_first(g);
    //int r2i_cur = trilat_r2_coord(coord);
    //if(r2i_cur > maxr2i)
    int r2i3_cur = trilat_3r2_coord_s(coord, g->hexshift);
    if(r2i3_cur > maxr2i3)
      break;
    g->ps.rs[g->ps.nrs] = sqrt(/*r2i_cur*/ r2i3_cur/3.) * g->a;
    g->ps.r_offsets[g->ps.nrs+1] = g->ps.r_offsets[g->ps.nrs]; // the difference is the number of points on the circle
    while(1) {
      coord = trilatgen_pointlist_first(g);
      //if(r2i_cur != trilat_r2_coord(coord))
      if (r2i3_cur != trilat_3r2_coord_s(coord, g->hexshift))
        break;
      else {
        trilatgen_pointlist_deletefirst(g);
        point2d thepoint;
        switch (g->orientation) {
          case TRIANGULAR_HORIZONTAL:
            thepoint = point2d_fromxy((coord.i+.5*coord.j)*g->a, (M_SQRT3_2*coord.j + g->hexshift*M_1_SQRT3)*g->a);
            break;
          case TRIANGULAR_VERTICAL:
            thepoint = point2d_fromxy(-(M_SQRT3_2*coord.j + g->hexshift*M_1_SQRT3)*g->a, (coord.i+.5*coord.j)*g->a);
            break;
          default:
            abort();
        }
        g->ps.base[g->ps.r_offsets[g->ps.nrs+1]] = thepoint;
        ++(g->ps.r_offsets[g->ps.nrs+1]);
      }
    }
    ++(g->ps.nrs);
  }
  g->priv->maxs = maxsteps;
  return 0;
}

honeycomb_lattice_gen_t *honeycomb_lattice_gen_init_h(double h, TriangularLatticeOrientation ori) {
  double a = M_SQRT3 * h;
  honeycomb_lattice_gen_t *g = honeycomb_lattice_gen_init_a(a, ori);
  g->h = h; // maybe it's not necessary as sqrt is "exact"
  return g;
}

honeycomb_lattice_gen_t *honeycomb_lattice_gen_init_a(double a, TriangularLatticeOrientation ori) {
  honeycomb_lattice_gen_t *g = calloc(1, sizeof(honeycomb_lattice_gen_t)); // this already inits g->ps to zeros
  g->a = a;
  g->h = a * M_1_SQRT3;
  g->tg = triangular_lattice_gen_init(a, ori, true, 1);
  return g;
}

void honeycomb_lattice_gen_free(honeycomb_lattice_gen_t *g) {
  free(g->ps.rs);
  free(g->ps.base);
  free(g->ps.r_offsets);
  triangular_lattice_gen_free(g->tg);
  free(g);
}

int honeycomb_lattice_gen_extend_to_r(honeycomb_lattice_gen_t *g, double maxr) {
  return honeycomb_lattice_gen_extend_to_steps(g, maxr/g->a); /*CHECKME whether g->a is the correct denom.*/
}

int honeycomb_lattice_gen_extend_to_steps(honeycomb_lattice_gen_t *g, const int maxsteps) {
  if (maxsteps <= g->tg->priv->maxs)  // nothing needed
    return 0;
  triangular_lattice_gen_extend_to_steps(g->tg, maxsteps);

  double *newmem = realloc(g->ps.rs, g->tg->ps.nrs * sizeof(double));
  if (NULL != newmem)
    g->ps.rs = newmem;
  else abort();
  ptrdiff_t *newmem2 = realloc(g->ps.r_offsets, (g->tg->ps.nrs+1) * sizeof(ptrdiff_t));
  if (NULL != newmem2)
    g->ps.r_offsets = newmem2;
  else abort();
  point2d *newmem3 = realloc(g->ps.base, 2 * (g->tg->ps.r_offsets[g->tg->ps.nrs]) * sizeof(point2d));
  if (NULL != newmem3)
    g->ps.base = newmem3;
  else abort();

  // Now copy (new) contents of g->tg->ps into g->ps, but with inverse copy of each point
  for (size_t ri = g->ps.nrs; ri <= g->tg->ps.nrs; ++ri) 
    g->ps.r_offsets[ri] = g->tg->ps.r_offsets[ri] * 2;
  for (ptrdiff_t i_orig = g->tg->ps.r_offsets[g->ps.nrs]; i_orig < g->tg->ps.r_offsets[g->tg->ps.nrs]; ++i_orig) {
    point2d p = g->tg->ps.base[i_orig];
    g->ps.base[2*i_orig] = p;
    p.x *= -1; p.y *= -1;
    g->ps.base[2*i_orig + 1] = p;
  }
  g->ps.nrs = g->tg->ps.nrs;
  return 0;
}



// THE NICE PART


/*
 * Lagrange-Gauss reduction of a 2D basis.
 * The output shall satisfy |out1| <= |out2| <= |out2 - out1|
 */
void l2d_reduceBasis(cart2_t b1, cart2_t b2, cart2_t *out1, cart2_t *out2){
  double B1 = cart2_dot(b1, b1);
  double mu = cart2_dot(b1, b2) / B1;
  b2 = cart2_substract(b2, cart2_scale(round(mu), b1));
  double B2 = cart2_dot(b2, b2);
  while(B2 < B1) {
    cart2_t b2t = b1;
    b1 = b2;
    b2 = b2t;
    B1 = B2;
    mu = cart2_dot(b1, b2) / B1;
    b2 = cart2_substract(b2, cart2_scale(round(mu), b1));
    B2 = cart2_dot(b2, b2);
  }
  *out1 = b1;
  *out2 = b2;
}

void l3d_reduceBasis(const cart3_t in[3], cart3_t out[3]) {
  memcpy(out, in, 3*sizeof(cart3_t));
  QPMS_ENSURE_SUCCESS(qpms_reduce_lattice_basis((double *)out, 3, 3, 1.));
}

/*
 * This gives the "ordered shortest triple" of base vectors (each pair from the triple
 * is a base) and there may not be obtuse angle between o1, o2 and between o2, o3
 */
void l2d_shortestBase3(cart2_t b1, cart2_t b2, cart2_t *o1, cart2_t *o2, cart2_t *o3){
  l2d_reduceBasis(b1, b2, &b1, &b2);
  *o1 = b1;
  if (l2d_is_obtuse_r(b1, b2, 0)) {
    *o3 = b2;
    *o2 = cart2_add(b2, b1);
  } else {
    *o2 = b2;
    *o3 = cart2_substract(b2, b1);
  }
}

// Determines whether angle between inputs is obtuse
bool l2d_is_obtuse_r(cart2_t b1, cart2_t b2, double rtol) {
  const double B1 = cart2_normsq(b1);
  const double B2 = cart2_normsq(b2);
  const cart2_t b3 = cart2_substract(b2, b1);
  const double B3 = cart2_normsq(b3);
  const double eps = rtol * (B1 + B2); // TODO check what kind of quantity this should be. Maybe rtol should relate to lengths, not lengths**2
  return (B3 - B2 - B1 > eps);
}


/*
 * TODO doc
 * return value is 4 or 6.
 */
int l2d_shortestBase46(const cart2_t i1, const cart2_t i2,  cart2_t *o1, cart2_t *o2, cart2_t *o3, cart2_t *o4, cart2_t *o5, cart2_t *o6, double rtol){
  cart2_t b1, b2, b3;
  l2d_reduceBasis(i1, i2, &b1, &b2);
  const double b1s = cart2_normsq(b1);
  const double b2s = cart2_normsq(b2);
  b3 = cart2_substract(b2, b1);
  const double b3s = cart2_normsq(b3);
  const double eps = rtol * (b1s + b2s); // TODO check the same as in l2d_is_obtuse_r
  if(fabs(b3s-b2s-b1s) < eps) {
    *o1 = b1; *o2 = b2; *o3 = cart2_scale(-1, b1); *o4 = cart2_scale(-1, b2);
    return 4;
  }
  else {
    if (b3s-b2s-b1s > eps) { //obtuse 
      b3 = b2;
      b2 = cart2_add(b2, b1);
    }
    *o1 = b1; *o2 = b2; *o3 = b3;
    *o4 = cart2_scale(-1, b1);
    *o5 = cart2_scale(-1, b2);
    *o6 = cart2_scale(-1, b3);
    return 6;
  }
}

/*
 * Given two basis vectors, returns 2D Bravais lattice type.
 */
LatticeType2 l2d_classifyLattice(cart2_t b1, cart2_t b2, double rtol)
{
  l2d_reduceBasis(b1, b2, &b1, &b2);
  cart2_t b3 = cart2_substract(b2, b1);
  double b1s = cart2_normsq(b1), b2s = cart2_normsq(b2), b3s = cart2_normsq(b3);
  double eps = rtol * (b2s + b1s); //FIXME what should eps be?
  // avoid obtuse angle between b1 and b2. TODO this should be yet tested
  // TODO use is_obtuse here?
  if (b3s - b2s - b1s > eps) {
    b3 = b2;
    b2 = cart2_add(b2, b1);
    // N.B. now the assumption |b3| >= |b2| is no longer valid
    // b3 = cart2_substract(b2, b1)
    b2s = cart2_normsq(b2); 
    b3s = cart2_normsq(b3);
  }
  if (fabs(b2s-b1s) < eps || fabs(b2s - b3s) < eps) { // isoscele
    if (fabs(b3s-b1s) <  eps)
      return EQUILATERAL_TRIANGULAR;
    else if (fabs(b3s - 2*b1s))
      return SQUARE;
    else
      return RHOMBIC;
  } else if (fabs(b3s-b2s-b1s) < eps)
    return RECTANGULAR;
  else
    return OBLIQUE;
}

LatticeFlags l2d_detectRightAngles(cart2_t b1, cart2_t b2, double rtol)
{
  l2d_reduceBasis(b1, b2, &b1, &b2);
  cart2_t ht = cart2_substract(b2, b1);
  double b1s = cart2_normsq(b1), b2s = cart2_normsq(b2), hts = cart2_normsq(ht);
  double eps = rtol * (b2s + b1s); //FIXME what should eps be?
  if (hts - b2s - b1s <= eps)
    return ORTHOGONAL_01;
  else 
    return NOT_ORTHOGONAL;
}

LatticeFlags l3d_detectRightAngles(const cart3_t basis_nr[3], double rtol) 
{
  cart3_t b[3];
  l3d_reduceBasis(basis_nr, b);
  LatticeFlags res = NOT_ORTHOGONAL;
  for (int i = 0; i < 3; ++i) {
    cart3_t ba = b[i], bb = b[(i+1) % 3];
    cart3_t ht = cart3_substract(ba, bb);
    double bas = cart3_normsq(ba), bbs = cart3_normsq(ba), hts = cart3_normsq(ht);
    double eps = rtol * (bas + bbs);
    if (hts - bbs - bas <= eps)
      res |= ((LatticeFlags[]){ORTHOGONAL_01, ORTHOGONAL_12, ORTHOGONAL_02})[i];
  }
  return res;
}

# if 0
// variant
int l2d_shortestBase46_arr(cart2_t i1, cart2_t i2,  cart2_t *oarr, double rtol);

// Determines whether angle between inputs is obtuse
bool l2d_is_obtuse_r(cart2_t i1, cart2_t i2, double rtol);


// Other functions in lattices2d.py: TODO?
// range2D()
// generateLattice()
// generateLatticeDisk()
// cutWS()
// filledWS()
// change_basis()

/*
 * Given basis vectors, returns the corners of the Wigner-Seits unit cell (W1, W2, -W1, W2)
 * for rectangular and square lattice or (w1, w2, w3, -w1, -w2, -w3) otherwise.
 */
int l2d_cellCornersWS(cart2_t i1, cart2_t i2,  cart2_t *o1, cart2_t *o2, cart2_t *o3, cart2_t *o4, cart2_t *o5, cart2_t *o6, double rtol);
// variant
int l2d_cellCornersWS_arr(cart2_t i1, cart2_t i2,  cart2_t *oarr, double rtol);

#endif

// Reciprocal bases; returns 0 on success, TODO non-zero if b1 and b2 are parallel
int l2d_reciprocalBasis1(cart2_t b1, cart2_t b2, cart2_t *rb1, cart2_t *rb2) {
  l2d_reduceBasis(b1, b2, &b1, &b2);
  const double det = b1.x * b2.y - b1.y * b2.x;
  if (!det) {
    rb1->x = rb1->y = rb2->x = rb2->y = NAN;
    return QPMS_ERROR; // TODO more specific error code
  } else {
    rb1->x =  b2.y / det;
    rb1->y = -b2.x / det;
    rb2->x = -b1.y / det;
    rb2->y =  b1.x / det;
    return QPMS_SUCCESS;
  }
}

int l2d_reciprocalBasis2pi(cart2_t b1, cart2_t b2, cart2_t *rb1, cart2_t *rb2) {
  int retval = l2d_reciprocalBasis1(b1, b2, rb1, rb2);
  if (retval == QPMS_SUCCESS) {
    *rb1 = cart2_scale(2 * M_PI, *rb1);
    *rb2 = cart2_scale(2 * M_PI, *rb2);
  }
  return retval;
};

// returns the radius of inscribed circle of a hexagon (or rectangle/square if applicable) created by the shortest base triple
double l2d_hexWebInCircleRadius(cart2_t i1, cart2_t i2) {
  cart2_t b1, b2, b3;
  l2d_shortestBase3(i1, i2, &b1, &b2, &b3);
  const double r1 = cart2norm(b1), r2 = cart2norm(b2), r3 = cart2norm(b3);
  const double p = (r1+r2+r3)*0.5;
  return 2*sqrt(p*(p-r1)*(p-r2)*(p-r3))/r3; // CHECK is r3 guaranteed to be longest?
}


double l2d_unitcell_area(cart2_t b1, cart2_t b2) {
  l2d_reduceBasis(b1, b2, &b1, &b2);
  const double det = b1.x * b2.y - b1.y * b2.x;
  return fabs(det);
}
  


