#ifndef LATTICES_H
#define LATTICES_H
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stddef.h>
#define M_SQRT3 1.7320508075688772935274463415058724
#define M_SQRT3_2 (M_SQRT3/2)
#define M_1_SQRT3 0.57735026918962576450914878050195746


// fuck, I already had had suitable type
#include "vectors.h"
typedef cart2_t point2d;

static inline point2d point2d_fromxy(const double x, const double y) {
	point2d p;
	p.x = x;
	p.y = y;
	return p;
}

/* 
 * General set of points ordered by the r-coordinate.
 * Typically, this will include all lattice inside a certain circle.
 * This structure is internally used by the "lattice generators" below.
 * It does not have its memory management of its own, as it is handled
 * by the "generators". For everything except the generators,
 * this structure shall be read-only.
 */
typedef struct {
	size_t nrs; // number of different radii
	double *rs; // the radii; of length nrs (largest contained radius == rs[nrs-1])
	point2d *base; 
	ptrdiff_t *r_offsets; // of length nrs+1 (using relative offsets due to possible realloc's)
	// the jth point of i-th radius is base[r_offsets[i]+j] or using the inline below..
	/* // redundand (therefore removed) members
         * point2d *points; // redundant as it is the same as points_at_r[0]
	 * size_t npoints;  // redundant as it is the same as points_at_r[nrs]-points_at_r[0]
	 */ 
} points2d_rordered_t;

// returns a copy but scaled by a factor
points2d_rordered_t *points2d_rordered_scale(const points2d_rordered_t *orig, double factor);

void points2d_rordered_free(points2d_rordered_t *); // use only for result of points2d_rordered_scale

static inline point2d points2d_rordered_get_point(const points2d_rordered_t *ps, int r_order, int i) {
	assert(i >= 0);
	assert(r_order < ps->nrs);
	assert(i < (ps->r_offsets[r_order+1] - ps->r_offsets[r_order]));
	return ps->base[ps->r_offsets[r_order] + i];
}

static inline double points2d_rordered_get_r(const points2d_rordered_t *ps, int r_order) {
	assert(r_order < ps->nrs);
	return ps->rs[r_order];
}

/* 
 * EQUILATERAL TRIANGULAR LATTICE
 */

typedef enum {
	TRIANGULAR_VERTICAL, // there is a lattice base vector parallel to the y-axis
  	TRIANGULAR_HORIZONTAL // there is a lattice base vector parallel to the x-axis
} TriangularLatticeOrientation;

// implementation data structures; not needed in the header file
typedef struct triangular_lattice_gen_privstuff_t triangular_lattice_gen_privstuff_t;

typedef struct {
	// public:
	points2d_rordered_t ps;
	TriangularLatticeOrientation orientation;
	double a; // lattice vector length
	
	// not sure if needed:
	bool includes_origin;

	// Denotes an offset of the "origin" point; meaning step hexshift * a / sqrt(2) upwards
	// or leftwards for the horizontal or vertical orientations, respectively.
	int hexshift; 

	// private:
	triangular_lattice_gen_privstuff_t *priv;

} triangular_lattice_gen_t;

triangular_lattice_gen_t *triangular_lattice_gen_init(double a, TriangularLatticeOrientation ori, bool include_origin,
		int halfoffset); 
const points2d_rordered_t * triangular_lattice_gen_getpoints(const triangular_lattice_gen_t *g);
int triangular_lattice_gen_extend_to_r(triangular_lattice_gen_t *g, double r);
int triangular_lattice_gen_extend_to_steps(triangular_lattice_gen_t *g, int maxsteps);
void triangular_lattice_gen_free(triangular_lattice_gen_t *g);

#if 0

/*
 * HONEYCOMB LATTICE
 */

typedef struct {

} honeycomb_lattice_generator_t;

#endif

#endif // LATTICES_H
