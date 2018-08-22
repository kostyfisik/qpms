#ifndef LATTICES_H
#define LATTICES_H
#include <math.h>
#include <stdbool.h>

#define M_SQRT3 1.7320508075688772935274463415058724

// This might be reduced to x, y only; not sure yet
typedef struct {
	double key; // distance key in a given lattice
	double x, y, r, phi;
} point2d;


static inline point2d point2d_fromxy(const double x, const double y) {
	point2d p;
	p.x = x;
	p.y = y;
	p.r = sqrt(x*x+y*y);
	p.phi = atan2(y, x);
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
	point2d **points_at_r; // of length nrs+1
	/* // redundand (therefore removed) members
         * point2d *points; // redundant as it is the same as points_at_r[0]
	 * size_t npoints;  // redundant as it is the same as points_at_r[nrs]-points_at_r[0]
	 */ 
} points2d_rordered_t;



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

	// private:
	triangular_lattice_gen_privstuff_t *priv;

} triangular_lattice_gen_t;

triangular_lattice_gen_t *triangular_lattice_gen_init(double a, TriangularLatticeOrientation ori, bool include_origin); 
const points2d_reordered_t * triangular_lattice_gen_getpoints(const triangular lattice_generator_t *g);
int triangular_lattice_gen_extend_to_r(triangular_lattice_generator_t *g, double r);
int triangular_lattice_gen_extend_to_steps(triangular_lattice_generator_t *g, int maxsteps);
void triangular_lattice_gen_free(triangular_lattice_generator_t *g);

#if 0

/*
 * HONEYCOMB LATTICE
 */

typedef struct {

} honeycomb_lattice_generator_t;

#endif

#endif // LATTICES_H
