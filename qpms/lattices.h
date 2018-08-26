#ifndef LATTICES_H
#define LATTICES_H


/* IMPORTANT TODO
 * ==============
 *
 * The current content of this part (and the implementation) is extremely ugly.
 * When I have some time, I have to rewrite this in the style of lattices2d.py
 *
 */



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
 * THE NICE PART (adaptation of lattices2d.py)
 * ===========================================
 *
 * all the functions are prefixed with l2d_
 * convention for argument order: inputs, *outputs, add. params
 */

#define BASIS_RTOL 1e-13

typedef enum {
	OBLIQUE = 1,
	RECTANGULAR = 2,
	SQUARE = 4,
	RHOMBIC = 5,
	EQUILATERAL_TRIANGULAR = 3,
	RIGHT_ISOSCELES=SQUARE,
	PARALLELOGRAMMIC=OBLIQUE,
	CENTERED_RHOMBIC=RECTANGULAR,
	RIGHT_TRIANGULAR=RECTANGULAR,
	CENTERED_RECTANGULAR=RHOMBIC,
	ISOSCELE_TRIANGULAR=RHOMBIC,
	RIGHT_ISOSCELE_TRIANGULAR=SQUARE,
	HEXAGONAL=EQUILATERAL_TRIANGULAR
} LatticeType;


/*
 * Lagrange-Gauss reduction of a 2D basis.
 * The output shall satisfy |out1| <= |out2| <= |out2 - out1|
 */
void l2d_reduceBasis(cart2_t in1, cart2_t in2, cart2_t *out1, cart2_t *out2);

/* 
 * This gives the "ordered shortest triple" of base vectors (each pair from the triple
 * is a base) and there may not be obtuse angle between o1, o2 and between o2, o3
 */
void l2d_shortestBase3(cart2_t i1, cart2_t i2, cart2_t *o1, cart2_t *o2, cart2_t *o3);

/* 
 * TODO doc
 * return value is 4 or 6.
 */
int l2d_shortestBase46(cart2_t i1, cart2_t i2,  cart2_t *o1, cart2_t *o2, cart2_t *o3, cart2_t *o4, cart2_t *o5, cart2_t *o6, double rtol);
// variant
int l2d_shortestBase46_arr(cart2_t i1, cart2_t i2,  cart2_t *oarr, double rtol);

// Determines whether angle between inputs is obtuse
bool l2d_is_obtuse(cart2_t i1, cart2_t i2, double rtol);

/* 
 * Given two basis vectors, returns 2D Bravais lattice type.
 */
LatticeType l2d_classifyLattice(cart2_t b1, cart2_t b2, double rtol);

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

// Reciprocal bases
void l2d_reciprocalBasis1(cart2_t b1, cart2_t b2, cart2_t *rb1, cart2_t *rb2);
void l2d_reciprocalBasis2pi(cart2_t b1, cart2_t b2, cart2_t *rb1, cart2_t *rb2);


/*
 * THE MORE OR LESS OK PART
 * ========================
 */

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


// sorts arbitrary points and creates points2d_rordered_t
points2d_rordered_t *points2d_rordered_frompoints(const point2d *orig_base,
	       	size_t nmemb, double rtol, double atol);

// returns a copy but shifted by a constant (actually in a stupid way, but whatever)
points2d_rordered_t *points2d_rordered_shift(const points2d_rordered_t *orig,
		point2d shift, double rtol, double atol);

// returns a copy but scaled by a factor
points2d_rordered_t *points2d_rordered_scale(const points2d_rordered_t *orig,
       	double factor);


/* The destructor: use only for results of
 *  - points2D_rordered_frompoints,
 *  - points2d_rordered_shift,
 *  - points2d_rordered_scale.
 */
void points2d_rordered_free(points2d_rordered_t *);

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


ptrdiff_t points2d_rordered_locate_r(const points2d_rordered_t *, double r);

// returns a "view" (does not copy any of the arrays)
// -- DO NOT FREE orig BEFORE THE END OF SCOPE OF THE RESULT
points2d_rordered_t points2d_rordered_annulus(const points2d_rordered_t *orig, double minr, bool minr_inc,
		double maxr, bool maxr_inc);



/*
 * THE UGLY PART
 * =============
 */

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


/*
 * HONEYCOMB LATTICE
 */

typedef struct {
	// public:
	points2d_rordered_t ps;
	TriangularLatticeOrientation orientation;
	double a;
	double h;

	// private:
	triangular_lattice_gen_t *tg;
} honeycomb_lattice_gen_t;

honeycomb_lattice_gen_t *honeycomb_lattice_gen_init_h(double h, TriangularLatticeOrientation ori);
honeycomb_lattice_gen_t *honeycomb_lattice_gen_init_a(double a, TriangularLatticeOrientation ori);
int honeycomb_lattice_gen_extend_to_steps(honeycomb_lattice_gen_t *g, int maxsteps);
int honeycomb_lattice_gen_extend_to_r(honeycomb_lattice_gen_t *g, double r);
void honeycomb_lattice_gen_free(honeycomb_lattice_gen_t *g);

#endif // LATTICES_H
