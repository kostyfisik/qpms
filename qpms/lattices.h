#ifndef LATTICES_H
#define LATTICES_H

#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#define M_SQRT3 1.7320508075688772935274463415058724
#define M_SQRT3_2 (M_SQRT3/2)
#define M_1_SQRT3 0.57735026918962576450914878050195746



/* IMPORTANT TODO
 * ==============
 *
 * The current content of this part (and the implementation) is extremely ugly.
 * When I have some time, I have to rewrite this in the style of lattices2d.py
 *
 */

typedef enum LatticeDimensionality {
	LAT1D = 1,
	LAT2D = 2,
	LAT3D = 4,
	SPACE1D = 8,
	SPACE2D = 16,
	SPACE3D = 32,
	LAT_1D_IN_3D = 33,
	LAT_2D_IN_3D = 34,
	LAT_3D_IN_3D = 40,
	// special coordinate arrangements (indicating possible optimisations)
	LAT_ZONLY = 64,
	LAT_XYONLY = 128,
	LAT_1D_IN_3D_ZONLY = 97, // LAT1D | SPACE3D | 64
	LAT_2D_IN_3D_XYONLY = 162 // LAT2D | SPACE3D | 128
} LatticeDimensionality;

inline static bool LatticeDimensionality_checkflags(
		LatticeDimensionality a, LatticeDimensionality flags_a_has_to_contain) {
	return ((a & flags_a_has_to_contain) == flags_a_has_to_contain);
}

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
 * GENERIC LATTICE POINT GENERATOR TYPE PGenSph
 * ============================================
 *
 * A bit of OOP-in-C brainfuck here.
 * 
 * The basic principle of operation is following:
 * Instead of a list (array) of points, an initialized PGenSph object 
 * is passed to a function that does something over a set of points.
 * Each time PGenSph-type object is "called", it returns PGenSphReturnData, 
 * which contains a point in spherical coordinates (sph_t) and some metadata.
 *
 * After the last generated point, the generator frees all internal memory
 * and returns PGenSphReturnData with PGEN_NOTDONE flag unset (the rest
 * shall be considered invalid data).
 * The caller can also decide not to use the rest and end getting the points
 * even when the PGEN_NOTDONE was set in the last returned data.
 * In such case, the caller shall call PGenSph_destroy() manually.
 *
 * MEMORY MANAGEMENT POLICY
 * ------------------------
 * The basic PGenSph structure shall be allocated on stack (it's only two pointers),
 * everything internal goes on heap.
 */

struct PGenSph;  // full definition below

typedef enum PGenPointFlags {
	PGEN_NOTDONE = 2, // The most important flag: when this is not set, the interation ended – other data returned should be considered nonsense and at this point, the generator should have de-allocated all internal memory.
	PGEN_NEWR = 1, // The r-coordinate is different than in the previous generated point (so radial parts of the calculation have to be redone);
	PGEN_AT_Z = 4, // This is set if we are at the z-axis (theta is either 0 or M_PI)
	PGEN_AT_XY = 8, // This is set if we are at the xy-plane (theta is M_PI2)
	PGEN_DONE = 0, // convenience value, not an actual flag
} PGenPointFlags;

typedef struct PGenSphReturnData {
  PGenPointFlags flags; // metatada
  sph_t point_sph; // the actual point data
} PGenSphReturnData;

static const PGenSphReturnData PGenSphDoneVal = {PGEN_DONE, {0,0,0}}; // convenience constant for use in the exctractor implementations

typedef struct PGenSphClassInfo { // static PGenSph info
	char * const name; // mainly for debugging purposes
	PGenSphReturnData (*next)(struct PGenSph *); // This contains the actual generator procedure (TODO shouldn't I rather point to stateData?)
	void (*destructor)(struct PGenSph *); // Destructor to be called by next() at iteration end, or by the caller if ending the generation prematurely
} PGenSphClassInfo;

// TOP DATA STRUCTURE DEFINITION HERE
typedef struct PGenSph {
	const PGenSphClassInfo * /*const*/ c;
	void *stateData; // shall be NULL if invalid (destroyed)
} PGenSph;

static inline void PGenSph_destroy(PGenSph *g) {
	g->c->destructor(g);
	assert(g->stateData == NULL); // this should be done by the destructor
}

static inline PGenSphReturnData PGenSph_next(PGenSph *g) {
	// TODO maybe some asserts around here
	return g->c->next(g);
}

static inline bool PGenSph_notDone(PGenSphReturnData data) {
	return data.flags & PGEN_NOTDONE ? true : false;
}

/*
 * Some basic lattice generators implementing the abstract interface above (implemented in latticegens.c).
 */

// This one simply iterates over an existing array of Point2d
extern const PGenSphClassInfo PGenSph_FromPoint2DArray; // TODO Do I even need this to be declared here?
PGenSph PGenSph_FromPoints2DArray_new(const point2d *points, size_t len);

extern const PGenSphClassInfo PGenSph_zAxis;
typedef enum PGenSph_zAxis_incrementDirection{
    //PGENSPH_ZAXIS_POSITIVE_INC, // not implemented
    //PGENSPH_ZAXIS_NEGATIVE_INC, // not implemented
    PGENSPH_ZAXIS_INC_FROM_ORIGIN,
    PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN
} PGenSph_zAxis_incrementDirection;
PGenSph PGenSph_zAxis_new_minMaxR(double period, double offset, double minR, bool inc_minR, double maxR, bool inc_maxR,
    PGenSph_zAxis_incrementDirection incdir);


/*
 * THE NICE PART (adaptation of lattices2d.py)
 * ===========================================
 *
 * all the functions are prefixed with l2d_
 * convention for argument order: inputs, *outputs, add. params
 */

#define BASIS_RTOL 1e-13

// Bravais lattice types
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
} LatticeType2;

#if 0
// Wallpaper groups
typedef enum {
	TODO
} SpaceGroup2;
#endif

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
bool l2d_is_obtuse_r(cart2_t i1, cart2_t i2, double rtol);
bool l2d_is_obtuse(cart2_t i1, cart2_t i2);

/* 
 * Given two basis vectors, returns 2D Bravais lattice type.
 */
LatticeType2 l2d_classifyLattice(cart2_t b1, cart2_t b2, double rtol);

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


// returns the radius of inscribed circle of a hexagon (or rectangle/square if applicable) created by the shortest base triple
double l2d_hexWebInCircleRadius(cart2_t b1, cart2_t b2);

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

static inline TriangularLatticeOrientation reverseTriangularLatticeOrientation(TriangularLatticeOrientation o){
	switch(o) {
		case TRIANGULAR_VERTICAL:
			return TRIANGULAR_HORIZONTAL;
			break;
		case TRIANGULAR_HORIZONTAL:
			return TRIANGULAR_VERTICAL;
			break;
		default:
			abort();
	}
	abort();
}

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
