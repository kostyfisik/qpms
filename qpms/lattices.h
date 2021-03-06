/*! \file lattices.h
 * \brief Lattice point generators and lattice vector analysis / transformation.
 *
 */
#ifndef LATTICES_H
#define LATTICES_H

#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#ifndef M_SQRT3
#define M_SQRT3 1.7320508075688772935274463415058724
#endif
#ifndef M_SQRT3_2
#define M_SQRT3_2 (M_SQRT3/2)
#endif
#ifndef M_1_SQRT3
#define M_1_SQRT3 0.57735026918962576450914878050195746
#endif



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

/// Lattice basis reduction.
/** This is currenty a bit naïve implementation of
 * Lenstra-Lenstra-Lovász algorithm.
 *
 * The reduction happens in-place, i.e. the basis vectors in \a b are
 * replaced with the reduced basis.
 */
int qpms_reduce_lattice_basis(double *b, ///< Array of dimension [bsize][ndim].
		const size_t bsize, ///< Number of the basis vectors (dimensionality of the lattice).
		const size_t ndim, ///< Dimension of the space into which the lattice is embedded.
		/// Lovász condition parameter \f$ \delta \f$.
		/** Polynomial time complexity guaranteed for \f$\delta \in (1/4,1)\f$.
		 */
		double delta
		);

/// Generic lattice point generator type.
/**
 * A bit of OOP-in-C brainfuck here.
 * 
 * The basic principle of operation is following:
 * Instead of a list (array) of points, an initialized PGen object 
 * is passed to a function that does something over a set of points.
 * Each time PGen-type object is "called" (more specifically, one of
 * the "methods" specified in the PGenClassInfo structure in @ref c, 
 * it returns PGenReturnData
 * which contains a point in given coordinates (depending on the generator 
 * class) and some metadata.
 *
 * After the last generated point, the generator frees all internal memory
 * and returns PGenSphReturnData with PGEN_NOTDONE flag unset (the rest
 * shall be considered invalid data).
 * The caller can also decide not to use the rest and end getting the points
 * even when the PGEN_NOTDONE was set in the last returned data.
 * In such case, the caller shall call PGen_destroy() manually.
 *
 * Methods
 * -------
 *
 *  The standard wrapper "methods" to generate a single point in a given
 *  coordinate system are
 *
 *  * PGen_next_z(),
 *  * PGen_next_cart2(),
 *  * PGen_next_cart3(),
 *  * PGen_next_pol(),
 *  * PGen_next_sph().
 *
 * Memory management policy
 * ------------------------
 *
 * The basic PGen structure shall be allocated on stack (it's only two pointers),
 * everything internal goes on heap.
 */
typedef struct PGen {
	/// Pointer to the "class" metadata defining the behaviour of the generator.
	const struct PGenClassInfo * /*const*/ c; 
	/// Pointer to internal state data; shall be NULL if invalid (destroyed);
	void *stateData;
} PGen;

typedef enum PGenPointFlags {
	/** The most important flag: when this is not set, the
	 *  interation ended – other data returned should be 
	 *  considered nonsense and at this point, the generator 
	 *  should have de-allocated all internal memory.
	 */
	PGEN_NOTDONE = 2, 
        /** Set if the r-coordinate is not different than in the 
	 *  previous generated point (so radial parts of the
	 *  calculation have to be redone).
	 *  Optional.
	 */
	PGEN_OLD_R = 1, 
	/** Set if the r-coordinate has not changed between the
	 *  first and the last point generated in the current
	 *  call.
	 *  Only for the bulk generator methods.
	 *  Optional.
	 */
	PGEN_SINGLE_R = 16,
	PGEN_AT_Z = 4, ///< Set if the point(s) lie(s) at the z-axis (theta is either 0 or M_PI).
	PGEN_AT_XY = 8, ///< Set if the point(s) lie(s) in the xy-plane (theta is M_PI2).
	PGEN_METHOD_UNAVAILABLE = 2048, ///< Set if no suitable method exists (no point generated).
	PGEN_DONE = 0, ///< Convenience identifier, not an actual flag.
	PGEN_COORDS_CART1 = QPMS_COORDS_CART1,
	PGEN_COORDS_CART2 = QPMS_COORDS_CART2,
	PGEN_COORDS_CART3 = QPMS_COORDS_CART3,
	PGEN_COORDS_POL = QPMS_COORDS_POL,
	PGEN_COORDS_SPH = QPMS_COORDS_SPH,
	PGEN_COORDS_BITRANGE = QPMS_COORDS_BITRANGE,
} PGenPointFlags;


/// Metadata generated by the fetch*() methods from PGenClassInfo
typedef struct PGenReturnDataBulk {
	/// Flags describing the returned data.
	PGenPointFlags flags;
	size_t generated; ///< Number of points really generated
} PGenReturnDataBulk;

/// Generic PGen return type that might contain point represented in any of the supported coordinate systems.
typedef struct PGenReturnData { 
	PGenPointFlags flags; ///< Metadata, must contain valid coordinate system defining flags.
	anycoord_point_t point; ///< Generated point in a coordinate system defined by flags.
} PGenReturnData;

/// PGen single-point return data type (1D).
typedef struct PGenZReturnData {
  PGenPointFlags flags; ///< Medatata.
  double point_z; ///< Generated point on a real axis.
} PGenZReturnData;

/// PGen single-point return data type (2D, polar coordinates).
typedef struct PGenPolReturnData {
  PGenPointFlags flags; ///< Metadata.
  pol_t point_pol; ///< Generated point in polar coordinates.
} PGenPolReturnData;

/// PGen single-point return data type (3D, spherical coordinates).
typedef struct PGenSphReturnData {
  PGenPointFlags flags; ///< Metadata.
  sph_t point_sph; ///< Generated point in spherical coordinates.
} PGenSphReturnData;

/// PGen single-point return data type (2D, cartesian coordinates).
typedef struct PGenCart2ReturnData {
  PGenPointFlags flags; ///< Metadata.
  cart2_t point_cart2; ///< Generated point in cartesian coordinates.
} PGenCart2ReturnData;

/// PGen single-point return data type (3D, cartesian coordinates).
typedef struct PGenCart3ReturnData {
  PGenPointFlags flags; ///< Metadata.
  cart3_t point_cart3; ///< Generated point in cartesian coordinates.
} PGenCart3ReturnData;

// convenience constants for use in the extractor implementations
static const PGenZReturnData PGenZDoneVal = {PGEN_DONE, 0}; 
static const PGenPolReturnData PGenPolDoneVal = {PGEN_DONE, {0,0}}; 
static const PGenSphReturnData PGenSphDoneVal = {PGEN_DONE, {0,0,0}}; 
static const PGenCart2ReturnData PGenCart2DoneVal = {PGEN_DONE, {0,0}};
static const PGenCart3ReturnData PGenCart3DoneVal = {PGEN_DONE, {0,0,0}};


/// PGen class metadata.
/**
 * This structure determines the behaviour of the PGen
 * instance pointing to it.
 *
 * For generating a single point, use the next() method.
 * For generating up to N points in a single call, use the 
 * fetch() method. 
 *
 * It is strongly recommended that at least the native-coordinate
 * fetch method and the native-coordinate next method are implemented.
 *
 * Usually, each generator uses internally one "native" coordinate
 * system (in lattice generators, this will typically be nD 
 * cartesian coordinates) in which the next() method gives its result.
 *
 * One does not have to explicitly implement every single method.
 *
 * TODO doc about the default transformations etc.
 */
typedef struct PGenClassInfo { // static PGenSph info
	char * const name; // mainly for debugging purposes
	int dimensionality; // lower-dimensional can be converted to higher-D, not vice versa; bit redundant with the following, whatever.
	/// Info about the generator native coordinate system.
	PGenPointFlags native_point_flags;
	/// Generate a single point in the native coordinates.
	PGenReturnData (*next)(struct PGen *);
	/// Generate a single 1D point.
	PGenZReturnData (*next_z)(struct PGen *);
	/// Generate a single 2D point in polar coordinates.
	PGenPolReturnData (*next_pol)(struct PGen *); 
	/// Generate a single 3D point in spherical coordinates.
	PGenSphReturnData (*next_sph)(struct PGen *);
	/// Generate a single 2D point in cartesian coordinates.
	PGenCart2ReturnData (*next_cart2)(struct PGen *);
	/// Generate a single 3D point in cartesian coordinates.
	PGenCart3ReturnData (*next_cart3)(struct PGen *);
	/// Generate up to \a n points in the native coordinates.
	PGenReturnDataBulk (*fetch)(struct PGen *, size_t, anycoord_point_t *);
	/// Generate up to \a n 1D points.
	PGenReturnDataBulk (*fetch_z)(struct PGen *, size_t, double *);
	/// Generate up to \a n 2D points in polar coordinates.
	PGenReturnDataBulk (*fetch_pol)(struct PGen *, size_t, pol_t *); 
	/// Generate up to \a n 3D points in spherical coordinates.
	PGenReturnDataBulk (*fetch_sph)(struct PGen *, size_t, sph_t *);
	/// Generate up to \a n 2D points in cartesian coordinates.
	PGenReturnDataBulk (*fetch_cart2)(struct PGen *, size_t, cart2_t *);
	/// Generate up to \a n 3D points in cartesian coordinates.
	PGenReturnDataBulk (*fetch_cart3)(struct PGen *, size_t, cart3_t *);
	/// Destructor.
	/** To be called by next() at iteration end, or by the caller 
	 * if ending the generation prematurely
	 */
	void (*destructor)(struct PGen *); 
} PGenClassInfo;

/// Generate a point with any of the next-methods.
static inline PGenReturnData PGen_next_nf(struct PGen *g) {
	if (g->c->next)
		return g->c->next(g);
	else {
		PGenReturnData r;
		if (g->c->next_z) {
			PGenZReturnData res = g->c->next_z(g);
			r.flags = res.flags;
			r.point.z = res.point_z;
		} else if (g->c->next_pol) {
			PGenPolReturnData res = g->c->next_pol(g);
			r.flags = res.flags;
			r.point.pol = res.point_pol;
		} else if (g->c->next_cart2) {
			PGenCart2ReturnData res = g->c->next_cart2(g);
			r.flags = res.flags;
			r.point.cart2 = res.point_cart2;
		} else if (g->c->next_sph) {
		        PGenSphReturnData res = g->c->next_sph(g);
			r.flags = res.flags;
			r.point.sph = res.point_sph;
		} else if (g->c->next_cart3) {
			PGenCart3ReturnData res = g->c->next_cart3(g);
			r.flags = res.flags;
			r.point.cart3 = res.point_cart3;
		} else 
			r.flags = PGEN_METHOD_UNAVAILABLE;
		return r;
	}
}

/// Generate multiple points with PGen in any coordinate system.
// Ultimate ugliness
static inline PGenReturnDataBulk PGen_fetch_any(struct PGen *g, size_t nmemb,
	       	anycoord_point_t *target) {
	if (g->c->fetch)
		return g->c->fetch(g, nmemb, target);
	else if (g->c->fetch_cart3) {
		cart3_t *t2 = (cart3_t*) ((char *) target 
			+ nmemb * (sizeof(anycoord_point_t)-sizeof(cart3_t)));
		PGenReturnDataBulk res = g->c->fetch_cart3(g, nmemb, t2);
		for (size_t i = 0; i < nmemb; ++i) 
			target[i].cart3 = t2[i];
		return res;
	} else if (g->c->fetch_sph) {
		sph_t *t2 = (sph_t*) ((char *) target 
			+ nmemb * (sizeof(anycoord_point_t)-sizeof(sph_t)));
		PGenReturnDataBulk res = g->c->fetch_sph(g, nmemb, t2);
		for (size_t i = 0; i < nmemb; ++i) 
			target[i].sph = t2[i];
		return res;
	} else if (g->c->fetch_cart2) {
		cart2_t *t2 = (cart2_t*) ((char *) target 
			+ nmemb * (sizeof(anycoord_point_t)-sizeof(cart2_t)));
		PGenReturnDataBulk res = g->c->fetch_cart2(g, nmemb, t2);
		for (size_t i = 0; i < nmemb; ++i) 
			target[i].cart2 = t2[i];
		return res;
	} else if (g->c->fetch_pol) {
		pol_t *t2 = (pol_t*) ((char *) target 
			+ nmemb * (sizeof(anycoord_point_t)-sizeof(pol_t)));
		PGenReturnDataBulk res = g->c->fetch_pol(g, nmemb, t2);
		for (size_t i = 0; i < nmemb; ++i) 
			target[i].pol = t2[i];
		return res;
	} else if (g->c->fetch_z) {
		double *t2 = (double*) ((char *) target 
			+ nmemb * (sizeof(anycoord_point_t)-sizeof(double)));
		PGenReturnDataBulk res = g->c->fetch_z(g, nmemb, t2);
		for (size_t i = 0; i < nmemb; ++i) 
			target[i].z = t2[i];
		return res;
	} else {
		// This is ridiculously inefficient
		PGenReturnDataBulk res = {PGEN_NOTDONE, 0};
		for (res.generated = 0; res.generated < nmemb;
			       	++res.generated) {
			PGenReturnData res1 = PGen_next_nf(g);
			QPMS_ENSURE(!(res1.flags & PGEN_METHOD_UNAVAILABLE),
				"No method found to generate points. The PGenClassInfo"
				" %s is apparently broken.", g->c->name);
			if (res1.flags & PGEN_NOTDONE) {
				target[res.generated] = res1.point;
				// The coordinate system generated by next() must be consistent:
				assert(!res.generated || ((res1.flags & PGEN_COORDS_BITRANGE) == (res.flags & PGEN_COORDS_BITRANGE)));
				res.flags |= res1.flags & PGEN_COORDS_BITRANGE;
			} else {
				res.flags &= ~PGEN_NOTDONE;
				break;
			}
		}
		// Do not guarantee anything for; low priority TODO
		res.flags &= ~(PGEN_OLD_R & PGEN_SINGLE_R);
		return res;
	}
}

/// Generate a point with any of the next-methods or fetch-methods.
static inline PGenReturnData PGen_next(struct PGen *g) {
	PGenReturnData res = PGen_next_nf(g);
	if (!(res.flags & PGEN_METHOD_UNAVAILABLE))
		return res;
	else { // Slow if implementation is stupid, but short!
		PGenReturnDataBulk resb = PGen_fetch_any(g, 1, &res.point);
		if (resb.generated) 
			// the | PGEN_NOTDONE may not be needed, but my brain melted
			res.flags = resb.flags | PGEN_NOTDONE;
		else
			res.flags = PGEN_DONE;
		return res;
	}
}

/// Generate multiple points in spherical coordinates.
static inline PGenReturnDataBulk PGen_fetch_sph(struct PGen *g, 
		size_t nmemb, sph_t *target) {
	if (g->c->fetch_sph)
		return g->c->fetch_sph(g, nmemb, target);
	else {
		anycoord_point_t *tmp;
		QPMS_CRASHING_MALLOC(tmp, sizeof(anycoord_point_t) * nmemb);
		PGenReturnDataBulk res = PGen_fetch_any(g, nmemb, tmp);
		anycoord_arr2something(target, QPMS_COORDS_SPH,
				tmp, res.flags, res.generated);
		free(tmp);
		res.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
			| QPMS_COORDS_SPH;
		return res;
	}
}

/// Generate multiple points in 3D cartesian coordinates.
static inline PGenReturnDataBulk PGen_fetch_cart3(struct PGen *g, 
		size_t nmemb, cart3_t *target) {
	if (g->c->fetch_cart3)
		return g->c->fetch_cart3(g, nmemb, target);
	else {
		anycoord_point_t *tmp;
		QPMS_CRASHING_MALLOC(tmp, sizeof(anycoord_point_t) * nmemb);
		PGenReturnDataBulk res = PGen_fetch_any(g, nmemb, tmp);
		anycoord_arr2something(target, QPMS_COORDS_CART3,
				tmp, res.flags, res.generated);
		free(tmp);
		res.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
			| QPMS_COORDS_CART3;
		return res;
	}
}

/// Generate multiple points in polar coordinates.
static inline PGenReturnDataBulk PGen_fetch_pol(struct PGen *g, 
		size_t nmemb, pol_t *target) {
	if (g->c->fetch_pol)
		return g->c->fetch_pol(g, nmemb, target);
	else {
		anycoord_point_t *tmp;
		QPMS_CRASHING_MALLOC(tmp, sizeof(anycoord_point_t) * nmemb);
		PGenReturnDataBulk res = PGen_fetch_any(g, nmemb, tmp);
		anycoord_arr2something(target, QPMS_COORDS_POL,
				tmp, res.flags, res.generated);
		free(tmp);
		res.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
			| QPMS_COORDS_POL;
		return res;
	}
}

/// Generate multiple points in 2D cartesian coordinates.
static inline PGenReturnDataBulk PGen_fetch_cart2(struct PGen *g, 
		size_t nmemb, cart2_t *target) {
	if (g->c->fetch_cart2)
		return g->c->fetch_cart2(g, nmemb, target);
	else {
		anycoord_point_t *tmp;
		QPMS_CRASHING_MALLOC(tmp, sizeof(anycoord_point_t) * nmemb);
		PGenReturnDataBulk res = PGen_fetch_any(g, nmemb, tmp);
		anycoord_arr2something(target, QPMS_COORDS_CART2,
				tmp, res.flags, res.generated);
		free(tmp);
		res.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
			| QPMS_COORDS_CART2;
		return res;
	}
}

/// Generate multiple points in 1D cartesian coordinates.
static inline PGenReturnDataBulk PGen_fetch_z(struct PGen *g, 
		size_t nmemb, double *target) {
	if (g->c->fetch_z)
		return g->c->fetch_z(g, nmemb, target);
	else {
		anycoord_point_t *tmp;
		QPMS_CRASHING_MALLOC(tmp, sizeof(anycoord_point_t) * nmemb);
		PGenReturnDataBulk res = PGen_fetch_any(g, nmemb, tmp);
		anycoord_arr2something(target, QPMS_COORDS_CART1,
				tmp, res.flags, res.generated);
		free(tmp);
		res.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
			| QPMS_COORDS_CART1;
		return res;
	}
}

/// Deallocate and invalidate a PGen point generator.
static inline void PGen_destroy(PGen *g) {
	g->c->destructor(g);
	assert(g->stateData == NULL); // this should be done by the destructor
}

/// Generate a point in a 1D real space.
static inline PGenZReturnData PGen_next_z(PGen *g) {
	if (g->c->next_z)
		return g->c->next_z(g);
	else { // Super-slow generic fallback.
		PGenReturnData res = PGen_next(g);
		if (res.flags & PGEN_NOTDONE) {
			PGenZReturnData r;
			r.point_z = anycoord2cart1(res.point, res.flags);
			r.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
				| QPMS_COORDS_CART1;
			return r;
		} else 
			return PGenZDoneVal;
	}
}

/// Generate a point in a 3D real space (spherical coordinates).
static inline PGenSphReturnData PGen_next_sph(PGen *g) {
	if (g->c->next_sph) 
		return g->c->next_sph(g);
	else { // Super-slow generic fallback.
		PGenReturnData res = PGen_next(g);
		if (res.flags & PGEN_NOTDONE) {
			PGenSphReturnData r;
			r.point_sph = anycoord2sph(res.point, res.flags);
			r.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
				| QPMS_COORDS_SPH;
			return r;
		} else 
			return PGenSphDoneVal;
	}
}

/// Generate a point in a 2D real space (polar coordinates).
static inline PGenPolReturnData PGen_next_pol(PGen *g) {
	if (g->c->next_pol) 
		return g->c->next_pol(g);
	else { // Super-slow generic fallback.
		PGenReturnData res = PGen_next(g);
		if (res.flags & PGEN_NOTDONE) {
			PGenPolReturnData r;
			r.point_pol = anycoord2pol(res.point, res.flags);
			r.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
				| QPMS_COORDS_POL;
			return r;
		} else 
			return PGenPolDoneVal;
	}
}

/// Generate a point in a 3D real space (cartesian coordinates).
static inline PGenCart3ReturnData PGen_next_cart3(PGen *g) {
	if (g->c->next_cart3) 
		return g->c->next_cart3(g);
	else { // Super-slow generic fallback.
		PGenReturnData res = PGen_next(g);
		if (res.flags & PGEN_NOTDONE) {
			PGenCart3ReturnData r;
			r.point_cart3 = anycoord2cart3(res.point, res.flags);
			r.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
				| QPMS_COORDS_CART3;
			return r;
		} else 
			return PGenCart3DoneVal;
	}
}

/// Generate a point in a 2D real space (cartesian coordinates).
static inline PGenCart2ReturnData PGen_next_cart2(PGen *g) {
	if (g->c->next_cart2) 
		return g->c->next_cart2(g);
	else { // Super-slow generic fallback.
		PGenReturnData res = PGen_next(g);
		if (res.flags & PGEN_NOTDONE) {
			PGenCart2ReturnData r;
			r.point_cart2 = anycoord2cart2(res.point, res.flags);
			r.flags = (res.flags & ~QPMS_COORDS_BITRANGE)
				| QPMS_COORDS_CART2;
			return r;
		} else 
			return PGenCart2DoneVal;
	}
}

#if 0
/// Generate up to \a n points.
static inline PGenReturnDataBulk PGen_fetch(PGen *g, size_t n, anycoord_point_t *arr){
	if (g->c->fetch)
		return g->c->fetch(g, n, arr);
	else abort();
}
#endif

static inline bool PGenSph_notDone(PGenSphReturnData data) {
	return data.flags & PGEN_NOTDONE ? true : false;
}
static inline bool PGenCart3_notDone(PGenCart3ReturnData data) {
	return data.flags & PGEN_NOTDONE ? true : false;
}

/// Standard PGen spherical coordinates -> 3d cartesian convertor.
static inline PGenCart3ReturnData PGenReturnDataConv_sph_cart3(PGenSphReturnData sphdata){
	PGenCart3ReturnData c3data;
	c3data.flags = sphdata.flags;
	c3data.point_cart3 = sph2cart(sphdata.point_sph);
	return c3data; 
}

/// Standard PGen 3d cartesian -> spherical coordinates convertor.
static inline PGenSphReturnData PGenReturnDataConv_cart3_sph(PGenCart3ReturnData c){
	PGenSphReturnData s;
	s.flags = c.flags;
	s.point_sph = cart2sph(c.point_cart3);
	return s; 
}

/*
 * Some basic lattice generators implementing the abstract interface above (implemented in latticegens.c).
 */

/// 2D point generator that simply iterates over an existing array of Point2d.
extern const PGenClassInfo PGen_FromPoint2DArray; // TODO Do I even need this to be declared here?
/// PGen_FromPoint2DArray constructor.
PGen PGen_FromPoints2DArray_new(const point2d *points, size_t len);

/// 1D equidistant point generator.
extern const PGenClassInfo PGen_1D;
typedef enum PGen_1D_incrementDirection{
    //PGEN_1D_POSITIVE_INC, // not implemented
    //PGEN_1D_NEGATIVE_INC, // not implemented
    PGEN_1D_INC_FROM_ORIGIN,
    PGEN_1D_INC_TOWARDS_ORIGIN
} PGen_1D_incrementDirection;
/// PGen_1D point generator constructor.
PGen PGen_1D_new_minMaxR(double period, ///< Distance between points.
	double offset, ///< Lattice offset from zero.
	double minR, ///< Lower bound of |z| of the generated points.
	bool inc_minR, ///< Include also |z| == minR (if applicable).
	double maxR, ///< Upper bound of |z| of the generated points.
	bool inc_maxR, ///< Include also |z| == maxR if applicable.
	PGen_1D_incrementDirection incdir ///< Order of generated points.
);

extern const PGenClassInfo PGen_xyWeb;
PGen PGen_xyWeb_new(cart2_t b1, cart2_t b2, double rtol, cart2_t offset, 
		double minR, bool inc_minR, double maxR, bool inc_maxR);

/// Returns a number larger or equal than the number of all the points generated by a PGen_xyWeb.
size_t PGen_xyWeb_sizecap(cart2_t b1, cart2_t b2, double rtol, cart2_t offset,
		double minR, bool inc_minR, double maxR, bool inc_maxR);


extern const PGenClassInfo PGen_LatticeRadialHeap2D;
extern const PGenClassInfo PGen_LatticeRadialHeap3D;
PGen PGen_LatticeRadialHeap2D_new(cart2_t b1, cart2_t b2, cart2_t offset, 
		double minR, bool inc_minR, double maxR, bool inc_maxR);
PGen PGen_LatticeRadialHeap3D_new(const cart3_t *b1, const cart3_t *b2, const cart3_t *b3,
	       	const cart3_t *offset, double minR, bool inc_minR, double maxR, bool inc_maxR);


/// A metagenerator generating points from another generator shifted by a constant.
extern const PGenClassInfo PGen_shifted;
PGen PGen_shifted_new(PGen orig, cart3_t shift);

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

// Just for detecting the right angles (needed for generators).
typedef enum {
	NOT_ORTHOGONAL = 0,
	ORTHOGONAL_01 = 1,
	ORTHOGONAL_12 = 2,
	ORTHOGONAL_02 = 4
} LatticeFlags;



/*
 * Lagrange-Gauss reduction of a 2D basis.
 * The output shall satisfy |out1| <= |out2| <= |out2 - out1|
 */
void l2d_reduceBasis(cart2_t in1, cart2_t in2, cart2_t *out1, cart2_t *out2);

// This one uses LLL reduction.
void l3d_reduceBasis(const cart3_t in[3], cart3_t out[3]);

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

// Detects right angles.
LatticeFlags l2d_detectRightAngles(cart2_t b1, cart2_t b2, double rtol);
LatticeFlags l3d_detectRightAngles(const cart3_t basis[3], double rtol);

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

// Reciprocal bases; returns 0 on success, possibly a non-zero if b1 and b2 are parallel
int l2d_reciprocalBasis1(cart2_t b1, cart2_t b2, cart2_t *rb1, cart2_t *rb2);
int l2d_reciprocalBasis2pi(cart2_t b1, cart2_t b2, cart2_t *rb1, cart2_t *rb2);
// 3D reciprocal bases; returns (direct) unit cell volume with possible sign. Assumes direct lattice basis already reduced.
double l3d_reciprocalBasis1(const cart3_t direct_basis[3], cart3_t reciprocal_basis[3]);
double l3d_reciprocalBasis2pi(const cart3_t direct_basis[3], cart3_t reciprocal_basis[3]);

double l2d_unitcell_area(cart2_t b1, cart2_t b2);

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


#if 0
/// Gives the frequency of \a n-th empty lattice mode at a given wave vector \a k.
double qpms_emptylattice2_mode_nth(
                cart2_t b1_rec, ///< First reciprocal lattice base vector
                cart2_t b2_rec, ///< Second reciprocal lattice base vector
                double rtol, ///< Relative tolerance to detect right angles
                cart2_t k, ///< The wave vector
                double wave_speed, ///< Wave speed in a given medium (i.e. vacuum speed / refractive index).
                size_t N ///< Index of the mode (note that degenerate modes are counted multiple times).
                );

/// Gives the first `maxindex` frequencies of empty lattice modes at a given wave vector \a k.
void qpms_emptylattice2_modes_maxindex(
                double target_freqs[], ///< Target array of size maxindex.
                cart2_t b1_rec, ///< First reciprocal lattice base vector
                cart2_t b2_rec, ///< Second reciprocal lattice base vector
                double rtol, ///< Relative tolerance to detect right angles
                cart2_t k, ///< The wave vector
                double wave_speed, ///< Wave speed in a given medium (i.e. vacuum speed / refractive index).
                size_t maxindex ///< Number of the frequencies generated.
                );
#endif

/// Gives the frequencies of empty lattice modes at a given wave vector \a k up to \a maxfreq and one more.
/**
 * The frequencies are saved to a newly allocated array *target_freqs (to be deallocated
 * using free() by the caller).
 *
 * \returns Number of found mode frequencies lower or equal than \a maxfreq plus one.
 */
size_t qpms_emptylattice2_modes_maxfreq(
                double **target_freqs,
                cart2_t b1_rec, ///< First reciprocal lattice base vector
                cart2_t b2_rec, ///< Second reciprocal lattice base vector
                double rtol, ///< Relative tolerance to detect right angles
                cart2_t k, ///< The wave vector
                double wave_speed, ///< Wave speed in a given medium (i.e. vacuum speed / refractive index).
                double maxfreq ///< The maximum frequency.
                );

/// Gives the frequencies of two empty lattice modes nearest to \a omega at a given wave vector \a k.
void qpms_emptylattice2_modes_nearest(
                double target[2], ///< Target array with lower ([0]) and upper ([1]) frequency.
                cart2_t b1_rec, ///< First reciprocal lattice base vector
                cart2_t b2_rec, ///< Second reciprocal lattice base vector
                double rtol, ///< Relative tolerance to detect right angles
                cart2_t k, ///< The wave vector
                double wave_speed, ///< Wave speed in a given medium (i.e. vacuum speed / refractive index).
                double omega ///< The frequency around which the frequencies are searched.
                );



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
