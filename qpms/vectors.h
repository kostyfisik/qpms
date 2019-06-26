/*! \file vectors.h
 * \brief Coordinate transforms and vector arithmetics.
 */
#ifndef VECTORS_H
#define VECTORS_H
#include <math.h>
#ifndef M_PI_2
#define M_PI_2 (1.570796326794896619231321691639751442098584699687552910487)
#endif
#include "qpms_types.h"
#include "qpms_error.h"

//static inline double vectors_h_sq(double x) {return x*x;}

static const cart2_t CART2_ZERO = {0, 0};
static const cart3_t CART3_ZERO = {0, 0, 0};

/// 2D vector addition.
static inline cart2_t cart2_add(const cart2_t a, const cart2_t b) {
	cart2_t res = {a.x+b.x, a.y+b.y};
	return res;
}

/// 2D vector substraction.
static inline cart2_t cart2_substract(const cart2_t a, const cart2_t b) {
	cart2_t res = {a.x-b.x, a.y-b.y};
	return res;
}

/// 2D vector scaling.
static inline cart2_t cart2_scale(const double c, const cart2_t v) {
	cart2_t res = {c * v.x, c * v.y};
	return res;
}

/// 2D vector dot product.
static inline double cart2_dot(const cart2_t a, const cart2_t b) {
	return a.x * b.x + a.y * b.y;
}

static inline double cart2_normsq(const cart2_t a) {
	return cart2_dot(a, a);
}

/// 2D vector euclidian norm.
static inline double cart2norm(const cart2_t v) {
	return hypot(v.x, v.y); //sqrt(v.x*v.x + v.y*v.y);
}

/// 2D cartesian to polar coordinates conversion. See @ref coord_conversions.
static inline pol_t cart2pol(const cart2_t cart) {
	pol_t pol;
	pol.r = cart2norm(cart);
	pol.phi = atan2(cart.y, cart.x);
	return pol;
}

/// Polar to spherical coordinates conversion. See @ref coord_conversions.
static inline sph_t pol2sph_equator(const pol_t pol) {
	sph_t sph;
	sph.r = pol.r;
	sph.phi = pol.phi;
	sph.theta = M_PI_2;
	return sph;
}

/// 2D cartesian to spherical coordinates conversion. See @ref coord_conversions.
static inline sph_t cart22sph(const cart2_t cart) {
	sph_t sph;
	sph.r = cart2norm(cart);
	sph.theta = M_PI_2;
	sph.phi = atan2(cart.y, cart.x);
	return sph;
}

/// 1D cartesian to spherical coordinates conversion. See @ref coord_conversions.
static inline sph_t cart12sph_zaxis(double z) {
	sph_t sph = {fabs(z), z < 0 ? M_PI : 0, 0};
	return sph;
}

/// 1D to 3D cartesian coordinates conversion. See @ref coord_conversions.
static inline cart3_t cart12cart3z(double z) {
	cart3_t c = {0, 0, z};
	return c;
}

/// 2D to 3D cartesian coordinates conversion. See @ref coord_conversions.
static inline cart3_t cart22cart3xy(const cart2_t a) {
	cart3_t c;
	c.x = a.x;
	c.y = a.y;
	c.z = 0;
	return c;
}

static inline cart2_t cart3xy2cart2(const cart3_t a) {
	cart2_t c = {a.x, a.y};
	return c;
}

/// 3D vector euclidian norm.
static inline double cart3norm(const cart3_t v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}


/// 3D cartesian to spherical coordinates conversion. See @ref coord_conversions.
static inline sph_t cart2sph(const cart3_t cart) {
	sph_t sph;
	sph.r = cart3norm(cart);
	sph.theta = sph.r ? acos(cart.z / sph.r) : M_PI_2;
	sph.phi = atan2(cart.y, cart.x);
	return sph;
}

/// Spherical to 3D cartesian coordinates conversion. See @ref coord_conversions.
static inline cart3_t sph2cart(const sph_t sph) {
	cart3_t cart;
	double sin_th = 
#ifdef QPMS_VECTORS_NICE_TRANSFORMATIONS
	   (sph.theta == M_PI) ? 0 :
#endif
		sin(sph.theta);
	cart.x = sph.r * sin_th * cos(sph.phi);
	cart.y = sph.r * sin_th * sin(sph.phi);
	cart.z = sph.r * cos(sph.theta);
	return cart;
}

/// Polar to 2D cartesian coordinates conversion. See @ref coord_conversions.
static inline cart2_t pol2cart(const pol_t pol) {
	cart2_t cart;
	cart.x = pol.r * cos(pol.phi);
	cart.y = pol.r * sin(pol.phi);
	return cart;
}

/// Polar to 3D cartesian coordinates conversion. See @ref coord_conversions.
static inline cart3_t pol2cart3_equator(const pol_t pol) {
	cart2_t c = pol2cart(pol);
	cart3_t cart3 = {c.x, c.y, 0};
	return cart3;
}

/// 3D vector addition.
static inline cart3_t cart3_add(const cart3_t a, const cart3_t b) {
	cart3_t res = {a.x+b.x, a.y+b.y, a.z+b.z};
	return res;
}

/// 3D vector substraction.
static inline cart3_t cart3_substract(const cart3_t a, const cart3_t b) {
	cart3_t res = {a.x-b.x, a.y-b.y, a.z-b.z};
	return res;
}

/// 3D vector scaling
static inline cart3_t cart3_scale(const double c, const cart3_t v) {
	cart3_t res = {c * v.x, c * v.y, c * v.z};
	return res;
}

/// Euclidian distance between two 3D points.
static inline double cart3_dist(const cart3_t a, const cart3_t b) {
	return cart3norm(cart3_substract(a,b));
}

static inline bool cart3_isclose(const cart3_t a, const cart3_t b, double rtol, double atol) {
	return cart3_dist(a,b) <= atol + rtol * (cart3norm(b) + cart3norm(a)) * .5;
}

/// Complex 3D vector scaling.
static inline ccart3_t ccart3_scale(const complex  double c, const ccart3_t v) {
	ccart3_t res = {c * v.x, c * v.y, c * v.z};
	return res;
}

/// Complex 3D vector adition.
static inline ccart3_t ccart3_add(const ccart3_t a, const ccart3_t b) {
	ccart3_t res = {a.x+b.x, a.y+b.y, a.z+b.z};
	return res;
}

/// Complex 3D vector substraction.
static inline ccart3_t ccart3_substract(const ccart3_t a, const ccart3_t b) {
	ccart3_t res = {a.x-b.x, a.y-b.y, a.z-b.z};
	return res;
}

/// Complex 3D vector (geographic coordinates) addition.
static inline csphvec_t csphvec_add(const csphvec_t a, const csphvec_t b) {
	csphvec_t res = {a.rc + b.rc, a.thetac + b.thetac, a.phic + b.phic};
	return res;
}

/// Complex 3D vector (geographic coordinates) substraction.
static inline csphvec_t csphvec_substract(const csphvec_t a, const csphvec_t b) {
	csphvec_t res = {a.rc - b.rc, a.thetac - b.thetac, a.phic - b.phic};
	return res;
}

/// Complex 3D vector (geographic coordinates) scaling.
static inline csphvec_t csphvec_scale(complex double c, const csphvec_t v) {
	csphvec_t res = {c * v.rc, c * v.thetac, c * v.phic};
	return res;
}

/// Complex 3D vector (geographic coordinates) "dot product" without conjugation.
static inline complex double csphvec_dotnc(const csphvec_t a, const csphvec_t b) {
	//N.B. no complex conjugation done here
	return a.rc * b.rc + a.thetac * b.thetac + a.phic * b.phic;
}

/// 3D vector dot product.
static inline double cart3_dot(const cart3_t a, const cart3_t b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// Spherical coordinate system scaling.
static inline sph_t sph_scale(double c, const sph_t s) {
	sph_t res = {c * s.r, s.theta, s.phi};
	return res;
}

/// "Complex spherical" coordinate system scaling.
static inline csph_t sph_cscale(complex double c, const sph_t s) {
	csph_t res = {c * s.r, s.theta, s.phi};
	return res;
}

/// Coordinate transform of a vector in local geographic to global cartesian system.
// equivalent to sph_loccart2cart in qpms_p.py
static inline ccart3_t csphvec2ccart(const csphvec_t sphvec, const sph_t at) {
	const double st = sin(at.theta);
	const double ct = cos(at.theta);
	const double sf = sin(at.phi);
	const double cf = cos(at.phi);
	const double rx = st * cf;
	const double ry = st * sf;
	const double rz = ct;
	const double tx = ct * cf;
	const double ty = ct * sf;
	const double tz = -st;
	const double fx = -sf;
	const double fy = cf;
	const double fz = 0.;
	ccart3_t res;
	res.x = rx * sphvec.rc + tx * sphvec.thetac + fx * sphvec.phic;
	res.y = ry * sphvec.rc + ty * sphvec.thetac + fy * sphvec.phic;
	res.z = rz * sphvec.rc + tz * sphvec.thetac + fz * sphvec.phic;
	return res;
}

/// Coordinate transform of a vector in local geographic to global cartesian system.
/**
 *  Same as csphvec2ccart, but with csph_t as second argument.
 *  (The radial part (which is the only complex part of csph_t) 
 *  of the second argument does not play role in the
 *  transformation, so this is completely legit
 */
static inline ccart3_t csphvec2ccart_csph(const csphvec_t sphvec, const csph_t at) {
	const sph_t atreal = {0 /*not used*/, at.theta, at.phi};
	return csphvec2ccart(sphvec, atreal);
}

/// Coordinate transform of a vector in global cartesian to local geographic system.
static inline csphvec_t ccart2csphvec(const ccart3_t cartvec, const sph_t at) {
	// this chunk is copy-pasted from csphvec2cart, so there should be a better way...
	const double st = sin(at.theta);
	const double ct = cos(at.theta);
	const double sf = sin(at.phi);
	const double cf = cos(at.phi);
	const double rx = st * cf;
	const double ry = st * sf;
	const double rz = ct;
	const double tx = ct * cf;
	const double ty = ct * sf;
	const double tz = -st;
	const double fx = -sf;
	const double fy = cf;
	const double fz = 0.;
	csphvec_t res;
	res.rc     = rx * cartvec.x + ry * cartvec.y + rz * cartvec.z;
	res.thetac = tx * cartvec.x + ty * cartvec.y + tz * cartvec.z;
	res.phic   = fx * cartvec.x + fy * cartvec.y + fz * cartvec.z;
	return res;
}

void print_csphvec(csphvec_t);
void print_ccart3(ccart3_t);
void print_cart3(cart3_t);
void print_sph(sph_t);

// kahan sums for various types... TODO make generic code using macros

/// Kanan sum initialisation for ccart3_t.
static inline void ccart3_kahaninit(ccart3_t *sum, ccart3_t *compensation) {
	sum->x = sum->y = sum->z = compensation->x = compensation->y = compensation->z = 0;
}
/// Kanan sum initialisation for csphvec_t.
static inline void csphvec_kahaninit(csphvec_t *sum, csphvec_t *compensation) {
	sum->rc = sum->thetac = sum->phic = compensation->rc = compensation->thetac = compensation->phic = 0;
}

/// Add element to Kahan sum (ccart3_t).
static inline void ccart3_kahanadd(ccart3_t *sum, ccart3_t *compensation, const ccart3_t input) {
	ccart3_t comped_input = ccart3_substract(input, *compensation);
	ccart3_t nsum = ccart3_add(*sum, comped_input);
	*compensation = ccart3_substract(ccart3_substract(nsum, *sum), comped_input);
	*sum = nsum;
}

/// Add element to Kahan sum (csphvec_t).
static inline void csphvec_kahanadd(csphvec_t *sum, csphvec_t *compensation, const csphvec_t input) {
	csphvec_t comped_input = csphvec_substract(input, *compensation);
	csphvec_t nsum = csphvec_add(*sum, comped_input);
	*compensation = csphvec_substract(csphvec_substract(nsum, *sum), comped_input);
	*sum = nsum;
}

/// Euclidian norm of a vector in geographic coordinates.
static inline double csphvec_norm(const csphvec_t a) {
	return  sqrt(creal(a.rc * conj(a.rc) + a.thetac * conj(a.thetac) + a.phic * conj(a.phic)));
}

static inline double csphvec_reldiff_abstol(const csphvec_t a, const csphvec_t b, double tolerance) {
	double anorm = csphvec_norm(a);
	double bnorm = csphvec_norm(b);
	if (anorm <= tolerance && bnorm <= tolerance) return 0;
	return csphvec_norm(csphvec_substract(a,b)) / (anorm + bnorm);
}

static inline double csphvec_reldiff(const csphvec_t a, const csphvec_t b) {
	return csphvec_reldiff_abstol(a, b, 0);
}


/*! \page coord_conversions Coordinate systems and default conversions
 *
 * The coordinate system transformations are defined as following:
 * 
 * \section coordtf_same_d Equal-dimension coordinate tranforms
 * \subsection sph_cart3 Spherical and 3D cartesian coordinates
 * * \f$ x = r \sin \theta \cos \phi \f$,
 * * \f$ y = r \sin \theta \sin \phi \f$,
 * * \f$ z = r \cos \theta \f$.
 * \subsection pol_cart2 Polar and 2D cartesian coordinates
 * * \f$ x = r \cos \phi \f$,
 * * \f$ y = r \sin \phi \f$.
 *
 * \section coordtf_123 Lower to higher dimension conversions.
 * * The 1D coordinate is identified with the \a z 3D cartesian coordinate.
 * * The 2D cartesian coordinates \a x, \a y are identified with the \a x, \a y
 *   3D cartesian coordinates.
 * * For the sake of consistency, default conversion between
 *   1D and 2D coordinates is not allowed and yields NAN values.
 *
 * \section coordtf_321 Higher to lower dimension conversions.
 * Default conversions from higher to lower-dimensional coordinate
 * systems are not allowed. Any projections have to be done explicitly.
 */

/// Conversion from anycoord_point_t to explicitly spherical coordinates.
/** See @ref coord_conversions for the conversion definitions.
 */ 
static inline sph_t anycoord2sph(anycoord_point_t p, qpms_coord_system_t t) {
	switch(t & QPMS_COORDS_BITRANGE) {
		case QPMS_COORDS_SPH:
			return p.sph;
			break;
		case QPMS_COORDS_POL:
			return pol2sph_equator(p.pol);
			break;
		case QPMS_COORDS_CART3:
			return cart2sph(p.cart3);
			break;
		case QPMS_COORDS_CART2:
			return cart22sph(p.cart2);
			break;
		case QPMS_COORDS_CART1:
			return cart12sph_zaxis(p.z);
			break;
	}
	QPMS_WTF;
}


/// Conversion from anycoord_point_t to explicitly 3D cartesian coordinates.
/** See @ref coord_conversions for the conversion definitions.
 */ 
static inline cart3_t anycoord2cart3(anycoord_point_t p, qpms_coord_system_t t) {
	switch(t & QPMS_COORDS_BITRANGE) {
		case QPMS_COORDS_SPH:
			return sph2cart(p.sph);
			break;
		case QPMS_COORDS_POL:
			return pol2cart3_equator(p.pol);
			break;
		case QPMS_COORDS_CART3:
			return p.cart3;
			break;
		case QPMS_COORDS_CART2:
			return cart22cart3xy(p.cart2);
			break;
		case QPMS_COORDS_CART1:
			return cart12cart3z(p.z);
			break;
	}
	QPMS_WTF;
}

#if 0
// Convenience identifiers for return values.
static const cart3_t CART3_INVALID = {NAN, NAN, NAN};
static const cart2_t CART2_INVALID = {NAN, NAN};
static const double CART1_INVALID = NAN;
static const sph_t SPH_INVALID = {NAN, NAN, NAN};
static const pol_t POL_INVALID = {NAN, NAN};
#endif

/// Conversion from anycoord_point_t to explicitly polar coordinates.
/** See @ref coord_conversions for the conversion definitions.
 */ 
static inline pol_t anycoord2pol(anycoord_point_t p, qpms_coord_system_t t) {
	switch(t & QPMS_COORDS_BITRANGE) {
		case QPMS_COORDS_SPH:
		case QPMS_COORDS_CART3:
			QPMS_PR_ERROR("Implicit conversion from 3D to 2D"
				       	" coordinates not allowed");
			break;
		case QPMS_COORDS_POL:
			return p.pol;
			break;
		case QPMS_COORDS_CART2:
			return cart2pol(p.cart2);
			break;
		case QPMS_COORDS_CART1:
			QPMS_PR_ERROR("Implicit conversion from 1D to 2D"
					" coordinates not allowed");
			break;
	}
	QPMS_WTF;
}


/// Conversion from anycoord_point_t to explicitly 2D cartesian coordinates.
/** See @ref coord_conversions for the conversion definitions.
 */ 
static inline cart2_t anycoord2cart2(anycoord_point_t p, qpms_coord_system_t t) {
	switch(t & QPMS_COORDS_BITRANGE) {
		case QPMS_COORDS_SPH:
		case QPMS_COORDS_CART3:
			QPMS_PR_ERROR("Implicit conversion from 3D to 2D"
				       	" coordinates not allowed");
			break;
		case QPMS_COORDS_POL:
			return pol2cart(p.pol);
			break;
		case QPMS_COORDS_CART2:
			return p.cart2;
			break;
		case QPMS_COORDS_CART1:
			QPMS_PR_ERROR("Implicit conversion from 1D to 2D"
					" coordinates not allowed");
			break;
	}
	QPMS_WTF;
}


/// Conversion from anycoord_point_t to explicitly 1D cartesian coordinates.
/** See @ref coord_conversions for the conversion definitions.
 */ 
static inline double anycoord2cart1(anycoord_point_t p, qpms_coord_system_t t) {
	if (t & QPMS_COORDS_BITRANGE == QPMS_COORDS_CART1) 
			return p.z;
	else
		QPMS_PR_ERROR("Implicit conversion from nD (n > 1)"
				" to 1D not allowed.");
}

typedef double matrix3d[3][3];
typedef double matrix2d[2][2];
typedef complex double cmatrix3d[3][3];
typedef complex double cmatrix2d[2][2];
#endif //VECTORS_H
