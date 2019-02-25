/*! \file qpms_types.h
 * \brief Common qpms types.
 */
#ifndef QPMS_TYPES_H
#define QPMS_TYPES_H
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

#ifndef M_PI_2
#define M_PI_2 (1.570796326794896619231321691639751442098584699687552910487)
#endif
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// integer index types
typedef int qpms_lm_t;
/// Type for spherical harmonic degree l.
typedef int qpms_l_t; /* can't be unsigned because of the behaviour under - operator;
			 also -1 needed as an invalid value for scalar waves. */

/// Type for spherical harmonic order m.
typedef qpms_lm_t qpms_m_t;

/// Type for the (l,m) multiindex of transversal (M or N-type) VSWFs.
/** This corresponds to the typical memory layout for various coefficient etc.
 *  Corresponds to the l-primary, m-secondary ordering, i.e.
 *  \f[ y = 0: l = 1, m = -1, \f]
 *  \f[ y = 1: l = 1, m =  0, \f]
 *  \f[ y = 2: l = 1, m = +1, \f]
 *  \f[ y = 3: l = 2, m = -2, \f]
 *  ...
 */
typedef size_t qpms_y_t;

/// Type for the (l,m) multiindex of spherical harmonics, including (0,0).
/** This differs from qpms_y_t by being shifted by one and including
 *  the l = 0 option. Suitable also for scalar and longitudinal waves.
 *  Corresponds to the l-primary, m-secondary ordering, i.e.
 *  \f[ y = 0: l = 0, m = 0,  \f]
 *  \f[ y = 1: l = 1, m = -1, \f]
 *  \f[ y = 2: l = 1, m = 0,  \f]
 *  \f[ y = 3: l = 1, m = +1, \f]
 *  \f[ y = 4: l = 2, m = -2, \f]
 *  ...
 */
typedef size_t qpms_y_sc_t;

/// Codes of the VSWF types (electric/N, magnetic/M, longitudinal/L).
typedef enum {
	QPMS_VSWF_ELECTRIC = 2, ///< "Electric" ($N$-type) transversal wave.
	QPMS_VSWF_MAGNETIC = 1, ///< "Magnetic" ($M$-type) transversal wave.
	QPMS_VSWF_LONGITUDINAL = 0 ///< Longitudinal ($L$-type) wave (not relevant for radiation).
} qpms_vswf_type_t;


/// Exhaustive index type for VSWF basis functions.
/** Carries information about the wave being of M/N/L (magnetic, electric,
 *  or longitudinal) type, as well as the wave's degree and order (l, m).
 *
 *  The formula is 4 * (qpms_y_sc_t) y_sc + (qmps_vswf_type_t) type_code,
 *  but don't rely on this and use the functions
 *  qpms_tmn2uvswfi() and qpms_uvswfi2tmn()
 *  from qpms_types.h instead
 *  as the formula might change in future versions.
 */
typedef size_t qpms_uvswfi_t; 

/// Error codes / return values for certain numerical functions.
/** These are de facto a subset of the GSL error codes. */
typedef enum {
	QPMS_SUCCESS = 0, ///< Success.
	QPMS_ERROR = 1, ///< Unspecified error.
	QPMS_ENOMEM = 8 ///< Out of memory.
} qpms_errno_t;

/// Vector spherical wavefuction normalisation (and sign) convention codes.
/** Throughout the literature, various conventions for VSWF bases are used.
 *  The meaningful ones are the "power" and "spherical harmonic" normalisation
 *  conventions, as the (l,m) and (l,-m) waves of the same type have the same
 *  intensities.
 *  One might also encounter a very inconvenient and messy "antinormalisation"
 *  used in Xu (TODO reference).
 *
 *  Moreover, VSWFs might use various sign convention. Usually they either
 *  carry the Condon-Shortley phase $(-1)^m$ or not, which is also saved here.
 *
 *  TODO references and exact definitions.
 */
//const int QPMS_NORMALISATION_T_CSBIT = 128;
#define QPMS_NORMALISATION_T_CSBIT 128
typedef enum {
#ifdef USE_XU_ANTINORMALISATION
	// As in TODO
	QPMS_NORMALISATION_XU = 4, ///< such that the numerical values in Xu's tables match, not recommended to use otherwise
	QPMS_NORMALISATION_XU_CS = QPMS_NORMALISATION_XU | QPMS_NORMALISATION_T_CSBIT, 
#endif
	QPMS_NORMALISATION_NONE = 3, ///< genuine unnormalised waves (with unnormalised Legendre polynomials)
	QPMS_NORMALISATION_KRISTENSSON = 2, ///< As in http://www.eit.lth.se/fileadmin/eit/courses/eit080f/Literature/book.pdf, power-normalised
	QPMS_NORMALISATION_POWER = QPMS_NORMALISATION_KRISTENSSON, 
	// as in TODO
	QPMS_NORMALISATION_TAYLOR = 1,
	QPMS_NORMALISATION_SPHARM = QPMS_NORMALISATION_TAYLOR,
	// Variants with Condon-Shortley phase
	QPMS_NORMALISATION_NONE_CS = QPMS_NORMALISATION_NONE | QPMS_NORMALISATION_T_CSBIT,
	QPMS_NORMALISATION_KRISTENSSON_CS = QPMS_NORMALISATION_KRISTENSSON | QPMS_NORMALISATION_T_CSBIT, 
	QPMS_NORMALISATION_POWER_CS = QPMS_NORMALISATION_KRISTENSSON_CS,
	QPMS_NORMALISATION_TAYLOR_CS = QPMS_NORMALISATION_TAYLOR | QPMS_NORMALISATION_T_CSBIT,
	QPMS_NORMALISATION_SPHARM_CS = QPMS_NORMALISATION_TAYLOR_CS,
	QPMS_NORMALISATION_UNDEF = 0
} qpms_normalisation_t;

/// Determine whether the convention includes Condon-Shortley phase (-1) or not (+1).
static inline int qpms_normalisation_t_csphase(qpms_normalisation_t norm) {
	return (norm & QPMS_NORMALISATION_T_CSBIT)? -1 : 1;
}

/// Returns the normalisation convention code without the Condon-Shortley phase.
static inline int qpms_normalisation_t_normonly(qpms_normalisation_t norm) {
	return norm & (~QPMS_NORMALISATION_T_CSBIT);
}

// TODO move the inlines elsewhere
/* Normalisation of the spherical waves is now scattered in at least three different files:
 * here, we have the norm in terms of radiated power of outgoing wave.
 * In file legendre.c, function qpms_pitau_get determines the norm used in the vswf.c
 * spherical vector wave norms. The "dual" waves in vswf.c use the ..._abssquare function below.
 * In file translations.c, the normalisations are again set by hand using the normfac and lognormfac
 * functions.
 */
#include <math.h>
#include <assert.h>
// relative to QPMS_NORMALISATION_KRISTENSSON_CS, i.e.
// P_l^m[normtype] = P_l^m[Kristensson]
static inline double qpms_normalisation_t_factor(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	int csphase = qpms_normalisation_t_csphase(norm);
	norm = qpms_normalisation_t_normonly(norm);
	double factor;
	switch (norm) {
		case QPMS_NORMALISATION_KRISTENSSON:
			factor = 1.;
			break;
		case QPMS_NORMALISATION_TAYLOR:
			factor = sqrt(l*(l+1));
			break;
		case QPMS_NORMALISATION_NONE:
			factor = sqrt(l*(l+1) * 4 * M_PI / (2*l+1) * exp(lgamma(l+m+1)-lgamma(l-m+1)));
			break;
#ifdef USE_XU_ANTINORMALISATION // broken probably in legendre.c
		case QPMS_NORMALISATION_XU:
			factor = sqrt(4 * M_PI) / (2*l+1) * exp(lgamma(l+m+1)-lgamma(l-m+1));
			break;
#endif
		default:
			assert(0);
	}
	factor *= (m%2)?(-csphase):1;
	return factor;
}


// TODO move elsewhere
static inline double qpms_normalisation_t_factor_abssquare(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	norm = qpms_normalisation_t_normonly(norm);
	switch (norm) {
		case QPMS_NORMALISATION_KRISTENSSON:
			return 1.;
			break;
		case QPMS_NORMALISATION_TAYLOR:
			return l*(l+1);
			break;
		case QPMS_NORMALISATION_NONE:
			return l*(l+1) * 4 * M_PI / (2*l+1) * exp(lgamma(l+m+1)-lgamma(l-m+1));
			break;
#ifdef USE_XU_ANTINORMALISATION // broken probably in legendre.c
		case QPMS_NORMALISATION_XU:
			{
			  double fac = sqrt(4 * M_PI) / (2*l+1) * exp(lgamma(l+m+1)-lgamma(l-m+1));
			  return fac * fac;
			}
			break;
#endif 
		default:
			assert(0);
			return NAN;
	}
}


/// Bessel function kinds.
typedef enum {
	QPMS_BESSEL_REGULAR = 1, ///< regular (spherical) Bessel function \a j (Bessel function of the first kind)
	QPMS_BESSEL_SINGULAR = 2, ///< singular (spherical)  Bessel function \a y (Bessel function of the second kind)
	QPMS_HANKEL_PLUS = 3, ///< (spherical) Hankel function \f$ h_1 = j + iy \f$
	QPMS_HANKEL_MINUS = 4, ///< (spherical) Hankel function \f$ h_2 = j - iy \f$
	QPMS_BESSEL_UNDEF = 0 ///< invalid / unspecified kind
} qpms_bessel_t;

// coordinate system types
/// 3D cartesian coordinates.
typedef struct cart3_t {
	double x, y, z;
} cart3_t;

/// 3D complex (actually 6D) coordinates.
typedef struct ccart3_t {
	complex double x, y, z;
} ccart3_t;

/// 2D cartesian coordinates.
typedef struct cart2_t {
	double x, y;
} cart2_t;

/// Spherical coordinates.
typedef struct sph_t {
	double r, theta, phi;
} sph_t;

/// Spherical coordinates with complex radial component.
typedef struct csph_t { // Do I really need this???
	complex double r;
	double	theta, phi;
} csph_t;

/// 3D complex vector components in local spherical basis.
typedef struct csphvec_t {
	complex double rc, thetac, phic; 
} csphvec_t;

/// 2D polar coordinates.
typedef struct pol_t {
	double r, phi;
} pol_t;

/// Union type capable to contain various 1D, 2D and 3D coordinates.
typedef union anycoord_point_t {
	double z; /// 1D cartesian coordinate.
	cart3_t cart3;
	cart2_t cart2;
	sph_t sph;
	pol_t pol;
} anycoord_point_t;

/// Enum codes for common coordinate systems.
typedef enum { 
	// IF EVER CHANGING THE CONSTANT VALUES HERE, 
	// CHECK THAT THEY DO NOT CLASH WITH THOSE IN PGenPointFlags!
	QPMS_COORDS_CART1 = 64, ///< 1D cartesian (= double).
	QPMS_COORDS_POL = 128, ///< 2D polar.
	QPMS_COORDS_SPH = 256, ///< 3D spherical.
	QPMS_COORDS_CART2 = 512, ///< 2D cartesian.
	QPMS_COORDS_CART3 = 1024, ///< 3D cartesian.
} qpms_coord_system_t;

/// Quaternion type.
/**
 * Internaly represented as a pair of complex numbers,
 * \f$ Q_a = Q_1 + iQ_z, Q_b = Q_y + i Q_x\f$.
 *
 * Check wigner.h for "methods".
 */
typedef struct qpms_quat_t {
        complex double a, b;
} qpms_quat_t;

/// Quaternion type as four doubles.
/** Check wigner.h for "methods".
 */
typedef struct qpms_quat4d_t {
        double c1, ci, cj, ck;
} qpms_quat4d_t;

/// 3D improper rotations represented as a quaternion and a sign of the determinant.
/** Check wigner.h for "methods".
 */
typedef struct qpms_irot3_t {
        qpms_quat_t rot; ///< Quaternion representing the rotation part.
        short det; ///< Determinant of the transformation (valid values are 1 (rotation) or -1 (improper rotation)
} qpms_irot3_t;


#define lmcheck(l,m) assert((l) >= 1 && abs(m) <= (l))
#endif // QPMS_TYPES
