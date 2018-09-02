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
typedef int qpms_l_t; // can't be unsigned because of the behaviour under - operator
typedef qpms_lm_t qpms_m_t;
typedef size_t qpms_y_t;

typedef enum {
	QPMS_SUCCESS = 0,
	QPMS_ERROR = 1,
	QPMS_ENOMEM = 8
} qpms_errno_t;

// Normalisations
//const int QPMS_NORMALISATION_T_CSBIT = 128;
#define QPMS_NORMALISATION_T_CSBIT 128
typedef enum {
#ifdef USE_XU_ANTINORMALISATION
	// As in TODO
	QPMS_NORMALISATION_XU = 4, // such that the numerical values in Xu's tables match, not recommended to use otherwise
	QPMS_NORMALISATION_XU_CS = QPMS_NORMALISATION_XU | QPMS_NORMALISATION_T_CSBIT, 
#endif
	QPMS_NORMALISATION_NONE = 3, // genuine unnormalised waves (with unnormalised Legendre polynomials)
	// As in http://www.eit.lth.se/fileadmin/eit/courses/eit080f/Literature/book.pdf, power-normalised
	QPMS_NORMALISATION_KRISTENSSON = 2,
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




static inline int qpms_normalisation_t_csphase(qpms_normalisation_t norm) {
	return (norm & QPMS_NORMALISATION_T_CSBIT)? -1 : 1;
}

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


typedef enum {
	QPMS_BESSEL_REGULAR = 1, // regular function j
	QPMS_BESSEL_SINGULAR = 2, // singular function y
	QPMS_HANKEL_PLUS = 3, // hankel function h1 = j + I*y
	QPMS_HANKEL_MINUS = 4, // hankel function h2 = j - I*y
	QPMS_BESSEL_UNDEF = 0
} qpms_bessel_t;

// coordinate system types
typedef struct {
	double x, y, z;
} cart3_t;

typedef struct {
	complex double x, y, z;
} ccart3_t;

typedef struct {
	double x, y;
} cart2_t;

typedef struct {
	double r, theta, phi;
} sph_t;

typedef struct { // Do I really need this???
	complex double r;
	double	theta, phi;
} csph_t;

// complex vector components in local spherical basis
typedef struct {
	complex double rc, thetac, phic; 
} csphvec_t;

typedef struct {
	double r, phi;
} pol_t;


#define lmcheck(l,m) assert((l) >= 1 && abs(m) <= (l))
#endif // QPMS_TYPES
