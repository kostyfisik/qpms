#ifndef QPMS_TYPES_H
#define QPMS_TYPES_H
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#define M_PI_2 (1.570796326794896619231321691639751442098584699687552910487)

// integer index types
typedef int qpms_lm_t;
typedef int qpms_l_t; // can't be unsigned because of the behaviour under - operator
typedef qpms_lm_t qpms_m_t;
typedef size_t qpms_y_t;

typedef enum {
	QPMS_SUCCESS = 0,
	QPMS_ERROR = 1
} qpms_errno_t;

// Normalisations
typedef enum {
	// As in TODO
	QPMS_NORMALISATION_XU = 3, // NI!
	// As in http://www.eit.lth.se/fileadmin/eit/courses/eit080f/Literature/book.pdf, power-normalised
	QPMS_NORMALISATION_KRISTENSSON = 2, // NI!
	QPMS_NORMALISATION_POWER = QPMS_NORMALISATION_KRISTENSSON, // NI!
	QPMS_NORMALISATION_TAYLOR = 1,
	QPMS_NORMALISATION_UNDEF = 0
} qpms_normalisation_t;

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
	double x, y;
} cart2_t;

typedef struct {
	double r, theta, phi;
} sph_t;

// complex vector components in local spherical basis
typedef struct {
	complex double rc, thetac, phic; 
} csphvec_t;

typedef struct {
	double r, phi;
} pol_t;

#endif // QPMS_TYPES
