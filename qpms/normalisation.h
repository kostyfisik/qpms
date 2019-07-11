/*! \file normalisation.h
 * \brief Convention-dependent coefficients for VSWFs.
 *
 * See also @ref qpms_normalisation_t and @ref vswf_conventions.
 */
#ifndef NORMALISATION_H
#define NORMALISATION_H

#include "qpms_types.h"
#include <math.h>
#include <complex.h>


/// Returns the (real positive) common norm factor of a given normalisation compared to the reference convention.
/** Does NOT perform the inversion if QPMS_NORMALISATION_INVERSE is set. */
static inline double qpms_normalisation_normfactor(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	switch (norm & QPMS_NORMALISATION_NORM_BITS) {
		case QPMS_NORMALISATION_NORM_POWER:
			return 1;
		case QPMS_NORMALISATION_NORM_SPHARM:
			return sqrt(l*(l+1));
		case QPMS_NORMALISATION_NORM_NONE: // TODO more precision
			return sqrt(l*(l+1) * 4*M_PI / (2*l+1)) *
			       	exp(0.5*(lgamma(l+m+1) - lgamma(l-m+1)));
		default:
			QPMS_WTF;
	}
}


/// Returns the factors of a magnetic VSWF of a given convention compared to the reference convention.
/**
 * This version ignores the Condon-Shortley phase bit (perhaps because the Condon-Shortley
 * phase is already taken into account in a `gsl_sf_legendre_*_e()` call.)
 */
static inline complex double qpms_normalisation_factor_M_noCS(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qpms_normalisation_normfactor(norm, l, m);
	if (norm & QPMS_NORMALISATION_M_MINUS) fac *= -1;
	if (norm & QPMS_NORMALISATION_M_I) fac *= I;
	if (norm & QPMS_NORMALISATION_INVERSE) fac = 1/fac;
	return fac;
}


/// Returns the factors of a magnetic VSWF of a given convention compared to the reference convention.
/**
 * This version takes into account the Condon-Shortley phase bit. 
 * Do not use if the C.-S. has already been taken into account e.g. in
 * a `gsl_sf_legendre_*_e()` call.
 */
static inline complex double qpms_normalisation_factor_M(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qn;
	return ((norm & QPMS_NORMALISATION_CSPHASE) && (m % 2)) ? -fac : fac;
}


/// Returns the factors of a electric VSWF of a given convention compared to the reference convention.
/**
 * This version ignores the Condon-Shortley phase bit (perhaps because the Condon-Shortley
 * phase is already taken into account in a `gsl_sf_legendre_*_e()` call.)
 */
static inline complex double qpms_normalisation_factor_N_noCS(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qpms_normalisation_normfactor(norm, l, m);
	if (norm & QPMS_NORMALISATION_N_MINUS) fac *= -1;
	if (norm & QPMS_NORMALISATION_N_I) fac *= I;
	if (norm & QPMS_NORMALISATION_INVERSE) fac = 1/fac;
	return fac;
}


/// Returns the factors of a electric VSWF of a given convention compared to the reference convention.
/**
 * This version takes into account the Condon-Shortley phase bit. 
 * Do not use if the C.-S. has already been taken into account e.g. in
 * a `gsl_sf_legendre_*_e()` call.
 */
static inline complex double qpms_normalisation_factor_N(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qn;
	return ((norm & QPMS_NORMALISATION_CSPHASE) && (m % 2)) ? -fac : fac;
}


#if 0
/// Returns the factors of a longitudinal VSWF of a given convention compared to the reference convention.
/**
 * This version ignores the Condon-Shortley phase bit (perhaps because the Condon-Shortley
 * phase is already taken into account in a `gsl_sf_legendre_*_e()` call.)
 */
static inline complex double qpms_normalisation_factor_L_noCS(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qpms_normalisation_normfactor(norm, l, m);
	if (norm & QPMS_NORMALISATION_L_MINUS) fac *= -1;
	if (norm & QPMS_NORMALISATION_L_I) fac *= I;
	if (norm & QPMS_NORMALISATION_INVERSE) fac = 1/fac;
	return fac;
}

/// Returns the factors of a longitudinal VSWF of a given convention compared to the reference convention.
/**
 * This version takes into account the Condon-Shortley phase bit. 
 * Do not use if the C.-S. has already been taken into account e.g. in
 * a `gsl_sf_legendre_*_e()` call.
 */
static inline complex double qpms_normalisation_factor_L(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qn;
	return ((norm & QPMS_NORMALISATION_CSPHASE) && (m % 2)) ? -fac : fac;
}
#endif

/// Returns normalisation flags corresponding to the dual spherical harmonics / waves.
/**
 * This reverses the normalisation factors returned by qpms_normalisation_factor_*
 * and conjugates the asimuthal part for complex spherical harmonics, 
 * \f$ e^{\pm im\phi} \leftrightarrow e^{\mp im\phi} \f$.
 */
static inline qpms_normalisation_t qpms_normalisation_dual(qpms_normalisation_t norm) {
	norm ^= QPMS_NORMALISATION_INVERSE;
	if (!(norm & QPMS_NORMALISATION_SPHARM_REAL))
		norm ^= QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE;
	return norm;
}

/// Returns the asimuthal part of a spherical harmonic.
/** Returns \f[ e^{im\phi} \f] for standard complex spherical harmonics,
 * \f[ e^{-im\phi \f] for complex spherical harmonics 
 * and QPMS_NORMALISATION_REVERSE_ASIMUTHAL_PHASE set.
 *
 * For real spherical harmonics, this gives
 * \f[ 
 * 	\sqrt{2}\cos{m \phi} \quad \mbox{if } m>0, \\
 * 	\sqrt{2}\sin{m \phi} \quad \mbox{if } m<0, \\
 * 	1 \quad \mbox{if } m>0. \\
 * \f]
 */
static inline complex double qpms_spharm_azimuthal_part(qpms_normalisation_t norm, qpms_m_t m, double phi) {
	switch(norm & (QPMS_NORMALISATION_REVERSE_ASIMUTHAL_PHASE | QPMS_NORMALISATION_SPHARM_REAL)) {
		case 0:
			return cexp(I*m*phi);
		case QPMS_NORMALISATION_REVERSE_ASIMUTHAL_PHASE:
			return cexp(-I*m*phi);
		case QPMS_NORMALISATION_SPHARM_REAL:
			if (m > 0) return M_SQRT2 * cos(m*phi);
			else if (m < 0) return M_SQRT2 * sin(m*phi);
			else return 1.;
		default:
			QPMS_WTF;
	}
}

#endif //NORMALISATION_H
