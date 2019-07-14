/*! \file normalisation.h
 * \brief Convention-dependent coefficients for VSWFs.
 *
 * See also @ref qpms_normalisation_t and @ref vswf_conventions.
 */
#ifndef NORMALISATION_H
#define NORMALISATION_H

#include "qpms_types.h"
#include "qpms_error.h"
#include <math.h>
#include <complex.h>
#include "indexing.h"

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



/// Returns the factors of a magnetic basis VSWF of a given convention compared to the reference convention.
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


/// Returns the factors of a magnetic basis VSWF of a given convention compared to the reference convention.
/**
 * This version takes into account the Condon-Shortley phase bit. 
 * Do not use if the C.-S. has already been taken into account e.g. in
 * a `gsl_sf_legendre_*_e()` call.
 */
static inline complex double qpms_normalisation_factor_M(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qpms_normalisation_factor_M_noCS(norm, l, m);
	return ((norm & QPMS_NORMALISATION_CSPHASE) && (m % 2)) ? -fac : fac;
}


/// Returns the factors of a electric basis VSWF of a given convention compared to the reference convention.
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


/// Returns the factors of a electric basis VSWF of a given convention compared to the reference convention.
/**
 * This version takes into account the Condon-Shortley phase bit. 
 * Do not use if the C.-S. has already been taken into account e.g. in
 * a `gsl_sf_legendre_*_e()` call.
 */
static inline complex double qpms_normalisation_factor_N(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qpms_normalisation_factor_N_noCS(norm, l, m);
	return ((norm & QPMS_NORMALISATION_CSPHASE) && (m % 2)) ? -fac : fac;
}


/// Returns the factors of a electric basis VSWF divided by the factor of a magnetic VWFS of a given convention, compared to the reference one.
static inline complex double qpms_normalisation_factor_N_M(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	return qpms_normalisation_factor_N_noCS(norm, l, m) 
		/ qpms_normalisation_factor_M_noCS(norm, l, m);
}


/// Returns the factors of a longitudinal basis VSWF of a given convention compared to the reference convention.
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

/// Returns the factors of a longitudinal basis VSWF of a given convention compared to the reference convention.
/**
 * This version takes into account the Condon-Shortley phase bit. 
 * Do not use if the C.-S. has already been taken into account e.g. in
 * a `gsl_sf_legendre_*_e()` call.
 */
static inline complex double qpms_normalisation_factor_L(qpms_normalisation_t norm, qpms_l_t l, qpms_m_t m) {
	complex double fac = qpms_normalisation_factor_L_noCS(norm, l, m);
	return ((norm & QPMS_NORMALISATION_CSPHASE) && (m % 2)) ? -fac : fac;
}

/// Returns the factors of a basis VSWF of a given convention compared to the reference convention.
static inline complex double qpms_normalisation_factor_uvswfi(const qpms_normalisation_t norm, qpms_uvswfi_t ui) {
	qpms_vswf_type_t t; qpms_m_t m; qpms_l_t l;
	qpms_uvswfi2tmn(ui, &t, &m, &l);
	switch(t) {
		case QPMS_VSWF_MAGNETIC:
			return qpms_normalisation_factor_M(norm, l, m);
		case QPMS_VSWF_ELECTRIC:
			return qpms_normalisation_factor_N(norm, l, m);
		case QPMS_VSWF_LONGITUDINAL:
			return qpms_normalisation_factor_L(norm, l, m);
		default:
			QPMS_WTF;
	}
}


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
 * and QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE set.
 *
 * For real spherical harmonics, this gives
 * \f[ 
 * 	\sqrt{2}\cos{m \phi} \quad \mbox{if } m>0, \\
 * 	\sqrt{2}\sin{m \phi} \quad \mbox{if } m<0, \\
 * 	0 \quad \mbox{if } m>0. \\
 * \f]
 */
static inline complex double qpms_spharm_azimuthal_part(qpms_normalisation_t norm, qpms_m_t m, double phi) {
	switch(norm & (QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE | QPMS_NORMALISATION_SPHARM_REAL)) {
		case 0:
			return cexp(I*m*phi);
		case QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE:
			return cexp(-I*m*phi);
		case QPMS_NORMALISATION_SPHARM_REAL:
			if (m > 0) return M_SQRT2 * cos(m*phi);
			else if (m < 0) return M_SQRT2 * sin(m*phi);
			else return 1.;
		default:
			QPMS_WTF;
	}
}

/// Returns derivative of the asimuthal part of a spherical harmonic divided by \a m.
/**
 *
 * This is used to evaluate the VSWFs together with the \a pi member array of the
 * qpms_pitau_t structure.
 *
 * Returns \f[ i e^{im\phi} \f] for standard complex spherical harmonics,
 * \f[-i e^{-i\phi \f] for complex spherical harmonics 
 * and QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE set.
 *
 * For real spherical harmonics, this gives
 * \f[ 
 * 	-\sqrt{2}\sin{m \phi} \quad \mbox{if } m>0, \\
 * 	\sqrt{2}\cos{m \phi} \quad \mbox{if } m<0, \\
 * 	-1 \quad \mbox{if } \mbox{if }m=0. \\
 * \f]
 *
 * (The value returned for \f$ m = 0 \f$ should not actually be used for
 * anything except for multiplying by zero.)
 *
 *
 */
static inline complex double qpms_spharm_azimuthal_part_derivative_div_m(qpms_normalisation_t norm, qpms_m_t m, double phi) {
	if(m==0) return 0;
	switch(norm & (QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE | QPMS_NORMALISATION_SPHARM_REAL)) {
		case 0:
			return I*cexp(I*m*phi);
		case QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE:
			return -I*cexp(-I*m*phi);
		case QPMS_NORMALISATION_SPHARM_REAL:
			if (m > 0) return -M_SQRT2 * sin(m*phi);
			else if (m < 0) return M_SQRT2 * cos(m*phi);
			else return -1;
		default:
			QPMS_WTF;
	}
}

#if 0 // legacy code moved from qpms_types.h. TODO cleanup
/// Returns the normalisation convention code without the Condon-Shortley phase.
static inline int qpms_normalisation_t_normonly(qpms_normalisation_t norm) {
	return norm & (~QPMS_NORMALISATION_T_CSBIT);
}

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
#endif

#endif //NORMALISATION_H
