/*! \file indexing.h
 * \brief Various index conversion functions.
 */
#ifndef QPMS_INDEXING_H
#define QPMS_INDEXING_H

#include "qpms_types.h"
#include <math.h>

static inline qpms_y_t qpms_mn2y(qpms_m_t m, qpms_l_t n) {
        return n * (n + 1) + m - 1;
}

static inline qpms_lm_t qpms_y2n(qpms_y_t y) {
        //return (sqrt(5+y)-2)/2; // the cast will truncate the fractional part, which is what we want
        return sqrt(y+1);
}

static inline qpms_m_t qpms_yn2m(qpms_y_t y, qpms_l_t n) {
        return y-qpms_mn2y(0,n);
}

static inline void qpms_y2mn_p(qpms_y_t y, qpms_m_t *m, qpms_l_t *n){
        *m=qpms_yn2m(y,*n=qpms_y2n(y));
}

static inline qpms_y_t qpms_lMax2nelem(qpms_l_t lmax){
	return lmax * ((qpms_y_t)lmax + 2);
}

// Scalar versions: they have a place for the 0, 0 term in the beginning

static inline qpms_y_t qpms_mn2y_sc(qpms_m_t m, qpms_l_t n) {
        return n * (n + 1) + m;
}

static inline qpms_lm_t qpms_y2n_sc(qpms_y_t y) {
        //return (sqrt(5+y)-2)/2; // the cast will truncate the fractional part, which is what we want
        return sqrt(y);
}

static inline qpms_m_t qpms_yn2m_sc(qpms_y_t y, qpms_l_t n) {
        return y-qpms_mn2y_sc(0,n);
}

static inline void qpms_y2mn_sc_p(qpms_y_t y, qpms_m_t *m, qpms_l_t *n){
        *m=qpms_yn2m_sc(y,*n=qpms_y2n_sc(y));
}

static inline qpms_y_t qpms_lMax2nelem_sc(qpms_l_t lmax){
	return lmax * ((qpms_y_t)lmax + 2) + 1;
}

// TODO maybe enable crashing / validity control by macro definitions...

/// Conversion from VSWF type, order and degree to universal index.
static inline qpms_uvswfi_t qpms_tmn2uvswfi(
		qpms_vswf_type_t t, qpms_m_t m, qpms_l_t n) {
	return t + 4 * qpms_mn2y_sc(m, n);
}

/// Conversion from universal VSWF index u to type, order and degree.
/** Returns a non-zero value if the u value is invalid. */
static inline qpms_errno_t qpms_uvswfi2tmn(qpms_uvswfi_t u,
		qpms_vswf_type_t *t, qpms_m_t *m, qpms_l_t *n) {
	*t = u & 3;
	qpms_y_sc_t y_sc = u / 4;
	qpms_y2mn_sc_p(y_sc, m, n);
	// Test validity
	if (*t == 3) return QPMS_ERROR; // VSWF type code invalid, TODO WARN
	if (*t && !y_sc) return QPMS_ERROR; // l == 0 for transversal wave, TODO WARN
	return QPMS_SUCCESS;
}

/// Extract degree \a m from an universal VSWF index \a u.
static inline qpms_m_t qpms_uvswfi2m(qpms_uvswfi_t u) {
	qpms_vswf_type_t t; qpms_m_t m; qpms_l_t n;
	qpms_uvswfi2tmn(u, &t,&m,&n);
	return m;
}


#endif //QPMS_INDEXING_H
