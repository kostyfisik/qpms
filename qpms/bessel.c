#include <assert.h>
#include "qpms_specfunc.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "kahansum.h"
#include <gsl/gsl_sf_bessel.h>
#include <complex.h>

#ifndef M_LN2
#define M_LN2   0.69314718055994530942  /* log_e 2 */
#endif

static inline complex double ipow(int x) {
	return cpow(I,x);
}

// There is a big issue with gsl's precision of spherical bessel function; these have to be implemented differently
qpms_errno_t qpms_sph_bessel_fill(qpms_bessel_t typ, qpms_l_t lmax, double x, complex double *result_array) {
	int retval;
	double tmparr[lmax+1];
	switch(typ) {
		case QPMS_BESSEL_REGULAR:
			retval = gsl_sf_bessel_jl_steed_array(lmax, x, tmparr);
			for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
			return retval;
			break;
		case QPMS_BESSEL_SINGULAR: //FIXME: is this precise enough? Would it be better to do it one-by-one?
			retval = gsl_sf_bessel_yl_array(lmax,x,tmparr);
			for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
			return retval;
			break;
		case QPMS_HANKEL_PLUS:
		case QPMS_HANKEL_MINUS:
			retval = gsl_sf_bessel_jl_steed_array(lmax, x, tmparr);
			for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
			if(retval) return retval;
			retval = gsl_sf_bessel_yl_array(lmax, x, tmparr);
			if (typ==QPMS_HANKEL_PLUS)
				for (int l = 0; l <= lmax; ++l) result_array[l] += I * tmparr[l];
			else 
				for (int l = 0; l <= lmax; ++l) result_array[l] +=-I * tmparr[l];
			return retval;
			break;
		default:
			abort();
			//return GSL_EDOM;
	}
	assert(0);
}

static inline ptrdiff_t akn_index(qpms_l_t n, qpms_l_t k) {
	assert(k <= n);
	return ((ptrdiff_t) n + 1) * n / 2 + k;
}
static inline ptrdiff_t bkn_index(qpms_l_t n, qpms_l_t k) {
	assert(k <= n+1);
	return ((ptrdiff_t) n + 2) * (n + 1) / 2 - 1 + k;
}

static inline qpms_errno_t qpms_sbessel_calculator_ensure_lMax(qpms_sbessel_calculator_t *c, qpms_l_t lMax) {
	if (lMax <= c->lMax)
		return QPMS_SUCCESS;
	else {
		if ( NULL == (c->akn = realloc(c->akn, sizeof(double) * akn_index(lMax + 2, 0))))
			abort();
		//if ( NULL == (c->bkn = realloc(c->bkn, sizeof(complex double) * bkn_index(lMax + 1, 0))))
		//	abort();
		for(qpms_l_t n = c->lMax+1; n <= lMax + 1; ++n)
			for(qpms_l_t k = 0; k <= n; ++k)
				c->akn[akn_index(n,k)] = exp(lgamma(n + k + 1) - k*M_LN2 - lgamma(k + 1) - lgamma(n - k + 1));
		// ... TODO derivace
		c->lMax = lMax;
		return QPMS_SUCCESS;
	}
}

complex double qpms_sbessel_calc_h1(qpms_sbessel_calculator_t *c, qpms_l_t n, double x) {
	if(QPMS_SUCCESS != qpms_sbessel_calculator_ensure_lMax(c, n))
		abort();
	complex double z = I/x; // FIXME this should be imaginary double, but gcc is broken?
	complex double result = 0;
	for (qpms_l_t k = n; k >= 0; --k) 
		// can we use fma for complex?
		//result = fma(result, z, c->akn(n, k));
		result = result * z + c->akn[akn_index(n,k)];
	result *= z * ipow(-n-2) * cexp(I * x);
	return result;
}

qpms_errno_t qpms_sbessel_calc_h1_fill(qpms_sbessel_calculator_t * const c,
	       const qpms_l_t lMax, const double x, complex double * const target) {
	if(QPMS_SUCCESS != qpms_sbessel_calculator_ensure_lMax(c, lMax))
		abort();
	memset(target, 0, sizeof(complex double) * lMax);
	complex double kahancomp[lMax];
	memset(kahancomp, 0, sizeof(complex double) * lMax);
	for(qpms_l_t k = 0; k <= lMax; ++k){
		double xp = pow(x, -k-1);
		for(qpms_l_t l = k; l <= lMax; ++l)
			ckahanadd(target + l, kahancomp + l, c->akn[akn_index(l,k)] * xp * ipow(k-l-1));
	}
	complex double eix = cexp(I * x);
	for (qpms_l_t l = 0; l <= lMax; ++l)
		target[l] *= eix; 
	return QPMS_SUCCESS;
}

qpms_sbessel_calculator_t *qpms_sbessel_calculator_init() {
	qpms_sbessel_calculator_t *c = malloc(sizeof(qpms_sbessel_calculator_t));
	c->akn = NULL;
	//c->bkn = NULL;
	c->lMax = -1;
	return c;
}

void qpms_sbessel_calculator_pfree(qpms_sbessel_calculator_t *c) {
	if(c->akn) free(c->akn);
	//if(c->bkn) free(c->bkn);
	free(c);
}
