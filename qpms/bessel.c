#include <assert.h>
#include "qpms_specfunc.h"
#include <stdlib.h>
#include <stddef.h>
#include <kahansum.h>
#include <gsl/gsl_sf_bessel.h>

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


qpms_bessel_calculator_t *qpms_bessel_calculator_init(...){
	TODO dudom;
}
