#include "translations.h"
#include "qpms_types.h"
#include "indexing.h"
#include <math.h>
#include <stdlib.h> //abort();

/* 
 * References:
 * [1] Yu-Lin Xu, Journal of Computational Physics 127, 285–298 (1996)
 */

/* 
 * The following value delimits the interval for which the computation
 * of the pi and tau functions [1,(37)] actually takes place.
 * see also DLMF §14.8
 * For x in [-1, -1 + PITAU_0THRESHOLD] the x->-1 limit is taken,
 * for x in [-1 + PITAU_0THRESHOLD, 1 - PITAU_0THRESHOLD] value for x is calculated,
 * for x in [1 - PITAU_THRESHOLD, 1] the x->1 limit is taken.
 *
 * low-priority TODO: take more than 0th order expansion at x->±1
 */
#define PITAU_0THRESHOLD 1e-7

// Legendre functions for negative m; see DLMF 14.9.3
// 


// tau, pi at singularities and zero; see also DLMF 14.8
int taumncos_Xu_zerolim_fill(qpms_l_t lmax, double *where)
{	

	return 0;
}

double * taumncos_Xu_zerolim_get(qpms_l_t lmax) {
	double *ar = malloc(qpms_lMax2nelem(lmax) * sizeof(double));
	taumncos_Xu_zerolim_fill(lmax,ar);
	return ar;
}

int taumncos_Xu_pilim_fill(qpms_l_t lmax, double *where);

double * taumncos_Xu_pilim_get(qpms_l_t lmax) {
	double *ar = malloc(qpms_lMax2nelem(lmax) * sizeof(double));
	taumncos_Xu_pilim_fill(lmax,ar);
	return ar;
}

int taumncos_Xu_pihalflim_fill(qpms_l_t lmax, double *where);

double * taumncos_Xu_pihalflim_get(qpms_l_t lmax) {
	double *ar = malloc(qpms_lMax2nelem(lmax) * sizeof(double));
	taumncos_Xu_pihalflim_fill(lmax,ar);
	return ar;
}


int pimncos_Xu_zerolim_fill(qpms_l_t lmax, double *where);
double * pimncos_Xu_zerolim_get(qpms_l_t lmax){
	double *ar = malloc(qpms_lMax2nelem(lmax) * sizeof(double));
	pimncos_Xu_pihalflim_fill(lmax,ar);
	return ar;
}
int pimncos_Xu_pilim_fill(qpms_l_t lmax, double *where);
double * pimncos_Xu_pilim_get(qpms_l_t lmax);


// [1] (37)
static inline double pimncos_Xu(double theta) {
	abort();//NI
}
static inline double taumncos_Xu(double theta) {
	abort();//NI
}

// [1] (36)
complex double qpms_vswf_single_mg_Xu(qpms_m_t m, qpms_l_t n, sph_t kdlj,
		qpms_bessel_t btyp)
{
	abort();//NI	
}
