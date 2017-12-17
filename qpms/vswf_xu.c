#include "translations.h"
#include <math.h>
#include <stdlib.h> //abort();

/* 
 * References:
 * [1] Yu-Lin Xu, Journal of Computational Physics 127, 285–298 (1996)
 */

/* 
 * The following value delimits the interval for which the computation
 * of the pi and tau functions [1,(37)] actually takes place.
 * For x in [-1, -1 + PITAU_0THRESHOLD] the x->-1 limit is taken,
 * for x in [-1 + PITAU_0THRESHOLD, 1 - PITAU_0THRESHOLD] value for x is calculated,
 * for x in [1 - PITAU_THRESHOLD, 1] the x->1 limit is taken.
 *
 * low-priority TODO: take more than 0th order expansion at x->±1
 */
#define PITAU_0THRESHOLD 1e-7

int taumncos_Xu_zerolim_fill(int lmax, double *where);
double * taumncos_Xu_zerolim_get(int lmax);
int taumncos_Xu_pilim_fill(int lmax, double *where);
double * taumncos_Xu_pilim_get(int lmax);

int pimncos_Xu_zerolim_fill(int lmax, double *where);
double * pimncos_Xu_zerolim_get(int lmax);
int pimncos_Xu_pilim_fill(int lmax, double *where);
double * pimncos_Xu_pilim_get(int lmax);


// [1] (37)
static inline double pimncos_Xu(double theta) {
	abort();
}
static inline double taumncos_Xu(double theta) {
	abort();
}

// [1] (36)
complex double qpms_vswf_single_mg_Xu(int m, int n, sph_t kdlj,
		qpms_bessel_t btyp)
{
	
}
