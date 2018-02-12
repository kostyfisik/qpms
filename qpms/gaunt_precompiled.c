#ifndef GAUNT_PRECOMPILED
#define GAUNT_PRECOMPILED
#endif
#include "gaunt.h"
#include <assert.h>
#include <stddef.h>
#include <math.h>

const int gaunt_table_lMax = 18;

const double gaunt_table[] = {
#include "data/gaunt18baredouble"
};

const size_t gaunt_table_qmaxcumsum[] = {
#include "data/gaunt18qmaxcumsum"
};

/* Pořadí indexů:
 * for (n = 0; n <= lMax; n++)
 *   for (m = -n; m <= n; m++)
 *     for (nu = 0; nu <= lMax; nu++)
 *       for (mu = -nu; mu <= nu; mu++)
 *         for (q = 0; q < min(n, nu, (n + nu - |m + mu|)/2); q++)
 */


double const * gaunt_table_retrieve_allq(int m, int n, int mu, int nu) {
	assert(abs(m) <= n);
	assert(abs(mu) <= nu);
	if (n > gaunt_table_lMax || nu > gaunt_table_lMax)
		return NULL;
	size_t x = n * (size_t) (n + 1) + (ptrdiff_t) m;
	size_t xu = nu * (size_t) (nu + 1) + (ptrdiff_t) mu;
	size_t xcount = gaunt_table_lMax * (size_t) (gaunt_table_lMax + 2) + 1;
	size_t idx = x * xcount + xu;
	return gaunt_table + gaunt_table_qmaxcumsum[idx];
}


double gaunt_table_retrieve_single(int m, int n, int mu, int nu, int q) {
	double const * gauntptr = gaunt_table_retrieve_allq(m, n, mu, nu);
	if (NULL == gauntptr)
		return NAN;
	else return gauntptr[q];
}

int gaunt_table_or_xu_fill(double *target, int m, int n, int mu, int nu) {
	double const * gauntptr = gaunt_table_retrieve_allq(m, n, mu, nu);
	int qmax = gaunt_q_max(m,n,mu,nu);
	if (gauntptr) { // table hit
		for(int q = 0; q <= qmax; ++qmax) target[q] = gauntptr[q];
		return 0;
	} else {
		int err;
		gaunt_xu(m, n, mu, nu, qmax, target, &err);
		return err;
	}
}


