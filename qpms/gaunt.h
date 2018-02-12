#ifndef GAUNT_H
#define GAUNT_H
#include <stdlib.h>

#define _GAUNT_H_MIN(x,y) (((x) > (y)) ? (y) : (x))
static inline int gaunt_q_max(int m, int n, int mu, int nu) {
	return _GAUNT_H_MIN(n, _GAUNT_H_MIN(nu, (n+nu-abs(m+mu))/2));
}
#undef _GAUNT_H_MIN

void gaunt_xu(int m, int n, int mu, int nu, int qmax, double *v_aq, int *err);
//int gaunt(int m, int n, int mu, int nu, double *v_aq);


#ifdef GAUNT_PRECOMPILED
extern const double  gaunt_table[];
extern const int gaunt_table_lMax;

// Returns given gaunt coeff. If not in range, return nan
// TODO distinguish between invalid input and limitations given by gaunt_table_lMax
double gaunt_table_retrieve_single(int m, int n, int mu, int nu, int q);

// returns pointer to the memory where gaunt(m, n, mu, nu, q=0) is placed
// returns NULL if invalid input or exceeding gaunt_table_lMax
double const * gaunt_table_retrieve_allq(int m, int n, int mu, int nu);

int gaunt_table_or_xu_fill(double *target, int m, int n, int mu, int nu);
#endif //GAUNT_PRECOMPILED

#endif //GAUNT_H
