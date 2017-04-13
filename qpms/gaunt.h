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
#endif //GAUNT_H
