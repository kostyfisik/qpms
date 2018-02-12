#define GAUNT_PRECOMPILED
#include "gaunt.h"
#include <stdio.h>
#include <math.h>

const double rerrth = 1e-12;

int main()
{
	int lMax = gaunt_table_lMax;
	for (int n = 1; n <= lMax; n++)
		for (int m = -n; m <= n; m++)
			for (int nu = 1; nu <= lMax; nu++)
				for (int mu = -nu; mu <= nu; mu++) {
					int err;
					int qmax = gaunt_q_max(m,n,mu,nu);
					double gc_xu[qmax+1];
					gaunt_xu(m,n,mu,nu,qmax,gc_xu,&err);
					double const * gc_table = gaunt_table_retrieve_allq(m, n, mu, nu);
					for (int q = 0; q < qmax; ++q) {
						double rerr = (gc_xu[q] || gc_table[q]) ? 2 * (gc_xu[q] - gc_table[q]) / fabs(gc_xu[q] + gc_table[q]) : 0;
						printf("%.5e %s %d %d %d %d %d %.16e %.16e\n",
								rerr, (fabs(rerr) > rerrth) ? "!" : " ", m, n, mu, nu, q, gc_xu[q], gc_table[q]);
					}
				}
	return 0;
}




