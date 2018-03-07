#include "translations.h"
#include "gaunt.h"
#include <stdio.h>
//#include <math.h>
#include <complex.h>


void qpms_trans_calculator_multipliers_B_general(
                qpms_normalisation_t norm,
                complex double *dest, int m, int n, int mu, int nu, int Qmax) ;

int lMax=13;

#define MIN(x,y) ((x)<(y)?(x):(y))

// Python test: Qmax(M, n, mu, nu) = floor(min(n,nu,(n+nu+1-abs(M+mu))/2))
// q in IntegerRange(1, Qmax(-m,n,mu,nu))

int main() {
	qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, QPMS_NORMALISATION_XU);
	complex double dest[lMax + 2];

	for(int n = 1; n <= lMax; ++n)
		for(int nu = 1; nu <= lMax; ++nu)
			for(int m = -n; m <= n; ++m)
				for(int mu = -nu; mu <= nu; ++mu){
					int Qmax = gaunt_q_max(-m, n+1, mu, nu);
					int Qmax_alt = MIN(n,MIN(nu,(n+nu+1-abs(mu-m))));
					qpms_trans_calculator_multipliers_B_general(QPMS_NORMALISATION_XU_CS,
							dest, m, n, mu, nu, Qmax);
					for(int q = 0; q <= Qmax; ++q) {
						// int p = n + nu - 2*q;
						int tubig = cabs(dest[q]) > 1e-8;
						printf("%.16g + %.16g*I, // %d, %d, %d, %d, %d,%s\n",
								creal(dest[q]), cimag(dest[q]), m, n, mu, nu, q,
								q > Qmax_alt ? (tubig?" //tubig":" //tu") : "");
					}
					fflush(stdout);
				}
					


	qpms_trans_calculator_free(c);
}


