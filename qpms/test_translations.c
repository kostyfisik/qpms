#include "translations.h"
#include <stdio.h>
//#include <math.h>
#include <complex.h>

typedef struct {
	int m, n, mu, nu;
	sph_t kdlj;
	qpms_bessel_t J;
	complex double result_A, result_B;
} testcase_single_trans_t;

testcase_single_trans_t testcases_Taylor[] = {
#include "testcases_taylor"
};

int main() {
	for(testcase_single_trans_t *tc = testcases_Taylor; tc->J != QPMS_BESSEL_UNDEF; tc++) {
		if (tc->n > 12 || tc->nu > 12 || !tc->n || !tc->nu ) continue;

		printf("m=%d, n=%d, mu=%d, nu=%d,\n", tc->m,tc->n,tc->mu,tc->nu);
		complex double A = qpms_trans_single_A_Taylor(tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
		complex double B = qpms_trans_single_B_Taylor(tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
		printf("A = %.16f+%.16fj, relerr=%.16f, J=%d\n",
			creal(A), cimag(A),
			cabs(tc->result_A - A)/((cabs(A) < cabs(tc->result_A)) ? cabs(A) : cabs(tc->result_A)),
			tc->J);
		printf("B = %.16f+%.16fj, relerr=%.16f, J=%d\n",
			creal(B), cimag(B),
			cabs(tc->result_B - B)/((cabs(B) < cabs(tc->result_B)) ? cabs(B) : cabs(tc->result_B)),
			tc->J);
	}
}


