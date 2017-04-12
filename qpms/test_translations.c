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
	complex double A = qpms_trans_single_A_Taylor(tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
	printf("m=%d, n=%d, mu=%d, nu=%d, relerr=%.16f\n", tc->m,tc->n,tc->mu,tc->nu,
			cabs(tc->result_A - A)/((cabs(A) < cabs(tc->result_A)) ? cabs(A) : cabs(tc->result_A)));
	}
}


