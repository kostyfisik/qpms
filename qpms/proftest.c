
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
	int repete = 500;
	int lMax = 3;
	qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, QPMS_NORMALIZATION_TAYLOR);
	for( int rr = 0; rr < repete; rr++)
	for(testcase_single_trans_t *tc = testcases_Taylor; tc->J != QPMS_BESSEL_UNDEF; tc++) {
		//if (tc->n > 40 || tc->nu > 40 ) continue;
		complex double A_array[c->nelem * c->nelem];
		complex double B_array[c->nelem * c->nelem];
		qpms_trans_calculator_get_AB_arrays(c, A_array, B_array, c->nelem, 1, tc->kdlj, true, tc->J);
#if 0
		complex double A = qpms_trans_single_A_Taylor(tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
		complex double B = qpms_trans_single_B_Taylor(tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
		printf("A = %.16f+%.16fj, relerr=%.16f, J=%d\n",
			creal(A), cimag(A),  (0 == cabs(tc->result_A - A)) ? 0 :
			cabs(tc->result_A - A)/((cabs(A) < cabs(tc->result_A)) ? cabs(A) : cabs(tc->result_A)),
			tc->J);
		printf("B = %.16f+%.16fj, relerr=%.16f, J=%d\n",
			creal(B), cimag(B), (0 == cabs(tc->result_B - B)) ? 0 :
			cabs(tc->result_B - B)/((cabs(B) < cabs(tc->result_B)) ? cabs(B) : cabs(tc->result_B)),
			tc->J);
#endif
	}
}


