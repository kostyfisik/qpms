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

int lMax=10;

int main() {
	qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, QPMS_NORMALIZATION_TAYLOR);

	for(testcase_single_trans_t *tc = testcases_Taylor; tc->J != QPMS_BESSEL_UNDEF; tc++) {
		if (!tc->n || !tc->nu || tc->n > lMax || tc->nu > lMax ) continue;

		printf("m=%d, n=%d, mu=%d, nu=%d,\n", tc->m,tc->n,tc->mu,tc->nu);
		complex double A = qpms_trans_single_A_Taylor(tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
		complex double B = qpms_trans_single_B_Taylor(tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
		complex double A2 = qpms_trans_calculator_get_A(c, tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
		complex double B2 = qpms_trans_calculator_get_B(c, tc->m, tc->n, tc->mu, tc->nu, tc->kdlj, true, tc->J);
		printf("A  = %.16f+%.16fj, relerr=%.16f, J=%d\n",
				creal(A), cimag(A),  (0 == cabs(tc->result_A - A)) ? 0 :
				cabs(tc->result_A - A)/((cabs(A) < cabs(tc->result_A)) ? cabs(A) : cabs(tc->result_A)),
				tc->J);
		printf("A' = %.16f+%.16fj, relerr=%.16f, relerr2=%.3e\n",
				creal(A2), cimag(A2),  (0 == cabs(tc->result_A - A2)) ? 0 :
				cabs(tc->result_A - A2)/((cabs(A2) < cabs(tc->result_A)) ? cabs(A2) : cabs(tc->result_A)),
				(0 == cabs(A - A2)) ? 0 :
				cabs(A - A2)/((cabs(A2) < cabs(A)) ? cabs(A2) : cabs(A))
		      );
		printf("B  = %.16f+%.16fj, relerr=%.16f, J=%d\n",
				creal(B), cimag(B),  (0 == cabs(tc->result_B - B)) ? 0 :
				cabs(tc->result_B - B)/((cabs(B) < cabs(tc->result_B)) ? cabs(B) : cabs(tc->result_B)),
				tc->J);
		printf("B' = %.16f+%.16fj, relerr=%.16f, relerr2=%.3e\n",
				creal(B2), cimag(B2),  (0 == cabs(tc->result_B - B2)) ? 0 :
				cabs(tc->result_B - B2)/((cabs(B2) < cabs(tc->result_B)) ? cabs(B2) : cabs(tc->result_B)),
				(0 == cabs(B - B2)) ? 0 :
				cabs(B - B2)/((cabs(B2) < cabs(B)) ? cabs(B2) : cabs(B))
		      );
	}
	qpms_trans_calculator_free(c);
}


