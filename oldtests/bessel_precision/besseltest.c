#include <gsl/gsl_sf_bessel.h>
#include <stdio.h>
#include <math.h>

#if defined JTEST_QPMS || defined DJTEST_QPMS || defined YTEST_QPMS || defined DYTEST_QPMS
#include "../../qpms/qpms_specfunc.h"
#endif


int main() {
	int lMax;
#if defined JTEST_QPMS || defined DJTEST_QPMS || defined YTEST_QPMS || defined DYTEST_QPMS
	qpms_sbessel_calculator_t *c = qpms_sbessel_calculator_init();
#endif
	while (1 == scanf("%d", &lMax)) {
		double x;
		if (1 != scanf("%lf", &x))
			abort();
		double gsl[lMax+2], relerr[lMax+1], orig[lMax+1];

		for (int l = 0; l <= lMax; l++)
			if (1 != scanf("%lf", orig+l))
				abort();
#if defined JTEST_QPMS || defined DJTEST_QPMS || defined YTEST_QPMS || defined DYTEST_QPMS
		complex double hankel[lMax+2];
		qpms_sbessel_calc_h1_fill(c, lMax+1, x, hankel);
		for(int l = 0; l <= lMax+1; l++)
	#if defined JTEST_QPMS
			gsl[l] = creal(hankel[l]);
	#endif
#elif defined JTEST || defined DJTEST
		if (gsl_sf_bessel_jl_array(lMax+1, x, gsl)) abort();
#elif defined YTEST || defined DYTEST
		if (gsl_sf_bessel_yl_array(lMax+1, x, gsl)) abort();
#elif defined JTEST_STEED || DJTEST_STEED
		if (gsl_sf_bessel_jl_steed_array(lMax+1, x, gsl)) abort();
#endif
#if defined DJTEST || defined DYTEST || defined DJTEST_STEED
		for (int l = 0; l <= lMax; l++)
			gsl[l] = -gsl[l+1] + ((double)l/x) * gsl[l];
#endif
		printf("x = %.16g, lMax = %d:\nsage: ", x, lMax);
		for (int l = 0; l <= lMax; l++)
			printf("%.16g ", orig[l]);
		printf("\ngsl:  ");
		for (int l = 0; l <= lMax; l++)
			printf("%.16g ", gsl[l]);
		printf("\nrell: ");
		for (int l = 0; l <= lMax; l++)
			printf("%.16g ", 2 * fabs(gsl[l] - orig[l]) / (fabs(gsl[l]) + fabs(gsl[l])));
		putchar('\n');
	}
#if defined JTEST_QPMS || defined DJTEST_QPMS || defined YTEST_QPMS || defined DYTEST_QPMS
	qpms_sbessel_calculator_pfree(c);
#endif
	return 0;
}

