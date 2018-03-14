#include <gsl/gsl_sf_bessel.h>
#include <stdio.h>
#include <math.h>


int main() {
	int lMax;
	while (1 == scanf("%d", &lMax)) {
		double x;
		if (1 != scanf("%lf", &x))
			abort();
		double gsl[lMax+2], relerr[lMax+1], orig[lMax+1];
		for (int l = 0; l <= lMax; l++)
			if (1 != scanf("%lf", orig+l))
				abort();
#if defined JTEST || defined DJTEST
		if (gsl_sf_bessel_jl_array(lMax+1, x, gsl))
#elif defined YTEST || defined DYTEST
		if (gsl_sf_bessel_yl_array(lMax+1, x, gsl))
#else
		if (gsl_sf_bessel_jl_steed_array(lMax+1, x, gsl))
#endif
			abort();
#if defined DJTEST || defined DYTEST || defined DJTEST_STEED
		for (int l = 0; l <= lMax; l++)
			gsl[l] = -gsl[l+1] + (l/x) * gsl[l];
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
	return 0;
}

