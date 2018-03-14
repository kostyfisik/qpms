#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

int main() {
	int lMax;
	while (1 == scanf("%d", &lMax)){
		double x;
		if (1 != scanf(" %lf", &x))
			abort();
		double orig[lMax+1], gsl[lMax+2], relerr[lMax+1];
		for(int l=0; l <= lMax; ++l) 
			if (1 != scanf(" %lf", orig+l))
				abort();
#if defined JTEST || defined DJTEST
		if(gsl_sf_bessel_jl_array(lMax+1,x,gsl))
#elif defined YTEST || defined DYTEST
		if(gsl_sf_bessel_yl_array(lMax+1,x,gsl))
#else
		if(gsl_sf_bessel_jl_steed_array(lMax+1,x,gsl))
#endif
			abort();
#if defined DJTEST || defined DYTEST || defined DJTEST_STEED
#if 1
		for (int l = 0; l <= lMax; ++l)
			gsl[l] = -gsl[l+1] + (l/x) * gsl[l];
#else
	// todo varianta 10.51.2
#endif			
#endif
		for (int l = 0; l <= lMax; ++l)
			relerr[l] = fabs(gsl[l] - orig[l]) / (fabs(gsl[l]) + fabs(orig[l])) * 2;
		printf("x = %.16g\n", x);
		printf("orig: ");
		for (int l = 0; l <= lMax; ++l)
			printf("%.16g ", orig[l]);
		printf("\ngsl:  ");
		for (int l = 0; l <= lMax; ++l)
			printf("%.16g ", gsl[l]);
		printf("\nrerr:  ");
		for (int l = 0; l <= lMax; ++l)
			printf("%.16g ", relerr[l]);
		putchar('\n');
	}
}


