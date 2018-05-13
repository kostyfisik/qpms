#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>


int main(int argc, char **argv) {
	int lmax;
	int scanned;
	double x;
	while (EOF != (scanned = scanf("%d %lf", &lmax, &x))) {
		if (scanned != 2) continue;
		size_t as = gsl_sf_legendre_array_n(lmax);
		double *r_none = calloc(as, sizeof(double));
		double *d_none = calloc(as, sizeof(double));
		double *r_schmidt = calloc(as, sizeof(double));
		double *d_schmidt = calloc(as, sizeof(double));
		double *r_spharm = calloc(as, sizeof(double));
		double *d_spharm = calloc(as, sizeof(double));
		double *r_full = calloc(as, sizeof(double));
		double *d_full = calloc(as, sizeof(double));
		if( 0
				|| gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_NONE,lmax,x,-1,r_none,d_none) 
				|| gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SCHMIDT,lmax,x,-1,r_schmidt,d_schmidt) 
				|| gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,lmax,x,-1,r_spharm,d_spharm) 
				|| gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_FULL,lmax,x,-1,r_full,d_full) 
		  ) fprintf(stderr, "Something gone wrong for lmax = %d, x = %.15e!\n", lmax, x);
		for (int l = 0; l <= lmax; ++l)
			for (int m = 0; m <= l; ++m){
				size_t i = gsl_sf_legendre_array_index(l,m);
				printf("P(%d,%d)\t%.16e\t%.16e\t%.16e\t%.16e\n", l, m,
						r_none[i], r_schmidt[i], r_spharm[i], r_full[i]);
				printf("dP(%d,%d)\t%.16e\t%.16e\t%.16e\t%.16e\n", l, m,
						d_none[i], d_schmidt[i], d_spharm[i], d_full[i]);

			}




		free(r_none);
		free(d_none);
		free(r_schmidt);
		free(d_schmidt);
		free(r_spharm);
		free(d_spharm);
		free(r_full);
		free(d_full);
	}
	return 0;
}



