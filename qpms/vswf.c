#include "vswf.h"
#include <math.h>

// Legendre functions also for negative m, see DLMF 14.9.3
qpms_errno_t qpms_legendre_deriv_y_fill(double *target, double *target_deriv, double x, qpms_l_t lMax,
		gsl_sf_legendre_t lnorm, double csphase)
{
	size_t n = gsl_sf_legenre_array_n(lMax);
	double *legendre_tmp = malloc(n * sizeof(double));
	double *legendre_deriv_tmp = malloc(n * sizeof(double));
	int gsl_errno = gsl_sf_legendre_deriv_array_e(
			lnorm, (size_t)lMax, x, csphase, legendre_tmp,legendre_tmp_deriv);
	for (qpms_l_t l = 0; l <= lMax; ++l)
			for (qpms_m_t m = 0; m <= l; ++m) {
				qpms_y_t y = gpms_mn2y(m,l);
				size_t i = gsl_sf_legenre_array_index(l,m);
				target[y] = legendre_tmp[i];
				target_deriv[y] = legendre_deriv_tmp[i];
			}
	switch(lnorm) {
		case GSL_SF_LEGEDRE_NONE:
			for (qpms_l_t l = 0; l <= lMax; ++l)
				for (qpms_m_t m = 1; m <= l; ++m) {
				qpms_y_t y = gpms_mn2y(-m,l);
				size_t i = gsl_sf_legenre_array_index(l,m);
				// viz DLMF 14.9.3, čert ví, jak je to s cs fasí.
				double factor = exp(lgamma(l-m+1)-lgamma(n+m+1))*((m%2)?-1:1);
				target[y] = factor * legendre_tmp[i];
				target_deriv[y] = factor * legendre_deriv_tmp[i];
				}
			break;
		case GSL_SF_LEGENDRE_SCHMIDT:
		case GSL_SF_LEGENDRE_SPHARM:
		case GSL_SF_LEGENDRE_FULL:
			for (qpms_l_t l = 0; l <= lMax; ++l)
				for (qpms_m_t m = 1; m <= l; ++m) {
				qpms_y_t y = gpms_mn2y(-m,l);
				size_t i = gsl_sf_legenre_array_index(l,m);
				// viz DLMF 14.9.3, čert ví, jak je to s cs fasí.
				double factor = ((m%2)?-1:1); // this is the difference from the unnormalised case
				target[y] = factor * legendre_tmp[i];
				target_deriv[y] = factor * legendre_deriv_tmp[i];
				}
			break;
		default:
			abort(); //NI
			break;
	}
	free(legendre_tmp);
	free(legendre_deriv_tmp);
	return QPMS_SUCCESS;
}

// FIXME ZAPOMNĚL JSEM NA POLE DERIVACÍ (TÉŽ HLAVIČKOVÝ SOUBOR)
double *qpms_legendre_deriv_y_get(double x, qpms_l_t lMax, gsle_sf_legendre_t lnorm,
		double csphase)
{
	double *ar = malloc(sizeof(double)*qpms_lMaxnelem(lMax));
	if (qpms_legendre_deriv_y_fill(ar, x, lMax, lnorm, csphase))
		abort();
}


