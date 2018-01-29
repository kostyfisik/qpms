#include "vswf.h"
#include "indexing.h"
#include "translations.h" // TODO move qpms_sph_bessel_fill elsewhere
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <stdlib.h>
#include <string.h>
#include "assert_cython_workaround.h"

// Legendre functions also for negative m, see DLMF 14.9.3
qpms_errno_t qpms_legendre_deriv_y_fill(double *target, double *target_deriv, double x, qpms_l_t lMax,
		gsl_sf_legendre_t lnorm, double csphase)
{
	size_t n = gsl_sf_legendre_array_n(lMax);
	double *legendre_tmp = malloc(n * sizeof(double));
	double *legendre_deriv_tmp = malloc(n * sizeof(double));
	int gsl_errno = gsl_sf_legendre_deriv_array_e(
			lnorm, (size_t)lMax, x, csphase, legendre_tmp,legendre_deriv_tmp);
	for (qpms_l_t l = 1; l <= lMax; ++l)
			for (qpms_m_t m = 0; m <= l; ++m) {
				qpms_y_t y = qpms_mn2y(m,l);
				size_t i = gsl_sf_legendre_array_index(l,m);
				target[y] = legendre_tmp[i];
				target_deriv[y] = legendre_deriv_tmp[i];
			}
	switch(lnorm) {
		case GSL_SF_LEGENDRE_NONE:
			for (qpms_l_t l = 1; l <= lMax; ++l)
				for (qpms_m_t m = 1; m <= l; ++m) {
				qpms_y_t y = qpms_mn2y(-m,l);
				size_t i = gsl_sf_legendre_array_index(l,m);
				// viz DLMF 14.9.3, čert ví, jak je to s cs fasí.
				double factor = exp(lgamma(l-m+1)-lgamma(l+m+1))*((m%2)?-1:1);
				target[y] = factor * legendre_tmp[i];
				target_deriv[y] = factor * legendre_deriv_tmp[i];
				}
			break;
		case GSL_SF_LEGENDRE_SCHMIDT:
		case GSL_SF_LEGENDRE_SPHARM:
		case GSL_SF_LEGENDRE_FULL:
			for (qpms_l_t l = 1; l <= lMax; ++l)
				for (qpms_m_t m = 1; m <= l; ++m) {
				qpms_y_t y = qpms_mn2y(-m,l);
				size_t i = gsl_sf_legendre_array_index(l,m);
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

qpms_errno_t qpms_legendre_deriv_y_get(double **target, double **dtarget, double x, qpms_l_t lMax, gsl_sf_legendre_t lnorm,
		double csphase)
{

	*target = malloc(sizeof(double)*qpms_lMax2nelem(lMax));
	*dtarget = malloc(sizeof(double)*qpms_lMax2nelem(lMax));
	return qpms_legendre_deriv_y_fill(*target, *dtarget, x, lMax, lnorm, csphase);
}


qpms_pitau_t qpms_pitau_get(double theta, qpms_l_t lMax, qpms_normalisation_t norm) 
{
	const double csphase = qpms_normalisation_t_csphase(norm);
	norm = qpms_normalisation_t_normonly(norm);
	qpms_pitau_t res;
	qpms_y_t nelem = qpms_lMax2nelem(lMax);
	res.pi = malloc(nelem * sizeof(double));
	res.tau = malloc(nelem * sizeof(double));
	double ct = cos(theta), st = sin(theta);
	if (1 == fabs(ct)) { // singular case, use DLMF 14.8.2
		memset(res.pi, 0, nelem*sizeof(double));
		memset(res.tau, 0, nelem*sizeof(double));
		res.leg = calloc(nelem, sizeof(double));
		switch(norm) {
			case QPMS_NORMALISATION_XU:
				for (qpms_l_t l = 1; l <= lMax; ++l) {
					res.leg[qpms_mn2y(0, l)] = (l%2)?ct:1.;
					double p = l*(l+1)/2;
					const double n = 0.5;
					int lpar = (l%2)?-1:1;
					res.pi [qpms_mn2y(+1, l)] = -((ct>0) ? -1 : lpar) * p * csphase;
					res.pi [qpms_mn2y(-1, l)] = -((ct>0) ? -1 : lpar) * n * csphase;
					res.tau[qpms_mn2y(+1, l)] = ((ct>0) ? +1 : lpar) * p * csphase;
					res.tau[qpms_mn2y(-1, l)] = -((ct>0) ? +1 : lpar) * n * csphase;
				}
				break;
			case QPMS_NORMALISATION_TAYLOR:
				for (qpms_l_t l = 1; l <= lMax; ++l) {
					res.leg[qpms_mn2y(0, l)] = ((l%2)?ct:1.)*sqrt((2*l+1)*0.25*M_1_PI);
					int lpar = (l%2)?-1:1;
					double fl = 0.25 * sqrt((2*l+1)*l*(l+1)*M_1_PI);
					res.pi [qpms_mn2y(+1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
					res.pi [qpms_mn2y(-1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
					res.tau[qpms_mn2y(+1, l)] = ((ct>0) ? +1 : lpar) * fl * csphase;
					res.tau[qpms_mn2y(-1, l)] = -((ct>0) ? +1 : lpar) * fl * csphase;
				}
				break;
			case QPMS_NORMALISATION_POWER:
				for (qpms_l_t l = 1; l <= lMax; ++l) {
					res.leg[qpms_mn2y(0, l)] = ((l%2)?ct:1.)*sqrt((2*l+1)/(4*M_PI *l*(l+1)));
					int lpar = (l%2)?-1:1;
					double fl = 0.25 * sqrt((2*l+1)*M_1_PI);
					res.pi [qpms_mn2y(+1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
					res.pi [qpms_mn2y(-1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
					res.tau[qpms_mn2y(+1, l)] = ((ct>0) ? +1 : lpar) * fl * csphase;
					res.tau[qpms_mn2y(-1, l)] = -((ct>0) ? +1 : lpar) * fl * csphase;
	
				}
				break;
			default:
				abort();
		}
	}
	else { // cos(theta) in (-1,1), use normal calculation
		double *legder = malloc(sizeof(double)*qpms_lMax2nelem(lMax));
		res.leg = malloc(sizeof(double)*qpms_lMax2nelem(lMax));
		if (qpms_legendre_deriv_y_fill(res.leg, legder, ct, lMax, 
					norm == QPMS_NORMALISATION_XU ? GSL_SF_LEGENDRE_NONE
						: GSL_SF_LEGENDRE_SPHARM, csphase)) 
			abort();
		if (norm == QPMS_NORMALISATION_POWER) 
			/* for Xu (=non-normalized) and Taylor (=sph. harm. normalized) 
			 * the correct normalisation is already obtained from gsl_sf_legendre_deriv_array_e().
			 * However, Kristensson ("power") normalisation differs from Taylor
			 * by 1/sqrt(l*(l+1)) factor.
			 */
			for (qpms_l_t l = 1; l <= lMax; ++l) {
				double prefac = 1./sqrt(l*(l+1));
				for (qpms_m_t m = -l; m <= l; ++m) {
					res.leg[qpms_mn2y(m,l)] *= prefac;
					legder[qpms_mn2y(m,l)] *= prefac;
				}
			}
		for (qpms_l_t l = 1; l <= lMax; ++l) {
			for (qpms_m_t m = -l; m <= l; ++m) {
						res.pi [qpms_mn2y(m,l)] = m / st * res.leg[qpms_mn2y(m,l)];
						res.tau[qpms_mn2y(m,l)] = - st * legder[qpms_mn2y(m,l)];
			}
		}
		free(legder);
	}
	res.lMax = lMax;
	return res;
}

void qpms_pitau_free(qpms_pitau_t x) {
	free(x.leg);
	free(x.pi);
	free(x.tau);
}


csphvec_t qpms_vswf_single_el(qpms_m_t m, qpms_l_t l, sph_t kdlj,
		qpms_bessel_t btyp, qpms_normalisation_t norm) {
	lmcheck(l,m);
	csphvec_t N;
	complex double *bessel = malloc((l+1)*sizeof(complex double));
	if(qpms_sph_bessel_fill(btyp, l, kdlj.r, bessel)) abort();
	qpms_pitau_t pt = qpms_pitau_get(kdlj.theta, l, norm);
	complex double eimf = cexp(m * kdlj.phi * I);
	qpms_y_t y = qpms_mn2y(m,l);

	N.rc = l*(l+1) * pt.leg[y] * bessel[l] / kdlj.r * eimf;
	complex double besselfac = bessel[l-1] - l * bessel[l] / kdlj.r;
	N.thetac = pt.tau[y] * besselfac * eimf;
	N.phic = pt.pi[y] * besselfac * I * eimf;

	qpms_pitau_free(pt);
	free(bessel);
	return N;
}
csphvec_t qpms_vswf_single_mg(qpms_m_t m, qpms_l_t l, sph_t kdlj,
		qpms_bessel_t btyp, qpms_normalisation_t norm) {
	lmcheck(l,m);
	csphvec_t M;
	complex double *bessel = malloc((l+1)*sizeof(complex double));
	if(qpms_sph_bessel_fill(btyp, l, kdlj.r, bessel)) abort();
	qpms_pitau_t pt = qpms_pitau_get(kdlj.theta, l, norm);
	complex double eimf = cexp(m * kdlj.phi * I);
	qpms_y_t y = qpms_mn2y(m,l);

	M.rc = 0.;
	M.thetac = pt.pi[y] * bessel[l] * I * eimf;
	M.phic = -pt.tau[y] * bessel[l] * eimf;

	qpms_pitau_free(pt);
	free(bessel);
	return M;
}

qpms_vswfset_sph_t *qpms_vswfset_make(qpms_l_t lMax, sph_t kdlj, 
		qpms_bessel_t btyp, qpms_normalisation_t norm) {
	qpms_vswfset_sph_t *res = malloc(sizeof(qpms_vswfset_sph_t));
	res->lMax = lMax;
	qpms_y_t nelem = qpms_lMax2nelem(lMax);
	res->el = malloc(sizeof(csphvec_t)*nelem);
	res->mg = malloc(sizeof(csphvec_t)*nelem);
	if(QPMS_SUCCESS != qpms_vswf_fill(res->mg, res->el, lMax, kdlj, btyp, norm))
		abort(); // or return NULL? or rather assert?
	return res;
}

void qpms_vswfset_sph_pfree(qpms_vswfset_sph_t *w) {
	assert(NULL != w && NULL != w->el && NULL != w->mg);
	free(w->el);
	free(w->mg);
	free(w);
}

qpms_errno_t qpms_vswf_fill(csphvec_t *mgtarget, csphvec_t *eltarget,
		qpms_l_t lMax, sph_t kr,
		qpms_bessel_t btyp, qpms_normalisation_t norm) {
	assert(lMax >= 1);
	complex double *bessel = malloc((lMax+1)*sizeof(complex double));
	if(qpms_sph_bessel_fill(btyp, lMax, kr.r, bessel)) abort();
	qpms_pitau_t pt = qpms_pitau_get(kr.theta, lMax, norm);
	complex double const *pbes = bessel + 1; // starting from l = 1
	double const *pleg = pt.leg;
	double const *ppi = pt.pi;
	double const *ptau = pt.tau;
	csphvec_t *pmg = mgtarget, *pel = eltarget;
	for(qpms_l_t l = 1; l <= lMax; ++l) {
		complex double besfac = *pbes / kr.r;
		complex double besderfac = *(pbes-1) - l * besfac;
		for(qpms_m_t m = -l; m <= l; ++m) {
			complex double eimf = cexp(m * kr.phi * I);
			pel->rc = l*(l+1) * (*pleg) * besfac * eimf;
			pel->thetac = *ptau * besderfac * eimf;
			pel->phic = *ppi * besderfac * I * eimf;
			pmg->rc = 0.;
			pmg->thetac = *ppi * (*pbes) * I * eimf;
			pmg->phic = - *ptau * (*pbes) * eimf;
			++pleg; ++ppi; ++ptau; ++pel; ++pmg;
		}
		++pbes;
	}
	free(bessel);
	qpms_pitau_free(pt);
	return QPMS_SUCCESS;
}
