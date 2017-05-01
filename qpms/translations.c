#include <math.h>
#include "gaunt.h"
#include "translations.h"
#include <stdbool.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include "assert_cython_workaround.h"


static const double sqrtpi = 1.7724538509055160272981674833411451827975494561223871;
//static const double ln2 = 0.693147180559945309417232121458176568075500134360255254120;

// Associated Legendre polynomial at zero argument (DLMF 14.5.1)
double qpms_legendre0(int m, int n) {
	return pow(2,m) * sqrtpi / tgamma(.5*n - .5*m + .5) / tgamma(.5*n-.5*m);
}

static inline int min1pow(int x) {
	return (x % 2) ? -1 : 1;
}

static inline complex double ipow(int x) {
	return cpow(I, x);
}

// Derivative of associated Legendre polynomial at zero argument (DLMF 14.5.2)
double qpms_legendreD0(int m, int n) {
	return -2 * qpms_legendre0(m, n);
}

int qpms_sph_bessel_array(qpms_bessel_t typ, int lmax, double x, complex double *result_array) {
	int retval;
	double tmparr[lmax+1];
	switch(typ) {
		case QPMS_BESSEL_REGULAR:
			retval = gsl_sf_bessel_jl_steed_array(lmax, x, tmparr);
			for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
			return retval;
			break;
		case QPMS_BESSEL_SINGULAR: //FIXME: is this precise enough? Would it be better to do it one-by-one?
			retval = gsl_sf_bessel_yl_array(lmax,x,tmparr);
			for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
			return retval;
			break;
		case QPMS_HANKEL_PLUS:
		case QPMS_HANKEL_MINUS:
			retval = gsl_sf_bessel_jl_steed_array(lmax, x, tmparr);
			for (int l = 0; l <= lmax; ++l) result_array[l] = tmparr[l];
			if(retval) return retval;
			retval = gsl_sf_bessel_yl_array(lmax, x, tmparr);
			if (typ==QPMS_HANKEL_PLUS)
				for (int l = 0; l <= lmax; ++l) result_array[l] += I * tmparr[l];
			else 
				for (int l = 0; l <= lmax; ++l) result_array[l] +=-I * tmparr[l];
			return retval;
			break;
		default:
			abort();
			//return GSL_EDOM;
	}
	assert(0);
}


complex double qpms_trans_single_A_Taylor(int m, int n, int mu, int nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J) { // TODO make J enum
	if(r_ge_d) J = QPMS_BESSEL_REGULAR;

	double costheta = cos(kdlj.theta);

	int qmax = gaunt_q_max(-m,n,mu,nu); // nemá tu být +m?
	// N.B. -m !!!!!!
	double a1q[qmax+1];
	int err;
	gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
	double a1q0 = a1q[0];
	if (err) abort();

	double leg[gsl_sf_legendre_array_n(n+nu)];
	if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu,costheta,-1,leg)) abort();
	complex double bes[n+nu+1];
	if (qpms_sph_bessel_array(J, n+nu, kdlj.r, bes)) abort();
	complex double sum = 0;
	for(int q = 0; q <= qmax; ++q) {
		int p = n+nu-2*q;
		int Pp_order = mu-m;
		//if(p < abs(Pp_order)) continue; // FIXME raději nastav lépe meze
		assert(p >= abs(Pp_order));
		double a1q_n = a1q[q] / a1q0;
		double Pp = leg[gsl_sf_legendre_array_index(p, abs(Pp_order))];
		if (Pp_order < 0) Pp *= min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
		complex double zp = bes[p];
		complex double summandq = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n * zp * Pp;
		sum += summandq;
	}

	double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
			+lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
			+lgamma(n+nu+1) - lgamma(2*(n+nu)+1));
	complex double presum = exp(exponent);
	presum *= cexp(I*(mu-m)*kdlj.phi) * min1pow(m) * ipow(nu+n) / (4*n);

	complex double prenormratio = ipow(nu-n) *  sqrt(((2.*nu+1)/(2.*n+1))* exp(
				lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1)));
	return (presum / prenormratio) * sum;
}

complex double qpms_trans_single_B_Taylor(int m, int n, int mu, int nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J) { 
	if(r_ge_d) J = QPMS_BESSEL_REGULAR;
	double costheta = cos(kdlj.theta);

	int q2max = gaunt_q_max(-m-1,n+1,mu+1,nu);
	int Qmax = gaunt_q_max(-m,n+1,mu,nu);
	double a2q[q2max+1], a3q[Qmax+1], a2q0, a3q0;
	int err;
	if (mu == nu) {
		for (int q = 0; q <= q2max; ++q) 
			a2q[q] = 0;
		a2q0 = 1;
	}
	else { 
		gaunt_xu(-m-1,n+1,mu+1,nu,q2max,a2q,&err); if (err) abort();
		a2q0 = a2q[0];
	}
	gaunt_xu(-m,n+1,mu,nu,Qmax,a3q,&err); if (err) abort();
	a3q0 = a3q[0];

	double leg[gsl_sf_legendre_array_n(n+nu+1)];
	if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,costheta,-1,leg)) abort();
	complex double bes[n+nu+2];
	if (qpms_sph_bessel_array(J, n+nu+1, kdlj.r, bes)) abort();

	complex double sum = 0;
	for (int q = 0; q <= Qmax; ++q) {
		int p = n+nu-2*q;
		double a2q_n = a2q[q]/a2q0;
		double a3q_n = a3q[q]/a3q0;
		complex double zp_ = bes[p+1];
		int Pp_order_ = mu-m;
		//if(p+1 < abs(Pp_order_)) continue; // FIXME raději nastav lépe meze
		assert(p+1 >= abs(Pp_order_));
		double Pp_ = leg[gsl_sf_legendre_array_index(p+1, abs(Pp_order_))];
		if (Pp_order_ < 0) Pp_ *= min1pow(mu-m) * exp(lgamma(1+1+p+Pp_order_)-lgamma(1+1+p-Pp_order_));
		complex double summandq = ((2*(n+1)*(nu-mu)*a2q_n
					-(-nu*(nu+1) - n*(n+3) - 2*mu*(n+1)+p*(p+3))* a3q_n)
				*min1pow(q) * zp_ * Pp_);
		sum += summandq;
	}

	double exponent=(lgamma(2*n+3)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
			+lgamma(n+nu+m-mu+2)-lgamma(n-m+1)-lgamma(nu+mu+1)
			+lgamma(n+nu+2) - lgamma(2*(n+nu)+3));
	complex double presum = exp(exponent);
	presum *= cexp(I*(mu-m)*kdlj.phi) * min1pow(m) * ipow(nu+n+1) / (
			(4*n)*(n+1)*(n+m+1));

	// Taylor normalisation v2, proven to be equivalent
	complex double prenormratio = ipow(nu-n) * sqrt(((2.*nu+1)/(2.*n+1))* exp(
				lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1)));

	return (presum / prenormratio) * sum;
}

complex double qpms_trans_single_A_Taylor_ext(int m, int n, int mu, int nu, 
		double kdlj_r, double kdlj_theta, double kdlj_phi, int r_ge_d, int J) {
	sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
	return qpms_trans_single_A_Taylor(m,n,mu,nu,kdlj,r_ge_d,J);
}

complex double qpms_trans_single_B_Taylor_ext(int m, int n, int mu, int nu, 
		double kdlj_r, double kdlj_theta, double kdlj_phi, int r_ge_d, int J) {
	sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
	return qpms_trans_single_B_Taylor(m,n,mu,nu,kdlj,r_ge_d,J);
}

void qpms_trans_calculator_free(qpms_trans_calculator *c) {
	free(c->A_multipliers[0]);
	free(c->A_multipliers);
	free(c->B_multipliers[0]);
	free(c->B_multipliers);
	free(c);
}

static inline size_t qpms_mn2y(int m, int n) {
	return (size_t) n * (n + 1) + m - 1;
}

static inline int qpms_y2n(size_t y) {
	//return (sqrt(5+y)-2)/2; // the cast will truncate the fractional part, which is what we want
	return sqrt(y+1);
}

static inline int qpms_yn2m(size_t y, int n) {
	return y-qpms_mn2y(0,n);
}

static inline void qpms_y2mn_p(size_t y, int *m, int *n){
	*m=qpms_yn2m(y,*n=qpms_y2n(y));
}

static inline size_t qpms_trans_calculator_index_mnmunu(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu){
	return c->nelem * qpms_mn2y(m,n) + qpms_mn2y(mu,nu);
}

static inline size_t qpms_trans_calculator_index_yyu(const qpms_trans_calculator *c,
		size_t y, size_t yu) {
	return c->nelem * y + yu;
}


#define SQ(x) ((x)*(x))

//#if 0
static void qpms_trans_calculator_multipliers_A_Taylor(
		complex double *dest, int m, int n, int mu, int nu, int qmax) {
	assert(qmax == gaunt_q_max(-m,n,mu,nu));
	double a1q[qmax+1];
	int err;
	gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
	if (err) abort();
	double a1q0 = a1q[0];

	double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
			+lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
			+lgamma(n+nu+1) - lgamma(2*(n+nu)+1)) - 0.5*( // ex-prenormratio
			lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1));
	double presum = exp(exponent);
	presum *=  min1pow(m+n) * sqrt((2.*n+1)/(2.*nu+1)) / (4*n);

	for(int q = 0; q <= qmax; q++) {
		int p = n+nu-2*q;
		int Pp_order = mu - m;
		assert(p >= abs(Pp_order));
		double a1q_n = a1q[q] / a1q0;
		// Assuming non_normalized legendre polynomials!
		double Ppfac = (Pp_order >= 0) ? 1 :
			min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
		double summandfac = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n;
		dest[q] = presum * summandfac * Ppfac;
		// FIXME I might not need complex here
	}
}
//#endif
#if 0
static void qpms_trans_calculator_multipliers_A_Taylor(
		complex double *dest, int m, int n, int mu, int nu, int qmax) {
	assert(qmax == gaunt_q_max(-m,n,mu,nu));
	double a1q[qmax+1];
	int err;
	gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
	if (err) abort();
	double a1q0 = a1q[0];
	for(int q = 0; q <= qmax; ++q) {
		int p = n+nu-2*q;
		int Pp_order = mu-m;
		//if(p < abs(Pp_order)) continue; // FIXME raději nastav lépe meze
		assert(p >= abs(Pp_order));
		double a1q_n = a1q[q] / a1q0;
		//double Pp = leg[gsl_sf_legendre_array_index(p, abs(Pp_order))];
		//complex double zp = bes[p];
		dest[q] = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n /* * zp * Pp*/;
		if (Pp_order < 0) dest[q] *= min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
		//sum += summandq;
	}

	double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
			+lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
			+lgamma(n+nu+1) - lgamma(2*(n+nu)+1));
	complex double presum = exp(exponent);
	presum *=/* cexp(I*(mu-m)*kdlj.phi) * */  min1pow(m) * ipow(nu+n) / (4*n);

	complex double prenormratio = ipow(nu-n) *  sqrt(((2.*nu+1)/(2.*n+1))* exp(
				lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1)));
	//return (presum / prenormratio) * sum;
	for(int q=0;q<=qmax;++q) dest[q] *= presum / prenormratio;
}
#endif	



static void qpms_trans_calculator_multipliers_B_Taylor(
		complex double *dest, int m, int n, int mu, int nu, int Qmax) {
	assert(Qmax == gaunt_q_max(-m,n+1,mu,nu));
	int q2max = gaunt_q_max(-m-1,n+1,mu+1,nu);
	assert(Qmax == q2max);
	// FIXME remove the q2max variable altogether, as it is probably equal
	// to Qmax
	double a2q[q2max+1], a3q[Qmax+1], a2q0, a3q0;
	int err;
	if (mu == nu) {
		for (int q = 0; q <= q2max; ++q)
			a2q[q] = 0;
		a2q0 = 1;
	}
	else {
		gaunt_xu(-m-1,n+1,mu+1,nu,q2max,a2q,&err); if (err) abort();
		a2q0 = a2q[0];
	}
	gaunt_xu(-m,n+1,mu,nu,Qmax,a3q,&err); if (err) abort();
	a3q0 = a3q[0];

	double exponent=(lgamma(2*n+3)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
			+lgamma(n+nu+m-mu+2)-lgamma(n-m+1)-lgamma(nu+mu+1)
			+lgamma(n+nu+2) - lgamma(2*(n+nu)+3)) - 0.5 * (
			lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)
			-lgamma(nu+mu+1));
	complex double presum = exp(exponent);
	presum *= I * (min1pow(m+n) *sqrt((2.*n+1)/(2.*nu+1)) / (
				(4*n)*(n+1)*(n+m+1)));

	for (int q = 0; q <= Qmax; ++q) {
		int p = n+nu-2*q;
		double a2q_n = a2q[q]/a2q0;
		double a3q_n = a3q[q]/a3q0;
		int Pp_order_ = mu-m;
		//if(p+1 < abs(Pp_order_)) continue; // FIXME raději nastav lépe meze
		assert(p+1 >= abs(Pp_order_));
		double Ppfac = (Pp_order_ >= 0) ? 1 :		

			min1pow(mu-m) * exp(lgamma(1+1+p+Pp_order_)-lgamma(1+1+p-Pp_order_));
		double summandq = ((2*(n+1)*(nu-mu)*a2q_n
					-(-nu*(nu+1) - n*(n+3) - 2*mu*(n+1)+p*(p+3))* a3q_n)
				*min1pow(q));
		dest[q] = Ppfac * summandq * presum;
	}
}

int qpms_trans_calculator_multipliers_A(qpms_normalization_t norm, complex double *dest, int m, int n, int mu, int nu, int qmax) {
	switch (norm) {
		case QPMS_NORMALIZATION_TAYLOR:
			qpms_trans_calculator_multipliers_A_Taylor(dest,m,n,mu,nu,qmax);
			return 0;
			break;
		default:
			abort();
	}
	assert(0);
}

int qpms_trans_calculator_multipliers_B(qpms_normalization_t norm, complex double *dest, int m, int n, int mu, int nu, int qmax) {
	switch (norm) {
		case QPMS_NORMALIZATION_TAYLOR:
			qpms_trans_calculator_multipliers_B_Taylor(dest,m,n,mu,nu,qmax);
			return 0;
			break;
		default:
			abort();
	}
	assert(0);
}

qpms_trans_calculator
*qpms_trans_calculator_init (int lMax, qpms_normalization_t normalization) {
	qpms_trans_calculator *c = malloc(sizeof(qpms_trans_calculator));
	c->lMax = lMax;
	c->nelem = lMax * (lMax+2);
	c->A_multipliers = malloc((1+SQ(c->nelem)) * sizeof(complex double *));
	c->B_multipliers = malloc((1+SQ(c->nelem)) * sizeof(complex double *));
	c->normalization = normalization;
	size_t *qmaxes = malloc(SQ(c->nelem) * sizeof(size_t));
	size_t qmaxsum = 0;
	for(size_t y = 0; y < c->nelem; y++)
		for(size_t yu = 0; yu < c->nelem; yu++) {
			int m,n, mu, nu;
			qpms_y2mn_p(y,&m,&n);
			qpms_y2mn_p(yu,&mu,&nu);
			qmaxsum += 1 + (
					qmaxes[qpms_trans_calculator_index_yyu(c,y,yu)] 
					= gaunt_q_max(-m,n,mu,nu));
		}
	c->A_multipliers[0] = malloc(qmaxsum * sizeof(complex double));
	// calculate multiplier beginnings
	for(size_t i = 0; i < SQ(c->nelem); ++i) 
		c->A_multipliers[i+1] = c->A_multipliers[i] + qmaxes[i] + 1;
	// calculate the multipliers
	for(size_t y = 0; y < c->nelem; ++y)
		for(size_t yu = 0; yu < c->nelem; ++yu) {
			size_t i = y * c->nelem + yu;
			int m, n, mu, nu;
			qpms_y2mn_p(y, &m, &n);
			qpms_y2mn_p(yu, &mu, &nu);
			qpms_trans_calculator_multipliers_A(normalization,
					c->A_multipliers[i], m, n, mu, nu, qmaxes[i]);
		}

	qmaxsum = 0;
	for(size_t y=0; y < c->nelem; y++)
		for(size_t yu = 0; yu < c->nelem; yu++) {
			int m, n, mu, nu;
			qpms_y2mn_p(y,&m,&n);
			qpms_y2mn_p(yu,&mu,&nu);
			qmaxsum += 1 + (
					qmaxes[qpms_trans_calculator_index_yyu(c,y,yu)] 
					= gaunt_q_max(-m,n+1,mu,nu));
		}
	c->B_multipliers[0] = malloc(qmaxsum * sizeof(complex double));
	// calculate multiplier beginnings
	for(size_t i = 0; i < SQ(c->nelem); ++i) 
		c->B_multipliers[i+1] = c->B_multipliers[i] + qmaxes[i] + 1;
	// calculate the multipliers
	for(size_t y = 0; y < c->nelem; ++y)
		for(size_t yu = 0; yu < c->nelem; ++yu) {
			size_t i = y * c->nelem + yu;
			int m, n, mu, nu;
			qpms_y2mn_p(y, &m, &n);
			qpms_y2mn_p(yu, &mu, &nu);
			qpms_trans_calculator_multipliers_B(normalization,
					c->B_multipliers[i], m, n, mu, nu, qmaxes[i]);
		}

	free(qmaxes);
	return c;
}

static inline complex double qpms_trans_calculator_get_A_precalcbuf(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J,
		const complex double *bessel_buf, const double *legendre_buf) {
	size_t i = qpms_trans_calculator_index_mnmunu(c, m, n, mu, nu);
	size_t qmax = c->A_multipliers[i+1] - c->A_multipliers[i] - 1;
	assert(qmax == gaunt_q_max(-m,n,mu,nu));
	complex double sum = 0;
	for(size_t q = 0; q <= qmax; ++q) {
		int p = n+nu-2*q;
		double Pp = legendre_buf[gsl_sf_legendre_array_index(p, abs(mu-m))];
		complex double zp = bessel_buf[p];
		complex double multiplier = c->A_multipliers[i][q];
		sum += Pp * zp *  multiplier;
	}
	complex double eimf =  cexp(I*(mu-m)*kdlj.phi);
	return sum * eimf;
}

complex double qpms_trans_calculator_get_A_buf(const qpms_trans_calculator *c,
								int m, int n, int mu, int nu, sph_t kdlj,
								bool r_ge_d, qpms_bessel_t J,
								complex double *bessel_buf, double *legendre_buf) {
	// This functions gets preallocated memory for bessel and legendre functions, but computes them itself
	if (r_ge_d) J = QPMS_BESSEL_REGULAR;
	if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) 
		// TODO warn? 
		return NAN+I*NAN;
	switch(c->normalization) {
		case QPMS_NORMALIZATION_TAYLOR:
			{
				double costheta = cos(kdlj.theta);
				if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu,
							costheta,-1,legendre_buf)) abort();
				if (qpms_sph_bessel_array(J, n+nu+1, kdlj.r, bessel_buf)) abort();
				return qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
						kdlj,r_ge_d,J,bessel_buf,legendre_buf);
			}
			break;
		default:
			abort();
	}
	assert(0);
}

static inline complex double qpms_trans_calculator_get_B_precalcbuf(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J,
		const complex double *bessel_buf, const double *legendre_buf) {
	size_t i = qpms_trans_calculator_index_mnmunu(c, m, n, mu, nu);
	size_t qmax = c->B_multipliers[i+1] - c->B_multipliers[i] - 1;
	assert(qmax == gaunt_q_max(-m,n+1,mu,nu));
	complex double sum = 0;
	for(int q = 0; q <= qmax; ++q) {
		int p = n+nu-2*q;
		double Pp_ = legendre_buf[gsl_sf_legendre_array_index(p+1, abs(mu-m))];
		complex double zp_ = bessel_buf[p+1];
		complex double multiplier = c->B_multipliers[i][q];
		sum += Pp_ * zp_ * multiplier;
	}
	complex double eimf =  cexp(I*(mu-m)*kdlj.phi);
	return sum * eimf;
}

complex double qpms_trans_calculator_get_B_buf(const qpms_trans_calculator *c,
								int m, int n, int mu, int nu, sph_t kdlj,
								bool r_ge_d, qpms_bessel_t J,
								complex double *bessel_buf, double *legendre_buf) {
	// This functions gets preallocated memory for bessel and legendre functions, but computes them itself
	if (r_ge_d) J = QPMS_BESSEL_REGULAR;
	if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) 
		// TODO warn? 
		return NAN+I*NAN;
	switch(c->normalization) {
		case QPMS_NORMALIZATION_TAYLOR:
			{
				double costheta = cos(kdlj.theta);
				if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,
							costheta,-1,legendre_buf)) abort();
				if (qpms_sph_bessel_array(J, n+nu+2, kdlj.r, bessel_buf)) abort();
				return qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
						kdlj,r_ge_d,J,bessel_buf,legendre_buf);
			}
			break;
		default:
			abort();
	}
	assert(0);
}

int qpms_trans_calculator_get_AB_buf_p(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		int m, int n, int mu, int nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J,
		complex double *bessel_buf, double *legendre_buf) {
	if (r_ge_d) J = QPMS_BESSEL_REGULAR;
	if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) {
		*Adest = NAN+I*NAN;
		*Bdest = NAN+I*NAN;
		// TODO warn? different return value?
		return 0;
	}
	switch(c->normalization) {
		case QPMS_NORMALIZATION_TAYLOR:
			{
				double costheta = cos(kdlj.theta);
				if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu+1,
							costheta,-1,legendre_buf)) abort();
				if (qpms_sph_bessel_array(J, n+nu+2, kdlj.r, bessel_buf)) abort();
				*Adest = qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
						kdlj,r_ge_d,J,bessel_buf,legendre_buf);
				*Bdest = qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
						kdlj,r_ge_d,J,bessel_buf,legendre_buf);
				return 0;
			}
			break;
		default:
			abort();
	}
	assert(0);
}




int qpms_trans_calculator_get_AB_arrays_buf(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		size_t deststride, size_t srcstride,
		sph_t kdlj, bool r_ge_d, qpms_bessel_t J,
		complex double *bessel_buf, double *legendre_buf) {
	if (r_ge_d) J = QPMS_BESSEL_REGULAR;
	if (0 == kdlj.r && J != QPMS_BESSEL_REGULAR) {
		for (size_t i = 0; i < c->nelem; ++i)
			for (size_t j = 0; j < c->nelem; ++j) {
				*(Adest + i*srcstride + j*deststride) = NAN+I*NAN;
				*(Bdest + i*srcstride + j*deststride) = NAN+I*NAN;
			}
		// TODO warn? different return value?
		return 0;
	}
	switch(c->normalization) {
		case QPMS_NORMALIZATION_TAYLOR:
			{
				double costheta = cos(kdlj.theta);
				if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,2*c->lMax+1,
							costheta,-1,legendre_buf)) abort();
				if (qpms_sph_bessel_array(J, 2*c->lMax+2, kdlj.r, bessel_buf)) abort();
				size_t desti = 0, srci = 0;
				for (int n = 1; n <= c->nelem; ++n) for (int m = -n; m <= n; ++m) {
					for (int nu = 1; nu <= c->nelem; ++nu) for (int mu = -nu; mu <= nu; ++mu) {
						assert(qpms_trans_calculator_index_mnmunu(c,m,n,mu,nu) == desti*c->nelem + srci);
						*(Adest + deststride * desti + srcstride * srci) = 
							qpms_trans_calculator_get_A_precalcbuf(c,m,n,mu,nu,
									kdlj,r_ge_d,J,bessel_buf,legendre_buf);
						*(Bdest + deststride * desti + srcstride * srci) = 
							qpms_trans_calculator_get_B_precalcbuf(c,m,n,mu,nu,
									kdlj,r_ge_d,J,bessel_buf,legendre_buf);
						++srci;
					}
					++desti;
				}
				return 0;
			}
			break;
		default:
			abort();
	}
	assert(0);
}

complex double qpms_trans_calculator_get_A(const qpms_trans_calculator *c,
                int m, int n, int mu, int nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J) {
	double leg[gsl_sf_legendre_array_n(n+nu)];
	complex double bes[n+nu+1];
	return qpms_trans_calculator_get_A_buf(c,m,n,mu,nu,kdlj,r_ge_d,J,
			bes,leg);
}

complex double qpms_trans_calculator_get_B(const qpms_trans_calculator *c,
                int m, int n, int mu, int nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J) {
	double leg[gsl_sf_legendre_array_n(n+nu+1)];
	complex double bes[n+nu+2];
	return qpms_trans_calculator_get_B_buf(c,m,n,mu,nu,kdlj,r_ge_d,J,
			bes,leg);
}

int qpms_trans_calculator_get_AB_p(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
                int m, int n, int mu, int nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J) {
	double leg[gsl_sf_legendre_array_n(2*c->lMax+1)];
	complex double bes[2*c->lMax+2];
	return qpms_trans_calculator_get_AB_buf_p(c,Adest, Bdest,m,n,mu,nu,kdlj,r_ge_d,J,
			bes,leg);
}

int qpms_trans_calculator_get_AB_arrays(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		size_t deststride, size_t srcstride,
		sph_t kdlj, bool r_ge_d, qpms_bessel_t J) {
	double leg[gsl_sf_legendre_array_n(c->lMax+c->lMax+1)];
	complex double bes[c->lMax+c->lMax+2];
	return qpms_trans_calculator_get_AB_arrays_buf(c, 
			Adest, Bdest, deststride, srcstride,
			kdlj, r_ge_d, J, 
			bes, leg);
}



complex double qpms_trans_calculator_get_A_ext(const qpms_trans_calculator *c,
                int m, int n, int mu, int nu,
	       	double kdlj_r, double kdlj_theta, double kdlj_phi,
                int r_ge_d, int J) {
	sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
	return qpms_trans_calculator_get_A(c,m,n,mu,nu,kdlj,r_ge_d,J);
}

complex double qpms_trans_calculator_get_B_ext(const qpms_trans_calculator *c,
                int m, int n, int mu, int nu,
	       	double kdlj_r, double kdlj_theta, double kdlj_phi,
                int r_ge_d, int J) {
	sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
	return qpms_trans_calculator_get_B(c,m,n,mu,nu,kdlj,r_ge_d,J);
}

int qpms_trans_calculator_get_AB_p_ext(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
                int m, int n, int mu, int nu,
	       	double kdlj_r, double kdlj_theta, double kdlj_phi,
                int r_ge_d, int J) {
	sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
	return qpms_trans_calculator_get_AB_p(c,Adest,Bdest,m,n,mu,nu,kdlj,r_ge_d,J);
}

int qpms_trans_calculator_get_AB_arrays_ext(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		size_t deststride, size_t srcstride,
		double kdlj_r, double kdlj_theta, double kdlj_phi,
	 	int r_ge_d, int J) {
	sph_t kdlj = {kdlj_r, kdlj_theta, kdlj_phi};
	return qpms_trans_calculator_get_AB_arrays(c,Adest,Bdest,deststride,srcstride,
						kdlj, r_ge_d, J);
}

