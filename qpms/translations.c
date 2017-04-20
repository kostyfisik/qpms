#include <math.h>
#include "gaunt.h"
#include "translations.h"
#include <stdbool.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <assert.h>


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
	    if(p < abs(Pp_order)) continue; // FIXME raději nastav lépe meze
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
		bool r_ge_d, qpms_bessel_t J) { // TODO make J enum
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
	return (sqrt(5+y)-2)/2; // the cast will truncate the fractional part, which is what we want
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
			qmaxsum += (qmaxes[qpms_trans_calculator_index_yyu(c,y,yu)] = gaunt_q_max(-m,n,mu,nu));
		}
	c->A_multipliers[0] = malloc(qmaxsum * sizeof(complex double));
	for(size_t i = 0; i < SQ(c->nelem); ++i)
		c->A_multipliers[i+1] = c->A_multipliers[i] + qmaxes[i];
	// TODO here comes the filling of A_multipliers
	


	qmaxsum = 0;
	for(size_t y=0; y < c->nelem; y++)
		for(size_t yu = 0; yu < c->nelem; yu++) {
			int m, n, mu, nu;
			qpms_y2mn_p(y,&m,&n);
			qpms_y2mn_p(yu,&mu,&nu);
			qmaxsum += (qmaxes[qpms_trans_calculator_index_yyu(c,y,yu)] = gaunt_q_max(-m,n+1,mu,nu));
		}
	c->B_multipliers[0] = malloc(qmaxsum * sizeof(complex double));
	for(size_t i = 0; i < SQ(c->nelem); ++i)
		c->B_multipliers[i+1] = c->B_multipliers[i] + qmaxes[i];

	// TODO here comes the filling of B_multipliers

	free(qmaxes);
	return c;
}


