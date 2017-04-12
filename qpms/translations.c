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
				for (int l = 0; l <= lmax; ++l) result_array[l] *= I * tmparr[l];
			else 
				for (int l = 0; l <= lmax; ++l) result_array[l] *=-I * tmparr[l];
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
    double exponent=(lgamma(2*n+1)-lgamma(n+2)+lgamma(2*nu+3)-lgamma(nu+2)
                +lgamma(n+nu+m-mu+1)-lgamma(n-m+1)-lgamma(nu+mu+1)
                +lgamma(n+nu+1) - lgamma(2*(n+nu)+1));
    double costheta = cos(kdlj.theta);
    int qmax = gaunt_q_max(-m,n,mu,nu); // nemá tu být +m?
    // N.B. -m !!!!!!
    double a1q[qmax+1];
    int err;
    gaunt_xu(-m,n,mu,nu,qmax,a1q,&err);
    double a1q0 = a1q[0];
    if (err) abort();
    //double *leg = malloc(sizeof(double)*gsl_sf_legendre_array_n(n+nu));
    //if (!leg) abort();
    double leg[gsl_sf_legendre_array_n(n+nu)];
    if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,n+nu,costheta,-1,leg)) abort();
    complex double bes[n+nu+1];
    if (qpms_sph_bessel_array(J, n+nu, kdlj.r, bes)) abort();
    complex double sum = 0;
    for(int q = 0; q <= qmax; ++q) {
	    int p = n+nu-2*q;
	    double a1q_n = a1q[q] / a1q0;
	    int Pp_order = mu-m;
	    double Pp = leg[gsl_sf_legendre_array_index(p, abs(Pp_order))];
	    if (Pp_order < 0) Pp *= min1pow(mu-m) * exp(lgamma(1+p+Pp_order)-lgamma(1+p-Pp_order));
	    complex double zp = bes[p];
	    complex double summandq = (n*(n+1) + nu*(nu+1) - p*(p+1)) * min1pow(q) * a1q_n * zp * Pp;
	    sum += summandq;
    }
    //free(leg);
    //free(bes);
    complex double presum = exp(exponent);
    presum *= cexp(I*(mu-m)*kdlj.phi) * min1pow(m) * ipow(nu+n) / (4*n);
 
    complex double prenormratio = ipow(nu-n) *  sqrt(((2.*nu+1)/(2.*n+1))* exp(
	        lgamma(n+m+1)-lgamma(n-m+1)+lgamma(nu-mu+1)-lgamma(nu+mu+1)));
    return (presum / prenormratio) * sum;
}

