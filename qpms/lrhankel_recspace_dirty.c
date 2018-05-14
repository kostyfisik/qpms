/*
 * This is a dirty implementation of lrhankel_recpart_fill() that calculates
 * the (cylindrical) Hankel transforms of the regularised part of the spherical Hankel
 * functions that are to be summed in the reciprocal space.
 *
 * For now, only the regularisation with κ == 5 && q <= 2 && n <= 5 is implemented
 * by writing down the explicit formula for each q,n pair and k>k0 vs k<k0 case,
 * only with the help of some macros to make the whole thing shorter.
 * 
 * N.B. the results for very small k/k0 differ significantly (sometimes even in the first
 * digit for n >= 3, probably due to catastrophic cancellation (hopefully not due
 * to an error in the formula!). On the other hand,
 * these numbers are tiny in their absolute value, so their contribution to the
 * lattice sum should be negligible.
 *
 * Therefore TODO use kahan summation.
 *
 */

#define MAXQM 1
#define MAXN 5
#define MAXKAPPA 5

#include "bessels.h"
//#include "mdefs.h"
#include <complex.h>
#include <string.h>
#include <assert.h>
#define SQ(x) ((x)*(x))

#define P4(x) (((x)*(x))*((x)*(x)))


/* 
 * General form of the κ == 5 transforms. One usually has to put a (-1) factor
 * to some part of the zeroth term of one of the cases k < k0 or k > k0 in order
 * to stay on the correct branch of complex square root...
 */
#define KAPPA5SUM(form) (\
                (form(0, 1))  \
            - 5*(form(1, 1)) \
            +10*(form(2, 1)) \
            -10*(form(3, 1)) \
            + 5*(form(4, 1)) \
            -   (form(5, 1)) \
)

#define KAPPA5SUMFF(form) (\
                (form(0, (-1)))  \
            - 5*(form(1, 1)) \
            +10*(form(2, 1)) \
            -10*(form(3, 1)) \
            + 5*(form(4, 1)) \
            -   (form(5, 1)) \
)

/* 
 * Prototype for the individual (per q,n) Bessel transform calculating functions.
 * a, b, d, e, ash are recurring pre-calculated intermediate results, see the definition
 * of lrhankel_recpart_fill() below to see their meaning
 */
#define LRHANKELDEF(fname) complex double fname(const double c, const double k0, const double k, \
		const complex double *a, const complex double *b, const complex double *d, \
		const complex double *e, const complex double *ash)

typedef complex double (*lrhankelspec)(const double, const double, const double, 
		const complex double *,
		const complex double *,
		const complex double *,
		const complex double *,
		const complex double *);
// complex double fun(double c, double k0, double k, ccd *a, ccd *b, ccd *d, ccd *e)

#define FORMK5Q1N0(i,ff) (ff*e[i])
LRHANKELDEF(fk5q1n0l){
	return (KAPPA5SUMFF(FORMK5Q1N0))/k0;
}
LRHANKELDEF(fk5q1n0s){
	return (KAPPA5SUM(FORMK5Q1N0))/k0;
}
#undef FORMK5Q1N0

#define FORMK5Q1N1(i,ff) (-ff*d[i])
LRHANKELDEF(fk5q1n1l){
	return (KAPPA5SUMFF(FORMK5Q1N1))/(k0*k);
}
LRHANKELDEF(fk5q1n1s){
	return (KAPPA5SUM(FORMK5Q1N1))/(k0*k);
}
#undef FORMK5Q1N1

#define FORMK5Q1N2(i,ff) (ff*e[i] - t*a[i] + ff*t*d[i]*a[i])
LRHANKELDEF(fk5q1n2l){
	double t = 2/(k*k);
	return (KAPPA5SUMFF(FORMK5Q1N2))/k0;
}
LRHANKELDEF(fk5q1n2s){
	double t = 2/(k*k);
	return (KAPPA5SUM(FORMK5Q1N2))/k0;
}
#undef FORMK5Q1N2

#define FORMK5Q1N3(i,ff) (-ff*d[i] * (kk3 + 4*a[i]*a[i]))
LRHANKELDEF(fk5q1n3l){
	double kk3 = 3*k*k;
	return (KAPPA5SUMFF(FORMK5Q1N3))/(k0*k*k*k);
}
LRHANKELDEF(fk5q1n3s){
	double kk3 = 3*k*k;
	return (KAPPA5SUM(FORMK5Q1N3))/(k0*k*k*k);	
}
#undef FORMK5Q1N3

#define FORMK5Q1N4(i,ff) (ff*e[i] * (kkkk + kk8*a[i]*a[i] + 8*P4(a[i])))
LRHANKELDEF(fk5q1n4l){
	double kk8 = k*k*8, kkkk = P4(k); 
	return (KAPPA5SUMFF(FORMK5Q1N4))/(k0*kkkk);
}
LRHANKELDEF(fk5q1n4s){
	double kk8 = k*k*8, kkkk = P4(k); 
	return (KAPPA5SUM(FORMK5Q1N4))/(k0*kkkk);
}
#undef FORMK5Q1N4

#define FORMK5Q1N5(i,ff) (d[i]*(kkkk*(-5*ff+b[i])-ff*kk20*a[i]*a[i]-ff*16*P4(a[i])))
LRHANKELDEF(fk5q1n5l){
	double kk20 = k*k*20, kkkk = P4(k); 
	return (KAPPA5SUMFF(FORMK5Q1N5))/(k0*kkkk*k);
}
LRHANKELDEF(fk5q1n5s){
	double kk20 = k*k*20, kkkk = P4(k); 
	return (KAPPA5SUM(FORMK5Q1N5))/(k0*kkkk*k);
}
#undef FORMK5Q1N5

#define FORMK5Q2N0(i,ff) (-ash[i])
LRHANKELDEF(fk5q2n0){
	return (KAPPA5SUM(FORMK5Q2N0)) / (k0*k0);
}
const lrhankelspec fk5q2n0s = fk5q2n0, fk5q2n0l = fk5q2n0;
#undef FORMK5Q2N0

#define FORMK5Q2N1(i,ff) (ff*b[i]*a[i])
LRHANKELDEF(fk5q2n1l){
	return (KAPPA5SUMFF(FORMK5Q2N1))/(k*k0*k0);
}
LRHANKELDEF(fk5q2n1s){
	return (KAPPA5SUM(FORMK5Q2N1))/(k*k0*k0);
}
#undef FORMK5Q2N1

#define FORMK5Q2N2(i,ff) (-ff*b[i]*a[i]*a[i])
LRHANKELDEF(fk5q2n2l){
	return  (KAPPA5SUMFF(FORMK5Q2N2)) / (k*k*k0*k0);
}
LRHANKELDEF(fk5q2n2s){
	return  (KAPPA5SUM(FORMK5Q2N2)) / (k*k*k0*k0);
}

#if 0
complex double fk5q3n0l(double c, double k0, double k, 
		const complex double *a, const complex double *b, const complex double *d, const complex double *e, const complex double *ash) { // FIXME	
	return  ( /* FIXME */
		  -     k*b[0] + a[0] * ash[0]
		  + 5 * k*b[1] + a[1] * ash[1]
		  -10 * k*b[2] + a[2] * ash[2]
		  +10 * k*b[3] + a[3] * ash[3]
		  - 5 * k*b[4] + a[4] * ash[4]
		  +     k*b[5] + a[5] * ash[5]
		)/(k0*k0*k0);
}
#endif

#define FORMK5Q2N3(i,ff) (ff*a[i]*b[i]*(kk + 4*a[i]*a[i]))
LRHANKELDEF(fk5q2n3l){
	double kk = k*k; 
	return (KAPPA5SUMFF(FORMK5Q2N3))/(3*k0*k0*kk*k);
}
LRHANKELDEF(fk5q2n3s){
	double kk = k*k; 
	return (KAPPA5SUM(FORMK5Q2N3))/(3*k0*k0*kk*k);
}
#undef FORMK5Q2N3

#define FORMK5Q2N4(i,ff) (-ff*b[i]*a[i]*a[i]*(kk+2*a[i]*a[i]))
LRHANKELDEF(fk5q2n4l){
	double kk = k*k; 
	return (KAPPA5SUMFF(FORMK5Q2N4))/(k0*k0*kk*kk);
}
LRHANKELDEF(fk5q2n4s){
	double kk = k*k; 
	return (KAPPA5SUM(FORMK5Q2N4))/(k0*k0*kk*kk);
}
#undef FORMK5Q2N4

#define FORMK5Q2N5(i,ff) (ff*a[i]*b[i]*(kkkk+12*kk*(a[i]*a[i])+16*P4(a[i]))  )
LRHANKELDEF(fk5q2n5l){
	double kk = k*k; 
        double kkkk = kk * kk;	
	return (
		KAPPA5SUMFF(FORMK5Q2N5)
	        +16*120*P4(c)*c // Stirling S2(5,5) is no longer zero
	       )/(5*k0*k0*kkkk*k);
}
LRHANKELDEF(fk5q2n5s){
	double kk = k*k; 
        double kkkk = kk * kk;	
	return (
		KAPPA5SUM(FORMK5Q2N5)
	        +16*120*P4(c)*c // Stirling S2(5,5) is no longer zero
	       )/(5*k0*k0*kkkk*k);
}
#undef FORMK5Q2N5


static const lrhankelspec transfuns_f[MAXKAPPA+1][MAXQM+1][MAXN+1] = {
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{fk5q1n0l,fk5q1n1l,fk5q1n2l,fk5q1n3l,fk5q1n4l,fk5q1n5l},{fk5q2n0,fk5q2n1l,fk5q2n2l,fk5q2n3l,fk5q2n4l,fk5q2n5l}}
};

static const lrhankelspec transfuns_n[MAXKAPPA+1][MAXQM+1][MAXN+1] = {
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{NULL,NULL,NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL,NULL,NULL}},
	{{fk5q1n0s,fk5q1n1s,fk5q1n2s,fk5q1n3s,fk5q1n4s,fk5q1n5s},{fk5q2n0,fk5q2n1s,fk5q2n2s,fk5q2n3s,fk5q2n4s,fk5q2n5s}}
};

void lrhankel_recpart_fill(complex  double *target,
               size_t maxp /* max. order of the Hankel transform */, 
	       size_t lrk_cutoff,
               complex double const *const hct,
               unsigned kappa, double c, double k0, double k)
{
	assert(5 == kappa); // Only kappa == 5 implemented so far
	assert(maxp <= MAXN); // only n <= 5 implemented so far
	assert(lrk_cutoff <= MAXQM + 1); // only q <= 2 implemented so far; TODO shouldn't it be only MAXQM ???
	const lrhankelspec (*funarr)[MAXQM+1][MAXN+1] = (k>k0) ? transfuns_f : transfuns_n;
	memset(target, 0, maxp*(maxp+1)/2*sizeof(complex double));
	complex double a[kappa+1], b[kappa+1], d[kappa+1], e[kappa+1], ash[kappa+1];
	for (size_t sigma = 0; sigma <= kappa; ++sigma) {
		a[sigma] = (sigma * c - I * k0);
		b[sigma] = csqrt(1+k*k/(a[sigma]*a[sigma]));
		d[sigma] = 1/b[sigma];
		e[sigma] = d[sigma] / a[sigma];
		ash[sigma] = casinh(a[sigma]/k);
	}
	for (size_t ql = 0; (ql <= maxp) && (ql < lrk_cutoff); ++ql)  // ql is q-1, i.e. corresponds to the hankel term power
		for (size_t deltam = 0; deltam <= maxp; ++deltam){
			complex double result = funarr[kappa][ql][deltam](c,k0,k,a,b,d,e,ash);
			for (size_t p = 0; p <= maxp; ++p)
				trindex_cd(target,p)[deltam] += result * hankelcoeffs_get(hct,p)[ql];
		}
}

#ifdef TESTING
#include <stdio.h>
int main() {
	double k0 = 0.7;
	double c = 0.1324;
	double kmin = 0.000;
	double kmax = 20;
	double kstep = 0.001;
	size_t kappa = 5;

	for (double k = kmin; k <= kmax; k += kstep) {
		printf("%f ", k);
		complex double a[kappa+1], b[kappa+1], d[kappa+1], e[kappa+1], ash[kappa+1];
		for (size_t sigma = 0; sigma <= kappa; ++sigma) {
			a[sigma] = (sigma * c - I * k0);
			b[sigma] = csqrt(1+k*k/(a[sigma]*a[sigma]));
			d[sigma] = 1/b[sigma];
			e[sigma] = d[sigma] / a[sigma];
			ash[sigma] = casinh(a[sigma]/k);
		}
		for (size_t qm = 0; qm <= MAXQM; ++qm)
			for (size_t n = 0; n <= MAXN; ++n) {
			//if (/*!*/((qm==1)&&(n==0))){ //  not skip q==2, n=0 for now
// complex double fun(double c, double k0, double k, ccd *a, ccd *b, ccd *d, ccd *e)
				complex double result = 
					(k < k0 ? transfuns_n : transfuns_f)[kappa][qm][n](c,k0,k,a,b,d,e,ash);
				printf("%.16e %.16e ", creal(result), cimag(result));
			}
		printf("\n");
	}
	return 0;
}
#endif
