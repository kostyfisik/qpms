// c99 -o test_planewave_decomposition -ggdb -I .. test_planewave_decomposition.c ../vswf.c -lgsl -lm -lblas ../legendre.c ../bessel.c
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <assert.h>
#include "vswf.h"
#include "vectors.h"
#include "string.h"
#include "indexing.h"

char *normstr(qpms_normalisation_t norm) {
	//int csphase = qpms_normalisation_t_csphase(norm);
	norm = qpms_normalisation_t_normonly(norm);
	switch (norm) {
		case QPMS_NORMALISATION_NONE:
			return "none";
		case QPMS_NORMALISATION_SPHARM:
			return "spharm";
		case QPMS_NORMALISATION_POWER:
			return "power";
		default:
			return "!!!undef!!!";
	}
}
int test_planewave_decomposition(cart3_t k, ccart3_t E, qpms_l_t lMax, qpms_normalisation_t norm, int npoints, cart3_t *points);
int test_planewave_decomposition_silent(cart3_t k, ccart3_t E, qpms_l_t lMax, qpms_normalisation_t norm, int npoints, cart3_t *points, double relerrthreshold, double *relerrs);

int main() {
	gsl_rng *rng =  gsl_rng_alloc(gsl_rng_ranlxs0);
	gsl_rng_set(rng, 666);

	qpms_l_t lMax = 50;
	int npoints = 10;
	double sigma = 1.3;
	cart3_t points[npoints];
	double relerrs[npoints];
	memset(points, 0, npoints * sizeof(cart3_t));
	points[1].x = points[2].y = points[3].z = 1.;
	double relerrthreshold = 1e-11;
	for (unsigned i = 4; i < npoints; ++i) {
		cart3_t *w = points+i;
		w->x = gsl_ran_gaussian(rng, sigma);
		w->y = gsl_ran_gaussian(rng, sigma);
		w->z = gsl_ran_gaussian(rng, sigma);
	}
	
	double ksigma = 15;

	for(int i = 1; i < 2; ++i) {
		cart3_t k = {0, 0, 2};
		ccart3_t E = {0., 1.1+I, 0.};
		if(i){
			k.x = gsl_ran_gaussian(rng, ksigma);
			k.y = gsl_ran_gaussian(rng, ksigma);
			k.z = gsl_ran_gaussian(rng, ksigma);
			E.x = gsl_ran_gaussian(rng, ksigma);
			E.x += I* gsl_ran_gaussian(rng, ksigma);	
			E.y = gsl_ran_gaussian(rng, ksigma);
			E.y += I* gsl_ran_gaussian(rng, ksigma);
			E.z = gsl_ran_gaussian(rng, ksigma);
			E.z += I* gsl_ran_gaussian(rng, ksigma);
		}
#if 0
		printf("Test wave k = %gx̂ + %gŷ + %gẑ", k.x, k.y, k.z);
		printf(", E_0 = (%g+%gj)x̂ + (%g+%gj)ŷ + (%g+%gj)ẑ\n",
			creal(E.x),cimag(E.x),
			creal(E.y),cimag(E.y),
			creal(E.z),cimag(E.z));
		printf("%d points, sigma = %g, rel. err. threshold = %g\n", npoints, sigma, relerrthreshold);
		for(int i = 1; i <= 3; ++i){
			int res = test_planewave_decomposition_silent(k, E, lMax, i, npoints, points, relerrthreshold, relerrs);
			printf("\n%s %d/%d %s %s\n", res ? "!!" : "OK", npoints-res, npoints, normstr(i), "");
			for(int p = 0; p < npoints; ++p) printf("%.3g ", relerrs[p]);
			res = test_planewave_decomposition_silent(k, E, lMax, i |QPMS_NORMALISATION_T_CSBIT, npoints, points, relerrthreshold, relerrs);
			printf("\n%s %d/%d %s %s\n", res ? "!!" : "OK", npoints-res, npoints, normstr(i), "with C.-S. phase");
			for(int p = 0; p < npoints; ++p) printf("%.3g ", relerrs[p]);
		}
#endif 
		for(int i = 1; i <= 3; ++i){
			test_planewave_decomposition(k, E, lMax, i, npoints, points);
			test_planewave_decomposition(k, E, lMax, i | QPMS_NORMALISATION_T_CSBIT, npoints, points);
		}

	}
	gsl_rng_free(rng);
}

int test_planewave_decomposition(cart3_t k, ccart3_t E, qpms_l_t lMax, qpms_normalisation_t norm, int npoints, cart3_t *points){
	qpms_y_t nelem = qpms_lMax2nelem(lMax);
	complex double lc[nelem], mc[nelem], ec[nelem];
	if (QPMS_SUCCESS != 
		qpms_planewave2vswf_fill_cart(k, E, lc, mc, ec, lMax, norm)) {
		printf("Error\n");
		return -1;
	}
	printf("==============================================================\n");
	printf("Test wave k = %fx̂ + %fŷ + %fẑ", k.x, k.y, k.z);
	sph_t k_sph = cart2sph(k);
	printf(" = %fr̂ @ θ = %f, φ = %f\n", k_sph.r, k_sph.theta, k_sph.phi);
	printf("        E_0 = (%f+%fj)x̂ + (%f+%fj)ŷ + (%f+%fj)ẑ", 
		creal(E.x),cimag(E.x),
		creal(E.y),cimag(E.y),
		creal(E.z),cimag(E.z));
	csphvec_t E_s = ccart2csphvec(E, k_sph);
	printf(" = (%f+%fj)r̂ + (%f+%fj)θ̂ + (%f+%fj)φ̂  @ k\n",
		creal(E_s.rc), cimag(E_s.rc), creal(E_s.thetac), cimag(E_s.thetac),
		creal(E_s.phic), cimag(E_s.phic));
	printf("        lMax = %d, norm: %s, csphase = %d\n",
			(int)lMax, normstr(norm), qpms_normalisation_t_csphase(norm));
	printf("a_L: ");
	for(qpms_y_t y = 0; y < nelem; ++y) printf("%g+%gj ", creal(lc[y]), cimag(lc[y]));
	printf("\na_M: ");
	for(qpms_y_t y = 0; y < nelem; ++y) printf("%g+%gj ", creal(mc[y]), cimag(mc[y]));
	printf("\na_N: ");
	for(qpms_y_t y = 0; y < nelem; ++y) printf("%g+%gj ", creal(ec[y]), cimag(ec[y]));
	printf("\n");
	for (int i = 0; i < npoints; i++) {
		cart3_t w = points[i];
		sph_t w_sph = cart2sph(w);
		printf("Point %d: x = %f, y = %f, z = %f,\n", i, w.x, w.y, w.z);
		printf("        |r| = %f, θ = %f, φ = %f:\n", w_sph.r, w_sph.theta, w_sph.phi);
		double phase = cart3_dot(k,w);
		printf("  k.r = %f\n", phase);
		complex double phfac = cexp(phase * I);
		ccart3_t Ew = ccart3_scale(phfac, E);
		printf("  pw  E(r) = (%f+%fj)x̂ + (%f+%fj)ŷ + (%f+%fj)ẑ", 
			creal(Ew.x),cimag(Ew.x),
			creal(Ew.y),cimag(Ew.y),
			creal(Ew.z),cimag(Ew.z));
		csphvec_t Ew_s = ccart2csphvec(Ew, w_sph);
		printf(" = (%f+%fj)r̂ + (%f+%fj)θ̂ + (%f+%fj)φ̂  @ r\n",
			creal(Ew_s.rc), cimag(Ew_s.rc), 
			creal(Ew_s.thetac), cimag(Ew_s.thetac),
			creal(Ew_s.phic), cimag(Ew_s.phic));
		w_sph.r *= k_sph.r; /// NEVER FORGET THIS!!!
		csphvec_t Ew_s_recomp = qpms_eval_vswf(w_sph,
			lc, mc, ec, lMax, QPMS_BESSEL_REGULAR, norm);
		ccart3_t Ew_recomp = csphvec2ccart(Ew_s_recomp, w_sph);
		printf(" rec E(r) = (%f+%fj)x̂ + (%f+%fj)ŷ + (%f+%fj)ẑ", 
			creal(Ew_recomp.x),cimag(Ew_recomp.x),
			creal(Ew_recomp.y),cimag(Ew_recomp.y),
			creal(Ew_recomp.z),cimag(Ew_recomp.z));
		printf(" = (%f+%fj)r̂ + (%f+%fj)θ̂ + (%f+%fj)φ̂  @ r\n",
			creal(Ew_s_recomp.rc), cimag(Ew_s_recomp.rc), 
			creal(Ew_s_recomp.thetac), cimag(Ew_s_recomp.thetac),
			creal(Ew_s_recomp.phic), cimag(Ew_s_recomp.phic));
		double relerrfac = 2/(cabs(Ew_s_recomp.rc) + cabs(Ew_s.rc)
				+cabs(Ew_s_recomp.thetac) + cabs(Ew_s.thetac)
				+cabs(Ew_s_recomp.phic) + cabs(Ew_s.phic));
		printf(" rel. err. magnitude: %g @ r̂, %g @ θ̂, %g @ φ̂\n",
			cabs(Ew_s_recomp.rc - Ew_s.rc) * relerrfac,
			cabs(Ew_s_recomp.thetac - Ew_s.thetac) * relerrfac,
			cabs(Ew_s_recomp.phic - Ew_s.phic) * relerrfac
		      );
	}
	return 0;
}

int test_planewave_decomposition_silent(cart3_t k, ccart3_t E, qpms_l_t lMax, qpms_normalisation_t norm, int npoints, cart3_t *points, double relerrthreshold, double *relerrs) {
	qpms_y_t nelem = qpms_lMax2nelem(lMax);
	int failcount = 0;
	complex double lc[nelem], mc[nelem], ec[nelem];
	if (QPMS_SUCCESS != 
		qpms_planewave2vswf_fill_cart(k, E, lc, mc, ec, lMax, norm)) {
		printf("Error\n");
		return -1;
	}
	sph_t k_sph = cart2sph(k);
	csphvec_t E_s = ccart2csphvec(E, k_sph);
	for (int i = 0; i < npoints; i++) {
		cart3_t w = points[i];
		sph_t w_sph = cart2sph(w);
		w_sph.r *= k_sph.r;
		double phase = cart3_dot(k,w);
		complex double phfac = cexp(phase * I);
		ccart3_t Ew = ccart3_scale(phfac, E);
		csphvec_t Ew_s = ccart2csphvec(Ew, w_sph);
		csphvec_t Ew_s_recomp = qpms_eval_vswf(w_sph,
			lc, mc, ec, lMax, QPMS_BESSEL_REGULAR, norm);
		ccart3_t Ew_recomp = csphvec2ccart(Ew_s_recomp, w_sph);
		double relerrfac = 2/(cabs(Ew_s_recomp.rc) + cabs(Ew_s.rc)
				+cabs(Ew_s_recomp.thetac) + cabs(Ew_s.thetac)
				+cabs(Ew_s_recomp.phic) + cabs(Ew_s.phic));
			
		double relerr = (cabs(Ew_s_recomp.rc - Ew_s.rc) 
			+	cabs(Ew_s_recomp.thetac - Ew_s.thetac)
			+ 	cabs(Ew_s_recomp.phic - Ew_s.phic)
			) * relerrfac;
		if(relerrs) relerrs[i] = relerr;
		if(relerr > relerrthreshold) ++failcount;
	}
	return failcount;
}
