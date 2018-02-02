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

int main() {
	gsl_rng *rng =  gsl_rng_alloc(gsl_rng_ranlxs0);
	gsl_rng_set(rng, 666);

	qpms_l_t lMax = 20;
	int npoints = 10;
	double sigma = 1.3;
	cart3_t points[npoints];
	memset(points, 0, npoints * sizeof(cart3_t));
	points[0].x = points[1].y = points[2].z = 1.;
	for (unsigned i = 3; i < npoints; ++i) {
		cart3_t *w = points+i;
		w->x = gsl_ran_gaussian(rng, sigma);
		w->y = gsl_ran_gaussian(rng, sigma);
		w->z = gsl_ran_gaussian(rng, sigma);
	}
	

	cart3_t k = {0, 0, 1};
	ccart3_t E = {0., 1., 0.};
	for(int i = 1; i <= 3; ++i){
		test_planewave_decomposition(k, E, lMax, i, npoints, points);
		test_planewave_decomposition(k, E, lMax, i | QPMS_NORMALISATION_T_CSBIT, npoints, points);
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
		double ph_2pi = cart3_dot(k,w);
		printf("  k.r = %f\n", ph_2pi);
		double phfac = cexp(2 * M_PI * ph_2pi * I);
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
		printf(" rel. err. magnitude: %e @ r̂, %e @ θ̂, %e @ φ̂\n",
			2. * cabs(Ew_s_recomp.rc - Ew_s.rc)
				/(cabs(Ew_s_recomp.rc) + cabs(Ew_s.rc)),
			2. * cabs(Ew_s_recomp.thetac - Ew_s.thetac)
				/(cabs(Ew_s_recomp.thetac) + cabs(Ew_s.thetac)),
			2. * cabs(Ew_s_recomp.phic - Ew_s.phic)
				/(cabs(Ew_s_recomp.phic) + cabs(Ew_s.phic)));
	}
	return 0;
}
