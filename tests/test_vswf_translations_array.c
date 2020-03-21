//c99 -o test_vswf_translations_array -ggdb -I .. test_vswf_translations_array.c -lqpms -lgsl -lm -lblas 
#include <qpms/translations.h>
#include <qpms/vswf.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <assert.h>
#include <qpms/vectors.h>
#include <string.h>
#include <qpms/indexing.h>

char *normstr(qpms_normalisation_t norm) {
	//int csphase = qpms_normalisation_t_csphase(norm);
	norm = norm & QPMS_NORMALISATION_NORM_BITS;
	switch (norm) {
		case QPMS_NORMALISATION_NORM_NONE:
			return "none";
		case QPMS_NORMALISATION_NORM_SPHARM:
			return "spharm";
		case QPMS_NORMALISATION_NORM_POWER:
			return "power";
		default:
			return "!!!undef!!!";
	}
}

int test_sphwave_translation(const qpms_trans_calculator *c, qpms_bessel_t wavetype, 
		cart3_t o2minuso1, int npoints, cart3_t *o1points);
//int test_planewave_decomposition(cart3_t k, ccart3_t E, qpms_l_t lMax, qpms_normalisation_t norm, int npoints, cart3_t *points);
//int test_planewave_decomposition_silent(cart3_t k, ccart3_t E, qpms_l_t lMax, qpms_normalisation_t norm, int npoints, cart3_t *points, double relerrthreshold, double *relerrs);

int main() {
	gsl_rng *rng =  gsl_rng_alloc(gsl_rng_ranlxs0);
	gsl_rng_set(rng, 666);

	qpms_l_t lMax = 8;
	//qpms_l_t viewlMax = 2;
	int npoints = 10;
	double sigma = 4;
  //double shiftsigma = 2.;

	cart3_t o2minuso1;
	o2minuso1.x = 1; //gsl_ran_gaussian(rng, shiftsigma);
	o2minuso1.y = 2; //gsl_ran_gaussian(rng, shiftsigma);
	o2minuso1.z = 5; //gsl_ran_gaussian(rng, shiftsigma);
	
	cart3_t points[npoints];
	double relerrs[npoints];
	memset(points, 0, npoints * sizeof(cart3_t));
	points[0].x = points[1].y = points[2].z = sigma;
	double relerrthreshold = 1e-11;
	for (unsigned i = 3; i < npoints; ++i) {
		cart3_t *w = points+i;
		w->x = gsl_ran_gaussian(rng, sigma);
		w->y = gsl_ran_gaussian(rng, sigma);
		w->z = gsl_ran_gaussian(rng, sigma);
	}

  for(int use_csbit = 0; use_csbit <= 1; ++use_csbit) {
    for(int i = 0; i < 3; ++i){
      qpms_normalisation_t norm = ((qpms_normalisation_t[])
          { QPMS_NORMALISATION_NORM_SPHARM,
            QPMS_NORMALISATION_NORM_POWER,
            QPMS_NORMALISATION_NORM_NONE
            })[i]
        | (use_csbit ? QPMS_NORMALISATION_CSPHASE : 0);
      qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, norm);
      for(int J = 1; J <= 4; ++J)
        test_sphwave_translation(c, J, o2minuso1, npoints, points);
      qpms_trans_calculator_free(c);
    }

  }	
	gsl_rng_free(rng);
}

int test_sphwave_translation(const qpms_trans_calculator *c, qpms_bessel_t wavetype, 
		cart3_t sc, int npoints, cart3_t *points) {
	puts("==============================================================");
	printf("Test translation o2-o1 = %fx̂ + %fŷ + %fẑ", sc.x, sc.y, sc.z);
	sph_t ss = cart2sph(sc);
	printf(" = %fr̂ @ o1, θ = %f, φ = %f\n", ss.r, ss.theta, ss.phi);
	printf("lMax = %d, norm: %s, csphase = %d\n", 
			(int)c->lMax, normstr(c->normalisation), qpms_normalisation_t_csphase(c->normalisation));
	printf("wave type J = %d\n", wavetype);
	
	qpms_l_t lMax = c->lMax;
	qpms_y_t nelem = c->nelem;
	csphvec_t N1[nelem], /* N2[nelem], */ M1[nelem] /*, M2[nelem]*/;

	for (int i = 0; i < npoints; i++) {
		printf("-------- Point %d --------\n", i);
		cart3_t w1c = points[i];
		cart3_t w2c = cart3_add(w1c, cart3_scale(-1, sc));
		sph_t w1s = cart2sph(w1c);
		sph_t w2s = cart2sph(w2c);
		printf(" = %fx̂ + %fŷ + %fẑ @o1\n", w1c.x, w1c.y, w1c.z);
		printf(" = %fx̂ + %fŷ + %fẑ @o2\n", w2c.x, w2c.y, w2c.z);
		printf(" = %fr̂ @ o1, θ = %f, φ = %f\n", w1s.r, w1s.theta, w1s.phi);
		printf(" = %fr̂ @ o2, θ = %f, φ = %f\n", w2s.r, w2s.theta, w2s.phi);
		printf("Outside the sphere centered in o1 intersecting o2: %s; by %f\n", (w1s.r > ss.r) ? "true" : "false", 
				w1s.r - ss.r);
		if(QPMS_SUCCESS != qpms_vswf_fill(NULL, M1, N1, lMax, w1s, wavetype, c->normalisation))
			abort(); // original wave set

    complex double A_whole[nelem][nelem], B_whole[nelem][nelem];
    if (qpms_trans_calculator_get_AB_arrays(c,(complex double *) A_whole, (complex double *) B_whole,
          1, nelem, sph2csph(ss), (w1s.r > ss.r), wavetype)) abort();


		for(qpms_y_t y1 = 0; y1 < nelem; ++y1) { //index of the wave originating in o1 that will be reconstructed in o2
			qpms_m_t m1;
			qpms_l_t l1;
			qpms_y2mn_p(y1, &m1, &l1);
			printf("*** wave l = %d, m = %d ***\n", l1, m1);

			//complex double A[nelem], B[nelem];
      /*
			for(qpms_y_t y2 = 0; y2 < nelem; ++y2){
				qpms_m_t m2; qpms_l_t l2;
				qpms_y2mn_p(y2, &m2, &l2);
				if(qpms_trans_calculator_get_AB_p(c, &(A[y2]), &(B[y2]), m2, l2, m1, l1, ss, (w1s.r > ss.r) , wavetype))
					abort();
			}
      */
      // array function instead above
      complex double *A = A_whole[y1]; // CHECKME right index ??
      complex double *B = B_whole[y1];


			printf("M = ");
			print_csphvec(M1[y1]);
			printf(" @ o1\n  = ");
			ccart3_t M1c = csphvec2ccart(M1[y1], w1s);
			print_ccart3(M1c);
			printf("\n  = ");
			csphvec_t M1s2 = ccart2csphvec(M1c, w2s);
			print_csphvec(M1s2);
			printf(" @ o2\n");
			csphvec_t M2s2 = qpms_eval_vswf(w2s, NULL, A, B, lMax, wavetype, c->normalisation);
			printf("Mr= ");
			print_csphvec(M2s2);
			printf(" @ o2\n");

			printf("N = ");
			print_csphvec(N1[y1]);
			printf(" @ o1\n  = ");
			ccart3_t N1c = csphvec2ccart(N1[y1], w1s);
			print_ccart3(N1c);
			printf("\n  = ");
			csphvec_t N1s2 = ccart2csphvec(N1c, w2s);
			print_csphvec(N1s2);
			printf(" @o2\nNr= ");
			csphvec_t N2s2 = qpms_eval_vswf(w2s, NULL, B, A, lMax, wavetype, c->normalisation);
			print_csphvec(N2s2);
			printf(" @o2\n");
		}
	}

	return 0; // FIXME something more meaningful here...
}




	
#if 0

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
#endif
