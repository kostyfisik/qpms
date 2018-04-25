#include "translations.h"
#include "vswf.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <assert.h>
#include "vectors.h"
#include "string.h"
#include "indexing.h"

#define MIN(x,y) (((x)<(y))?(x):(y))

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

int test_sphwave_translation(const qpms_trans_calculator *c, qpms_bessel_t wavetype, 
    cart3_t o2minuso1, int npoints, cart3_t *o1points);
//int test_planewave_decomposition(cart3_t k, ccart3_t E, qpms_l_t lMax, qpms_normalisation_t norm, int npoints, cart3_t *points);
//int test_planewave_decomposition_silent(cart3_t k, ccart3_t E, qpms_l_t lMax, qpms_normalisation_t norm, int npoints, cart3_t *points, double relerrthreshold, double *relerrs);

int main(int argc, char **argv) {
  gsl_rng *rng =  gsl_rng_alloc(gsl_rng_ranlxs0);
  gsl_rng_set(rng, 666);

  // parametry: prefix, n, m, norma, lmax, J
  // defaults
  qpms_m_t m1 = 0;
  qpms_l_t l1 = 1;
  qpms_normalisation_t norm = QPMS_NORMALISATION_NONE_CS;
  qpms_l_t lMax = 10;
  qpms_bessel_t J = 1;

  switch(MIN(argc,7)) {
    case 7: J = atoi(argv[6]);
    case 6: lMax = atoi(argv[5]);
    case 5: norm = atoi(argv[4]);
    case 4: m1 = atoi(argv[3]);
    case 3: l1 = atoi(argv[2]);
    case 2: break; // first argument is the filename prefix
    case 1: fputs("At least one argument needed (filename prefix)\n", stderr);
  }
  qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, norm);
  qpms_y_t nelem = c->nelem;
  qpms_y_t y1 = qpms_mn2y(m1, l1);

  //qpms_l_t viewlMax = 2;
  int npoints = 10;
  double sigma = 4.;
  //double shiftsigma = 0.01;

  cart3_t o2minuso1_max;
  // o2minuso1.x = gsl_ran_gaussian(rng, shiftsigma);
  // o2minuso1.y = gsl_ran_gaussian(rng, shiftsigma);
  // o2minuso1.z = gsl_ran_gaussian(rng, shiftsigma);
  o2minuso1_max.x = 6;
  o2minuso1_max.y = 0;
  o2minuso1_max.z = 0;
  int nraysteps = 512;
  cart3_t o2mo1_step = cart3_scale(1./nraysteps, o2minuso1_max);

  // measurement points
  cart3_t points[npoints]; 
  double relerrs[npoints];
  memset(points, 0, npoints * sizeof(cart3_t));
  points[0].x = points[1].y = points[2].z = sigma;
  points[3].x = points[3].y = 1./M_SQRT2;
  double relerrthreshold = 1e-11;
  for (unsigned i = 3; i < npoints; ++i) {
    cart3_t *w = points+i;
    w->x = gsl_ran_gaussian(rng, sigma);
    w->y = gsl_ran_gaussian(rng, sigma);
    w->z = gsl_ran_gaussian(rng, sigma);
  }
  // each point will have its file
  FILE *pfile[npoints];
  {
    char fname[strlen(argv[1])+20];
    for (int p = 0; p < npoints; ++p) {
      sprintf(fname, "%s.%d", argv[1], p);
      pfile[p] = fopen(fname, "w");
      fputs("##\tnu\tmu\t%d\tnorm\tlMax\tnelem\tJ\n", pfile[p]);
      fprintf(pfile[p], "#\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
          (int)l1, (int)m1, (int)y1, (int)norm, (int)lMax, (int)nelem, (int)J);
      fprintf(pfile[p], "#\tThis point cart:\t%g\t%g\t%g\n", 
          points[p].x, points[p].y, points[p].z);
      sph_t ps = cart2sph(points[p]);
      fprintf(pfile[p], "#\tThis point sph:\t%g\t%g\t%g\n",
          ps.r, ps.theta, ps.phi);
      csphvec_t M1[nelem], N1[nelem]; // mohl bych to počítat jednotlivě...
      if(QPMS_SUCCESS != qpms_vswf_fill(NULL, M1, N1, lMax, ps, J, norm)) abort();
      fprintf(pfile[p],"##\tM(%d,%d) at this point:\n", l1, m1);
      fputs("##\tRe(rc)\tIm(rc)\tRe(tc)\tIm(tc)\tRe(fc)\tIm(fc)\n",pfile[p]);
      fprintf(pfile[p],"#\t%g\t%g\t%g\t%g\t%g\t%g\n",
          creal(M1[y1].rc),cimag(M1[y1].rc),
          creal(M1[y1].thetac),cimag(M1[y1].thetac),
          creal(M1[y1].phic),cimag(M1[y1].phic)
      );
      fprintf(pfile[p],"##\tN(%d,%d) at this point:\n", l1, m1);
      fputs("##\tRe(rc)\tIm(rc)\tRe(tc)\tIm(tc)\tRe(fc)\tIm(fc)\n",pfile[p]);
      fprintf(pfile[p],"#\t%g\t%g\t%g\t%g\t%g\t%g\n",
          creal(N1[y1].rc),cimag(N1[y1].rc),
          creal(N1[y1].thetac),cimag(N1[y1].thetac),
          creal(N1[y1].phic),cimag(N1[y1].phic)
      );
      // TODO print column headers here
      fputs("## TODO print column headers here\n", pfile[p]);
      fputs("#step\t"
        "(o2-o1).x\t(o2-o1).y\t(o2-o1).z\t(o2-o1).r\t(o2-o1).θ\t(o2-o1).φ\t"
        "(x-o2).x\t(x-o2).y\t(x-o2).z\t(x-o2).r\t(x-o2).θ\t(x-o2).φ\t",
        pfile[p]
      );

      fprintf(pfile[p], //original and reconstructed waves at new origin
          "R(M1@2(%d,%d).r)\tI(M1@2(%d,%d).r)\t"
          "R(M1@2(%d,%d).θ)\tI(M1@2(%d,%d).θ)\t"
          "R(M1@2(%d,%d).φ)\tI(M1@2(%d,%d).φ)\t"
          "R(N1@2(%d,%d).r)\tI(N1@2(%d,%d).r)\t"
          "R(N1@2(%d,%d).θ)\tI(N1@2(%d,%d).θ)\t"
          "R(N1@2(%d,%d).φ)\tI(N1@2(%d,%d).φ)\t"
          "R(Mr1@2(%d,%d).r)\tI(Mr1@2(%d,%d).r)\t"
          "R(Mr1@2(%d,%d).θ)\tI(Mr1@2(%d,%d).θ)\t"
          "R(Mr1@2(%d,%d).φ)\tI(Mr1@2(%d,%d).φ)\t"
          "R(Nr1@2(%d,%d).r)\tI(Nr1@2(%d,%d).r)\t"
          "R(Nr1@2(%d,%d).θ)\tI(Nr1@2(%d,%d).θ)\t"
          "R(Nr1@2(%d,%d).φ)\tI(Nr1@2(%d,%d).φ)\n",
          l1,m1,l1,m1, l1,m1,l1,m1, l1,m1,l1,m1,
          l1,m1,l1,m1, l1,m1,l1,m1, l1,m1,l1,m1,
          l1,m1,l1,m1, l1,m1,l1,m1, l1,m1,l1,m1,
          l1,m1,l1,m1, l1,m1,l1,m1, l1,m1,l1,m1
      );

      for(qpms_y_t y = 0; y < nelem; ++y) {
        qpms_l_t l2; qpms_m_t m2;
        qpms_y2mn_p(y, &m2, &l2);
        fprintf(pfile[p], "R(A(%d,%d))\tI(A(%d,%d))\t"
          "R(B(%d,%d))\tI(B(%d,%d))\t"
          "R(M2(%d,%d).r)\tI(M2(%d,%d).r)\t"
          "R(M2(%d,%d).θ)\tI(M2(%d,%d).θ)\t"
          "R(M2(%d,%d).φ)\tI(M2(%d,%d).φ)\t"
          "R(N2(%d,%d).r)\tI(N2(%d,%d).r)\t"
          "R(N2(%d,%d).θ)\tI(N2(%d,%d).θ)\t"
          "R(N2(%d,%d).φ)\tI(N2(%d,%d).φ)\t",
          l2,m2,l2,m2, l2,m2,l2,m2, l2,m2,l2,m2, l2,m2,l2,m2, l2,m2,l2,m2, 
          l2,m2,l2,m2, l2,m2,l2,m2, l2,m2,l2,m2
        );
      }
    }
  }

  for (int r = 0; r <= nraysteps; ++r) {
    cart3_t sc = cart3_scale(r, o2mo1_step);
    sph_t ss = cart2sph(sc);
    csphvec_t N1[nelem], M1[nelem];
    for(int p = 0; p < npoints; p++){
      FILE *f = pfile[p];
      cart3_t w1c = points[p];
      cart3_t w2c = cart3_add(w1c, cart3_scale(-1.,sc));
      sph_t w1s = cart2sph(w1c);
      sph_t w2s = cart2sph(w2c);
      fprintf(f, "%d\t" "%g\t%g\t%g\t%g\t%g\t%g\t" "%g\t%g\t%g\t%g\t%g\t%g\t", 
          r, sc.x, sc.y, sc.z, ss.r, ss.theta, ss.phi,
          w2c.x, w2c.y, w2c.z, w2s.r, w2s.theta, w2s.phi);
      complex double A[nelem], B[nelem]; // translation ceofficients
      csphvec_t M2at2[nelem], N2at2[nelem], M1at2, N1at2, M1Rat2, N1Rat2;
      if(QPMS_SUCCESS != qpms_vswf_fill(NULL, M2at2, N2at2, lMax, w2s, J, norm)) abort();
      for(qpms_y_t y2 = 0; y2 < nelem; ++y2) {
        qpms_l_t l2; qpms_m_t m2;
        qpms_y2mn_p(y2, &m2, &l2);
        if (QPMS_SUCCESS != qpms_trans_calculator_get_AB_p(c, &(A[y2]), &(B[y2]),
            m2, l2, m1, l1, // !!! FIXME mám správné pořadí??? !!!
            ss, true /* FIXME Pro J != 1 */,  J)) abort();
      }

      M1at2 = ccart2csphvec(csphvec2ccart(M1[y1], w1s), w2s);
      N1at2 = ccart2csphvec(csphvec2ccart(N1[y1], w1s), w2s);
      fprintf(f, "%g\t%g\t%g\t%g\t%g\t%g\t" "%g\t%g\t%g\t%g\t%g\t%g\t", 
        creal(M1at2.rc), cimag(M1at2.rc),
        creal(M1at2.thetac), cimag(M1at2.thetac),
        creal(M1at2.phic), cimag(M1at2.phic),
        creal(N1at2.rc), cimag(N1at2.rc),
        creal(N1at2.thetac), cimag(N1at2.thetac),
        creal(N1at2.phic), cimag(N1at2.phic)
      );
      M1Rat2 = qpms_eval_vswf(w2s, NULL, A, B, lMax, J, norm);
      N1Rat2 = qpms_eval_vswf(w2s, NULL, B, A, lMax, J, norm);
      fprintf(f, "%g\t%g\t%g\t%g\t%g\t%g\t" "%g\t%g\t%g\t%g\t%g\t%g\t", 
        creal(M1Rat2.rc), cimag(M1Rat2.rc),
        creal(M1Rat2.thetac), cimag(M1Rat2.thetac),
        creal(M1Rat2.phic), cimag(M1Rat2.phic),
        creal(N1Rat2.rc), cimag(N1Rat2.rc),
        creal(N1Rat2.thetac), cimag(N1Rat2.thetac),
        creal(N1Rat2.phic), cimag(N1Rat2.phic)
      );

      for(qpms_y_t y2 = 0; y2 < nelem; ++y2){
        fprintf(f, "%g\t%g\t" "%g\t%g\t" "%g\t%g\t%g\t%g\t%g\t%g\t"
          "%g\t%g\t%g\t%g\t%g\t%g\t", 
          creal(A[y2]), cimag(A[y2]),
          creal(B[y2]), cimag(B[y2]), 
          creal(M2at2[y2].rc), cimag(M2at2[y2].rc),
          creal(M2at2[y2].thetac), cimag(M2at2[y2].thetac),
          creal(M2at2[y2].phic), cimag(M2at2[y2].phic),
          creal(N2at2[y2].rc), cimag(N2at2[y2].rc),
          creal(N2at2[y2].thetac), cimag(N2at2[y2].thetac),
          creal(N2at2[y2].phic), cimag(N2at2[y2].phic)
        );
      }
      fputc('\n', f);
    }
  }

  for (int p = 0; p < npoints; ++p)
    fclose(pfile[p]);
  gsl_rng_free(rng);
  return 0;
}

#if 0
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
    // cart3_t w2c = cart3_add(w1c, cart3_scale(-1, sc));
    cart3_t w2c = cart3_add(w1c, sc);
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

    for(qpms_y_t y1 = 0; y1 < nelem; ++y1) { //index of the wave originating in o1 that will be reconstructed in o2
      qpms_m_t m1;
      qpms_l_t l1;
      qpms_y2mn_p(y1, &m1, &l1);
      printf("*** wave l = %d, m = %d ***\n", l1, m1);

      complex double A[nelem], B[nelem];
      for(qpms_y_t y2 = 0; y2 < nelem; ++y2){
        qpms_m_t m2; qpms_l_t l2;
        qpms_y2mn_p(y2, &m2, &l2);
        if(qpms_trans_calculator_get_AB_p(c, A+y2, B+y2, m2, l2, m1, l1, ss, (w1s.r > ss.r), wavetype))
          abort();
      }

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
#endif

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
