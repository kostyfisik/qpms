// c99 -o ../../tests/raytests/n_nu0 -ggdb -I .. vswf_translation_test_rays.c ../translations.c ../vswf.c ../gaunt.c ../legendre.c -lgsl -lm -lblas

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
  int npoints = 3;
  double sigma = 4.;
  //double shiftsigma = 0.01;

  cart3_t o2minuso1_max;
  // o2minuso1.x = gsl_ran_gaussian(rng, shiftsigma);
  // o2minuso1.y = gsl_ran_gaussian(rng, shiftsigma);
  // o2minuso1.z = gsl_ran_gaussian(rng, shiftsigma);
  o2minuso1_max.x = 1;
  o2minuso1_max.y = 2;
  o2minuso1_max.z = 5;
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
          "R(Nr1@2(%d,%d).φ)\tI(Nr1@2(%d,%d).φ)\t",
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
      fputc('\n', pfile[p]);
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
      if(QPMS_SUCCESS != qpms_vswf_fill(NULL, M1, N1, lMax, w1s, J, norm)) abort();
      csphvec_t M2at2[nelem], N2at2[nelem], M1at2, N1at2, M1Rat2, N1Rat2;
      if(QPMS_SUCCESS != qpms_vswf_fill(NULL, M2at2, N2at2, lMax, w2s, J, norm)) abort();
      for(qpms_y_t y2 = 0; y2 < nelem; ++y2) {
        qpms_l_t l2; qpms_m_t m2;
        qpms_y2mn_p(y2, &m2, &l2);
        if (QPMS_SUCCESS != qpms_trans_calculator_get_AB_p(c, &(A[y2]), &(B[y2]),
#ifdef REVERSE
              m1, l1, m2, l2,
#else
              m2, l2, m1, l1, // !!! FIXME mám správné pořadí??? !!!
#endif
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

