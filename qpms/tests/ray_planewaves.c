// c99 -o raypwtest1 -I .. ray_planewaves.c -lgsl -lm ../legendre.c ../vswf.c ../bessel.c -lblas
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
#define MIN(x,y) (((x)<(y))?(x):(y))

int main(int argc, char **argv) {
	gsl_rng *rng =  gsl_rng_alloc(gsl_rng_ranlxs0);
	gsl_rng_set(rng, 666);


  int nraysteps = 512;


  // parametry: prefix, kx, ky, kz, E0x, E0y, E0z, norma, lmax
  // defaults
  double kx = 1, ky = 0, kz = 0, E0x = 0, E0y = 0, E0z = 1;
	qpms_l_t lMax = 10;
  qpms_normalisation_t norm = QPMS_NORMALISATION_NONE_CS;

  switch(MIN(argc, 10)) {
    case 10: lMax = atoi(argv[9]);
    case 9: norm = atoi(argv[8]);
    case 8: E0z = strtod(argv[7], NULL);
    case 7: E0y = strtod(argv[6], NULL);
    case 6: E0x = strtod(argv[5], NULL);
    case 5: kz = strtod(argv[4], NULL);
    case 4: ky = strtod(argv[3], NULL);
    case 3: kx = strtod(argv[2], NULL);
    case 2: break;
    case 1: fputs("At least one argument needed (filename prefix)\n", stderr);
  }

  qpms_y_t nelem = qpms_lMax2nelem(lMax);
	int nrays = 10;
	double sigma = 4;
	cart3_t rays[nrays];
	double relerrs[nrays];
	memset(rays, 0, nrays * sizeof(cart3_t));
	rays[1].x = rays[2].y = rays[3].z = sigma;
	double relerrthreshold = 1e-11;
	for (unsigned i = 4; i < nrays; ++i) {
		cart3_t *w = rays+i;
		w->x = gsl_ran_gaussian(rng, sigma);
		w->y = gsl_ran_gaussian(rng, sigma);
		w->z = gsl_ran_gaussian(rng, sigma);
	}
	


  cart3_t k = {kx, ky, kz};
  k = cart3_scale(1./cart3norm(k),k); // automatically normalise

	ccart3_t E0 = {E0x, E0y, E0z};
  sph_t ks = cart2sph(k);
  csphvec_t E0s = ccart2csphvec(E0, ks);
  complex double Lc[nelem], Mc[nelem], Nc[nelem];
  if (qpms_planewave2vswf_fill_sph(ks, E0s, Lc, Mc, Nc, lMax, norm)) abort();

  
  // each ray will have its file
  FILE *rfile[nrays];
  {
    char fname[strlen(argv[1])+20];
    for (int p = 0; p < nrays; ++p) {
      sprintf(fname, "%s.%d", argv[1], p);
      rfile[p] = fopen(fname, "w");
      fputs("##\tnorm\tlMax\tnelem\n", rfile[p]);
      fprintf(rfile[p], "#\t%d\t%d\t%d\n",
          (int)norm, (int)lMax, (int)nelem);
      fputs("##\tkx\tky\tkz\tR(E0x)\tI(E0x)\tR(E0y)\tI(E0y)\tR(E0z)\tI(E0z)\n", rfile[p]);
      fprintf(rfile[p], "#\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
          k.x, k.y, k.z,
          creal(E0.x), cimag(E0.x),
          creal(E0.y), cimag(E0.y),
          creal(E0.z), cimag(E0.z)
      );
      sph_t ks = cart2sph(k);
      csphvec_t E0s = ccart2csphvec(E0, ks);
      fputs("##\tkr\tkθ\tkφ\tR(E0r)\tI(E0r)\tR(E0θ)\tI(E0θ)\tR(E0φ)\tI(E0φ)\n", rfile[p]);
      fprintf(rfile[p], "#\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
          ks.r, ks.theta, ks.phi,
          creal(E0s.rc), cimag(E0s.rc),
          creal(E0s.thetac), cimag(E0s.thetac),
          creal(E0s.phic), cimag(E0s.phic)
      );
      fprintf(rfile[p], "#\tThis ray cart:\t%g\t%g\t%g\n",
          rays[p].x, rays[p].y, rays[p].z);
      sph_t ps = cart2sph(rays[p]);
      fprintf(rfile[p], "#\tThis ray sph:\t%g\t%g\t%g\n",
          ps.r, ps.theta, ps.phi);
      fputs("## TODO Coefficients of the wave decomposition\n", rfile[p]); // TODO print Lc, Mc, Nc
      fputs("#step\t"
        "x\ty\tz\tr\tθ\tφ\tphase\t"
        "R(Ex)\tI(Ex)\tR(Ey)\tI(Ey)\tR(Ez)\tI(Ez)\t"
        "R(Er)\tI(Er)\tR(Eθ)\tI(Eθ)\tR(Eφ)\tI(Eφ)\t"
        "R(rEx)\tI(rEx)\tR(rEy)\tI(rEy)\tR(rEz)\tI(rEz)\t"
        "R(rEr)\tI(rEr)\tR(rEθ)\tI(rEθ)\tR(rEφ)\tI(rEφ)\t",
        rfile[p]
      );


      for(qpms_y_t y = 0; y < nelem; ++y) {
        qpms_l_t l2; qpms_m_t m2;
        qpms_y2mn_p(y, &m2, &l2);
        fprintf(rfile[p], // TODO print also the longitudinal waves? 
          "R(M2(%d,%d).r)\tI(M2(%d,%d).r)\t"
          "R(M2(%d,%d).θ)\tI(M2(%d,%d).θ)\t"
          "R(M2(%d,%d).φ)\tI(M2(%d,%d).φ)\t"
          "R(N2(%d,%d).r)\tI(N2(%d,%d).r)\t"
          "R(N2(%d,%d).θ)\tI(N2(%d,%d).θ)\t"
          "R(N2(%d,%d).φ)\tI(N2(%d,%d).φ)\t",
          l2,m2,l2,m2, l2,m2,l2,m2, l2,m2,l2,m2,
          l2,m2,l2,m2, l2,m2,l2,m2, l2,m2,l2,m2
        );
      }
      fputc('\n', rfile[p]);
    }
  }



  for (int step = 0; step <= nraysteps; ++step)  {
    double stepratio = (double) step / (double) nraysteps;
    for(int p = 0; p < nrays; ++p) {
      cart3_t cart = cart3_scale(stepratio, rays[p]);
      sph_t sph = cart2sph(cart);
      double phase = cart3_dot(k, cart);
      ccart3_t E = ccart3_scale(cexp(I*phase), E0);
      csphvec_t L[nelem], M[nelem], N[nelem], Es, rEs;
      Es = ccart2csphvec(E, sph);
      rEs = qpms_eval_vswf(sph, Lc, Mc, Nc, lMax, QPMS_BESSEL_REGULAR, norm); // TODO maybe also check this by summing manually Lc * L etc.
      if (qpms_vswf_fill(L, M, N, lMax, sph, QPMS_BESSEL_REGULAR, norm)) abort();

      fprintf(rfile[p],"%d\t"
          "%g\t%g\t%g\t%g\t%g\t%g\t" "%g\t"
          "%g\t%g\t%g\t%g\t%g\t%g\t"
          "%g\t%g\t%g\t%g\t%g\t%g\t"
          "%g\t%g\t%g\t%g\t%g\t%g\t",
          step, cart.x, cart.y, cart.z, sph.r, sph.theta, sph.phi, phase,
          creal(E.x), cimag(E.x), creal(E.y), cimag(E.y), creal(E.z), cimag(E.z), 
          creal(Es.rc), cimag(Es.rc), creal(Es.thetac), cimag(Es.thetac), creal(Es.phic), cimag(Es.phic), 
          creal(rEs.rc), cimag(rEs.rc), creal(rEs.thetac), cimag(rEs.thetac), creal(rEs.phic), cimag(rEs.phic)
      );
      
      for(qpms_y_t y = 0; y < nelem; ++y) {
        csphvec_t /*Lsc,*/ Msc, Nsc; // TODO print also the longitudinal waves?
        //Lsc = csphvec_scale(Lc[y], L);
        Msc = csphvec_scale(Mc[y], M[y]);
        Nsc = csphvec_scale(Nc[y], N[y]);
        fprintf(rfile[p], "%g\t%g\t%g\t%g\t%g\t%g\t""%g\t%g\t%g\t%g\t%g\t%g\t",
            creal(Msc.rc), cimag(Msc.rc),
            creal(Msc.thetac), cimag(Msc.thetac),
            creal(Msc.phic), cimag(Msc.phic),
            creal(Nsc.rc), cimag(Nsc.rc),
            creal(Nsc.thetac), cimag(Nsc.thetac),
            creal(Nsc.phic), cimag(Nsc.phic)
        );
      }
      fputc('\n', rfile[p]);
    }
  }

  for (int p = 0; p < nrays; ++p)
    fclose(rfile[p]);
	gsl_rng_free(rng);
}

