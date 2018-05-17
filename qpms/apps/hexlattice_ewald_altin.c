// c99 -ggdb -O2 -DLATTICESUMS -I .. hexlattice_ewald.c ../translations.c ../bessels.c ../lrhankel_recspace_dirty.c ../gaunt.c -lm  -lgsl -lblas  
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "kahansum.h"
#include "vectors.h"
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include "qpms_types.h"
#include "translations.h"

const double s3 = 1.732050807568877293527446341505872366942805253810380628055;

// IMPORTANT: lattice properties here
const qpms_y_t lMax = 2;
const double REFINDEX = 1.52;
const double LATTICE_H = 576e-9;
static const double SCUFF_OMEGAUNIT = 3e14;
static const double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
static const double eV = GSL_CONST_MKSA_ELECTRON_CHARGE;
static const double c0 = GSL_CONST_MKSA_SPEED_OF_LIGHT;
static const double CC = 0.1;

// For sorting the points by distance from origin / radius
int cart2_cmpr (const void *p1, const void *p2) {
  const  cart2_t *p1t = (const  cart2_t *)p1;
  const  cart2_t *p2t = (const  cart2_t *)p2;
  double r21 = cart2norm(*p1t);
  double r22 = cart2norm(*p2t);
  if (r21 < r22) return -1;
  else if (r21 > r22) return 1;
  else return 0;
}

typedef struct {
  ptrdiff_t npoints; // number of lattice points.
  ptrdiff_t capacity; // for how much points memory is allocated
  double maxR; // circle radius, points <= R
   cart2_t *points;
} latticepoints_circle_t;

void sort_cart2points_by_r(cart2_t *points, size_t nmemb) {
  qsort(points, nmemb, sizeof(cart2_t), cart2_cmpr);
}

void latticepoints_circle_free(latticepoints_circle_t *c) {
  free(c->points);
  c->capacity = 0;
}

// "horizontal" orientation of the adjacent A, B points
latticepoints_circle_t generate_hexpoints_hor(double h, double R, cart2_t offset /* if zero, an A is in the origin */) {
  latticepoints_circle_t lat;
  lat.maxR = R;
  lat.npoints = 0;
  int nmax = R / (1.5 * h) + 2; // max no of lattice shifts in each direction (with reserve)
  double unitcellS = s3 * 3 / 2 * h * h; // unit cell area
  double flcapacity =  5 + 2 * (R + 5*h) * (R + 5*h) * M_PI / unitcellS; // should be enough with some random reserve
  lat.capacity = flcapacity;

  lat.points = malloc(lat.capacity *sizeof(cart2_t));

  cart2_t BAoffset = {h, 0};

  cart2_t a1 = {-1.5*h, s3/2 *h};
  cart2_t a2 = {1.5*h, s3/2 *h};

  for (ptrdiff_t i1 = -nmax; i1 <= nmax; ++i1)
    for (ptrdiff_t i2 = -nmax; i2 <= nmax; ++i2) {
      cart2_t Apoint = cart2_add(offset, cart2_add(cart2_scale(i1, a1), cart2_scale(i2, a2)));
         if (lat.npoints >= lat.capacity) 
          printf("%zd %zd %g %g %g %g\n", lat.npoints, lat.capacity, flcapacity, R, h, unitcellS);     if (cart2norm(Apoint) <= R) {
        assert(lat.npoints < lat.capacity);
        lat.points[lat.npoints] = Apoint;
        lat.npoints++;
      }
      cart2_t Bpoint = cart2_add(Apoint, BAoffset);
      if (cart2norm(Bpoint) <= R) {
        assert(lat.npoints < lat.capacity);
        lat.points[lat.npoints] = Bpoint;
        lat.npoints++;
      }
    }
  sort_cart2points_by_r(lat.points, lat.npoints);
  return lat;
}

latticepoints_circle_t generate_tripoints_ver(double a, double R, cart2_t offset /* if zero, an A is in the origin */) {
  double h = a / s3;
  latticepoints_circle_t lat;
  lat.maxR = R;
  lat.npoints = 0;
  int nmax = R / (1.5 * h) + 2; // max no of lattice shifts in each direction (with reserve)
  double unitcellS = (s3 * 3) / 2 * h * h; // unit cell area
  double flcapacity =  5 +  (R + 3*a) * (R + 3*a) * M_PI / unitcellS; // should be enough with some random reserve
  lat.capacity = flcapacity; // should be enough with some random reserve
  lat.points = malloc(lat.capacity *sizeof(cart2_t));

  cart2_t a1 = {-1.5*h, s3/2 *h};
  cart2_t a2 = {1.5*h, s3/2 *h};

  for (ptrdiff_t i1 = -nmax; i1 <= nmax; ++i1)
    for (ptrdiff_t i2 = -nmax; i2 <= nmax; ++i2) {
      cart2_t Apoint = cart2_add(offset, cart2_add(cart2_scale(i1, a1), cart2_scale(i2, a2)));
      if (cart2norm(Apoint) <= R) {
        if (lat.npoints >= lat.capacity) 
          printf("%zd %zd %g %g %g %g\n", lat.npoints, lat.capacity, flcapacity, R, a, unitcellS);
        assert(lat.npoints < lat.capacity);
        lat.points[lat.npoints] = Apoint;
        lat.npoints++;
      }
    }
  sort_cart2points_by_r(lat.points, lat.npoints);
  return lat;
}

latticepoints_circle_t generate_tripoints_hor(double a, double R, cart2_t offset /* if zero, an A is in the origin */) {
  double h = a / s3;
  latticepoints_circle_t lat;
  lat.maxR = R;
  lat.npoints = 0;
  int nmax = R / (1.5 * h) + 2; // max no of lattice shifts in each direction (with reserve)
  double unitcellS = s3 * 3 / 2 * h * h; // unit cell area
  double flcapacity =  5 +  (R + 3*a) * (R + 3*a) * M_PI / unitcellS; // should be enough with some random reserve
  lat.capacity = flcapacity; // should be enough with some random reserve
  lat.points = malloc(lat.capacity *sizeof(cart2_t));

  cart2_t a1 = {s3/2 *h, -1.5*h};
  cart2_t a2 = {s3/2 *h, 1.5 * h};

  for (int i1 = -nmax; i1 <= nmax; ++i1)
    for (int i2 = -nmax; i2 <= nmax; ++i2) {
      if (lat.npoints >= lat.capacity) 
          printf("%zd %zd %.12g %g %g %g\n", lat.npoints, lat.capacity, flcapacity, R, a, unitcellS);
      cart2_t Apoint = cart2_add(offset, cart2_add(cart2_scale(i1, a1), cart2_scale(i2, a2)));
      if (cart2norm(Apoint) <= R) {
        assert(lat.npoints < lat.capacity);
        lat.points[lat.npoints] = Apoint;
        lat.npoints++;
      }
    }
  sort_cart2points_by_r(lat.points, lat.npoints);
  return lat;
}

int main (int argc, char **argv) {
  const double LATTICE_A = s3*LATTICE_H;
  const double INVLATTICE_A = 4 * M_PI / s3 / LATTICE_A;
  const double MAXR_REAL = 100 * LATTICE_H;
  const double MAXR_K = 100 * INVLATTICE_A;


  char *omegastr = argv[1];
  // char *kfile = argv[2]; // not used;, will be read from stdin
  char *outprefix = argv[2];

  double scuffomega = strtod(omegastr, NULL);
  char *outfile = outprefix;
  char outlongfile[strlen(outprefix)+10];
  char outshortfile[strlen(outprefix)+10];
  sprintf(outlongfile, "%s.long", outprefix);
  sprintf(outshortfile, "%s.long", outprefix);



  //cart2_t klist[MAXKCOUNT];
  /*f = fopen(kfile, "r");
  int kcount = 100;
  while (fscanf(f, "%lf %lf", &(klist[kcount].x), &(klist[kcount].y)) == 2) {
    assert(kcount < MAXKCOUNT);
    ++kcount;
  }
  fclose(f);
  */

  const double refindex = REFINDEX;
  const double h = LATTICE_H;
  const double a = h * s3;
  const double rec_a = 4*M_PI/s3/a;
  
  // generation of the real-space lattices
  const cart2_t cart2_0 = {0, 0};
  const cart2_t ABoffset = {h, 0};
  const cart2_t BAoffset = {-h, 0};
  //const cart2_t ab_particle_offsets[2][2] = {{ {0, 0}, {h, 0} }, {-h, 0}, {0, 0}};

  // THIS IS THE LATTICE OF r_b
  latticepoints_circle_t lattice_0offset = generate_tripoints_ver(a, MAXR_REAL, cart2_0);
  // these have to have the same point order, therefore we must make the offset verision manually to avoid sorting;
  latticepoints_circle_t lattice_ABoffset, lattice_BAoffset;
  lattice_ABoffset.points = malloc(lattice_0offset.npoints * sizeof(cart2_t));
  lattice_ABoffset.capacity = lattice_0offset.npoints * sizeof(cart2_t);
  lattice_ABoffset.npoints = lattice_ABoffset.capacity;
  lattice_BAoffset.points = malloc(lattice_0offset.npoints * sizeof(cart2_t));
  lattice_BAoffset.capacity = lattice_0offset.npoints * sizeof(cart2_t);
  lattice_BAoffset.npoints = lattice_BAoffset.capacity;
  for (int i = 0; i < lattice_0offset.npoints; ++i) {
    lattice_ABoffset.points[i] = cart2_add(lattice_0offset.points[i], ABoffset);
    lattice_BAoffset.points[i] = cart2_add(lattice_0offset.points[i], BAoffset);
  }

  // reciprocal lattice, without offset – DON'T I NEED REFINDEX HERE? (I DON'T THINK SO.)
  latticepoints_circle_t reclattice = generate_tripoints_hor(rec_a, MAXR_K, cart2_0);

  qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, QPMS_NORMALISATION_POWER_CS);

  FILE *out = fopen(outfile, "w");
  FILE *outlong = fopen(outlongfile, "w");
  FILE *outshort = fopen(outshortfile, "w");

  // as in eq. (5) in my notes
  double WL_prefactor = 4*M_PI/(a*a)/s3 / /*??*/ (4*M_PI*M_PI);

    //double scuffomega = scuffomegas[omegai];
    double omega = scuffomega * SCUFF_OMEGAUNIT;
    double EeV = omega * hbar / eV;
    double k0_vac = omega / c0;
    double k0_eff = k0_vac * refindex; // this one will be used with the real x geometries
    double cv = CC * k0_eff;

    complex double Abuf[c->nelem][c->nelem], Bbuf[c->nelem][c->nelem];
    // indices : destpart (A/B-particle), srcpart (A/B-particle), coeff type (A/B- type), desty, srcy
    complex double WS[2][2][2][c->nelem][c->nelem];
    complex double WS_comp[2][2][2][c->nelem][c->nelem];
    complex double WL[2][2][2][c->nelem][c->nelem];
    complex double WL_comp[2][2][2][c->nelem][c->nelem];

    //for (int ki = 0; ki < kcount; ++ki) {
    //  cart2_t k = klist[ki];
    cart2_t k;
    while(scanf("%lf %lf", &(k.x), &(k.y)) == 2) {
      memset(WS, 0, sizeof(WS));
      memset(WS_comp, 0, sizeof(WS_comp));
      memset(WL, 0, sizeof(WL));
      memset(WL_comp, 0, sizeof(WL_comp));

      for (int bi = 0; bi < lattice_0offset.npoints; ++bi) {
        cart2_t point0 = lattice_0offset.points[bi];
        double phase = cart2_dot(k,point0);
        complex double phasefac = cexp(I*phase);

        if (point0.x || point0.y) { // skip the singular point
          qpms_trans_calculator_get_shortrange_AB_arrays(c, (complex double *) Abuf, (complex double *) Bbuf, c->nelem, 1,
              cart22sph(cart2_scale(k0_eff,lattice_0offset.points[bi])), 3, 2, 5, CC);
          for (int desty = 0; desty < c->nelem; ++desty)
            for (int srcy = 0; srcy < c->nelem; ++srcy) {
              ckahanadd(&(WS[0][0][0][desty][srcy]),&(WS_comp[0][0][0][desty][srcy]),Abuf[desty][srcy] * phasefac);
              ckahanadd(&(WS[0][0][1][desty][srcy]),&(WS_comp[0][0][1][desty][srcy]),Bbuf[desty][srcy] * phasefac);
            }
        }
        qpms_trans_calculator_get_shortrange_AB_arrays(c, (complex double *) Abuf, (complex double *) Bbuf, c->nelem, 1,
            cart22sph(cart2_scale(k0_eff,lattice_ABoffset.points[bi])), 3, 2, 5, CC);
        for (int desty = 0; desty < c->nelem; ++desty)
            for (int srcy = 0; srcy < c->nelem; ++srcy) {
              ckahanadd(&(WS[0][1][0][desty][srcy]),&(WS_comp[0][1][0][desty][srcy]),Abuf[desty][srcy] * phasefac);
              ckahanadd(&(WS[0][1][1][desty][srcy]),&(WS_comp[0][1][1][desty][srcy]),Bbuf[desty][srcy] * phasefac);
            }

        qpms_trans_calculator_get_shortrange_AB_arrays(c, (complex double *) Abuf, (complex double *) Bbuf, c->nelem, 1,
            cart22sph(cart2_scale(k0_eff,lattice_BAoffset.points[bi])), 3, 2, 5, CC);
        for (int desty = 0; desty < c->nelem; ++desty)
            for (int srcy = 0; srcy < c->nelem; ++srcy) {
              ckahanadd(&(WS[1][0][0][desty][srcy]),&(WS_comp[1][0][0][desty][srcy]),Abuf[desty][srcy] * phasefac);
              ckahanadd(&(WS[1][0][1][desty][srcy]),&(WS_comp[1][0][1][desty][srcy]),Bbuf[desty][srcy] * phasefac);
            }
        // WS[1][1] is the same as WS[0][0], so copy in the end rather than double-summing
      }
      for (int desty = 0; desty < c->nelem; ++desty)
        for (int srcy = 0; srcy < c->nelem; ++srcy)
          for (int ctype = 0; ctype < 2; ctype++)
            WS[1][1][ctype][desty][srcy] = WS[0][0][ctype][desty][srcy];
      // WS DONE
      for (int Ki = 0; Ki < reclattice.npoints; ++Ki) {
        cart2_t K = reclattice.points[Ki];
        cart2_t k_K = cart2_substract(k, K);
        double phase_AB = 
#ifdef SWAPSIGN1
          -
#endif
          cart2_dot(k_K, ABoffset); // And maybe the sign is excactly opposite!!! FIXME TODO CHECK
        complex double phasefacs[2][2];
        phasefacs[0][0] = phasefacs[1][1] = 1;
        phasefacs[1][0] = cexp(I * phase_AB); // sign???
        phasefacs[0][1] = cexp(- I * phase_AB); // sign???

        // FIXME should I skip something (such as the origin?)
        qpms_trans_calculator_get_2DFT_longrange_AB_arrays(c, (complex double *) Abuf, (complex double *) Bbuf, c->nelem, 1,
            cart22sph(k_K), 3, 2, 5, cv, k0_eff);
        for (int dp = 0; dp < 2; dp++)
          for (int sp = 0; sp < 2; sp++)
            for (int dy = 0; dy < c->nelem; dy++)
              for (int sy = 0; sy < c->nelem; sy++) {
                ckahanadd(&(WL[dp][sp][0][dy][sy]), &(WL_comp[dp][sp][0][dy][sy]), phasefacs[dp][sp] * Abuf[dy][sy] * WL_prefactor);
                ckahanadd(&(WL[dp][sp][1][dy][sy]), &(WL_comp[dp][sp][1][dy][sy]), phasefacs[dp][sp] * Bbuf[dy][sy] * WL_prefactor);
              }
      }
      fprintf(outshort, "%.16g\t%.16g\t%16g\t%.16g\t%.16g\t",
          scuffomega, EeV, k0_eff, k.x, k.y);
      fprintf(outlong, "%.16g\t%.16g\t%16g\t%.16g\t%.16g\t",
          scuffomega, EeV, k0_eff, k.x, k.y);
      fprintf(out, "%.16g\t%.16g\t%16g\t%.16g\t%.16g\t",
          scuffomega, EeV, k0_eff, k.x, k.y);
      size_t totalelems = sizeof(WL) / sizeof(complex double);
      for (int i = 0; i < totalelems; ++i) {
        complex double ws = ((complex double *)WS)[i];
        complex double wl = ((complex double *)WL)[i];
        complex double w = ws+wl;
        fprintf(outshort, "%.16g\t%.16g\t", creal(ws), cimag(ws));
        fprintf(outlong, "%.16g\t%.16g\t", creal(wl), cimag(wl));
        fprintf(out, "%.16g\t%.16g\t", creal(w), cimag(w));
      }
      fputc('\n', outshort);
      fputc('\n', outlong);
      fputc('\n', out);
      fflush(outshort);
      fflush(outlong);
      fflush(out);
    }
  fclose(out);
  fclose(outlong);
  fclose(outshort);

}

  


#if 0
int main (int argc, char **argv) {
  cart2_t offset = {0,0};

  latticepoints_circle_t lat = generate_tripoints_ver(1, 200, offset);
  for (int i = 0; i < lat.npoints; ++i) 
    printf("%g %g %g\n", lat.points[i].x, lat.points[i].y, cart2norm(lat.points[i]));
  latticepoints_circle_free(&lat);
}
#endif
