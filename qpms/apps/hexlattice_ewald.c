// c99 -o ew_altin -DALTIN -Wall -I ../.. -O2 -ggdb -DLATTICESUMS32 hexlattice_ewald.c ../translations.c ../ewald.c ../ewaldsf.c ../gaunt.c ../lattices2d.c -lgsl -lm -lblas 
// c99 -o ew_kin -DALTIN -DKSTDIN -Wall -I ../.. -O2 -ggdb -DLATTICESUMS32 hexlattice_ewald.c ../translations.c ../ewald.c ../ewaldsf.c ../gaunt.c ../lattices2d.c -lgsl -lm -lblas 
// c99 -o ew -Wall -I ../.. -O2 -ggdb -DLATTICESUMS32 hexlattice_ewald.c ../translations.c ../ewald.c ../ewaldsf.c ../gaunt.c ../lattices2d.c -lgsl -lm -lblas 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <qpms/translations.h>
#include <qpms/lattices.h>
#include <gsl/gsl_const_mksa.h>



#define MAXOMEGACOUNT 1000
#define MAXKCOUNT 100 // serves as klist default buffer size if KSTDIN is defined
#define KMINCOEFF 0.998 // not used if KSTDIN defined
#define KMAXCOEFF 1.002 // not used if KSTDIN defined
#define KLAYERS 20
#define RLAYERS 20

const double s3 = 1.732050807568877293527446341505872366942805253810380628055;

// IMPORTANT: lattice properties here
const qpms_y_t lMax = 3;
const double REFINDEX = 1.52;
const double LATTICE_H = 576e-9;
static const double SCUFF_OMEGAUNIT = 3e14;
static const double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
static const double eV = GSL_CONST_MKSA_ELECTRON_CHARGE;
static const double c0 = GSL_CONST_MKSA_SPEED_OF_LIGHT;
static const TriangularLatticeOrientation rs_orientation = TRIANGULAR_VERTICAL;

int main (int argc, char **argv) {
  const double LATTICE_A = s3*LATTICE_H;
  const double INVLATTICE_A = 4*M_PI / s3 / LATTICE_A;

  if (argc < 3) abort();

  // char *kfile = argv[2]; // not used
  char *outfile = argv[2];
  char *errfile = NULL;
  if (argc > 3)
    errfile = argv[3];

#ifdef ALTIN // omega is provided on command line
  char *omegastr = argv[1];
  const double scuffomega = strtod(omegastr, NULL);
#else
  char *omegafile = argv[1];
  double scuffomegas[MAXOMEGACOUNT];
  FILE *f = fopen(omegafile, "r");
  size_t omegacount = 0;
  while (fscanf(f, "%lf", scuffomegas + omegacount) == 1){
    assert(omegacount < MAXOMEGACOUNT);
    ++omegacount;
  }
  fclose(f);
#endif

  /*f = fopen(kfile, "r");
  int kcount = 100;
  while (fscanf(f, "%lf %lf", &(klist[kcount].x), &(klist[kcount].y)) == 2) {
    assert(kcount < MAXKCOUNT);
    ++kcount;
  }
  fclose(f);
  */

#ifdef KSTDIN
  size_t kcount = 0;
  size_t klist_capacity = MAXKCOUNT;
  cart2_t *klist = malloc(sizeof(cart2_t) * klist_capacity);
  while (scanf("%lf %lf", &(klist[kcount].x), &(klist[kcount].y)) == 2) {
    ++kcount;
    if(kcount >= klist_capacity) {
      klist_capacity *= 2;
      klist = realloc(klist, sizeof(cart2_t) * klist_capacity);
      if (klist == NULL) abort();
    }
  }
#else
  cart2_t klist[MAXKCOUNT];
  int kcount = MAXKCOUNT;
  for (int i = 0; i <  kcount; ++i) { // TODO this should depend on orientation...
    klist[i].x = 0;
    klist[i].y = (4.* M_PI / 3. / LATTICE_A) * (KMINCOEFF + (KMAXCOEFF-KMINCOEFF)/kcount*i);
  }
#endif

  const double refindex = REFINDEX;
  const double h = LATTICE_H;
  const double a = h * s3;
  const double unitcell_area = s3*a*a/2.;
  const double rec_a = 4*M_PI/s3/a;

  //fprintf(stderr, "%.16g\n", rec_a);

  triangular_lattice_gen_t *Rlg = triangular_lattice_gen_init(a, rs_orientation, true, 0);
  triangular_lattice_gen_extend_to_steps(Rlg, RLAYERS);
  triangular_lattice_gen_t *Klg = triangular_lattice_gen_init(rec_a, reverseTriangularLatticeOrientation(rs_orientation), true, 0);
  triangular_lattice_gen_extend_to_steps(Klg, KLAYERS);

  points2d_rordered_t Rselwo0 = points2d_rordered_annulus(&(Rlg->ps), 0, false, INFINITY, false);

  const point2d *Rpoints = Rlg->ps.base;
  size_t nR = Rlg->ps.r_offsets[Rlg->ps.nrs];
  // TODO automatic skip of the zero point directly in translations.c
  const point2d *Rpoints_wo0 = Rselwo0.base + Rselwo0.r_offsets[0];
  size_t nR_wo0 = Rselwo0.r_offsets[Rselwo0.nrs] - Rselwo0.r_offsets[0];
  assert(nR - nR_wo0 == 1);
  const point2d *Kpoints = Klg->ps.base;
  size_t nK = Klg->ps.r_offsets[Klg->ps.nrs];

  const point2d pshift0 = {0, 0};
  point2d pshiftAB = {0, 0}, pshiftBA = {0,0};
  if(rs_orientation == TRIANGULAR_VERTICAL) { // CHECKSIGN
    pshiftAB.x = h/2;
    pshiftBA.x = -h/2;
  } else { // CHECKSIGN
    pshiftAB.y = -h/2;
    pshiftBA.y = h/2;
  }

  qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, QPMS_NORMALISATION_POWER); // vai POWER_CS?

  FILE *out = fopen(outfile, "w");
  FILE *err = NULL;
  if (errfile)
    err = fopen(errfile, "w");

#ifndef ALTIN
  for (size_t omegai = 0; omegai < omegacount; ++omegai) {
    const double scuffomega = scuffomegas[omegai];
#else
  {
#endif
    const double omega = scuffomega * SCUFF_OMEGAUNIT;
    const double EeV = omega * hbar / eV;
    const double k0_vac = omega / c0;
    const double k0_eff  = k0_vac * refindex;
    const double eta = 4*rec_a; // FIXME quite arbitrary

    // indices : destpart (A/B-particle), srcpart (A/B-particle), coeff type (A/B- type), desty, srcy
    complex double W[2][2][2][c->nelem][c->nelem];
    double Werr[2][2][2][c->nelem][c->nelem];

    for (size_t ki = 0; ki < kcount; ++ki) {
      cart2_t beta = klist[ki];
      memset(W, 0, sizeof(W));
      if(err)
        memset(Werr, 0, sizeof(Werr));

      const ptrdiff_t deststride = &(W[0][0][0][1][0]) - &(W[0][0][0][0][0]);
      const ptrdiff_t srcstride = &(W[0][0][0][0][1]) - &(W[0][0][0][0][0]);
      assert (srcstride == 1 && deststride == c->nelem);
      
      // A<-A
      qpms_trans_calculator_get_AB_arrays_e32_both_points_and_shift(c,
          &(W[0][0][0][0][0]), err ? &(Werr[0][0][0][0][0]) : NULL, // Adest, Aerr,
          &(W[0][0][1][0][0]), err ? &(Werr[0][0][1][0][0]) : NULL, // Bdest, Berr,
          deststride, srcstride,
          eta, k0_eff, unitcell_area, nR_wo0, Rpoints_wo0, nK, Kpoints, beta, pshift0
          );
      // B<-B
      for(qpms_y_t desty = 0; desty < c->nelem; ++desty) 
        for (qpms_y_t srcy = 0; srcy < c->nelem; ++srcy)
          for (int t = 0; t < 2; ++t) {
          W[1][1][t][desty][srcy] = W[0][0][t][desty][srcy];
          if (err)
            Werr[1][1][t][desty][srcy] = Werr[0][0][t][desty][srcy];
          }

      // A<-B
      qpms_trans_calculator_get_AB_arrays_e32_both_points_and_shift(c,
          &(W[0][1][0][0][0]), err ? &(Werr[0][1][0][0][0]) : NULL, // Adest, Aerr,
          &(W[0][1][1][0][0]), err ? &(Werr[0][1][1][0][0]) : NULL, // Bdest, Berr,
          deststride, srcstride,
          eta, k0_eff, unitcell_area, nR, Rpoints, nK, Kpoints, beta, pshiftAB
          );
      // B<-A
      qpms_trans_calculator_get_AB_arrays_e32_both_points_and_shift(c,
          &(W[1][0][0][0][0]), err ? &(Werr[1][0][0][0][0]) : NULL, // Adest, Aerr,
          &(W[1][0][1][0][0]), err ? &(Werr[1][0][1][0][0]) : NULL, // Bdest, Berr,
          deststride, srcstride,
          eta, k0_eff, unitcell_area, nR, Rpoints, nK, Kpoints, beta, pshiftBA
          );
      // TODO CHECK B<-A vs. A<-B relation

      fprintf(out, "%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t",
          scuffomega, EeV, k0_eff, beta.x, beta.y);
      if(err) fprintf(err, "%.16g\t%.16g\t%16g\t%.16g\t%.16g\t",
          scuffomega, EeV, k0_eff, beta.x, beta.y);
      size_t totalelems = sizeof(W) / sizeof(complex double);
      for (size_t i = 0; i < totalelems; ++i) {
        complex double w = ((complex double *)W)[i];
        fprintf(out, "%.16g\t%.16g\t", creal(w), cimag(w));
        if (err) 
          fprintf(err, "%.3g\t", ((double *)Werr)[i]);
      }
      fputc('\n', out);
      if(err) fputc('\n', err);
    }
  }

  fclose(out);
  if(err) fclose(err);

#ifdef KSTDIN
  free(klist);
#endif

  triangular_lattice_gen_free(Klg);
  triangular_lattice_gen_free(Rlg);
}
