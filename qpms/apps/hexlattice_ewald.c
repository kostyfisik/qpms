// c99 -o ew -Wall -I ../.. -O2 -ggdb -DLATTICESUMS32 hexlattice_ewald.c ../translations.c ../ewald.c ../ewaldsf.c ../gaunt.c ../lattices2d.c -lgsl -lm -lblas 
#include <stdio.h>
#include <math.h>
#include <qpms/translations.h>
#include <qpms/lattices.h>
#include <gsl/gsl_const_mksa.h>



#define MAXOMEGACOUNT 1000
#define MAXKCOUNT 100
#define KLAYERS 20
#define RLAYERS 20

const double s3 = 1.732050807568877293527446341505872366942805253810380628055;

// IMPORTANT: lattice properties here
const qpms_y_t lMax = 2;
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
  
  char *omegafile = argv[1];
  // char *kfile = argv[2]; // not used
  char *outfile = argv[3];
  char *errfile = NULL;
  if (argc > 4)
    errfile = argv[4];
  //char *outlongfile = argv[4];
  //char *outshortfile = argv[5];
  double scuffomegas[MAXOMEGACOUNT];
  cart2_t klist[MAXKCOUNT];
  FILE *f = fopen(omegafile, "r");
  size_t omegacount = 0;
  while (fscanf(f, "%lf", scuffomegas + omegacount) == 1){
    assert(omegacount < MAXOMEGACOUNT);
    ++omegacount;
  }
  fclose(f);
  /*f = fopen(kfile, "r");
  int kcount = 100;
  while (fscanf(f, "%lf %lf", &(klist[kcount].x), &(klist[kcount].y)) == 2) {
    assert(kcount < MAXKCOUNT);
    ++kcount;
  }
  fclose(f);
  */

  int kcount = MAXKCOUNT;
  for (int i = 0; i <  kcount; ++i) { // TODO this should depend on orientation...
    klist[i].x = 0;
    klist[i].y = 2. * 4. * M_PI / 3. / LATTICE_A / kcount * i;
  }

  const double refindex = REFINDEX;
  const double h = LATTICE_H;
  const double a = h * s3;
  const double unitcell_area = s3*a*a/2.;
  const double rec_a = 4*M_PI/s3/a;

  triangular_lattice_gen_t *Rlg = triangular_lattice_gen_init(a, rs_orientation, true, 0);
  triangular_lattice_gen_extend_to_steps(Rlg, RLAYERS);
  triangular_lattice_gen_t *Klg = triangular_lattice_gen_init(rec_a, reverseTriangularLatticeOrientation(rs_orientation), true, 0);
  triangular_lattice_gen_extend_to_steps(Klg, KLAYERS);

  const point2d *Rpoints = Rlg->ps.base;
  size_t nR = Rlg->ps.r_offsets[Rlg->ps.nrs];
  const point2d *Kpoints = Klg->ps.base;
  size_t nK = Klg->ps.r_offsets[Klg->ps.nrs];

  const point2d pshift0 = {0, 0};
  point2d pshiftAB = {0, 0}, pshiftBA = {0,0};
  if(rs_orientation == TRIANGULAR_VERTICAL) { // CHECKSIGN
    pshiftAB.x = h;
    pshiftBA.x = -h;
  } else { // CHECKSIGN
    pshiftAB.y = -h;
    pshiftBA.y = h;
  }

  qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, QPMS_NORMALISATION_POWER); // vai POWER_CS?

  FILE *out = fopen(outfile, "w");
  FILE *err = NULL;
  if (errfile)
    err = fopen(errfile, "w");

  for (size_t omegai = 0; omegai < omegacount; ++omegai) {
    const double scuffomega = scuffomegas[omegai];
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
          eta, k0_eff, unitcell_area, nR, Rpoints, nK, Kpoints, beta, pshift0
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

      fprintf(out, "%.16g\t%.16g\t%16g\t%.16g\t%.16g\t",
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

  triangular_lattice_gen_free(Klg);
  triangular_lattice_gen_free(Rlg);
}
