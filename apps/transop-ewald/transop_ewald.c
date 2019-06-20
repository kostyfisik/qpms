// c99 -o ew_gen_kin -Wall -I ../.. -I ../../amos/ -O2 -ggdb -DQPMS_VECTORS_NICE_TRANSFORMATIONS -DLATTICESUMS32 2dlattice_ewald.c ../translations.c ../ewald.c ../ewaldsf.c ../gaunt.c ../lattices2d.c ../latticegens.c ../bessel.c -lgsl -lm -lblas ../../amos/libamos.a -lgfortran ../error.c

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "transop_ewald_cmdline.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#define LATTICESUMS32
#include <qpms/translations.h>
#include <qpms/lattices.h>
#include <qpms/qpms_error.h>
#include <gsl/gsl_const_mksa.h>
#include <qpms/parsing.h>

// Command line order:
// outfile b1.x b1.y b2.x b2.y lMax scuffomega refindex npart part0.x part0.y [part1.x part1.y [...]]
//
// Standard input (per line):
// k.x k.y
//
// Output data format (line):
//

#define MAXKCOUNT 200 // 200 // serves as klist default buffer size 
//#define KMINCOEFF 0.783 //0.9783 // 0.783 // not used if KSTDIN defined
//#define KMAXCOEFF 1.217 //1.0217 // 1.217 // not used if KSTDIN defined
#define KLAYERS 20
#define RLAYERS 20

const double s3 = 1.732050807568877293527446341505872366942805253810380628055;

//const qpms_y_t lMax = 3;
//const double REFINDEX = 1.52;
static const double SCUFF_OMEGAUNIT = 3e14;
static const double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
static const double eV = GSL_CONST_MKSA_ELECTRON_CHARGE;
static const double c0 = GSL_CONST_MKSA_SPEED_OF_LIGHT;

int main (int argc, char **argv) {
  struct gengetopt_args_info args_info;

  int retval = cmdline_parser(argc, argv, *args_info);
  if (retval) return retval;

  // Parse lattice vectors
  const int latdim = args_info.base_vector_given;
  QPMS_ENSURE(latdim == 2,
      "Sorry, only 2d lattices supported, but %d basis vectors were given\n",
      latdim);
  cart2_t b[latdim];
  for (int i = 0; i < latdim; ++i) {
    const int gotnumbers = qpms_parse_ndoubles(
          (*double) &(b[i].x), latdim,
          args_info.base_vector_arg[i]);
    QPMS_ENSURE(latdim == gotnumbers,
        "%d. base vector contained %d numbers, expected %d\n",
        i, gotnumbers, latdim);
  }

  // N.B. this is 2D specific, TODO generalize when Nd sum supported
  const double unitcell_area = l2d_unitcell_area(b[0], b[1]);
  l2d_reduceBasis(b[0], b[1], b, b+1);
  
  const qpms_l_t lMax = args_info.lMax_arg;
  QPMS_ENSURE(lMax > 0, "invalid value of lMax: %d", (int)lMax);

  const double refindex = args_info.refractive_index_arg;

  // Parse all particle positions
  const int npart = args_info.particle_given;
  if(!npart) ++npart;
  cart2_t part_positions[npart];
  if(!args_info.particle_given)
    part_positions[0].x = part_positions[0].y = 0;
  else for (int i = 0; i < npart; ++i) {
    const int gotnumbers = qpms_parse_ndoubles(
        (*double) &(part_positions[i].x), latdim,
        args_info.particle_arg[i]);
    QPMS_ENSURE(latdim == gotnumbers,
        "%d. particle position contained %d coordinates, expected %d\n",
        i, gotnumbers, latdim);
  }

  QPMS_ENSURE(!args_info.k_omega_meshgrid_mode_counter != !args_info.k_omega_points_mode_counter,
      "THIS IS A BUG. Only one mode ((k, ω) tuples, or k, ω lists) allowed.");
  if (args_info.k_omega_meshgrid_mode_counter) { // grid mode
    size_t omegacount = 0;
    double *omegalist = NULL;
    for (int i = 0; i < args_info.omega_given; ++i) 
      omegacount = qpms_parse_doubles(&omegalist, omegacount,
          args_info.omega_arg[i]);

    for (int i = 0; i < args_info.omegafile_given; ++i) {
      FILE *f omegaf = fopen(args_info.omegafile_arg[i], "r");




  } else if (args_info.k_omega_points_mode_counter) { // explic. point mode
    TODO;
  }


  const double scuffomega = strtod(argv[7], NULL);
//#ifdef KSTDIN
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
//#else
#if 0
  cart2_t klist[MAXKCOUNT];
  int kcount = MAXKCOUNT;
  for (int i = 0; i <  kcount; ++i) { // TODO this should depend on orientation...
    klist[i].x = 0;
    klist[i].y = (4.* M_PI / 3. / LATTICE_A) * (KMINCOEFF + (KMAXCOEFF-KMINCOEFF)/kcount*i);
  }
#endif

  
  // TODO more clever way of determining the cutoff
  const double a = sqrt(unitcell_area); // N.B. different meaning than before
  const double maxR = 25 * a;
  const double maxK = 25 * 2*M_PI/a;

  qpms_trans_calculator *c = qpms_trans_calculator_init(lMax, QPMS_NORMALISATION_POWER_CS); // vai POWER_CS?

  FILE *out, *ferr = NULL;
  if (args_info.error_estimate_output_given) {
    if (!strcmp(args_info.error_estimate_output_arg, "-"))
      ferr = stdout;
    else 
      ferr = fopen(args_info.error_estimate_output_arg, "w");
    QPMS_ENSURE(ferr, "Could not open error output file %s", 
        args_info.error_estimate_output_arg);
  if (args_info.output_given && !strcmp(args_info.output_arg, "-")
      && args_info.output_arg[0]) {
    out = fopen(args_info.output_arg, "w");
    QPMS_ENSURE(out, "Could not open output file %s", args_info.output_arg);
  } else
    out = stdout;

  {
    const double omega = scuffomega * SCUFF_OMEGAUNIT;
    const double EeV = omega * hbar / eV;
    const double k0_vac = omega / c0;
    const double k0_eff  = k0_vac * refindex;
    const double eta = 5.224/a; // FIXME quite arbitrary, but this one should work

    // indices : destpart (A/B-particle), srcpart (A/B-particle), coeff type (A/B- type), desty, srcy
    complex double W[npart][npart][2][c->nelem][c->nelem];
    double Werr[npart][npart][npart][c->nelem][c->nelem];

    for (size_t ki = 0; ki < kcount; ++ki) {
      cart2_t beta = klist[ki];
      memset(W, 0, sizeof(W));
      if(ferr)
        memset(Werr, 0, sizeof(Werr));

      const ptrdiff_t deststride = &(W[0][0][0][1][0]) - &(W[0][0][0][0][0]);
      const ptrdiff_t srcstride = &(W[0][0][0][0][1]) - &(W[0][0][0][0][0]);
      assert (srcstride == 1 && deststride == c->nelem);

      for (size_t ps = 0; ps < npart; ++ps) for (size_t pd = 0; pd < npart; ++pd) 
        // TODO optimize (calculate only once for each particle shift; especially if pd == ps)
        qpms_trans_calculator_get_AB_arrays_e32(c,
            &(W[pd][ps][0][0][0]), ferr ? &(Werr[pd][ps][0][0][0]) : NULL, // Adest, Aerr,
            &(W[pd][ps][1][0][0]), ferr ? &(Werr[pd][ps][1][0][0]) : NULL, // Bdest, Berr,
            deststride, srcstride,
            eta, k0_eff, b1, b2,
            beta,
            cart2_substract(part_positions[pd], part_positions[ps]),  // CHECKSIGN
            maxR, maxK 
            );
        // TODO CHECK B<-A vs. A<-B relation

      fprintf(out, "%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t",
          scuffomega, EeV, k0_eff, beta.x, beta.y);
      if(ferr) fprintf(ferr, "%.16g\t%.16g\t%16g\t%.16g\t%.16g\t",
          scuffomega, EeV, k0_eff, beta.x, beta.y);
      size_t totalelems = sizeof(W) / sizeof(complex double);
      for (size_t i = 0; i < totalelems; ++i) {
        complex double w = ((complex double *)W)[i];
        fprintf(out, "%.16g\t%.16g\t", creal(w), cimag(w));
        if (ferr) 
          fprintf(ferr, "%.3g\t", ((double *)Werr)[i]);
      }
      fputc('\n', out);
      if(ferr) fputc('\n', ferr);
    }
  }

  fclose(out);
  if(ferr) fclose(ferr);

//#ifdef KSTDIN
  free(klist);
//#endif

  qpms_trans_calculator_free(c);

}
