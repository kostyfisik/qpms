// TODO complex kr version

#include <stdbool.h>
#include <complex.h>
#include <string.h>
#include <qpms/translations.h>
#include <qpms/qpms_error.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <qpms/indexing.h>

#define RTOL (1e-8)
#define ATOL (1e-14)

int isclose_cmplx(complex double a, complex double b) {
  return cabs(a-b) <= ATOL + RTOL * .5 * (cabs(a) + cabs(b));
}

static inline size_t ssq(size_t s) { return s * s; }

int test_AB_single_vs_array(const qpms_trans_calculator *c, qpms_bessel_t wavetype,
    cart3_t kd_cart) 
{
  int fails = 0;
  csph_t kd_sph = sph2csph(cart2sph(kd_cart));

  complex double A[ssq(c->nelem)], B[ssq(c->nelem)];
  QPMS_ENSURE_SUCCESS(qpms_trans_calculator_get_AB_arrays(c, A, B, c->nelem, 1, kd_sph, false, wavetype));
  
  for (qpms_y_t ydest = 0; ydest < c->nelem; ++ydest) {
    qpms_l_t ldest; qpms_m_t mdest; qpms_y2mn_p(ydest, &mdest, &ldest);
    for (qpms_y_t ysrc = 0; ysrc < c->nelem; ++ysrc) {
      qpms_l_t lsrc; qpms_m_t msrc; qpms_y2mn_p(ysrc, &msrc, &lsrc);
      complex double Asingle = qpms_trans_single_A(c->normalisation, mdest, ldest, msrc, lsrc, kd_sph, false, wavetype);
      complex double Aarr = A[ysrc + c->nelem * ydest];
      if (!isclose_cmplx(Asingle, Aarr)) {
        ++fails;
        fprintf(stderr, "l=%d,m=%+d <- l=%d,m=%+d: A_single=%.16g%+.16gj,\tA_arr=%.16g%+.16gj,\tdiff=%.g\t(norm=%x)\n",
            (int)ldest, (int)mdest, (int)lsrc, (int)msrc, creal(Asingle), cimag(Asingle), 
            creal(Aarr), cimag(Aarr), cabs(Aarr-Asingle), (unsigned int)(c->normalisation));
      }
      complex double Bsingle = qpms_trans_single_B(c->normalisation, mdest, ldest, msrc, lsrc, kd_sph, false, wavetype);
      complex double Barr = B[ysrc + c->nelem * ydest];
      if (!isclose_cmplx(Bsingle, Barr)) {
        ++fails;
        fprintf(stderr, "l=%d,m=%+d <- l=%d,m=%+d: B_single=%.16g%+.16gj,\tB_arr=%.16g%+.16gj,\tdiff=%.g\t(norm=%x)\n",
            (int)ldest, (int)mdest, (int)lsrc, (int)msrc, creal(Bsingle), cimag(Bsingle), 
            creal(Barr), cimag(Barr), cabs(Barr-Bsingle), (unsigned int)(c->normalisation));
      }
    }
  }
  return fails;
}  

int main() {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxs0);
  gsl_rng_set(rng, 666);

	qpms_l_t lMax = 3;
	int npoints = 10;
	double sigma = 12;

  cart3_t points[npoints];
  double relerrs[npoints];
  memset(points, 0, npoints * sizeof(cart3_t));
  points[0].x = points[1].y = points[2].z = sigma;
  for (unsigned i = 3; i < npoints; ++i) {
    cart3_t *w = points+i;
    w->x = gsl_ran_gaussian(rng, sigma);
    w->y = gsl_ran_gaussian(rng, sigma);
    w->z = gsl_ran_gaussian(rng, sigma);
  }

  int fails = 0;

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
        for(int p = 0; p < npoints; ++p)
          fails += test_AB_single_vs_array(c, J, points[p]);
      qpms_trans_calculator_free(c);
    }
  }  

  gsl_rng_free(rng);
  
  return fails;
}

