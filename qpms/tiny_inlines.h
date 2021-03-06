#ifndef TINY_INLINES_H
#define TINY_INLINES_H
#include <stdlib.h>

static inline int min1pow(int pow) { return (pow % 2) ? -1 : 1; }


// This is useful for calculating spherical harmonics with negative m
// if spharm-normalised legendre functions for positive m are available.
// TODO: write a function that gets legendre buffer, m, n, and returns the correct spharm
// and use it in the code (mainly translations.c, ewald.c).
static inline int min1pow_m_neg(int m) {
	return (m < 0) ? min1pow(m) : 1;
}


#if 0
#ifdef __GSL_SF_LEGENDRE_H__
static inline complex double
spharm_eval(gsl_sf_legendre_t P_normconv, int P_csphase, qpms_l_t l, qpms_m_t m, double P_n_abs_m, complex double exp_imf) {

	return;
}
#endif
#endif

// this has shitty precision:
// static inline complex double ipow(int x) { return cpow(I, x); }

static inline complex double ipow(int x) {
  x = ((x % 4) + 4) % 4;
  switch(x) {
    case 0:
      return 1;
    case 1:
      return I;
    case 2:
      return -1;
    case 3:
      return -I;
    default:
      abort();
  }
}

static inline int isq(int x) {return x * x;}

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) (((x) >= (y)) ? (x) : (y))
#endif

#ifndef SQ
#define SQ(x) ((x) * (x))
#endif


#endif // TINY_INLINES_H
