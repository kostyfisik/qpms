#ifndef QPMS_TRANSLATIONS_H
#define QPMS_TRANSLATIONS_H
#include "vectors.h"
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
double qpms_legendre0(int m, int n); 
double qpms_legendred0(int m, int n); 

typedef enum {
	QPMS_NORMALIZATION_TAYLOR = 1,
	QPMS_NORMALIZATION_UNDEF = 0
} qpms_normalization_t;

typedef enum {
        QPMS_BESSEL_REGULAR = 1, // regular function j
        QPMS_BESSEL_SINGULAR = 2, // singular function y
        QPMS_HANKEL_PLUS = 3, // hankel function h1 = j + I*y
        QPMS_HANKEL_MINUS = 4, // hankel function h2 = j - I*y
	QPMS_BESSEL_UNDEF = 0
} qpms_bessel_t;

int qpms_sph_bessel_array(qpms_bessel_t typ, int lmax, double x, complex double *result_array);

complex double qpms_trans_single_A_Taylor(int m, int n, int mu, int nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_single_B_Taylor(int m, int n, int mu, int nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_single_A_Taylor_ext(int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

complex double qpms_trans_single_B_Taylor_ext(int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);


typedef struct qpms_trans_calculator {
        int lMax;
        size_t nelem;
        complex double **A_multipliers;
        complex double **B_multipliers;
        qpms_normalization_t normalization;
} qpms_trans_calculator;

qpms_trans_calculator *qpms_trans_calculator_init(int lMax, qpms_normalization_t nt);
void qpms_trans_calculator_free(qpms_trans_calculator *);

#endif // QPMS_TRANSLATIONS_H
