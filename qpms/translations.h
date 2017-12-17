#ifndef QPMS_TRANSLATIONS_H
#define QPMS_TRANSLATIONS_H
#include "vectors.h"
#include "qpms_types.h"
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
// Associated Legendre polynomial at zero argument (DLMF 14.5.1)
double qpms_legendre0(int m, int n); 
// Associated Legendre polynomial derivative at zero argument (DLMF 14.5.2)
double qpms_legendred0(int m, int n); 
int qpms_sph_bessel_array(qpms_bessel_t typ, int lmax, double x, complex double *result_array);

// TODO replace the xplicit "Taylor" functions with general,
// taking qpms_bessel_t argument.
complex double qpms_trans_single_A_Taylor(int m, int n, int mu, int nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_single_B_Taylor(int m, int n, int mu, int nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_single_A_Taylor_ext(int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

complex double qpms_trans_single_B_Taylor_ext(int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

// Electric wave N; NI
complex double qpms_vswf_single_el(int m, int n, sph_t kdlj, 
		qpms_bessel_t btyp, qpms_normalization_t norm);
// Magnetic wave M; NI
complex double qpms_vswf_single_mg(int m, int n, sph_t kdlj, 
		qpms_bessel_t btyp, qpms_normalization_t norm);


typedef struct qpms_trans_calculator {
        int lMax;
        size_t nelem;
        complex double **A_multipliers;
        complex double **B_multipliers;
        qpms_normalization_t normalization;
} qpms_trans_calculator;

qpms_trans_calculator *qpms_trans_calculator_init(int lMax, qpms_normalization_t nt);
void qpms_trans_calculator_free(qpms_trans_calculator *);

complex double qpms_trans_calculator_get_A(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J);
complex double qpms_trans_calculator_get_B(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_calculator_get_A_ext(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);
complex double qpms_trans_calculator_get_B_ext(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

int qpms_trans_calculator_get_AB_p(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		int m, int n, int mu, int nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J);
int qpms_trans_calculator_get_AB_p_ext(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

int qpms_trans_calculator_get_AB_arrays(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		size_t deststride, size_t srcstride,
		sph_t kdlj, bool r_ge_d, qpms_bessel_t J); 
int qpms_trans_calculator_get_AB_arrays_ext(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		size_t deststride, size_t srcstride,
		double kdlj_r, double kdlj_theta, double kdlj_phi,
		int r_ge_d, int J);

#ifdef QPMS_COMPILE_PYTHON_EXTENSIONS
#include <Python.h>
#include <numpy/npy_common.h>
int qpms_cython_trans_calculator_get_AB_arrays_loop(
                const qpms_trans_calculator *c, qpms_bessel_t J, const int resnd,
                int daxis, int saxis,
                char *A_data, const npy_intp *A_shape, const npy_intp *A_strides,
                char *B_data, const npy_intp *B_shape, const npy_intp *B_strides,
                const char *r_data, const npy_intp *r_shape, const npy_intp *r_strides,
                const char *theta_data, const npy_intp *theta_shape, const npy_intp *theta_strides,
                const char *phi_data, const npy_intp *phi_shape, const npy_intp *phi_strides,
                const char *r_ge_d_data, const npy_intp *r_ge_d_shape, const npy_intp *r_ge_d_strides);


#endif //QPMS_COMPILE_PYTHON_EXTENSIONS


#endif // QPMS_TRANSLATIONS_H
