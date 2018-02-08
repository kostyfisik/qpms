#ifndef QPMS_TRANSLATIONS_H
#define QPMS_TRANSLATIONS_H
#include "vectors.h"
#include "qpms_types.h"
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

// TODO replace the xplicit "Taylor" functions with general,
// taking qpms_normalisation_t argument.
complex double qpms_trans_single_A_Taylor(qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_single_B_Taylor(qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_single_A_Taylor_ext(qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

complex double qpms_trans_single_B_Taylor_ext(qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

complex double qpms_trans_single_A(qpms_normalisation_t norm, qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_single_B(qpms_normalisation_t norm, qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
                bool r_ge_d, qpms_bessel_t J);

typedef struct qpms_trans_calculator {
        qpms_normalisation_t normalisation;
        qpms_l_t lMax;
        qpms_y_t nelem;
        complex double **A_multipliers;
        complex double **B_multipliers;
#if 0
	// Normalised values of the Legendre functions and derivatives
	// for θ == π/2, i.e. for the 2D case.
	double *leg0; 
	double *pi0;
	double *tau0;
	// Spherical Bessel function coefficients:
	// TODO
#endif
} qpms_trans_calculator;

qpms_trans_calculator *qpms_trans_calculator_init(qpms_l_t lMax, qpms_normalisation_t nt);
void qpms_trans_calculator_free(qpms_trans_calculator *);

complex double qpms_trans_calculator_get_A(const qpms_trans_calculator *c,
		qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J);
complex double qpms_trans_calculator_get_B(const qpms_trans_calculator *c,
		qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J);
int qpms_trans_calculator_get_AB_p(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J);
int qpms_trans_calculator_get_AB_arrays(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		size_t deststride, size_t srcstride,
		sph_t kdlj, bool r_ge_d, qpms_bessel_t J); 


// TODO update the types later
complex double qpms_trans_calculator_get_A_ext(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

complex double qpms_trans_calculator_get_B_ext(const qpms_trans_calculator *c,
		int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

int qpms_trans_calculator_get_AB_p_ext(const qpms_trans_calculator *c,
		complex double *Adest, complex double *Bdest,
		int m, int n, int mu, int nu, double kdlj_r,
		double kdlj_th, double kdlj_phi, int r_ge_d, int J);

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
