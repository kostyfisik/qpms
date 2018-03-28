#ifndef QPMS_SPECFUNC_H
#define QPMS_SPECFUNC_H
#include "qpms_types.h"
#include <gsl/gsl_sf_legendre.h>

/******************************************************************************
 *                    Spherical Bessel functions                              *
 ******************************************************************************/

// TODO unify types
qpms_errno_t qpms_sph_bessel_fill(qpms_bessel_t typ, qpms_l_t lmax, double x, complex double *result_array);



typedef struct {
        qpms_l_t lMax;
        double *akn; // coefficients as in DLMF 10.49.1
	//complex double *bkn; // coefficients of the derivatives
} qpms_sbessel_calculator_t;

qpms_sbessel_calculator_t *qpms_sbessel_calculator_init();
void qpms_sbessel_calculator_pfree(qpms_sbessel_calculator_t *c);

qpms_errno_t qpms_sbessel_calc_fill(qpms_sbessel_calculator_t *c, qpms_bessel_t typ, qpms_l_t lmax,
	double x, complex double *result_array);

complex double qpms_sbessel_calc_h1(qpms_sbessel_calculator_t *c, qpms_l_t n, double x);
qpms_errno_t qpms_sbessel_calc_h1_fill(qpms_sbessel_calculator_t *c, qpms_l_t lmax,
		double x, complex double *result_array);


/******************************************************************************
 *         Legendre functions and their "angular derivatives"                 *
 ******************************************************************************/

/*
 * N.B. for the norm definitions, see
 * https://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html
 * ( gsl/specfunc/legendre_source.c and 7.24.2 of gsl docs
 */

qpms_errno_t qpms_legendre_deriv_y_get(double **result, double **result_deriv, double x, qpms_l_t lMax, 
		gsl_sf_legendre_t lnorm, double csphase); // free() result and result_deriv yourself!
qpms_errno_t qpms_legendre_deriv_y_fill(double *where, double *where_deriv, double x, 
		qpms_l_t lMax, gsl_sf_legendre_t lnorm, double csphase); 


double *qpms_legendre_y_get(double x, qpms_l_t lMax, qpms_normalisation_t norm);//NI
double *qpms_legendre0d_y_get(qpms_l_t lMax, qpms_normalisation_t norm); //NI
double *qpms_legendre_plus1d_y_get(qpms_l_t lMax, qpms_normalisation_t norm); //NI
double *qpms_legendre_minus1d_y_get(qpms_l_t lMax, qpms_normalisation_t norm); //NI



// array of Legendre and pi, tau auxillary functions (see [1,(37)])
// This should handle correct evaluation for theta -> 0 and theta -> pi
typedef struct {
	//qpms_normalisation_t norm;
	qpms_l_t lMax;
	//qpms_y_t nelem;
	double *leg, *pi, *tau;
} qpms_pitau_t;
qpms_pitau_t qpms_pitau_get(double theta, qpms_l_t lMax, qpms_normalisation_t norm);
void qpms_pitau_free(qpms_pitau_t);//NI
void qpms_pitau_pfree(qpms_pitau_t*);//NI

// Associated Legendre polynomial at zero argument (DLMF 14.5.1) DEPRECATED?
double qpms_legendre0(int m, int n);
// Associated Legendre polynomial derivative at zero argument (DLMF 14.5.2)
double qpms_legendred0(int m, int n);

#endif // QPMS_SPECFUNC_H
