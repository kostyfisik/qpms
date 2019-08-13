/*! \file qpms_specfunc.h
 * \brief Various special and auxillary functions.
 */
#ifndef QPMS_SPECFUNC_H
#define QPMS_SPECFUNC_H
#include "qpms_types.h"
#include <gsl/gsl_sf_legendre.h>

/******************************************************************************
 *                    Spherical Bessel functions                              *
 ******************************************************************************/

// TODO unify types
qpms_errno_t qpms_sph_bessel_fill(qpms_bessel_t typ, qpms_l_t lmax, complex double x, complex double *result_array);



typedef struct {
        qpms_l_t lMax;
        double *akn; // coefficients as in DLMF 10.49.1
	//complex double *bkn; // coefficients of the derivatives
} qpms_sbessel_calculator_t;

qpms_sbessel_calculator_t *qpms_sbessel_calculator_init(void);
void qpms_sbessel_calculator_pfree(qpms_sbessel_calculator_t *c);

qpms_errno_t qpms_sbessel_calc_fill(qpms_sbessel_calculator_t *c, qpms_bessel_t typ, qpms_l_t lmax,
	double x, complex double *result_array);

complex double qpms_sbessel_calc_h1(qpms_sbessel_calculator_t *c, qpms_l_t n, complex double x);
qpms_errno_t qpms_sbessel_calc_h1_fill(qpms_sbessel_calculator_t *c, qpms_l_t lmax,
		complex double x, complex double *result_array);


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



/// Array of Legendre and and auxillary \f$\pi_{lm}, \tau_{lm} \f$ functions.
/**
 * See qpms_pitau_get() for definitions.
 *
 * The leg, pi, tau arrays are indexed using the standard qpms_mn2y() VSWF indexing.
 */
typedef struct {
	//qpms_normalisation_t norm;
	qpms_l_t lMax;
	//qpms_y_t nelem;
	double *leg, *pi, *tau;
} qpms_pitau_t;

/// Returns an array of normalised Legendre and auxillary \f$\pi_{lm}, \tau_{lm} \f$ functions.
/**
 * The normalised Legendre function here is defined as
 * \f[ 
 * 	\Fer[norm.]{l}{m} = \csphase^{-1} 
 * 		\sqrt{\frac{1}{l(l+1)}\frac{(l-m)!(2l+1)}{4\pi(l+m)!}},
 * \f] i.e. obtained using `gsl_sf_legendre_array_e()` with 
 * `norm = GSL_SF_LEGENDRE_SPHARM` and divided by \f$ \sqrt{l(l+1)} \f$.
 *
 * The auxillary functions are defined as
 * \f[
 * 	\pi_{lm}(\cos \theta) = \frac{m}{\sin \theta} \Fer[norm.]{l}{m}(\cos\theta),\\
 * 	\tau_{lm}(\cos \theta) = \frac{\ud}{\ud \theta} \Fer[norm.]{l}{m}(\cos\theta)
 * \f]
 * with appropriate limit expression used if \f$ \abs{\cos\theta} = 1 \f$.
 *
 * When done, don't forget to deallocate the memory using qpms_pitau_free().
 *
 */
qpms_pitau_t qpms_pitau_get(double theta, qpms_l_t lMax, double csphase);

/// Directly fills (pre-allocated) arrays of normalised Legendre and auxillary \f$\pi_{lm}, \tau_{lm} \f$ functions.
/**
 * Arrays must be preallocated for `lMax * (lMax + 2)` elements. `NULL` targets are skipped. 
 * For details, see qpms_pitau_get().
 */
qpms_errno_t qpms_pitau_fill(double *target_leg, double *target_pi, double *target_tau,
		double theta, qpms_l_t lMax, double csphase);

/// Frees the dynamically allocated arrays from qpms_pitau_t.
void qpms_pitau_free(qpms_pitau_t);

//void qpms_pitau_pfree(qpms_pitau_t*);//NI

// Associated Legendre polynomial at zero argument (DLMF 14.5.1) DEPRECATED?
double qpms_legendre0(int m, int n);
// Associated Legendre polynomial derivative at zero argument (DLMF 14.5.2)
double qpms_legendred0(int m, int n);

#endif // QPMS_SPECFUNC_H
