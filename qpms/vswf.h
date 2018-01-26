#ifndef QPMS_VSWF_H
#define QPMS_VSWF_H
#include "qpms_types.h"
#include <gsl/gsl_sf_legendre.h>

// Electric wave N; NI
csphvec_t qpms_vswf_single_el(int m, int n, sph_t kdlj,
                qpms_bessel_t btyp, qpms_normalization_t norm);
// Magnetic wave M; NI
csphvec_t qpms_vswf_single_mg(int m, int n, sph_t kdlj,
                qpms_bessel_t btyp, qpms_normalization_t norm);

// Set of electric and magnetic VSWF in spherical coordinate basis
typedef struct {
	//qpms_normalisation_t norm;
	qpms_l_t lMax;
	//qpms_y_t nelem;
	//sph_t kdlj
	csphvec_t *el, *mg;
} qpms_vswfset_sph_t;


/*
 * N.B. for the norm definitions, see
 * https://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html
 * ( gsl/specfunc/legendre_source.c and 7.24.2 of gsl docs
 */

int qpms_legendre_deriv_y_get(double **result, double **result_deriv, double x, qpms_l_t lMax, 
		gsle_sf_legendre_t lnorm, double csphase);
int qpms_errno_t qpms_legendre_deriv_y_fill(double *where, double *where_deriv, double x, 
		qpms_l_t lMax, gsl_sf_legendre_t lnorm, double csphase); 

qpms_vswfset_sph_t *qpms_vswfset_make(qpms_l_t lMax, sph_t kdlj,
	qpms_bessel_t btyp, qpms_normalization_t norm);//NI
void qpms_vswfst_sph_pfree(qpms_vswfset_t *);//NI

double *qpms_legendre_y_get(double x, qpms_l_t lMax, qpms_normalisation_t norm);//NI
double *qpms_legendre0d_y_get(qpms_l_t lMax, qpms_normalization_t norm); //NI
double *qpms_legendre_plus1d_y_get(qpms_l_t lMax, qpms_normalization_t norm); //NI
double *qpms_legendre_minus1d_y_get(qpms_l_t lMax, qpms_normalization_t norm); //NI



// array of pi, tau auxillary function (see [1,(37)])
typedef struct {
	//qpms_normalization_t norm;
	qpms_l_t lMax;
	//qpms_y_t nelem;
	double *pi, *tau;
} qpms_pitau_t;
void qpms_pitau_free(qpms_pitau_t);//NI
void qpms_pitau_pfree(qpms_pitau_t*);//NI

#endif // QPMS_VSWF_H
