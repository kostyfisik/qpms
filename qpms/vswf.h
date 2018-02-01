#ifndef QPMS_VSWF_H
#define QPMS_VSWF_H
#include "qpms_types.h"
#include <gsl/gsl_sf_legendre.h>

// Electric wave N; NI
csphvec_t qpms_vswf_single_el(int m, int n, sph_t kdlj,
                qpms_bessel_t btyp, qpms_normalisation_t norm);
// Magnetic wave M; NI
csphvec_t qpms_vswf_single_mg(int m, int n, sph_t kdlj,
                qpms_bessel_t btyp, qpms_normalisation_t norm);

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

qpms_errno_t qpms_legendre_deriv_y_get(double **result, double **result_deriv, double x, qpms_l_t lMax, 
		gsl_sf_legendre_t lnorm, double csphase); // free() result and result_deriv yourself!
qpms_errno_t qpms_legendre_deriv_y_fill(double *where, double *where_deriv, double x, 
		qpms_l_t lMax, gsl_sf_legendre_t lnorm, double csphase); 

/* some of the result targets may be NULL */
qpms_errno_t qpms_vswf_fill(csphvec_t *resultL, csphvec_t *resultM, csphvec_t *resultN, qpms_l_t lMax, sph_t kdrj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);
// Should give the same results: for consistency checks
qpms_errno_t qpms_vswf_fill_alternative(csphvec_t *resultL, csphvec_t *resultM, csphvec_t *resultN, qpms_l_t lMax, sph_t kdrj,
		qpms_bessel_t btyp, qpms_normalisation_t norm);

qpms_errno_t qpms_vecspharm_fill(csphvec_t *const a1target, csphvec_t *const a2target, csphvec_t *const a3target,
		                qpms_l_t lMax, sph_t dir, qpms_normalisation_t norm);
qpms_errno_t qpms_vecspharm_dual_fill(csphvec_t *const a1target, csphvec_t *const a2target, csphvec_t *const a3target,
		                qpms_l_t lMax, sph_t dir, qpms_normalisation_t norm);

qpms_errno_t qpms_planewave2vswf_fill_cart(cart3_t wavedir, ccart3_t amplitude,
		complex double *targt_longcoeff, complex double *target_mgcoeff, complex double *target_elcoeff,
		qpms_l_t lMax, qpms_normalisation_t norm);
qpms_errno_t qpms_planewave2vswf_fill_sph(sph_t wavedir, csphvec_t amplitude,
		complex double *targt_longcoeff, complex double *target_mgcoeff, complex double *target_elcoeff,
		qpms_l_t lMax, qpms_normalisation_t norm);


csphvec_t qpms_eval_vswf(sph_t where,
		complex double *longcoeffs, complex double *mgcoeffs, complex double *elcoeffs,
		qpms_l_t lMax, qpms_bessel_t btyp, qpms_normalisation_t norm);


qpms_vswfset_sph_t *qpms_vswfset_make(qpms_l_t lMax, sph_t kdlj,
	qpms_bessel_t btyp, qpms_normalisation_t norm);//NI
void qpms_vswfset_sph_pfree(qpms_vswfset_sph_t *);//NI

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

#endif // QPMS_VSWF_H
