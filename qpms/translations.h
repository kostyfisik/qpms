/*! \file translations.h
 *  \brief VSWF translation operator.
 *
 * ### Argument conventions
 *
 * A single wave with indices mu, nu is re-expanded at kdlj into waves with indices m, n,
 * i.e. in the following functions, the first arguments over which one sums (multiplied
 * by the waves with new origin).
 *
 * HOWEVER, this means that if a field has an expansion with coeffs a(mu, nu)
 * at the original origin, with the basis at the new origin, the coeffs will be
 * a(m, n) = \sum_{mu,nu} A(m, n, mu, nu) a(mu, nu).
 *
 * With qpms_trans_calculator_get_AB_arrays_buf (and other functions from *AB_arrays*
 * family), one can choose the stride. And it seems that the former stride argument (now called
 * destride) and the latter (now called srcstride) are connected to (m,n) and (mu,nu) indices,
 * respectively. Seems consistent.
 *
 *
 * #### r_ge_d argument:
 *
 * If r_ge_d == true, the translation coefficients are calculated using regular bessel functions, 
 * regardless of what J argument is.
 *
 *
 */
#ifndef QPMS_TRANSLATIONS_H
#define QPMS_TRANSLATIONS_H
#include "vectors.h"
#include "qpms_types.h"
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

#if defined LATTICESUMS32 || defined LATTICESUMS31
#include "ewald.h"
#endif


complex double qpms_trans_single_A(qpms_normalisation_t norm, qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J);

complex double qpms_trans_single_B(qpms_normalisation_t norm, qpms_m_t m, qpms_l_t n, qpms_m_t mu, qpms_l_t nu, sph_t kdlj,
		bool r_ge_d, qpms_bessel_t J);

/// Structure holding the constant factors in normalisation operators.
/**
 * The VSWF translation operator elements are rather complicated linear
 * combinations of Bessel functions and spherical harmonics.
 * The preceding factors are rather complicated but need to be calculated
 * (for a single normalisation convention)
 * only once and then recycled during the operator evaluation for different
 * translations.
 *
 * This structure is initialised with qpms_trans_calculator_t() and
 * holds all the constant factors up to a truncation
 * degree \a lMax.
 *
 * The destructor function is qpms_trans_calculator_free().
 *
 * If Ewald sums are enabled at build, it also holds the constant
 * factors useful for lattice sums of translation operator.
 */
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

#if defined LATTICESUMS32 || defined LATTICESUMS31
	qpms_ewald3_constants_t *e3c;
#endif
	double *legendre0; // Zero-argument Legendre functions – this might go outside #ifdef in the end...
} qpms_trans_calculator;

/// Initialise a qpms_trans_calculator_t instance for a given convention and truncation degree.
qpms_trans_calculator *qpms_trans_calculator_init(qpms_l_t lMax, ///< Truncation degree.
	       	qpms_normalisation_t nt
		);
/// Destructor for qpms_trans_calculator_t.
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

// Convenience functions using VSWF base specs
qpms_errno_t qpms_trans_calculator_get_trans_array(const qpms_trans_calculator *c,
		complex double *target, 
		/// Must be destspec->lMax <= c-> lMax && destspec->norm == c->norm. 
		const qpms_vswf_set_spec_t *destspec, size_t deststride,
		/// Must be srcspec->lMax <= c-> lMax && srcspec->norm == c->norm. 
		const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
		sph_t kdlj, bool r_ge_d, qpms_bessel_t J);

/// Version with \a k and cartesian particle positions
/// and with automatic \a r_ge_d = `false`.
qpms_errno_t qpms_trans_calculator_get_trans_array_lc3p(
		const qpms_trans_calculator *c,
		complex double *target, 
		/// Must be destspec->lMax <= c-> lMax && destspec->norm == c->norm. 
		const qpms_vswf_set_spec_t *destspec, size_t deststride,
		/// Must be srcspec->lMax <= c-> lMax && srcspec->norm == c->norm. 
		const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
		double k, cart3_t destpos, cart3_t srcpos,
		qpms_bessel_t J
		/// Workspace has to be at least 2 * c->neleme**2 long
		);

#ifdef LATTICESUMS32
// for the time being there are only those relatively low-level quick&dirty functions
// according to what has been implemented from ewald.h; 
// TODO more high-level functions with more advanced lattice generators etc. (after
// the prerequisities from lattices2d.h are implememted)


int qpms_trans_calculator_get_AB_arrays_e32(const qpms_trans_calculator *c,
		complex double *Adest, double *Aerr,
		complex double *Bdest, double *Berr,
		const ptrdiff_t deststride, const ptrdiff_t srcstride,
		const double eta, const complex double k,
		cart2_t b1, cart2_t b2,
		const cart2_t beta,
		const cart2_t particle_shift,
		double maxR, double maxK
		);

int qpms_trans_calculator_get_AB_arrays_e32_e(const qpms_trans_calculator *c,
		complex double *Adest, double *Aerr,
		complex double *Bdest, double *Berr,
		const ptrdiff_t deststride, const ptrdiff_t srcstride,
		const double eta, const complex double k,
		cart2_t b1, cart2_t b2,
		const cart2_t beta,
		const cart2_t particle_shift,
		double maxR, double maxK,
		qpms_ewald_part parts
		);

// Convenience functions using VSWF base specs
qpms_errno_t qpms_trans_calculator_get_trans_array_e32(const qpms_trans_calculator *c,
		complex double *target, double *err,
		/// Must be destspec->lMax <= c-> lMax && destspec->norm == c->norm. 
		const qpms_vswf_set_spec_t *destspec, size_t deststride,
		/// Must be srcspec->lMax <= c-> lMax && srcspec->norm == c->norm. 
		const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
		const double eta, const complex double k,
		cart2_t b1, cart2_t b2,
		const cart2_t beta,
		const cart2_t particle_shift,
		double maxR, double maxK
		);

qpms_errno_t qpms_trans_calculator_get_trans_array_e32_e(const qpms_trans_calculator *c,
		complex double *target, double *err,
		/// Must be destspec->lMax <= c-> lMax && destspec->norm == c->norm. 
		const qpms_vswf_set_spec_t *destspec, size_t deststride,
		/// Must be srcspec->lMax <= c-> lMax && srcspec->norm == c->norm. 
		const qpms_vswf_set_spec_t *srcspec, size_t srcstride,
		const double eta, const complex double k,
		cart2_t b1, cart2_t b2,
		const cart2_t beta,
		const cart2_t particle_shift,
		double maxR, double maxK,
		qpms_ewald_part parts
		);

#endif //LATTICESUMS32

#ifdef LATTICESUMS31
// e31z means that the particles are positioned along the z-axis;
// their positions and K-values are then denoted by a single z-coordinate
int qpms_trans_calculator_get_AB_arrays_e31z_both_points_and_shift(const qpms_trans_calculator *c,
		complex double *Adest, double *Aerr,
		complex double *Bdest, double *Berr,
		const ptrdiff_t deststride, const ptrdiff_t srcstride,
		const double eta, const double k,
		const double unitcell_area, // just lattice period
		const size_t nRpoints, const cart2_t *Rpoints,
		const size_t nKpoints, const cart2_t *Kpoints,
		const double beta,
		const double particle_shift
		);
#endif




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
