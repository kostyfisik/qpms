/*! \file scatsystem.h
 * \brief Modern interface for finite lattice calculations, including symmetries.
 */
#ifndef QPMS_SCATSYSTEM_H
#define QPMS_SCATSYSTEM_H
#include "vswf.h"
#include <stdbool.h>
#include <gsl/gsl_spline.h>

/// A T-matrix.
/** In the future, I might rather use a more abstract approach in which T-matrix
 *  is a mapping (function) of the field expansion coefficients.
 *  So the interface might change.
 *  For now, let me stick to the square dense matrix representation.
 */
typedef struct qpms_tmatrix_t {
	/** \brief VSWF basis specification, NOT owned by qpms_tmatrix_t by default.
	 *
	 *  Usually not checked for meaningfulness by the functions (methods),
	 *  so the caller should take care that \a spec->ilist does not 
	 *  contain any duplicities and that for each wave with order \a m
	 *  there is also one with order \a −m.
	 */
	const qpms_vswf_set_spec_t *spec; 
	complex double *m; ///< Matrix elements in row-major order.
	bool owns_m; ///< Information wheter m shall be deallocated with qpms_tmatrix_free()
} qpms_tmatrix_t;

/// Returns a pointer to the beginning of the T-matrix row \a rowno.
static inline complex double *qpms_tmatrix_row(qpms_tmatrix_t *t, size_t rowno){
	return t->m + (t->spec->n * rowno);
}

/// Initialises a zero T-matrix.
qpms_tmatrix_t *qpms_tmatrix_init(const qpms_vswf_set_spec_t *bspec);

/// Copies a T-matrix, allocating new array for the T-matrix data.
qpms_tmatrix_t *qpms_tmatrix_copy(const qpms_tmatrix_t *T);

/// Destroys a T-matrix.
void qpms_tmatrix_free(qpms_tmatrix_t *t);

/// A particle, defined by its T-matrix and position.
typedef struct qpms_particle_t {
	// Does it make sense to ever use other than cartesian coords for this?
	cart3_t pos; ///< Particle position in cartesian coordinates.
	const qpms_tmatrix_t *tmatrix; ///< T-matrix; not owned by qpms_particle_t.
} qpms_particle_t;

/// Check T-matrix equality/similarity.
/** 
 *  This function actually checks for identical vswf specs.
 * TODO define constants with "default" atol, rtol for this function.
 */
bool qpms_tmatrix_isnear(const qpms_tmatrix_t *T1, const qpms_tmatrix_t *T2,
		const double rtol, const double atol);

/// Creates a T-matrix from another matrix and a symmetry operation.
/** The symmetry operation is expected to be a unitary (square) 
 *  matrix \a M with the same
 *  VSWF basis spec as the T-matrix (i.e. \a t->spec). The new T-matrix will then
 *  correspond to CHECKME \f[ T' = MTM^\dagger \f] 
 */
qpms_tmatrix_t *qpms_tmatrix_apply_symop(
		const qpms_tmatrix_t *T, ///< the original T-matrix
		const complex double *M ///< the symmetry op matrix in row-major format
		);

/// Symmetrizes a T-matrix with an involution symmetry operation.
/** The symmetry operation is expected to be a unitary (square) 
 *  matrix \a M with the same
 *  VSWF basis spec as the T-matrix (i.e. \a t->spec). The new T-matrix will then
 *  correspond to CHECKME \f[ T' = \frac{T + MTM^\dagger}{2} \f] 
 */
qpms_tmatrix_t *qpms_tmatrix_symmetrise_involution(
		const qpms_tmatrix_t *T, ///< the original T-matrix
		const complex double *M ///< the symmetry op matrix in row-major format
		);

/// Creates a \f$ C_\infty \f$ -symmetrized version of a T-matrix.
/**
 *  \f[ {T'}_{tlm}^{\tau\lambda\mu} = T_{tlm}^{\tau\lambda\mu} \delta_{m\mu} \f]
 */
qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_inf(
		const qpms_tmatrix_t *T ///< the original T-matrix
		);
/// Creates a \f$ C_N \f$ -symmetrized version of a T-matrix.
/**
 *  \f[ {T'}_{tlm}^{\tau\lambda\mu} = \begin{cases}
 *  T{}_{lm}^{\lambda\mu} & (m-\mu)\mod N=0\\
 *  0 & (m-\mu)\mod N\ne0
 *  \end{cases} . \f]
 */
qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_N(
		const qpms_tmatrix_t *T, ///< the original T-matrix
		int N ///< number of z-axis rotations in the group
		);

/// Symmetrizes a T-matrix with an involution symmetry operation, rewriting the original one.
/** The symmetry operation is expected to be a unitary (square) 
 *  matrix \a M with the same
 *  VSWF basis spec as the T-matrix (i.e. \a t->spec). The new T-matrix will then
 *  correspond to CHECKME \f[ T' = \frac{T + MTM^\dagger}{2} \f] 
 */
qpms_tmatrix_t *qpms_tmatrix_symmetrise_involution_inplace(
		qpms_tmatrix_t *T, ///< the original T-matrix
		const complex double *M ///< the symmetry op matrix in row-major format
		);

/// Creates a \f$ C_\infty \f$ -symmetrized version of a T-matrix, rewriting the original one.
/**
 *  \f[ {T'}_{tlm}^{\tau\lambda\mu} = T_{tlm}^{\tau\lambda\mu} \delta_{m\mu} \f]
 */
qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_inf_inplace(
		qpms_tmatrix_t *T ///< the original T-matrix
		);
/// Creates a \f$ C_N \f$ -symmetrized version of a T-matrix, rewriting the original one.
/**
 *  \f[ {T'}_{tlm}^{\tau\lambda\mu} = \begin{cases}
 *  T{}_{lm}^{\lambda\mu} & (m-\mu)\mod N=0\\
 *  0 & (m-\mu)\mod N\ne0
 *  \end{cases} . \f]
 */
qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_N_inplace(
		qpms_tmatrix_t *T, ///< the original T-matrix
		int N ///< number of z-axis rotations in the group
		);

/// NOT IMPLEMENTED Loads scuff-tmatrix generated files.
/** 
 * freqs, freqs_su, tmatrices_array and tmdata arrays are allocated by this function
 * and have to be freed by the caller after use.
 * The contents of tmatrices_array is NOT supposed to be freed element per element.
 */
qpms_errno_t qpms_load_scuff_tmatrix(
		const char *path, ///< Path to the TMatrix file
		const qpms_vswf_set_spec_t *bspec, ///< VSWF set spec
		size_t *n, ///< Number of successfully loaded t-matrices
		double **freqs, ///< Frequencies in SI units
		double **freqs_su, ///< Frequencies in SCUFF units (optional)
		qpms_tmatrix_t *tmatrices_array, ///< The resulting T-matrices.
		complex double **tmdata ///< The t-matrices raw contents
		);

/* Fuck this, include the whole <gsl/gsl_spline.h>
typedef struct gsl_spline gsl_spline; // Forward declaration for the interpolator struct
typedef struct gsl_interp_type gsl_interp_type;
extern const gsl_interp_type * gsl_interp_linear;
extern const gsl_interp_type * gsl_interp_polynomial;
extern const gsl_interp_type * gsl_interp_cspline;
extern const gsl_interp_type * gsl_interp_cspline_periodic;
extern const gsl_interp_type * gsl_interp_akima;
extern const gsl_interp_type * gsl_interp_akima_periodic;
extern const gsl_interp_type * gsl_interp_steffen;
*/

// struct gsl_interp_accel; // use if lookup shows to be too slow
typedef struct qpms_tmatrix_interpolator_t {
	const qpms_vswf_set_spec_t *bspec;
	//bool owns_bspec;
	gsl_spline **splines_real; ///< There will be a spline object for each nonzero element
	gsl_spline **splines_imag; ///< There will be a spline object for each nonzero element
	// gsl_interp_accel **accel_real;
	// gsl_interp_accel **accel_imag;
} qpms_tmatrix_interpolator_t;

/// Free a T-matrix interpolator.
void qpms_tmatrix_interpolator_free(qpms_tmatrix_interpolator_t *interp);

/// Evaluate a T-matrix interpolated value.
qpms_tmatrix_t *qpms_tmatrix_interpolator_eval(const qpms_tmatrix_interpolator_t *interp, double freq);

/// Create a T-matrix interpolator from frequency and T-matrix arrays.
qpms_tmatrix_interpolator_t *qpms_tmatrix_interpolator_create(size_t n, ///< Number of freqs and T-matrices provided.
		const double *freqs, const qpms_tmatrix_t *tmatrices_array, ///< N.B. array of qpms_tmatrix_t, not pointers!
		const gsl_interp_type *iptype
	       //, bool copy_bspec ///< if true, copies its own copy of basis spec from the first T-matrix.
       	       /*, ...? */);

#endif //QPMS_SCATSYSTEM_H
