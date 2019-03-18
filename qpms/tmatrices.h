/* \file tmatrices.h
 * \brief T-matrices for scattering systems.
 */
#ifndef TMATRICES_H
#define TMATRICES_H
#include "qpms_types.h"
#include <gsl/gsl_spline.h>

struct qpms_finite_group_t;
typedef struct qpms_finite_group_t qpms_finite_group_t;

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

/// Check T-matrix equality/similarity.
/** 
 *  This function actually checks for identical vswf specs.
 * TODO define constants with "default" atol, rtol for this function.
 */
bool qpms_tmatrix_isclose(const qpms_tmatrix_t *T1, const qpms_tmatrix_t *T2,
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

/// Applies a symmetry operation onto a T-matrix, rewriting the original T-matrix data.
/** The symmetry operation is expected to be a unitary (square) 
 *  matrix \a M with the same
 *  VSWF basis spec as the T-matrix (i.e. \a t->spec). The new T-matrix will then
 *  correspond to CHECKME \f[ T' = MTM^\dagger \f] 
 */
qpms_tmatrix_t *qpms_tmatrix_apply_symop_inplace(
		qpms_tmatrix_t *T, ///< the original T-matrix
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

/// Reads an open scuff-tmatrix generated file.
/** 
 * \a *freqs, \a *freqs_su, \a *tmatrices_array and \a *tmdata 
 * arrays are allocated by this function
 * and have to be freed by the caller after use.
 * \a freqs_su and \a tmatrices_array can be NULL, in that case
 * the respective arrays are not filled nor allocated.
 *
 * The contents of tmatrices_array is NOT 
 * supposed to be freed element per element.
 *
 * TODO more checks and options regarding NANs etc.
 *
 */
qpms_errno_t qpms_load_scuff_tmatrix(
		const char *path, ///< Path to the TMatrix file
		const qpms_vswf_set_spec_t *bspec, ///< VSWF set spec
		size_t *n, ///< Number of successfully loaded t-matrices
		double **freqs, ///< Frequencies in SI units..
		double **freqs_su, ///< Frequencies in SCUFF units (optional).
		/// The resulting T-matrices (optional).
		qpms_tmatrix_t **tmatrices_array, 
		complex double **tmdata ///< The T-matrices raw contents
		);
/// Loads a scuff-tmatrix generated file.
/** A simple wrapper over qpms_read_scuff_tmatrix() that needs a 
 * path instead of open FILE.
 */
qpms_errno_t qpms_read_scuff_tmatrix(
		FILE *f, ///< An open stream with the T-matrix data.
		const qpms_vswf_set_spec_t *bspec, ///< VSWF set spec
		size_t *n, ///< Number of successfully loaded t-matrices
		double **freqs, ///< Frequencies in SI units.
		double **freqs_su, ///< Frequencies in SCUFF units (optional).
		/// The resulting T-matrices (optional).
		qpms_tmatrix_t **tmatrices_array, 
		/// The T-matrices raw contents. 
		/** The coefficient of outgoing wave defined by 
		 * \a bspec->ilist[desti] as a result of incoming wave
		 * \a bspec->ilist[srci] at frequency \a (*freqs)[fi]
		 * is accessed via
		 * (*tmdata)[bspec->n*bspec->n*fi + desti*bspec->n + srci].
		 */
		complex double ** tmdata 
		);

/// In-place application of point group elements on raw T-matrix data.
/** \a tmdata can be e.g. obtained by qpms_load_scuff_tmatrix().
 * The \a symops array should always contain all elements of a finite
 * point (sub)group, including the identity operation.
 *
 * TODO more doc.
 */
qpms_errno_t qpms_symmetrise_tmdata_irot3arr(
		complex double *tmdata, const size_t tmcount,
		const qpms_vswf_set_spec_t *bspec,
		size_t n_symops,
		const qpms_irot3_t *symops
		);

/// In-place application of a point group on raw T-matrix data.
/** This does the same as qpms_symmetrise_tmdata_irot3arr(),
 * but takes a valid finite point group as an argument.
 *
 * TODO more doc.
 */
qpms_errno_t qpms_symmetrise_tmdata_finite_group(
		complex double *tmdata, const size_t tmcount,
		const qpms_vswf_set_spec_t *bspec,
		const qpms_finite_group_t *pointgroup
		);

/// In-place application of point group elements on a T-matrix.
/** The \a symops array should always contain all elements of a finite
 * point (sub)group, including the identity operation.
 *
 * TODO more doc.
 */
qpms_tmatrix_t *qpms_tmatrix_symmetrise_irot3arr_inplace(
		qpms_tmatrix_t *T,
		size_t n_symops,
		const qpms_irot3_t *symops
		);

/// In-place application of point group elements on a T-matrix.
/** This does the same as qpms_tmatrix_symmetrise_irot3arr(),
 * but takes a valid finite point group as an argument.
 *
 * TODO more doc.
 */
qpms_tmatrix_t *qpms_tmatrix_symmetrise_finite_group_inplace(
		qpms_tmatrix_t *T,
		const qpms_finite_group_t *pointgroup
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

// struct gsl_interp_accel; // use if lookup proves to be too slow
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
/** The result is to be freed using qpms_tmatrix_free().*/
qpms_tmatrix_t *qpms_tmatrix_interpolator_eval(const qpms_tmatrix_interpolator_t *interp, double freq);

/// Create a T-matrix interpolator from frequency and T-matrix arrays.
qpms_tmatrix_interpolator_t *qpms_tmatrix_interpolator_create(size_t n, ///< Number of freqs and T-matrices provided.
		const double *freqs, const qpms_tmatrix_t *tmatrices_array, ///< N.B. array of qpms_tmatrix_t, not pointers!
		const gsl_interp_type *iptype
	       //, bool copy_bspec ///< if true, copies its own copy of basis spec from the first T-matrix.
       	       /*, ...? */);


#if 0
// Abstract types that describe T-matrix/particle/scatsystem symmetries
// To be implemented later. See also the thoughts in the beginning of groups.h.

typedef *char qpms_tmatrix_id_t; ///< Maybe I want some usual integer type instead.

///Abstract T-matrix type draft.
/**
 * TODO.
 */
typedef struct qpms_abstract_tmatrix_t{
	qpms_tmatrix_id_t id;
	/// Generators of the discrete point group under which T-matrix is invariant.
	qpms_irot3_t *invar_gens;
	/// Length of invar_gens.
	qpms_gmi_t invar_gens_size;

} qpms_abstract_tmatrix_t;


typedef struct qpms_abstract_particle_t{
} qpms_abstract_particle_t;

/// An abstract particle, defined by its position and abstract T-matrix.
typedef struct qpms_abstract_particle_t {
	cart3_t pos; ///< Particle position in cartesian coordinates.
	const qpms_abstract_tmatrix_t *tmatrix; ///< T-matrix; not owned by this.
} qpms_abstract_particle_t;


/** This is just an alias, as the same index can be used for 
 * abstract T-matrices as well.
 */
typedef qpms_particle_tid_t qpms_abstract_particle_tid_t;
#endif // 0

#endif //TMATRICES_H
