/* \file tmatrices.h
 * \brief T-matrices for scattering systems.
 */
#ifndef TMATRICES_H
#define TMATRICES_H
// #include "qpms_types.h" // included via materials.h
// #include <gsl/gsl_spline.h> // included via materials.h
#include "materials.h"
#include <stdio.h>

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

/// Tells whether qpms_load_scuff_tmatrix should crash if fopen() fails.
/** If true (default), the function causes the program
 * die e.g. when the tmatrix file
 * does not exists.
 *
 * If false, it does nothing and returns an appropriate error value instead.
 * This is desirable e.g. when used in Python (so that proper exception can
 * be thrown).
 */
extern bool qpms_load_scuff_tmatrix_crash_on_failure;

/// Loads a scuff-tmatrix generated file.
/** A simple wrapper over qpms_read_scuff_tmatrix() that needs a 
 * path instead of open FILE.
 *
 * The T-matrix is transformed from the VSWF basis defined by
 * QPMS_NORMALISATION_CONVENTION_SCUFF into the basis defined
 * by convention bspec->norm.
 * 
 * Right now, bspec->norm with real or "reversed complex" spherical
 * harmonics are not supported.
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

/// Application of T-matrix on a vector of incident field coefficients, \f$ f = Ta \f$.
complex double *qpms_apply_tmatrix(
		complex double *f_target, ///< Scattered field coefficient array of size T->spec->n; if NULL, a new one is allocated.
		const complex double *a, ///< Incident field coefficient array of size T->spec->n.
		const qpms_tmatrix_t *T
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

/// Fills an existing T-matrix with new interpolated values.
qpms_errno_t qpms_tmatrix_interpolator_eval_fill(qpms_tmatrix_t *target, ///< T-matrix to be updated, not NULL.
		const qpms_tmatrix_interpolator_t *interp, double freq);

/// Evaluate a T-matrix interpolated value.
/** The result is to be freed using qpms_tmatrix_free().*/
qpms_tmatrix_t *qpms_tmatrix_interpolator_eval(const qpms_tmatrix_interpolator_t *interp, double freq);

/// Create a T-matrix interpolator from frequency and T-matrix arrays.
qpms_tmatrix_interpolator_t *qpms_tmatrix_interpolator_create(size_t n, ///< Number of freqs and T-matrices provided.
		const double *freqs, const qpms_tmatrix_t *tmatrices_array, ///< N.B. array of qpms_tmatrix_t, not pointers!
		const gsl_interp_type *iptype
	       //, bool copy_bspec ///< if true, copies its own copy of basis spec from the first T-matrix.
       	       /*, ...? */);

/// Calculates the reflection Mie-Lorentz coefficients for a spherical particle.
/**
 * This function is based on the previous python implementation mie_coefficients() from qpms_p.py,
 * so any bugs therein should affect this function as well and perhaps vice versa.
 *
 * Most importantly, the case of magnetic material, \a mu_i != 0 or \a mu_e != 0 has never been tested
 * and might give wrong results.
 *
 * \return Array with the Mie-Lorentz reflection coefficients in the order determined by bspec.
 * If \a target was not NULL, this is target, otherwise a newly allocated array.
 *
 * TODO better doc.
 */
complex double *qpms_mie_coefficients_reflection(
                complex double *target, ///< Target array of length bspec->n. If NULL, a new one will be allocated.
                const qpms_vswf_set_spec_t *bspec, ///< Defines which of the coefficients are calculated.
                double a, ///< Radius of the sphere.
                complex double k_i, ///< Wave number of the internal material of the sphere.
                complex double k_e, ///< Wave number of the surrounding medium.
                complex double mu_i, ///< Relative permeability of the sphere material.
                complex double mu_e, ///< Relative permeability of the surrounding medium.
                qpms_bessel_t J_ext, ///< Kind of the "incoming" waves. Most likely QPMS_BESSEL_REGULAR.
                qpms_bessel_t J_scat ///< Kind of the "scattered" waves. Most likely QPMS_HANKEL_PLUS.
                );

/// Replaces the contents of an existing T-matrix with that of a spherical nanoparticle calculated using the Lorentz-mie theory.
qpms_errno_t qpms_tmatrix_spherical_fill(qpms_tmatrix_t *t, ///< T-matrix whose contents are to be replaced. Not NULL.
		double a, ///< Radius of the sphere.
		complex double k_i, ///< Wave number of the internal material of the sphere.
		complex double k_e, ///< Wave number of the surrounding medium.
		complex double mu_i, ///< Relative permeability of the sphere material.
		complex double mu_e ///< Relative permeability of the surrounding medium.
		);

/// Creates a new T-matrix of a spherical particle using the Lorentz-Mie theory.
static inline qpms_tmatrix_t *qpms_tmatrix_spherical(
		const qpms_vswf_set_spec_t *bspec,
		double a, ///< Radius of the sphere.
		complex double k_i, ///< Wave number of the internal material of the sphere.
		complex double k_e, ///< Wave number of the surrounding medium.
		complex double mu_i, ///< Relative permeability of the sphere material.
		complex double mu_e ///< Relative permeability of the surrounding medium.
		) {
	qpms_tmatrix_t *t = qpms_tmatrix_init(bspec);
	qpms_tmatrix_spherical_fill(t, a, k_i, k_e, mu_i, mu_e);
	return t;
}

/// Convenience function to calculate T-matrix of a non-magnetic spherical \
particle using the permittivity values, replacing existing T-matrix data.
qpms_errno_t qpms_tmatrix_spherical_mu0_fill(
		qpms_tmatrix_t *t, ///< T-matrix whose contents are to be replaced. Not NULL.
		double a, ///< Radius of the sphere.
		double omega, ///< Angular frequency.
		complex double epsilon_fg, ///< Relative permittivity of the sphere.
		complex double epsilon_bg ///< Relative permittivity of the background medium.
		); 

/// Convenience function to calculate T-matrix of a non-magnetic spherical particle using the permittivity values.
static inline qpms_tmatrix_t *qpms_tmatrix_spherical_mu0(
		const qpms_vswf_set_spec_t *bspec,
		double a, ///< Radius of the sphere.
		double omega, ///< Angular frequency.
		complex double epsilon_fg, ///< Relative permittivity of the sphere.
		complex double epsilon_bg ///< Relative permittivity of the background medium.
		) {
	qpms_tmatrix_t *t = qpms_tmatrix_init(bspec);
	qpms_tmatrix_spherical_mu0_fill(t, a, omega, epsilon_fg, epsilon_bg);
	return t;
};

/// Return value type for qpms_arc_function_t.
typedef struct qpms_arc_function_retval_t {
	double r; ///< Distance from the origin.
	double beta; ///< Angle between surface normal and radial direction.
} qpms_arc_function_retval_t;

/// Prototype for general parametrisation of \f$ C_\infty \f$-symmetric particle's surface.
typedef struct qpms_arc_function_t {
	/// Arc parametrisation function.
	/** TODO link to notes.
	 *
	 * Implemented by:
	 * qpms_arc_cylinder(),
	 * qpms_arc_sphere().
	 */
	qpms_arc_function_retval_t (*function) (
			double theta, ///< Polar angle from interval \f$ [0, \pi] \f$. 
			const void *params ///< Pointer to implementation specific parameters.
			);
	const void *params;
} qpms_arc_function_t;

/// Parameter structure for qpms_arc_cylinder().
typedef struct qpms_arc_cylinder_params_t {
	double R; ///< Cylinder radius.
	double h; ///< Cylinder height.
} qpms_arc_cylinder_params_t;

/// Arc parametrisation of cylindrical particle; for qpms_arc_function_t.
qpms_arc_function_retval_t qpms_arc_cylinder(double theta,
		const void *params; ///< Points to qpms_arc_cylinder_params_t
		);

/// Arc parametrisation of spherical particle; for qpms_arc_function_t.
/** Useful mostly only for benchmarks, as one can use the Mie-Lorentz solution. */
qpms_arc_function_retval_t qpms_arc_sphere(double theta,
		const void *R; ///< Points to double containing particle's radius.
		);


/// Replaces T-matrix contents with those of a particle with \f$ C_\infty \f$ symmetry.
qpms_errno_t qpms_tmatrix_axialsym_fill(
		qpms_tmatrix_t *t, ///< T-matrix whose contents are to be replaced. Not NULL.
		complex double omega, ///< Angular frequency.
		qpms_epsmu_t outside, ///< Optical properties of the outside medium.
		qpms_epsmu_t inside, ///< Optical properties of the particle's material.
		qpms_arc_function_t shape, ///< Particle surface parametrisation.
		/** If `lMax_extend > t->bspec->lMax`, then the internal \a Q, \a R matrices will be 
		 * trimmed at this larger lMax; the final T-matrix will then be trimmed
		 * according to bspec.
		 */
		qpms_l_t lMax_extend 
		);

/// Creates a new T-matrix of a particle with \f$ C_\infty \f$ symmetry.
static inline qpms_tmatrix_t *qpms_tmatrix_axialsym(
		const qpms_vswf_set_spec_t *bspec,
		complex double omega, ///< Angular frequency.
		qpms_epsmu_t outside, ///< Optical properties of the outside medium.
		qpms_epsmu_t inside, ///< Optical properties of the particle's material.
		qpms_arc_function_t shape, ///< Particle surface parametrisation.
		/// Internal lMax to improve precision of the result.
		/** If `lMax_extend > bspec->lMax`, then the internal \a Q, \a R matrices will be 
		 * trimmed at this larger lMax; the final T-matrix will then be trimmed
		 * according to bspec.
		 */
		qpms_l_t lMax_extend 
		) {
	qpms_tmatrix_t *t = qpms_tmatrix_init(bspec);
	qpms_tmatrix_axialsym_fill(t, omega, outside, inside, shape, lMax_extend);
	return t;
}



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
