/*! \file scatsystem.h
 * \brief Modern interface for finite lattice calculations, including symmetries.
 */
#ifndef QPMS_SCATSYSTEM_H
#define QPMS_SCATSYSTEM_H
#include "qpms_types.h"
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
qpms_tmatrix_t *qpms_tmatrix_interpolator_eval(const qpms_tmatrix_interpolator_t *interp, double freq);

/// Create a T-matrix interpolator from frequency and T-matrix arrays.
qpms_tmatrix_interpolator_t *qpms_tmatrix_interpolator_create(size_t n, ///< Number of freqs and T-matrices provided.
		const double *freqs, const qpms_tmatrix_t *tmatrices_array, ///< N.B. array of qpms_tmatrix_t, not pointers!
		const gsl_interp_type *iptype
	       //, bool copy_bspec ///< if true, copies its own copy of basis spec from the first T-matrix.
       	       /*, ...? */);


/// A particle, defined by its T-matrix INDEX and position, to be used in qpms_scatsys_t.
typedef struct qpms_particle_tid_t {
	// Does it make sense to ever use other than cartesian coords for this?
	cart3_t pos; ///< Particle position in cartesian coordinates.
	qpms_ss_tmi_t tmatrix_id; ///< T-matrix index
} qpms_particle_tid_t;

struct qpms_finite_group_t;

typedef qpms_gmi_t qpms_ss_orbit_pi_t; ///< Auxilliary type used in qpms_ss_orbit_type_t for labeling particles inside orbits.
typedef qpms_ss_tmi_t qpms_ss_oti_t; ///< Auxilliary type used for labeling orbit types.

/// Structure describing a particle's "orbit type" under symmetry group actions in a system.
/**
 * Each particle has its orbit with respect to a symmetry group of the system in which the particle lies,
 * i.e. a set of particles onto which the group operations map the original particle.
 * 
 * (TODO DOC improve the following explanation:)
 * Typically, there will be only a small number of different (T-matrix, particle 
 * <a href="https://en.wikipedia.org/wiki/Group_action#Fixed_points_and_stabilizer_subgroups">stabiliser</a>)
 * pairs in the system. We can group the particles accordingly, into the same "orbit types"
 * that will allow to do certain operations only once for each "orbit type", saving memory and (perhaps) time.
 *
 * Each particle will then have assigned:
 *   1. an orbit type,
 *   2. an ID inside that orbit.
 *
 *
 * TODO DOC how the process of assigning the particle IDs actually work, orbit type (non-)uniqueness.
 *
 *
 * Memory is managed by qpms_scatspec_t; qpms_ss_orbit_type_t does not own anything.
 *
 */
typedef struct qpms_ss_orbit_type_t {
	qpms_gmi_t size; ///< Size of the orbit (a divisor of the group order).
	/// Action of the group elements onto the "first" element in this orbit.
	/** Its size is sym->order and its values lie between 0 and \a this.size − 1.
	 *  
	 *  The corresponding stabilizer {\a g} of the i-th particle on the orbit
	 *  is given by action[i] = g.
	 *
	 */
	qpms_ss_orbit_pi_t *action;
	/// T-matrix IDs of the particles on this orbit (type).
	/**
	 * We save all the tmi's of the particles on the orbit here to make the number of array lookups
	 * and pointer dereferences constant.
	 *
	 * The size of this array is \a size.
	 */
	qpms_ss_tmi_t *tmatrices;
} qpms_ss_orbit_type_t;

/// Auxillary type used in qpms_scatsys_t that identifies the particle's orbit and its id inside that orbit.
typedef struct qpms_ss_particle_orbitinfo {
	qpms_ss_oti_t t; ///< Orbit type.
#define QPMS_SS_P_ORBITINFO_UNDEF (-1) ///< This labels that the particle has not yet been assigned to an orbit.
	qpms_ss_orbit_pi_t p; ///< "Order" of the particle inside that orbit type.
} qpms_ss_particle_orbitinfo_t;


typedef struct qpms_scatsys_t {
	// TODO does bspec belong here?
	qpms_tmatrix_t **tm; ///< T-matrices in the system
	qpms_ss_tmi_t tm_count; ///< Number of all different T-matrices
	qpms_ss_tmi_t tm_capacity; ///< Capacity of tm[].
	qpms_particle_tid_t *p; ///< Particles.
	qpms_ss_pi_t p_count; ///< Size of particles array.
	qpms_ss_pi_t p_capacity; ///< Capacity of p[].

	//TODO the index types do not need to be so big.
	const struct qpms_finite_group_t *sym; ///< Symmetry group of the array
	qpms_ss_pi_t *p_sym_map; ///< Which particles map onto which by the symmetry ops.
	///< p_sym_map[idi + pi * sym->order] gives the index of pi-th particle under the idi'th sym op.
	qpms_ss_tmi_t *tm_sym_map; ///< Which t-matrices map onto which by the symmetry ops. Lookup by tm_sum_map[idi + tmi * sym->order].

	qpms_ss_oti_t orbit_type_count;
	qpms_ss_orbit_type_t *orbit_types; ///< (Array length is \a orbit_type_count.)

	qpms_ss_particle_orbitinfo_t *p_orbitinfo; ///< Orbit type identification of each particle. (Array length is \a p_count.)

	size_t fecv_size; ///< Number of elements of a full excitation coefficient vector size. 
	//size_t *saecv_sizes; ///< NI. Number of elements of symmetry-adjusted coefficient vector sizes (order as in sym->irreps). 

	size_t *fecv_pstarts; ///< Indices of where pi'th particle's excitation coeffs start in a full excitation coefficient vector.
	//size_t **saecv_pstarts; ///< NI. Indices of where pi'th particle's excitation coeff start in a symmetry-adjusted e.c.v.
	///**< First index is irrep index as in sym->irreps, second index is particle index. */

	// TODO shifted origin of the symmetry group etc.
	// TODO some indices for fast operations here.
	// private
	
	// We keep the p_orbitinfo arrays in this chunk in order to avoid memory fragmentation
	char *otspace;
	char *otspace_end;
	double lenscale; // radius of the array, used as a relative tolerance measure
} qpms_scatsys_t;

/// Creates a new scatsys by applying a symmetry group, copying particles if needed.
/** In fact, it copies everything except the vswf set specs, so keep them alive until scatsys is destroyed.
 */
qpms_scatsys_t *qpms_scatsys_apply_symmetry(const qpms_scatsys_t *orig, const struct qpms_finite_group_t *sym);
/// Destroys the result of qpms_scatsys_apply_symmetry or qpms_scatsys_load.
void qpms_scatsys_free(qpms_scatsys_t *s);

/// NOT IMPLEMENTED Dumps a qpms_scatsys_t structure to a file.
qpms_errno_t qpms_scatsys_dump(qpms_scatsys_t *ss, char *path);

/// NOT IMPLEMENTED Reads a qpms_scatsys_t structure from a file.
qpms_scatsys_t *qpms_scatsys_load(char *path);



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
#endif //QPMS_SCATSYSTEM_H
