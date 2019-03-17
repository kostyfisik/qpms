/*! \file scatsystem.h
 * \brief Modern interface for finite lattice calculations, including symmetries.
 *
 * N.B. Only "reasonably normalised" waves are supported now in most of the
 * functions defined here, i.e. those that can be rotated by the usual
 * Wigner matrices, i.e. the "power" or "spharm" -normalised ones.
 * 
 * TODO FIXME check whether Condon-Shortley phase can have some nasty influence
 * here; I fear that yes.
 */
#ifndef QPMS_SCATSYSTEM_H
#define QPMS_SCATSYSTEM_H
#include "qpms_types.h"
#include <stdbool.h>
#include <gsl/gsl_spline.h>
#include <stdio.h> // only because of qpms_read_scuff_tmatrix()

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


/// A particle, defined by its T-matrix INDEX and position, to be used in qpms_scatsys_t.
typedef struct qpms_particle_tid_t {
	// Does it make sense to ever use other than cartesian coords for this?
	cart3_t pos; ///< Particle position in cartesian coordinates.
	qpms_ss_tmi_t tmatrix_id; ///< T-matrix index
} qpms_particle_tid_t;


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
	qpms_ss_orbit_pi_t size; ///< Size of the orbit (a divisor of the group order), i.e. number of particles on the orbit.
	size_t bspecn; ///< Individual particle's coefficient vector length. The same as ss->tm[this.tmatrices[0]]->spec->n.
	/// Action of the group elements onto the elements in this orbit.
	/** Its size is sym->order * this.size
	 *  and its values lie between 0 and \a this.size − 1.
	 *
	 *  Action of the group element g onto the pi-th particle
	 *  is given by action[g + pi*sym->order].
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
	/// Sizes of the per-orbit irrep bases. 
	/**
	 * The order of the irreps corresponds to the order in \a ss->sym->irreps.
	 * The size of this array is (obviously) \a ss->sym->nirreps.
	 *
	 * TODO different type? 
	 * TODO doc.
	 */
	size_t *irbase_sizes;
	//The following are pretty redundant, TODO reduce them at some point.
	/// Cumulative sums of irbase_sizes.
	size_t *irbase_cumsizes;
	/// TODO doc.
	size_t *irbase_offsets;
	/// Per-orbit irreducible representation orthonormal bases.
	/** This also defines the unitary operator that transforms the orbital excitation coefficients
	 * in the symmetry-adapted basis.
	 *
	 * The size is (\a this->size * \a this->tmatrices[0].spec->n)**2.
	 *
	 * TODO doc.
	 */
	complex double *irbases;
	/// TODO doc.
	size_t instance_count;
	/// Cumulative sum of the preceding ot->siza * ot->instance_count;
	qpms_ss_pi_t p_offset;
} qpms_ss_orbit_type_t;

typedef ptrdiff_t qpms_ss_osn_t; ///< "serial number" of av orbit in a given type.

/// Auxillary type used in qpms_scatsys_t that identifies the particle's orbit and its id inside that orbit.
typedef struct qpms_ss_particle_orbitinfo {
	qpms_ss_oti_t t; ///< Orbit type.
#define QPMS_SS_P_ORBITINFO_UNDEF (-1) ///< This labels that the particle has not yet been assigned to an orbit.
	qpms_ss_osn_t osn; ///< "Serial number" of the orbit in the given type. TODO type and more doc.
	qpms_ss_orbit_pi_t p; ///< Order (sija, ei rankki) of the particle inside that orbit type.
} qpms_ss_particle_orbitinfo_t;

struct qpms_trans_calculator;

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
	size_t *saecv_sizes; ///< Number of elements of symmetry-adjusted coefficient vector sizes (order as in sym->irreps).

	size_t *fecv_pstarts; ///< Indices of where pi'th particle's excitation coeffs start in a full excitation coefficient vector.
	size_t *saecv_ot_offsets; ///< TODO DOC. In the packed vector, saecv_ot_offsets[iri * orbit_type_count + oti] indicates start of ot
	/**< TODO maybe move it to qpms_ss_orbit_type_t, ffs. */
	//size_t **saecv_pstarts; ///< NI. Indices of where pi'th particle's excitation coeff start in a symmetry-adjusted e.c.v.
	///**< First index is irrep index as in sym->irreps, second index is particle index. */

	// TODO shifted origin of the symmetry group etc.
	// TODO some indices for fast operations here.
	// private
	size_t max_bspecn; ///< Maximum tm[...]->spec->n. Mainly for workspace allocation.
	
	/// Particles grouped by orbit (in the order corresponding to the packed memory layout).
	qpms_ss_pi_t *p_by_orbit;

	// We keep the p_orbitinfo arrays in this chunk in order to avoid memory fragmentation
	char *otspace;
	char *otspace_end;
	double lenscale; // radius of the array, used as a relative tolerance measure
	struct qpms_trans_calculator *c;
} qpms_scatsys_t;

/// Creates a new scatsys by applying a symmetry group, copying particles if needed.
/** In fact, it copies everything except the vswf set specs, so keep them alive until scatsys is destroyed.
 * The following fields must be filled:
 * orig->tm
 * orig->tm_count
 * orig->p
 * orig->p_count
 */
qpms_scatsys_t *qpms_scatsys_apply_symmetry(const qpms_scatsys_t *orig, const struct qpms_finite_group_t *sym);
/// Destroys the result of qpms_scatsys_apply_symmetry or qpms_scatsys_load.
void qpms_scatsys_free(qpms_scatsys_t *s);

/// Creates a "full" transformation matrix U that takes a full vector and projects it onto an symmetry adapted basis.
/** Mostly as a reference and a debugging tool, as multiplicating these big matrices would be inefficient.
 * 
 * TODO doc about shape etc.
 */
complex double *qpms_scatsys_irrep_transform_matrix(complex double *target_U,
		const qpms_scatsys_t *ss, qpms_iri_t iri);

/// Projects a "big" matrix onto an irrep (slow reference implementation).
/** TODO doc */
complex double *qpms_scatsys_irrep_pack_matrix_stupid(complex double *target_packed,
		const complex double *orig_full, const qpms_scatsys_t *ss,
		qpms_iri_t iri);

/// Transforms a big "packed" matrix into the full basis (slow reference implementation).
/** TODO doc */
complex double *qpms_scatsys_irrep_unpack_matrix_stupid(complex double *target_full,
		const complex double *orig_packed, const qpms_scatsys_t *ss,
		qpms_iri_t iri, bool add);

/// Projects a "big" matrix onto an irrep.
/** TODO doc */
complex double *qpms_scatsys_irrep_pack_matrix(complex double *target_packed,
		const complex double *orig_full, const qpms_scatsys_t *ss,
		qpms_iri_t iri);

/// Transforms a big "packed" matrix into the full basis.
/** TODO doc */
complex double *qpms_scatsys_irrep_unpack_matrix(complex double *target_full,
		const complex double *orig_packed, const qpms_scatsys_t *ss,
		qpms_iri_t iri, bool add);

/// Projects a "big" vector onto an irrep.
/** TODO doc */
complex double *qpms_scatsys_irrep_pack_vector(complex double *target_packed,
		const complex double *orig_full, const qpms_scatsys_t *ss,
		qpms_iri_t iri);

/// Transforms a big "packed" vector into the full basis.
/** TODO doc */
complex double *qpms_scatsys_irrep_unpack_vector(complex double *target_full,
		const complex double *orig_packed, const qpms_scatsys_t *ss,
		qpms_iri_t iri, bool add);

complex double *qpms_scatsys_build_translation_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		double k ///< Wave number to use in the translation matrix.
		);
complex double *qpms_scatsys_build_modeproblem_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		double k ///< Wave number to use in the translation matrix.
		);
/// Creates the mode problem matrix directly in the irrep-packed form.
complex double *qpms_scatsys_build_modeproblem_matrix_irrep_packed(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss, qpms_iri_t iri,
		double k ///< Wave number to use in the translation matrix.
		);
/// Alternative implementation of qpms_scatsys_build_modeproblem_matrix_irrep_packed().
complex double *qpms_scatsys_build_modeproblem_matrix_irrep_packed_orbitorderR(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss, qpms_iri_t iri,
		double k ///< Wave number to use in the translation matrix.
		);

/// NOT IMPLEMENTED Dumps a qpms_scatsys_t structure to a file.
qpms_errno_t qpms_scatsys_dump(qpms_scatsys_t *ss, char *path);

/// NOT IMPLEMENTED Reads a qpms_scatsys_t structure from a file.
qpms_scatsys_t *qpms_scatsys_load(char *path);

struct qpms_finite_group_t;

/// Constructs a "full matrix action" of a point group element for an orbit type.
/** TODO detailed doc */
complex double *qpms_orbit_action_matrix(
		/// Target array. If NULL, a new one is allocated.
		/** The size of the array is (orbit->size * bspec->n)**2
		 * (it makes sense to assume all the T-matrices share their spec).
		 */
		complex double *target,
		/// The orbit (type).
		const qpms_ss_orbit_type_t *orbit,
		/// Base spec of the t-matrices (we don't know it from orbit, as it has 
		/// only T-matrix indices.
		const qpms_vswf_set_spec_t *bspec,
		/// The symmetry group used to generate the orbit (must have rep3d filled).
		const struct qpms_finite_group_t *sym,
		/// The index of the operation in sym to represent.
		const qpms_gmi_t g);

/// Constructs a dense matrix representation of a irrep projector for an orbit type.
/** TODO detailed doc */
complex double *qpms_orbit_irrep_projector_matrix(
		/// Target array. If NULL, a new one is allocated.
		/** The size of the array is (orbit->size * bspec->n)**2
		 * (it makes sense to assume all the T-matrices share their spec).
		 */
		complex double *target,
		/// The orbit (type).
		const qpms_ss_orbit_type_t *orbit,
		/// Base spec of the t-matrices (we don't know it from orbit, as it has 
		/// only T-matrix indices.
		const qpms_vswf_set_spec_t *bspec,
		/// The symmetry group used to generate the orbit (must have rep3d filled).
		const struct qpms_finite_group_t *sym,
		/// The index of the irreducible representation of sym.
		const qpms_iri_t iri);

/// TODO DOC!!!!!
complex double *qpms_orbit_irrep_basis(
		/// Here theh size of theh basis shall be saved,
		size_t *basis_size,
		/// Target array. If NULL, a new one is allocated.
		/** The size of the array is basis_size * (orbit->size * bspec->n)
		 * (it makes sense to assume all the T-matrices share their spec).
		 */
		complex double *target,
		/// The orbit (type).
		const qpms_ss_orbit_type_t *orbit,
		/// Base spec of the t-matrices (we don't know it from orbit, as it has 
		/// only T-matrix indices.
		const qpms_vswf_set_spec_t *bspec,
		/// The symmetry group used to generate the orbit (must have rep3d filled).
		const struct qpms_finite_group_t *sym,
		/// The index of the irreducible representation of sym.
		const qpms_iri_t iri);



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
