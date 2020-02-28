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
#include "vswf.h"
#include "tmatrices.h"
#include <stdbool.h>

/// Overrides the number of threads spawned by the paralellized functions.
/** TODO MORE DOC which are those? */
void qpms_scatsystem_set_nthreads(long n);

/// A particle, defined by its T-matrix and position.
/** This is rather only an auxillary intermediate structure to ultimately 
 * build an qpms_scatsys_t instance */
typedef struct qpms_particle_t {
	// Does it make sense to ever use other than cartesian coords for this?
	cart3_t pos; ///< Particle position in cartesian coordinates.
	const qpms_tmatrix_function_t *tmg; ///< T-matrix function; not owned by qpms_particle_t.
	qpms_tmatrix_operation_t op; ///< T-matrix transformation operation w.r.t. \a tmg.
} qpms_particle_t;

struct qpms_finite_group_t;
typedef struct qpms_finite_group_t qpms_finite_group_t;

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

/// Auxillary type used in qpms_scatsys_t: A recepy to create another T-matrices by symmetry operations.
typedef struct qpms_ss_derived_tmatrix_t {
	qpms_ss_tmgi_t tmgi; ///< Index of the corresponding qpms_scatsys_t::tm element.
	struct qpms_tmatrix_operation_t op; ///< Operation to derive this particular T-matrix.
} qpms_ss_derived_tmatrix_t;

typedef struct qpms_scatsys_periodic_info_t {	
	/// (Direct) lattice basis of the system (only \a lattice_dimension elements are used)
	/** This is mandatory for \a lattice_dimension != 0 */
	cart3_t lattice_basis[3];

	/// Reciprocal lattice basis. 
	/**(TODO specify whether it includes 2π or not) */
	cart3_t reciprocal_basis[3];

	/// Unitcell volume (irrelevant for non-periodic systems).
	/** The dimensionality of the volume corresponds to lattice_dimension, so
	 * for lattice_dimension == 1, this will actually be lenght and for
	 * lattice_dimension == 2, a 2D area.
	 */
	double unitcell_volume;

	/// Ewald parameter \f$ \eta \f$.
	double eta;
} qpms_scatsys_periodic_info_t;


struct qpms_trans_calculator;
struct qpms_epsmu_generator_t;

/// Common "class" for system of scatterers, both periodic and non-periodic.
/**
 * Infinite periodic structures (those with \a lattice_dimension > 0) 
 * have the \a per element allocated and filled.
 * These are ignored for finite systems (lattice_dimension == 0).
 */
typedef struct qpms_scatsys_t {
	/// Number of dimensions in which the system is periodic from the range 0–3.
	int lattice_dimension;
	struct qpms_epsmu_generator_t medium; ///< Optical properties of the background medium.
	
	/// (Template) T-matrix functions in the system.
	/** The qpms_abstract_tmatrix_t objects (onto which this array member point)
	 *  are NOT owned by this and must be kept alive for the whole lifetime
	 *  of all qpms_scatsys_t objects that are built upon them.
	 */
	qpms_tmatrix_function_t *tmg; 
	qpms_ss_tmgi_t tmg_count; ///< Number of all different original T-matrix generators in the system.

	/// All the different T-matrix functions in the system, including those derived from \a tmg elements by symmetries.
	qpms_ss_derived_tmatrix_t *tm; 
	qpms_ss_tmi_t tm_count; ///< Number of all different T-matrices in the system (length of tm[]).
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

	/// Periodic lattice metadata.
	qpms_scatsys_periodic_info_t per;
} qpms_scatsys_t;

/// Retrieve the bspec of \a tmi'th element of \a ss->tm.
static inline const qpms_vswf_set_spec_t *qpms_ss_bspec_tmi(const qpms_scatsys_t *ss, qpms_ss_tmi_t tmi) {
	return ss->tmg[ss->tm[tmi].tmgi].spec;
}

/// Retrieve the bspec of \a pi'th particle in \a ss->p.
static inline const qpms_vswf_set_spec_t *qpms_ss_bspec_pi(const qpms_scatsys_t *ss, qpms_ss_pi_t pi) {
	return ss->tmg[ss->tm[ss->p[pi].tmatrix_id].tmgi].spec;
}

typedef struct qpms_scatsys_at_omega_t {
	const qpms_scatsys_t *ss; ///< Parent scattering system.
	/// T-matrices from \a ss, evaluated at \a omega.
	/** The T-matrices are in the same order as in \a ss,
	 * i.e in the order corresponding to \a ss->tm.
	 */
	qpms_tmatrix_t **tm;
	complex double omega; ///< Angular frequency
	qpms_epsmu_t medium; ///< Background medium optical properties at the given frequency
	complex double wavenumber; ///< Background medium wavenumber
} qpms_scatsys_at_omega_t;


/// Creates a new scatsys by applying a symmetry group onto a "proto-scatsys", copying particles if needed.
/** In fact, it copies everything except the vswf set specs and qpms_abstract_tmatrix_t instances,
 * so keep them alive until scatsys is destroyed.
 *  
 *  The following fields must be filled in the "proto- scattering system" \a orig:
 *  * orig->lattice_dimension
 *  * orig->medium – The pointers are copied to the new qpms_scatsys_t instance; 
 *    the target qpms_abstract_tmatrix_t objects must be kept alive before all the resulting 
 *    qpms_scatsys_t instances are properly destroyed.
 *  * orig->tmg – The pointers are copied to the new qpms_scatsys_t instance; 
 *    the target qpms_abstract_tmatrix_t objects must be kept alive before all the resulting 
 *    qpms_scatsys_t instances are properly destroyed. The pointers from orig->tmg, however, are copied.
 *  * orig->tmg_count
 *  * orig->tm – Must be filled, although the operations will typically be identities
 *    (QPMS_TMATRIX_OPERATION_NOOP). N.B. these NOOPs might be replaced with some symmetrisation operation
 *    in the resulting "full" qpms_scatsys_t instance.
 *  * orig->tm_count
 *  * orig->p
 *  * orig->p_count
 *
 *  For periodic systems, the corresponding number of orig->per->lattice_basis[] elements
 *  must be filled as well.
 *
 *  For periodic systems, only trivial group is currently supported. Non-trivial 
 *  groups will cause undefined behaviour.
 *
 *  The resulting qpms_scatsys_t is obtained by actually evaluating the T-matrices
 *  at the given frequency \a omega and where applicable, these are compared
 *  by their values with given tolerances. The T-matrix generators are expected
 *  to preserve the point group symmetries for all frequencies.
 *
 *  This particular function will likely change in the future.
 *
 *  \returns An instance \a sso of qpms_scatsys_omega_t. Note that \a sso->ss
 *  must be saved by the caller before destroying \a sso 
 *  (with qpms_scatsys_at_omega_free(), and destroyed only afterwards with 
 *  qpms_scatsys_free() when not needed anymore.
 *  \a sso->ss can be reused for different frequency by a
 *  qpms_scatsys_at_omega() call.
 *
 */
qpms_scatsys_at_omega_t *qpms_scatsys_apply_symmetry(const qpms_scatsys_t *orig, const struct qpms_finite_group_t *sym,
		complex double omega, const struct qpms_tolerance_spec_t *tol);

/// Destroys the result of qpms_scatsys_apply_symmetry or qpms_scatsys_load.
void qpms_scatsys_free(qpms_scatsys_t *s);

/// Destroys a qpms_scatsys_at_omega_t.
/** Used on results of qpms_scatsys_apply_symmetry() and qpms_scatsys_at_omega(). */
void qpms_scatsys_at_omega_free(qpms_scatsys_at_omega_t *ssw);

/// Evaluates scattering system T-matrices at a given frequency.
/** Free the result using qpms_scatsys_at_omega_free() when done. */
qpms_scatsys_at_omega_t *qpms_scatsys_at_omega(const qpms_scatsys_t *ss,
		complex double omega);

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

/// Global translation matrix.
/**
 * The diagonal (particle self-) block are filled with zeros (even for regular Bessel waves).
 * This may change in the future.
 */
complex double *qpms_scatsys_build_translation_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		complex double k ///< Wave number to use in the translation matrix.
		);

/// Creates the full \f$ (I - WS) \f$ matrix of the periodic scattering system.
/**
 * \returns \a target on success, NULL on error.
 */
complex double *qpms_scatsysw_build_modeproblem_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_at_omega_t *ssw
		);

/// As qpms_scatsys_build_translation_full() but with choice of Bessel function type.
/** Might be useful for evaluation of cross sections and testing.
 */
complex double *qpms_scatsys_build_translation_matrix_e_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		complex double k, ///< Wave number to use in the translation matrix.
		qpms_bessel_t J
		);

/// Global translation matrix with selectable Bessel function, projected on an irrep.
/**
 * The diagonal (particle self-) blocks are currently filled with zeros.
 * This may change in the future.
 */
complex double *qpms_scatsys_build_translation_matrix_e_irrep_packed(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		qpms_iri_t iri,
		complex double k, ///< Wave number to use in the translation matrix.
		qpms_bessel_t J
		);

/// Creates the mode problem matrix \f$ (I - TS) \f$ directly in the irrep-packed form.
/**
 * \returns \a target on success, NULL on error.
 */
complex double *qpms_scatsysw_build_modeproblem_matrix_irrep_packed(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_at_omega_t *ssw,
		qpms_iri_t iri ///< Index of the irreducible representation in ssw->ss->sym
		);
/// Alternative implementation of qpms_scatsysw_build_modeproblem_matrix_irrep_packed().
/**
 * \returns \a target on success, NULL on error.
 */
complex double *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_orbitorderR(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_at_omega_t *ssw,
		qpms_iri_t iri ///< Index of the irreducible representation in ssw->ss->sym
		);
/// Alternative (serial reference) implementation of qpms_scatsysw_build_modeproblem_matrix_irrep_packed().
/**
 * \returns \a target on success, NULL on error.
 */
complex double *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_serial(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_at_omega_t *ssw,
		qpms_iri_t iri ///< Index of the irreducible representation in ssw->ss->sym
		);

struct qpms_scatsys_at_omega_k_t; // Defined below.
/// LU factorisation (LAPACKE_zgetrf) result holder.
typedef struct qpms_ss_LU {
	const qpms_scatsys_at_omega_t *ssw;
	const struct qpms_scatsys_at_omega_k_t *sswk; ///< Only for periodic systems, otherwise NULL.
	bool full; ///< true if full matrix; false if irrep-packed.
	qpms_iri_t iri; ///< Irrep index if `full == false`.
	/// LU decomposition array.
	complex double *a;
	/// Pivot index array, size at least max(1,min(m, n)).
	int *ipiv; 
} qpms_ss_LU;
void qpms_ss_LU_free(qpms_ss_LU);

/// Builds an LU-factorised mode/scattering problem \f$ (I - TS) \f$ matrix from scratch. Nonperiodic systems only.
qpms_ss_LU qpms_scatsysw_build_modeproblem_matrix_full_LU(
		complex double *target, ///< Pre-allocated target array. Optional (if NULL, new one is allocated).
		int *target_piv, ///< Pre-allocated pivot array. Optional (if NULL, new one is allocated).
		const qpms_scatsys_at_omega_t *ssw
		);

/// Builds an irrep-packed LU-factorised mode/scattering problem matrix from scratch.
qpms_ss_LU qpms_scatsysw_build_modeproblem_matrix_irrep_packed_LU(
		complex double *target, ///< Pre-allocated target array. Optional (if NULL, new one is allocated).
		int *target_piv, ///< Pre-allocated pivot array. Optional (if NULL, new one is allocated).
		const qpms_scatsys_at_omega_t *ssw,
		qpms_iri_t iri
		);

/// Computes LU factorisation of a pre-calculated mode/scattering problem matrix, replacing its contents.
qpms_ss_LU qpms_scatsysw_modeproblem_matrix_full_factorise(
		complex double *modeproblem_matrix_full, ///< Pre-calculated mode problem matrix (I-TS). Mandatory.
		int *target_piv, ///< Pre-allocated pivot array. Optional (if NULL, new one is allocated).
		const qpms_scatsys_at_omega_t *ssw, ///< Must be filled for non-periodic systems.
		const struct qpms_scatsys_at_omega_k_t *sswk ///< Must be filled for periodic systems, otherwise must be NULL.
		);

/// Computes LU factorisation of a pre-calculated irrep-packed mode/scattering problem matrix, replacing its contents.
qpms_ss_LU qpms_scatsysw_modeproblem_matrix_irrep_packed_factorise(
		complex double *modeproblem_matrix_irrep_packed, ///< Pre-calculated mode problem matrix (I-TS). Mandatory.
		int *target_piv, ///< Pre-allocated pivot array. Optional (if NULL, new one is allocated).
		const qpms_scatsys_at_omega_t *ssw,
		qpms_iri_t iri
		);

/// Solves a (possibly partial, irrep-packed) scattering problem \f$ (I-TS)f = Ta_\mathrm{inc} \f$ using a pre-factorised \f$ (I-TS) \f$.
complex double *qpms_scatsys_scatter_solve(
		complex double *target_f, ///< Target (full or irrep-packed, depending on `ludata.full`) array for \a f. If NULL, a new one is allocated.
		const complex double *a_inc, ///< Incident field expansion coefficient vector \a a (full or irrep-packed, depending on `ludata.full`).
		qpms_ss_LU ludata ///< Pre-factorised \f$ I - TS \f$ matrix data.
		);

// ======================= Periodic system -only related stuff =============================

/// Scattering system at a given frequency and k-vector. Used only with periodic systems.
/**
 * N.B. use as a stack variable now, but this might become heap-allocated in the future (with own con- and destructor)
 */
typedef struct qpms_scatsys_at_omega_k_t {
	const qpms_scatsys_at_omega_t *ssw;
	double k[3]; ///< The k-vector's cartesian coordinates.
} qpms_scatsys_at_omega_k_t;

/// Creates the full \f$ (I - WS) \f$ matrix of the periodic scattering system. 
/**
 * \returns \a target on success, NULL on error.
 */
complex double *qpms_scatsyswk_build_modeproblem_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_at_omega_k_t *sswk
		);

/// Global translation matrix.
complex double *qpms_scatsys_periodic_build_translation_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		complex double wavenumber, ///< Wave number to use in the translation matrix.
		const cart3_t *wavevector ///< Wavevector / pseudomomentum in cartesian coordinates.
		);

/// Global translation matrix. 
complex double *qpms_scatsyswk_build_translation_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_at_omega_k_t *sswk
		);


/// Builds an LU-factorised mode/scattering problem \f$ (I - TS) \f$ matrix from scratch. Periodic systems only.
qpms_ss_LU qpms_scatsyswk_build_modeproblem_matrix_full_LU(
		complex double *target, ///< Pre-allocated target array. Optional (if NULL, new one is allocated).
		int *target_piv, ///< Pre-allocated pivot array. Optional (if NULL, new one is allocated).
		const qpms_scatsys_at_omega_k_t *sswk
		);

/// Searches for periodic scattering system's eigenmodes using Beyn's algorithm.
/**
 * Currently, elliptical contour is used.
 *
 * TODO In the future, this will probably support irrep decomposition as well,
 * but for the case of periodic / small systems,
 * the bottleneck is the T-matrix and translation matrix evaluation
 * rather than the linear algebra.
 */
struct beyn_result_t *qpms_scatsys_periodic_find_eigenmodes(
		const qpms_scatsys_t *ss, 
		/// Wavevector in cartesian coordinates (must lie in the lattice plane).
		const double k[3],
		complex double omega_centre, ///< Center of the ellipse inside which the eigenfreqs are searched for.
		double omega_rr, ///< Real half-axis of the ellipse inside which the eigenfreqs are searched for.
		double omega_ri, ///< Imaginary half-axis of the ellipse inside which the eigenfreqs are searched for.
		size_t contour_npoints, ///< Number of elliptic contour discretisation points (preferably even number)
		double rank_tol, ///< (default: `1e-4`) TODO DOC.
		size_t rank_sel_min, ///< Minimum number of eigenvalue candidates, even if they don't pass \a rank_tol.
		double res_tol ///< (default: `0.0`) TODO DOC.
		);

// ======================= Periodic system -only related stuff end =========================


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


/// Creates an incident field vector in the full basis, given a function that evaluates the field expansions at points.
/** TODO detailed doc!
 * \returns target_full if target_full was not NULL, otherwise the newly allocated array. */
complex double *qpms_scatsys_incident_field_vector_full(
		/// Target array. If NULL, a new one is allocated.
		/** The length of the array is ss->fecv_size. */
		complex double *target_full,
		const qpms_scatsys_t *ss,
		qpms_incfield_t field_at_point,
		const void *args, ///< Pointer passed as the last argument to (*field_at_point)()
		bool add ///< If true, add to target_full; rewrite target_full if false.
		);

/// Applies T-matrices onto an incident field vector in the full basis.
complex double *qpms_scatsysw_apply_Tmatrices_full(
		complex double *target_full, /// Target vector array. If NULL, a new one is allocated.
		const complex double *inc_full, /// Incident field coefficient vector. Must not be NULL.
		const qpms_scatsys_at_omega_t *ssw
		);

struct beyn_result_t; // See beyn.h for full definition

/// Searches for finite scattering system's eigenmodes using Beyn's algorithm.
/**
 * Currently, elliptical contour is used.
 *
 * TODO In the future, this will probably support irrep decomposition as well,
 * but it does not make much sense for periodic / small systems, as in their
 * case the bottleneck is the T-matrix and translation matrix evaluation
 * rather than the linear algebra.
 */
struct beyn_result_t *qpms_scatsys_finite_find_eigenmodes(
		const qpms_scatsys_t *ss, 
		/// A valid irrep index to search only in one irrep, or QPMS_NO_IRREP for solving the full system.
		qpms_iri_t iri,
		complex double omega_centre, ///< Center of the ellipse inside which the eigenfreqs are searched for.
		double omega_rr, ///< Real half-axis of the ellipse inside which the eigenfreqs are searched for.
		double omega_ri, ///< Imaginary half-axis of the ellipse inside which the eigenfreqs are searched for.
		size_t contour_npoints, ///< Number of elliptic contour discretisation points (preferably even number)
		double rank_tol, ///< (default: `1e-4`) TODO DOC.
		size_t rank_sel_min, ///< Minimum number of eigenvalue candidates, even if they don't pass \a rank_tol.
		double res_tol ///< (default: `0.0`) TODO DOC.
		);

#if 0
/// Searches for scattering system's eigenmodes using Beyn's algorithm.
/**
 * Currently, elliptical contour is used.
 *
 * TODO In the future, this will probably support irrep decomposition as well,
 * but it does not make much sense for periodic / small systems, as in their
 * case the bottleneck is the T-matrix and translation matrix evaluation
 * rather than the linear algebra.
 */
struct beyn_result_t *qpms_scatsys_find_eigenmodes(
		const qpms_scatsys_t *ss,
		double eta, ///< Ewald sum parameter
		const double *beta_, ///< k-vector of corresponding dimensionality, NULL/ignored for finite system.
		complex double omega_centre, ///< Center of the ellipse inside which the eigenfreqs are searched for.
		double omega_rr, ///< Real half-axis of the ellipse inside which the eigenfreqs are searched for.
		double omega_ri, ///< Imaginary half-axis of the ellipse inside which the eigenfreqs are searched for.
		size_t contour_npoints, ///< Number of elliptic contour discretisation points (preferably even number)
		double rank_tol, ///< (default: `1e-4`) TODO DOC.
		size_t rank_sel_min, ///< Minimum number of eigenvalue candidates, even if they don't pass \a rank_tol.
		double res_tol ///< (default: `0.0`) TODO DOC.
		);
#endif


#if 0
/// Creates a (partial) incident field vector in the symmetry-adapted basis, given a function that evaluates the field expansions at points.
/** TODO detailed doc! */
complex double *qpms_scatsys_incident_field_vector_irrep_packed(
		/// Target array. If NULL, a new one is allocated.
		/** The length of the array is ss->fecv_size. */
		complex double *target_full,
		const qpms_scatsys_t *ss,
		const qpms_iri_t iri, ///< The index of given irreducible representation of ss->sym.
		qpms_incfield_t field_at_point,
		const void *args, ///< Pointer passed as the last argument to (*field_at_point)()
		bool add ///< If true, add to target_full; rewrite target_full if false.
		);
#endif

/// Evaluates scattered fields (corresponding to a given excitation vector) at a given point.
/**
 * By scattered field, one assumes a linear combination of positive-Hankel-type
 * spherical waves.
 *
 * \return Complex electric field at the point defined by \a where.
 */
ccart3_t qpms_scatsys_eval_E(const qpms_scatsys_t *ss,
		const complex double *coeff_vector, ///< A full-length excitation vector (outgoing wave coefficients).
		cart3_t where, ///< Evaluation point.
		complex double k ///< Wave number.
		);


#if 0
/** Evaluates partial scattered fields (corresponding to a given irrep-reduced excitation vector)
 * at a given point.
 *
 * \return Complex electric field at the point defined by \a where.
 */
ccart3_t qpms_scatsys_eval_E_irrep(const qpms_scatsys_t *ss,
		qpms_iri_t iri, ///< Irreducible representation
		const complex double *coeff_vector, ///< A reduced excitation vector, corresponding to \a iri.
		cart3_t where, ///< Evaluation point.
		complex double k ///< Wave number.
		);
#endif

#endif //QPMS_SCATSYSTEM_H
