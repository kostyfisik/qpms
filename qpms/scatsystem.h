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
#include <stdbool.h>


/// Overrides the number of threads spawned by the paralellized functions.
/** TODO MORE DOC which are those? */
void qpms_scatsystem_set_nthreads(long n);

/// A particle, defined by its T-matrix and position.
typedef struct qpms_particle_t {
	// Does it make sense to ever use other than cartesian coords for this?
	cart3_t pos; ///< Particle position in cartesian coordinates.
	const qpms_tmatrix_t *tmatrix; ///< T-matrix; not owned by qpms_particle_t.
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

/// Convenience function to access pi'th particle's bspec.
static inline const qpms_vswf_set_spec_t *qpms_ss_bspec_pi(const qpms_scatsys_t *ss, qpms_ss_pi_t pi) {
	return ss->tm[ss->p[pi].tmatrix_id]->spec;
}

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

/// Global translation matrix.
/**
 * The diagonal (particle self-) block are filled with zeros.
 * This may change in the future.
 */
complex double *qpms_scatsys_build_translation_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		double k ///< Wave number to use in the translation matrix.
		);

/// As qpms_scatsys_build_translation_full() but with choice of Bessel function type.
/** Might be useful for evaluation of cross sections and testing.
 */
complex double *qpms_scatsys_build_translation_matrix_e_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		double k, ///< Wave number to use in the translation matrix.
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
		double k, ///< Wave number to use in the translation matrix.
		qpms_bessel_t J
		);

complex double *qpms_scatsys_build_modeproblem_matrix_full(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss,
		/*COMPLEXIFY*/double k ///< Wave number to use in the translation matrix.
		);
/// Creates the mode problem matrix directly in the irrep-packed form.
complex double *qpms_scatsys_build_modeproblem_matrix_irrep_packed(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss, qpms_iri_t iri,
		/*COMPLEXIFY*/double k ///< Wave number to use in the translation matrix.
		);
/// Alternative implementation of qpms_scatsys_build_modeproblem_matrix_irrep_packed().
complex double *qpms_scatsys_build_modeproblem_matrix_irrep_packed_orbitorderR(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss, qpms_iri_t iri,
		/*COMPLEXIFY*/double k ///< Wave number to use in the translation matrix.
		);
/// Alternative implementation of qpms_scatsys_build_modeproblem_matrix_irrep_packed().
complex double *qpms_scatsys_build_modeproblem_matrix_irrep_packed_orbitorder_parallelR(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target,
		const qpms_scatsys_t *ss, qpms_iri_t iri,
		/*COMPLEXIFY*/double k ///< Wave number to use in the translation matrix.
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
