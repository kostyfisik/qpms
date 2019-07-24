/*! \file pointgroups.h
 * \brief Quaternion-represented 3D point groups.
 */

#ifndef POINTGROUPS_H
#define POINTGROUPS_H

#include "qpms_error.h"
#include "quaternions.h"


/// Returns true if the point group class belongs to one of the seven "axial" group series.
static inline _Bool qpms_pg_is_finite_axial(qpms_pointgroup_class cls) {
	switch(cls) {
        case QPMS_PGS_CN: 
        case QPMS_PGS_S2N: 
        case QPMS_PGS_CNH: 
        case QPMS_PGS_CNV: 
        case QPMS_PGS_DN: 
        case QPMS_PGS_DND: 
        case QPMS_PGS_DNH: 
		return true;
	default:
		return false;
	}
}

/// Absolute tolerance threshold used internally to consider two different `qpms_irot3_t` instances equal.
/**
 * Used by @ref qpms_pg_irot3_approx_cmp_v.
 * By default, set to @ref QPMS_QUAT_ATOL.
 * It should work fine if the point group orders stay reasonable (thousands or less).
 * Do not touch if unsure what you are doing.
 */
extern double qpms_pg_quat_cmp_atol;

/// An ordering of qpms_irot3_t.
int qpms_pg_irot3_cmp(const qpms_irot3_t *, const qpms_irot3_t *);
/// A `search.h` and `qsort()` compatible ordering of qpms_irot3_t.
int qpms_pg_irot3_cmp_v(const void *, const void *);
/// An ordering of qpms_irot3_t that considers close enough elements equal.
int qpms_pg_irot3_approx_cmp(const qpms_irot3_t *, const qpms_irot3_t *, 
		double atol ///< Absolute tolerance for the quaternion part difference.
		);
/// A `search.h` compatible ordering of qpms_irot3_t that considers close enough elements equal.
/** The tolerance is determined by global variable @ref qpms_pg_quat_cmp_atol.
 */
int qpms_pg_irot3_approx_cmp_v(const void *, const void *);

/// Returns the order of a given 3D point group type.
/** For infinite groups returns 0. */
qpms_gmi_t qpms_pg_order(qpms_pointgroup_class cls, ///< Point group class.
		qpms_gmi_t n ///< Number of rotations around main axis (only for finite axial groups).
                         );

/// Generates the canonical elements of a given 3D point group type.
/** Uses the canonical generators and DPS. */
qpms_irot3_t *qpms_pg_canonical_elems(
		qpms_irot3_t *target, ///< Target array (optional; if NULL, a new one is allocated)
		qpms_pointgroup_class cls, ///< Point group class.
		qpms_gmi_t ///< Number of rotations around \a z axis (only for axial group classes).
		);

/// Returns the number of canonical generators of a given 3D point group type.
/** TODO what does it do for infinite (Lie) groups? */
qpms_gmi_t qpms_pg_genset_size(qpms_pointgroup_class cls, ///< Point group class.
		qpms_gmi_t n ///< Number of rotations around main axis (only for axial groups).
		);

/// Fills an array of canonical generators of a given 3D point group type.
/** TODO what does it do for infinite (Lie) groups? */
qpms_gmi_t qpms_pg_genset(qpms_pointgroup_class cls, ///< Point group class.
		qpms_gmi_t n, ///< Number of rotations around main axis (only for axial groups).
		qpms_irot3_t gen[] ///< Target generator array
		);

/// Generates all elements of a given point group.
/** The order of elements corresponds to the one obtained from qpms_pg_canonical_elems().
 */
qpms_irot3_t *qpms_pg_elems(qpms_irot3_t *target, ///< Target array (optional; if NULL, a new one is allocated)
	       	qpms_pointgroup_t g ///< Specification of the point group.
		);

/// Checks whether \a a is subgroup of \a b (in a dirty general way).
_Bool qpms_pg_is_subgroup(qpms_pointgroup_t a, qpms_pointgroup_t b);



#endif //POINTGROUPS_H
