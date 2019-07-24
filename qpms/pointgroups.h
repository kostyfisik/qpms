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
static inline qpms_gmi_t qpms_pg_order(qpms_pointgroup_class cls, ///< Point group class.
		qpms_gmi_t n ///< Number of rotations around main axis (only for finite axial groups).
		) {
	if (qpms_pg_is_finite_axial(cls))
		QPMS_ENSURE(n > 0, "n must be at least 1 for finite axial groups");
	switch(cls) {
	// Axial groups
        case QPMS_PGS_CN: 
		return n;
        case QPMS_PGS_S2N: 
        case QPMS_PGS_CNH: 
        case QPMS_PGS_CNV: 
        case QPMS_PGS_DN: 
		return 2*n;
        case QPMS_PGS_DND: 
        case QPMS_PGS_DNH: 
		return 4*n;

        // Remaining polyhedral groups
        case QPMS_PGS_T: 
		return 12;
        case QPMS_PGS_TD: 
        case QPMS_PGS_TH: 
        case QPMS_PGS_O:
		return 24;
        case QPMS_PGS_OH: 
	        return 48;	
        case QPMS_PGS_I:  
		return 60;
        case QPMS_PGS_IH: 
		return 120;

        // Continuous axial groups
        case QPMS_PGS_CINF: 
        case QPMS_PGS_CINFH: 
        case QPMS_PGS_CINFV: 
        case QPMS_PGS_DINF: 
        case QPMS_PGS_DINFH: 

        // Remaining continuous groups
        case QPMS_PGS_SO3:   
        case QPMS_PGS_O3:    
		return 0; // 0 is infinity :-)
	default:
		QPMS_WTF;
	}
}

/// Generates the canonical elements of a given 3D point group type.
/** Uses the canonical generators and DPS. */
qpms_irot3_t *qpms_pg_canonical_elems(
		qpms_irot3_t *target, ///< Target array (optional; if NULL, a new one is allocated)
		qpms_pointgroup_class cls, ///< Point group class.
		qpms_gmi_t ///< Number of rotations around \a z axis (only for axial group classes).
		);

/// Returns the number of canonical generators of a given 3D point group type.
/** TODO what does it do for infinite (Lie) groups? */
static inline qpms_gmi_t qpms_pg_genset_size(qpms_pointgroup_class cls, ///< Point group class.
		qpms_gmi_t n ///< Number of rotations around main axis (only for axial groups).
		) {
	if (qpms_pg_is_finite_axial(cls)) {
		QPMS_ENSURE(n > 0, "n must be at least 1 for finite axial groups");
		// n = 1 needs special care:
		if (n==1)
			switch(cls) {
				case QPMS_PGS_CN:  return 0; // triv.
				case QPMS_PGS_S2N: return 1; // Z_2
				case QPMS_PGS_CNH: return 1; // Dih_1
				case QPMS_PGS_CNV: return 1; // Dih_1
				case QPMS_PGS_DN:  return 1; // Dih_1
				case QPMS_PGS_DND: return 2; // Dih_2
				case QPMS_PGS_DNH: return 2; // Dih_1 x Dih_1
				default: QPMS_WTF;
			}
	}

	switch(cls) {
	// Axial groups
        case QPMS_PGS_CN:  return 1; // Z_n
        case QPMS_PGS_S2N: return 1; // Z_{2n}
        case QPMS_PGS_CNH: return 2; // Z_n x Dih_1
        case QPMS_PGS_CNV: return 2; // Dih_n
        case QPMS_PGS_DN:  return 2; // Dih_n
        case QPMS_PGS_DND: return 2; // Dih_2n
        case QPMS_PGS_DNH: return 3; // Dih_n x Dih_1

        // Remaining polyhedral groups
        case QPMS_PGS_T:   // return ???; // A_4
        case QPMS_PGS_TD:  // return 2; // S_4
        case QPMS_PGS_TH:  // A_4 x Z_2
        case QPMS_PGS_O:   // S_4
        case QPMS_PGS_OH:  // return 3; //  S_4 x Z_2
        case QPMS_PGS_I:   // A_5
        case QPMS_PGS_IH:  // A_5 x Z_2

        // Continuous axial groups
        case QPMS_PGS_CINF: 
        case QPMS_PGS_CINFH: 
        case QPMS_PGS_CINFV: 
        case QPMS_PGS_DINF: 
        case QPMS_PGS_DINFH: 

        // Remaining continuous groups
        case QPMS_PGS_SO3:   
        case QPMS_PGS_O3:    
		QPMS_NOT_IMPLEMENTED("Not yet implemented for this point group class.");
	default:
		QPMS_WTF;
	}
}

/// Fills an array of canonical generators of a given 3D point group type.
/** TODO what does it do for infinite (Lie) groups? */
static inline qpms_gmi_t qpms_pg_genset(qpms_pointgroup_class cls, ///< Point group class.
		qpms_gmi_t n, ///< Number of rotations around main axis (only for axial groups).
		qpms_irot3_t gen[] ///< Target generator array
		) {
	if (qpms_pg_is_finite_axial(cls)) {
		QPMS_ENSURE(n > 0, "n must be at least 1 for finite axial groups");
		// n = 1 needs special care:
		if (n==1)
			switch(cls) {
				case QPMS_PGS_CN:  
					return 0; // triv.
				case QPMS_PGS_S2N: 
					gen[0] = QPMS_IROT3_INVERSION;
					return 1; // Z_2
				case QPMS_PGS_CNH: 
					gen[0] = QPMS_IROT3_ZFLIP;
					return 1; // Dih_1
				case QPMS_PGS_CNV: 
					gen[0] = QPMS_IROT3_XFLIP;
					return 1; // Dih_1
				case QPMS_PGS_DN:
					gen[0] = QPMS_IROT3_XROT_PI; // CHECKME
					return 1; // Dih_1
				case QPMS_PGS_DND: 
					gen[0] = QPMS_IROT3_INVERSION;
					gen[1] = QPMS_IROT3_XFLIP;
					return 2; // Dih_2
				case QPMS_PGS_DNH: // CHECKME
					gen[0] = QPMS_IROT3_ZFLIP;
					gen[1] = QPMS_IROT3_XFLIP;
					return 2; // Dih_1 x Dih_1
				default: QPMS_WTF;
			}
	}

	switch(cls) {
		// Axial groups
		case QPMS_PGS_CN:
			gen[0] = qpms_irot3_zrot_Nk(n, 1);
			return 1; // Z_n
		case QPMS_PGS_S2N:
			gen[0].rot = qpms_quat_zrot_Nk(2*n, 1);
			gen[0].det = -1;
			return 1; // Z_{2n}
		case QPMS_PGS_CNH:
			gen[0] = qpms_irot3_zrot_Nk(n, 1);
			gen[1] = QPMS_IROT3_ZFLIP;
			return 2; // Z_n x Dih_1
		case QPMS_PGS_CNV: 
			gen[0] = qpms_irot3_zrot_Nk(n, 1);
			gen[1] = QPMS_IROT3_XFLIP;
			return 2; // Dih_n
		case QPMS_PGS_DN:
			gen[0] = qpms_irot3_zrot_Nk(n, 1);
			gen[1] = QPMS_IROT3_XROT_PI;
			return 2; // Dih_n
		case QPMS_PGS_DND: 
			gen[0].rot = qpms_quat_zrot_Nk(2*n, 1);
			gen[0].det = -1;
			gen[1] = QPMS_IROT3_XFLIP;
			return 2; // Dih_2n
		case QPMS_PGS_DNH:
			gen[0] = qpms_irot3_zrot_Nk(n, 1);
			gen[1] = QPMS_IROT3_ZFLIP;
			gen[2] = QPMS_IROT3_XFLIP;
			return 3; // Dih_n x Dih_1

				   // Remaining polyhedral groups
		case QPMS_PGS_T:   // return ???; // A_4
		case QPMS_PGS_TD:  // return 2; // S_4
		case QPMS_PGS_TH:  // A_4 x Z_2
		case QPMS_PGS_O:   // S_4
		case QPMS_PGS_OH:  // return 3; //  S_4 x Z_2
		case QPMS_PGS_I:   // A_5
		case QPMS_PGS_IH:  // A_5 x Z_2

				   // Continuous axial groups
		case QPMS_PGS_CINF: 
		case QPMS_PGS_CINFH: 
		case QPMS_PGS_CINFV: 
		case QPMS_PGS_DINF: 
		case QPMS_PGS_DINFH: 

				   // Remaining continuous groups
		case QPMS_PGS_SO3:   
		case QPMS_PGS_O3:    
				   QPMS_NOT_IMPLEMENTED("Not yet implemented for this point group class.");
		default:
		QPMS_WTF;
	}
}

#endif //POINTGROUPS_H
