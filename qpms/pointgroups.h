#ifndef POINTGROUPS_H
#define POINTGROUPS_H

#include "qpms_error.h"
#include "wigner.h"


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

/// Returns the order of a given 3D point group type.
/** For infinite groups returns 0. */
static inline size_t qpms_pg_order(qpms_pointgroup_class cls, ///< Point group class.
		size_t n ///< Number of rotations around main axis (only for finite axial groups).
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

/// Returns the number of canonical generators of a given 3D point group type.
/** TODO what does it do for infinite (Lie) groups? */
static inline size_t qpms_pg_genset_size(qpms_pointgroup_class cls, ///< Point group class.
		size_t n ///< Number of rotations around main axis (only for axial groups).
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
static inline size_t qpms_pg_genset(qpms_pointgroup_class cls, ///< Point group class.
		size_t n, ///< Number of rotations around main axis (only for axial groups).
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
