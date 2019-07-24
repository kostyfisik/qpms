#include "pointgroups.h"
#include <search.h>
#include <stdlib.h>

#define PAIRCMP(a, b) {\
  if ((a) < (b)) return -1;\
  if ((a) > (b)) return 1;\
}

int qpms_pg_irot3_cmp(const qpms_irot3_t *a, const qpms_irot3_t *b) {
  PAIRCMP(a->det, b->det);
  PAIRCMP(creal(a->rot.a), creal(b->rot.a));
  PAIRCMP(cimag(a->rot.a), cimag(b->rot.a));
  PAIRCMP(creal(a->rot.b), creal(b->rot.b));
  PAIRCMP(cimag(a->rot.b), cimag(b->rot.b));
  return 0;
}

int qpms_pg_irot3_cmp_v(const void *av, const void *bv) {
  const qpms_irot3_t *a = av, *b = bv;
  return qpms_pg_irot3_cmp(a, b);
}

int qpms_pg_irot3_approx_cmp(const qpms_irot3_t *a, const qpms_irot3_t *b, double atol) {
  if (qpms_irot3_isclose(*a, *b, atol)) return 0;
  else return qpms_pg_irot3_cmp(a, b);
}

double qpms_pg_quat_cmp_atol = QPMS_QUAT_ATOL;

int qpms_pg_irot3_approx_cmp_v(const void *av, const void *bv) {
  const qpms_irot3_t *a = av, *b = bv;
  return qpms_pg_irot3_approx_cmp(a, b, qpms_pg_quat_cmp_atol);
}


/// Generates the canonical elements of a given 3D point group type.
qpms_irot3_t *qpms_pg_canonical_elems(qpms_irot3_t *target,
    qpms_pointgroup_class cls, const qpms_gmi_t then) {
  QPMS_UNTESTED;
  qpms_gmi_t order = qpms_pg_order(cls, then);
  QPMS_ENSURE(order, "Cannot generate an infinite group!");
  if (!target) QPMS_CRASHING_MALLOC(target, order * sizeof(qpms_irot3_t));
  target[0] = QPMS_IROT3_IDENTITY;
  qpms_gmi_t ngen = qpms_pg_genset_size(cls, then);
  qpms_irot3_t gens[ngen]; 
  (void) qpms_pg_genset(cls, then, gens);

  // Let's try it with a binary search tree, as an exercise :)
  qpms_irot3_t tree[order];
  void *root = NULL;
  // put the starting element (identity) to the tree
  (void) tsearch((void *) target, &root, qpms_pg_irot3_approx_cmp_v);
  qpms_gmi_t n = 1; // No. of generated elements.

  // And let's do the DFS without recursion; the "stack size" here (=order) might be excessive, but whatever
  qpms_gmi_t gistack[order], //< generator indices (related to gens[])
             srcstack[order], //< pre-image indices (related to target[])
             si = 0; //< stack index
  gistack[0] = 0;
  srcstack[0] = 0;

  while (si >= 0) { // DFS
    if (gistack[si] < ngen) { // are there generators left at this level? If so, try another elem
      if (n >= order) QPMS_WTF; // TODO some error message
      target[n] = qpms_irot3_mult(gens[gistack[si]], target[srcstack[si]]);
      if (tfind((void *) &(target[n]), &root, qpms_pg_irot3_approx_cmp_v))
        // elem found, try it with another gen in the next iteration
        gistack[si]++;
      else {
        // elem not found, add it to the tree and proceed to next level
        (void) tsearch( &(target[n]), &root, qpms_pg_irot3_approx_cmp_v);
        ++si;
        gistack[si] = 0;
        srcstack[si] = n;
        ++n;
      }
    } else { // no generators left at this level, get to the previous level
      --si; 
      if (si >= 0) ++gistack[si];
    }
  }

  QPMS_ENSURE(n == order, "Point group generation failure "
      "(assumed group order = %d, got %d; qpms_pg_quat_cmp_atol = %g)",
      order, n, qpms_pg_quat_cmp_atol);

  while(root) tdelete(root, &root, qpms_pg_irot3_approx_cmp_v); // I hope this is the correct way.

  return target;
}

qpms_irot3_t *qpms_pg_elems(qpms_irot3_t *target, qpms_pointgroup_t g) {
  QPMS_UNTESTED;
  target = qpms_pg_canonical_elems(target, g.c, g.n);
  qpms_gmi_t order = qpms_pg_order(g.c, g.n);
  const qpms_irot3_t o = g.orientation, o_inv = qpms_irot3_inv(o);
  for(qpms_gmi_t i = 0 ; i < order; ++i) 
    target[i] = qpms_irot3_mult(o_inv, qpms_irot3_mult(target[i], o));
  return target;
}

_Bool qpms_pg_is_subgroup(qpms_pointgroup_t small, qpms_pointgroup_t big) {
  QPMS_UNTESTED;
  qpms_gmi_t order_big = qpms_pg_order(big.c, big.n);
  qpms_gmi_t order_small = qpms_pg_order(small.c, small.n);
  if (!order_big || !order_small)
    QPMS_NOT_IMPLEMENTED("Subgroup testing for infinite groups not implemented");
  if (order_big < order_small) return false;
  // TODO maybe some other fast-track negative decisions

  qpms_irot3_t *elems_small = qpms_pg_elems(NULL, small);
  qpms_irot3_t *elems_big = qpms_pg_elems(NULL, small);
  qsort(elems_big, order_big, sizeof(qpms_irot3_t), qpms_pg_irot3_cmp_v);

  for(qpms_gmi_t smalli = 0; smalli < order_small; ++smalli) {
    qpms_irot3_t *match = bsearch(&elems_small[smalli], elems_big, order_big,
        sizeof(qpms_irot3_t), qpms_pg_irot3_approx_cmp_v);
    if (!match) return false;
  }

  return true;
}

/// Returns the order of a given 3D point group type.
/** For infinite groups returns 0. */
qpms_gmi_t qpms_pg_order(qpms_pointgroup_class cls, ///< Point group class.
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


/// Returns the number of canonical generators of a given 3D point group type.
/** TODO what does it do for infinite (Lie) groups? */
qpms_gmi_t qpms_pg_genset_size(qpms_pointgroup_class cls, ///< Point group class.
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
qpms_gmi_t qpms_pg_genset(qpms_pointgroup_class cls, ///< Point group class.
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
