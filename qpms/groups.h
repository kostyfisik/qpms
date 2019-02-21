/*! \file groups.h
 * \brief Point groups.
 *
 * Right now, the instances of qpms_finite_group_t are created at compilation time
 * from source code generated by Python script TODO (output groups.c)
 * and they are not to be constructed dynamically.
 */
#ifndef QPMS_GROUPS_H
#define QPMS_GROUPS_H

/// Group member index.
typedef int qpms_gmi_t;

/// Permutation representation of a group element.
/** For now, it's just a string of the form "(0,1)(3,4,5)"
 */
typedef const char * qpms_permutation_t;

/// To be used only in qpms_finite_group_t
struct qpms_finite_group_irrep_t {
	int dim; ///< Irrep dimension.
	char name[]; ///< Irrep label.
	/// Irrep matrix data.
	/** The r-th row, c-th column of the representation of the i'th element is retrieved as
	 *  m[i * dim * dim + r * dim + c]
	 */
	complex double *m; 
};

/// A point group with its irreducible representations and some metadata.
/** 
 *  The structure of the group is given by the multiplication table \a mt.
 *
 *  Each element of the group has its index from 0 to order.
 *  The metadata about some element are then accessed using that index.
 *
 *  All members are in principle optional except \a order and \a mt.
 *
 *  Note: after changing this struct, don't forget to update the Python method
 *  SVWFPointGroupInfo.generate_c_source().
 */
typedef struct qpms_finite_group_t {
	char *name;
	size_t order; ///< Group order (number of elements)
	qpms_gmi_t idi; ///< Identity element index
	qpms_gmi_t *mt; ///< Group multiplication table. If c = a*b, then ic = mt[order * ia + ib].
	qpms_gmi_t *gens; ///< A canonical set of group generators.
	int ngens; ///< Number of the generators in gens;
	qpms_permutation_t permrep[]; ///< Permutation representations of the elements.
	char **elemlabels; ///< Optional human readable labels for the group elements.
	int permrep_nelem; ///< Number of the elements over which the permutation representation acts.
	cmatrix3d rep3d[]; ///< The 'natural' 3D matrix representation of a 3D point group.
	int nirreps; ///< How many irreps does the group have
	struct qpms_finite_group_irrep_t irreps[]; ///< Irreducible representations of the group.
} qpms_finite_group_t;


#endif // QPMS_GROUPS_H
