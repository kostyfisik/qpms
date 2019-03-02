#include "groups.h"

/// Trivial group, with one (reduntant) generator.
/**
 * For the trivial group, zero generators are enough.
 * However, some functions might be not robust enough and require
 * a first generator to work properly.
 */
const qpms_finite_group_t QPMS_FINITE_GROUP_TRIVIAL_G = {
  "trivial_g", // name
  1, // order
  0, // idi
  (qpms_gmi_t[]) { // mt
    0,
  },
  (qpms_gmi_t[]) { // invi
    0
  },
  (qpms_gmi_t[]) {0}, // gens
  1, // ngens
  (qpms_permutation_t[]){ // permrep
    "()",
  },
  NULL, // elemlabels
  0, // permrep_nelem
  (qpms_irot3_t[]) { // rep3d
   {{1.0+0.0*I, 0.0+0.0*I}, 1},
  },
  1, // nirreps
  (struct qpms_finite_group_irrep_t[]) { // irreps
    {
      1, // dim
      "A", //name
      (complex double []) {1} // m
    },
  } // end of irreps
};

/// Trivial group.
const qpms_finite_group_t QPMS_FINITE_GROUP_TRIVIAL = {
  "trivial", // name
  1, // order
  0, // idi
  (qpms_gmi_t[]) { // mt
    0,
  },
  (qpms_gmi_t[]) { // invi
    0
  },
  NULL, // gens
  0, // ngens
  (qpms_permutation_t[]){ // permrep
    "()",
  },
  NULL, // elemlabels
  0, // permrep_nelem
  (qpms_irot3_t[]) { // rep3d
   {{1.0+0.0*I, 0.0+0.0*I}, 1},
  },
  1, // nirreps
  (struct qpms_finite_group_irrep_t[]) { // irreps
    {
      1, // dim
      "A", //name
      (complex double []) {1} // m
    },
  } // end of irreps
};


