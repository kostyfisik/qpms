#include "pointgroups.h"
#include <search.h>

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





