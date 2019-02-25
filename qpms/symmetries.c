#include "symmetries.h"
#include "tiny_inlines.h"
#include "indexing.h"
#include "wigner.h"

// TODO at some point, maybe support also other norms.
// (perhaps just use qpms_normalisation_t_factor() at the right places)
static inline void check_norm_compat(const qpms_vswf_set_spec_t *s)
{
  switch (qpms_normalisation_t_normonly(s->norm)) {
    case QPMS_NORMALISATION_POWER:
      break;
    case QPMS_NORMALISATION_SPHARM:
      break;
    default:
      abort(); // Only SPHARM and POWER norms are supported right now.
  }
}

// Used in the functions below to ensure memory allocation and checks for bspec validity
static inline complex double *ensure_alloc(complex double *target,
    const qpms_vswf_set_spec_t *bspec) {
  check_norm_compat(bspec);
  const size_t n = bspec->n;
  if (target == NULL)
    target = malloc(n * n * sizeof(complex double));
  if (target == NULL) abort();
  return target;
}



complex double *qpms_zflip_uvswi_dense(
    complex double *target,
    const qpms_vswf_set_spec_t *bspec) 
{
  target = ensure_alloc(target, bspec);
  const size_t n = bspec->n;

  for (size_t row = 0; row < n; row++) {
    qpms_vswf_type_t rt;
    qpms_l_t rl;
    qpms_m_t rm;
    qpms_uvswfi2tmn(bspec->ilist[row], &rt, &rm, &rl);
    for (size_t col = 0; col < n; col++) {
      qpms_vswf_type_t ct;
      qpms_l_t cl;
      qpms_m_t cm;
      if(qpms_uvswfi2tmn(bspec->ilist[col], &ct, &cm, &cl)) abort();
      if (rl == cl && rm == cm && rt == ct)
        switch(rt) {
          case QPMS_VSWF_ELECTRIC:
          case QPMS_VSWF_LONGITUDINAL:
            target[n*row + col] = min1pow(cm + cl);
            break;
          case QPMS_VSWF_MAGNETIC:
            target[n*row + col] = -min1pow(cm + cl);
            break;
          default:
            abort();
        }
      else target[n*row + col] = 0;
    }
  }
  return target;
}

complex double *qpms_yflip_uvswi_dense(
    complex double *target,
    const qpms_vswf_set_spec_t *bspec) 
{
  target = ensure_alloc(target, bspec);
  const size_t n = bspec->n;

  for (size_t row = 0; row < n; row++) {
    qpms_vswf_type_t rt;
    qpms_l_t rl;
    qpms_m_t rm;
    qpms_uvswfi2tmn(bspec->ilist[row], &rt, &rm, &rl);
    for (size_t col = 0; col < n; col++) {
      qpms_vswf_type_t ct;
      qpms_l_t cl;
      qpms_m_t cm;
      if(qpms_uvswfi2tmn(bspec->ilist[col], &ct, &cm, &cl)) abort();
      if (rl == cl && rm == -cm && rt == ct)
        switch(rt) {
          case QPMS_VSWF_ELECTRIC:
          case QPMS_VSWF_LONGITUDINAL:
            target[n*row + col] = min1pow(rm);
            break;
          case QPMS_VSWF_MAGNETIC:
            target[n*row + col] = -min1pow(rm);
            break;
          default:
            abort();
        }
      else target[n*row + col] = 0;
    }
  }
  return target;
}

complex double *qpms_xflip_uvswi_dense(
    complex double *target,
    const qpms_vswf_set_spec_t *bspec) 
{
  target = ensure_alloc(target, bspec);
  const size_t n = bspec->n;

  for (size_t row = 0; row < n; row++) {
    qpms_vswf_type_t rt;
    qpms_l_t rl;
    qpms_m_t rm;
    qpms_uvswfi2tmn(bspec->ilist[row], &rt, &rm, &rl);
    for (size_t col = 0; col < n; col++) {
      qpms_vswf_type_t ct;
      qpms_l_t cl;
      qpms_m_t cm;
      if(qpms_uvswfi2tmn(bspec->ilist[col], &ct, &cm, &cl)) abort();
      if (rl == cl && rm == -cm && rt == ct)
        switch(rt) {
          case QPMS_VSWF_ELECTRIC:
          case QPMS_VSWF_LONGITUDINAL:
            target[n*row + col] = 1;
            break;
          case QPMS_VSWF_MAGNETIC:
            target[n*row + col] = -1;
            break;
          default:
            abort();
        }
      else target[n*row + col] = 0;
    }
  }
  return target;
}

// Dense matrix representation of a rotation around the z-axis
complex double *qpms_zrot_uvswi_dense(
                complex double *target, ///< If NULL, a new array is allocated.
                const qpms_vswf_set_spec_t *bspec,
                double phi ///< Rotation angle
                )
{
  target = ensure_alloc(target, bspec);
  const size_t n = bspec->n;

  for (size_t row = 0; row < n; row++) {
    qpms_vswf_type_t rt;
    qpms_l_t rl;
    qpms_m_t rm;
    qpms_uvswfi2tmn(bspec->ilist[row], &rt, &rm, &rl);
    for (size_t col = 0; col < n; col++) {
      qpms_vswf_type_t ct;
      qpms_l_t cl;
      qpms_m_t cm;
      if(qpms_uvswfi2tmn(bspec->ilist[col], &ct, &cm, &cl)) abort();
      if (rl == cl && rm == cm && rt == ct) // TODO COMPARE WITH PYTHON
        target[n*row + col] = cexp(/* - ?*/I * rm * phi);
      else target[n*row + col] = 0;
    }
  }
  return target;
}

// Dense matrix representation of a "rational" rotation around the z-axis
/* Just for convenience. Corresponds to the angle \f$ \phi = 2\piw/N \f$.
 */
complex double *qpms_zrot_rational_uvswi_dense(
                complex double *target, ///< If NULL, a new array is allocated.
                const qpms_vswf_set_spec_t *bspec,
                int N,
                int w
) 
{
  double phi = 2 * M_PI * w / N;
  return qpms_zrot_uvswi_dense(target, bspec, phi);
}

complex double *qpms_irot3_uvswfi_dense(
    complex double *target,
    const qpms_vswf_set_spec_t *bspec,
    const qpms_irot3_t t)
{
  target = ensure_alloc(target, bspec);
  const size_t n = bspec->n;

  for (size_t row = 0; row < n; row++) {
    qpms_vswf_type_t rt;
    qpms_l_t rl;
    qpms_m_t rm;
    qpms_uvswfi2tmn(bspec->ilist[row], &rt, &rm, &rl);
    for (size_t col = 0; col < n; col++) {
      qpms_vswf_type_t ct;
      qpms_l_t cl;
      qpms_m_t cm;
      if(qpms_uvswfi2tmn(bspec->ilist[col], &ct, &cm, &cl)) abort();
      if (rl == cl && rm == cm && rt == ct)
        // TODO qpms_vswf_irot_elem_from_irot3 might be slow and not too accurate for large l
        target[n*row + col] = // Checkme rm and cm order
          qpms_vswf_irot_elem_from_irot3(t,
              rl, rm /* CHECKME here */, cm /* and here */,
              rt == QPMS_VSWF_MAGNETIC);
      else target[n*row + col] = 0;
    }
  }
  return target;
}

 
