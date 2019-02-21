#include "scatsystem.h"

#include <cblas.h>

qpms_tmatrix_t *qpms_tmatrix_init(const qpms_vswf_set_spec_t *bspec) {
  qpms_tmatrix_t *t = malloc(sizeof(qpms_tmatrix_t));
  if (!t) abort();
  else {
    t->spec = bspec;
    size_t n = bspec->n;
    t->m = calloc(n*n, sizeof(complex double));
    if (!t->m) abort();
  }
  return t;
}

void qpms_tmatrix_free(qpms_tmatrix_t *t){
  if(t) free(t->m);
  free(t);
}

qpms_tmatrix_t *qpms_tmatrix_apply_symop(
                const qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  qpms_tmatrix_t t = qpms_tmatrix_init(T->spec);
  const size_t n = T->spec->n;
  complex double tmp[n][n];
  // tmp = M T
  const complex double one = 1, zero = 0;
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasNoTrans,
      n, n, n, &one, M, n, T->m, n, &zero, tmp, n);
  // t->m = tmp M* = M T M*
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasConjTrans,
      n, n, n, &one, tmp, n, M, n, &zero, t->m, n);
  return t;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_involution(
                const qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  qpms_tmatrix_t t = qpms_tmatrix_init(T->spec);
  const size_t n = T->spec->n;
  complex double tmp[n][n];
  // tmp = M T
  const complex double one = 1, zero = 0;
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasNoTrans,
      n, n, n, &one, M, n, T->m, n, &zero, tmp, n);
  // t->m = tmp M* = M T M*
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasConjTrans,
      n, n, n, &one, tmp, n, M, n, &zero, t->m, n);
  for(size_t i = 0; i < n*n; ++i)
    t->m[i] = 0.5 * (t->m[i] + T->m[i]);
  return t;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_inf(const qpms_tmatrix_t *T) {
  qpms_tmatrix_t t = qpms_tmatrix_init(T->spec);
  const size_t n = T->spec->n;
  for (size_t row = 0; row < n; row++) {
    qpms_m_t rm = qpms_uvswfi2m(T->spec->ilist[row]);
    for (size_t col = 0; col < n; col++) {
      qpms_m_t cm = qpms_uvswfi2m(T->spec->ilist[col]);
      if (rm == cm)
        t->m[n*row + col] = T->m[n*row + col];
      else
        t->m[n*row + col] = 0;
    }
  }
  return t;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_N(const qpms_tmatrix_t *T, int N) {
  qpms_tmatrix_t t = qpms_tmatrix_init(T->spec);
  const size_t n = T->spec->n;
  for (size_t row = 0; row < n; row++) {
    qpms_m_t rm = qpms_uvswfi2m(T->spec->ilist[row]);
    for (size_t col = 0; col < n; col++) {
      qpms_m_t cm = qpms_uvswfi2m(T->spec->ilist[col]);
      if (((rm - cm) % N) == 0)
        t->m[n*row + col] = T->m[n*row + col];
      else
        t->m[n*row + col] = 0;
    }
  }
  return t;
}


