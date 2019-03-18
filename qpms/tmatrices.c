#define _POSIX_C_SOURCE 200809L // for getline()
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <unistd.h>
#include "scatsystem.h"
#include "indexing.h"
#include "vswf.h"
#include "groups.h"
#include "symmetries.h"
#include <gsl/gsl_spline.h>
#include <assert.h>
#include <unistd.h>
#include "vectors.h"
#include "wigner.h"
#include <string.h>
#include "qpms_error.h"
#include "tmatrices.h"

#define HBAR (1.05457162825e-34)
#define ELECTRONVOLT (1.602176487e-19)
#define SCUFF_OMEGAUNIT (3e14)

#define SQ(x) ((x)*(x))
qpms_tmatrix_t *qpms_tmatrix_init(const qpms_vswf_set_spec_t *bspec) {
  qpms_tmatrix_t *t = malloc(sizeof(qpms_tmatrix_t));
  if (!t) abort();
  else {
    t->spec = bspec;
    size_t n = bspec->n;
    t->m = calloc(n*n, sizeof(complex double));
    if (!t->m) abort();
    t->owns_m = true;
  }
  return t;
}

qpms_tmatrix_t *qpms_tmatrix_copy(const qpms_tmatrix_t *T) {
  qpms_tmatrix_t *t = qpms_tmatrix_init(T->spec);
  size_t n = T->spec->n;
  for(size_t i = 0; i < n*n; ++i)
    t->m = T->m;
  return t;
}

void qpms_tmatrix_free(qpms_tmatrix_t *t){
  if(t && t->owns_m) free(t->m);
  free(t);
}

qpms_tmatrix_t *qpms_tmatrix_apply_symop_inplace(
                qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  //qpms_tmatrix_t *t = qpms_tmatrix_init(T->spec);
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
      n, n, n, &one, tmp, n, M, n, &zero, T->m, n);
  return T;
}

qpms_tmatrix_t *qpms_tmatrix_apply_symop(
                const qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  qpms_tmatrix_t *t = qpms_tmatrix_init(T->spec);
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

qpms_errno_t qpms_symmetrise_tmdata_irot3arr(
    complex double *tmdata, const size_t tmcount,
    const qpms_vswf_set_spec_t *bspec,
    const size_t n_symops, const qpms_irot3_t *symops) {
  const size_t n = bspec->n;
  qpms_tmatrix_t *tmcopy = qpms_tmatrix_init(bspec);
  complex double *symop_matrices = malloc(n*n*sizeof(complex double) * n_symops);
  if(!symop_matrices) qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
      "malloc() failed.");
  for (size_t i = 0; i < n_symops; ++i) 
    qpms_irot3_uvswfi_dense(symop_matrices + i*n*n, bspec, symops[i]);
  complex double tmp[n][n];
  const complex double one = 1, zero = 0;
  for (size_t tmi = 0; tmi < tmcount; ++tmi) {
    // Move the data in tmcopy; we will then write the sum directly into tmdata.
    memcpy(tmcopy->m, tmdata+n*n*tmi, n*n*sizeof(complex double));
    memset(tmdata+n*n*tmi, 0, n*n*sizeof(complex double));
    for (size_t i = 0; i < n_symops; ++i) {
      const complex double *const M = symop_matrices + i*n*n;
      // tmp = M T
      cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
          n, n, n, &one, M, n, tmcopy->m, n, &zero, tmp, n);
      // tmdata[...] += tmp M* = M T M*
      cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
          n, n, n, &one, tmp, n, M, n, &one, tmdata + tmi*n*n, n);
    }
    for (size_t ii = 0; ii < n*n; ++ii)
      tmdata[n*n*tmi+ii] /= n_symops;
  }
  free(symop_matrices);
  qpms_tmatrix_free(tmcopy);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_symmetrise_tmdata_finite_group(
    complex double *tmdata, const size_t tmcount,
    const qpms_vswf_set_spec_t *bspec,
    const qpms_finite_group_t *pointgroup) {
  if (!(pointgroup->rep3d)) qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
      "This function requires pointgroup->rep3d to be set correctly!");
  return qpms_symmetrise_tmdata_irot3arr(tmdata, tmcount, bspec,
      pointgroup->order, pointgroup->rep3d);
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_irot3arr_inplace(
    qpms_tmatrix_t *T,
    size_t n_symops,
    const qpms_irot3_t *symops
    ) {
  if(qpms_symmetrise_tmdata_irot3arr(T->m, 1,
        T->spec, n_symops, symops) != QPMS_SUCCESS)
    return NULL;
  else return T;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_finite_group_inplace(
    qpms_tmatrix_t *T,
    const qpms_finite_group_t *pointgroup
    ) {
  if(qpms_symmetrise_tmdata_finite_group(T->m, 1,
        T->spec, pointgroup) != QPMS_SUCCESS)
    return NULL;
  else return T;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_involution_inplace(
                qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  qpms_tmatrix_t *t = qpms_tmatrix_apply_symop(T, M);
  const size_t n = T->spec->n;
  for(size_t i = 0; i < n*n; ++i)
    T->m[i] = 0.5 * (t->m[i] + T->m[i]);
  qpms_tmatrix_free(t);
  return T;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_involution(
                const qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  qpms_tmatrix_t *t = qpms_tmatrix_init(T->spec);
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
  qpms_tmatrix_t *t = qpms_tmatrix_copy(T);
  return qpms_tmatrix_symmetrise_C_inf_inplace(t);
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_inf_inplace(qpms_tmatrix_t *T) {
  const size_t n = T->spec->n;
  for (size_t row = 0; row < n; row++) {
    qpms_m_t rm = qpms_uvswfi2m(T->spec->ilist[row]);
    for (size_t col = 0; col < n; col++) {
      qpms_m_t cm = qpms_uvswfi2m(T->spec->ilist[col]);
      if (rm == cm)
        ;// No-op // t->m[n*row + col] = T->m[n*row + col];
      else
        T->m[n*row + col] = 0;
    }
  }
  return T;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_N(const qpms_tmatrix_t *T, int N) {
  qpms_tmatrix_t *t = qpms_tmatrix_copy(T);
  return qpms_tmatrix_symmetrise_C_N_inplace(t, N);
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_N_inplace(qpms_tmatrix_t *T, int N) {
  const size_t n = T->spec->n;
  for (size_t row = 0; row < n; row++) {
    qpms_m_t rm = qpms_uvswfi2m(T->spec->ilist[row]);
    for (size_t col = 0; col < n; col++) {
      qpms_m_t cm = qpms_uvswfi2m(T->spec->ilist[col]);
      if (((rm - cm) % N) == 0)
        ; // T->m[n*row + col] = T->m[n*row + col];
      else
        T->m[n*row + col] = 0;
    }
  }
  return T;
}

bool qpms_tmatrix_isclose(const qpms_tmatrix_t *A, const qpms_tmatrix_t *B,
    const double rtol, const double atol)
{
  if (!qpms_vswf_set_spec_isidentical(A->spec, B->spec))
    return false;
  if (A->m == B->m)
    return true;
  const size_t n = A->spec->n;
  for (size_t i = 0; i < n*n; ++i) {
    const double tol = atol + rtol * (cabs(B->m[i]));
    if ( cabs(B->m[i] - A->m[i]) > tol )
      return false;
  }
  return true;
}
    
qpms_tmatrix_interpolator_t *qpms_tmatrix_interpolator_create(const size_t incount,
    const double *freqs, const qpms_tmatrix_t *ta, const gsl_interp_type *iptype//, const bool copy_bspec
    ) {
  if (incount <= 0) return NULL;
  qpms_tmatrix_interpolator_t *ip = malloc(sizeof(qpms_tmatrix_interpolator_t));
  /*
  if (copy_bspec) {
    ip->bspec = qpms_vswf_set_spec_copy(ta[0].spec);
    ip->owns_bspec = true;
  }
  else {
  */
    ip->bspec = ta[0].spec;
  //  ip->owns_bspec = false;
  //}
  const size_t n = ip->bspec->n;

  // check if all matrices have the same bspec
  for (size_t i = 0; i < incount; ++i)
    if (!qpms_vswf_set_spec_isidentical(ip->bspec, ta[i].spec))
      abort();

  if (!(ip->splines_real = calloc(n*n,sizeof(gsl_spline *)))) abort();
  if (!(ip->splines_imag = calloc(n*n,sizeof(gsl_spline *)))) abort();
  for (size_t row = 0; row < n; ++row)
    for (size_t col = 0; col < n; ++col) {
      double y_real[incount], y_imag[incount];
      bool n0_real = false, n0_imag = false;
      for (size_t i = 0; i < incount; ++i) {
        complex double telem = ta[i].m[n * row + col];
        if ((y_real[i] = creal(telem))) n0_real = true;
        if ((y_imag[i] = cimag(telem))) n0_imag = true;
      }
      if (n0_real) {
        gsl_spline *s =
        ip->splines_real[n * row + col] = gsl_spline_alloc(iptype, incount);
        if (gsl_spline_init(s, freqs, y_real, incount) != 0 /*GSL_SUCCESS*/) abort();
      }
      else ip->splines_real[n * row + col] = NULL;
     if (n0_imag) {
        gsl_spline *s =
        ip->splines_imag[n * row + col] = gsl_spline_alloc(iptype, incount);
        if (gsl_spline_init(s, freqs, y_imag, incount) != 0 /*GSL_SUCCESS*/) abort();
      }
      else ip->splines_imag[n * row + col] = NULL;
    }
  return ip;
}

void qpms_tmatrix_interpolator_free(qpms_tmatrix_interpolator_t *ip) {
  if (ip) {
    const size_t n = ip->bspec->n;
    for (size_t i = 0; i < n*n; ++i) {
      if (ip->splines_real[i]) gsl_spline_free(ip->splines_real[i]);
      if (ip->splines_imag[i]) gsl_spline_free(ip->splines_imag[i]);
    }
    //if (ip->owns_bspec)
    //  qpms_vswf_set_spec_free(ip->bspec);
    free(ip);
  }
}

qpms_tmatrix_t *qpms_tmatrix_interpolator_eval(const qpms_tmatrix_interpolator_t *ip, double freq) {
  qpms_tmatrix_t *t = qpms_tmatrix_init(ip->bspec);
  const size_t n = ip->bspec->n;
  for (size_t i = 0; i < n*n; ++i){
    if (ip->splines_real[i]) t->m[i] = gsl_spline_eval(ip->splines_real[i], freq, NULL /*does this work?*/);
    if (ip->splines_imag[i]) t->m[i] += I* gsl_spline_eval(ip->splines_imag[i], freq, NULL /*does this work?*/);
  }
  return t;
}



double qpms_SU2eV(double e_SU) {
  return e_SU * SCUFF_OMEGAUNIT / (ELECTRONVOLT / HBAR);
}

double qpms_SU2SI(double e_SU) {
  return e_SU * SCUFF_OMEGAUNIT;
}

/// TODO doc and more checks
qpms_errno_t qpms_read_scuff_tmatrix(
    FILE *f, ///< file handle
    const qpms_vswf_set_spec_t * bs, ///< VSWF set spec
    size_t *const n, ///< Number of successfully loaded t-matrices
    double* *const freqs, ///< Frequencies in SI units
    double* *const freqs_su, ///< Frequencies in SCUFF units (optional)
    qpms_tmatrix_t* *const tmatrices_array, ///< The resulting T-matrices (optional).
    complex double* *const tmdata
    ) {
  if (!(freqs && n && tmdata)) 
    qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
      "freqs, n, and tmdata are mandatory arguments and must not be NULL.");
  if (bs->norm != QPMS_NORMALISATION_POWER_CS) // CHECKME CORRECT?
    qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
        "Not implemented; only QPMS_NORMALISATION_POWER_CS (CHECKME)"
        " norm supported right now.");
  int n_alloc = 128; // First chunk to allocate
  *n = 0;
  *freqs = malloc(n_alloc * sizeof(double));
  if (freqs_su) *freqs_su = malloc(n_alloc * sizeof(double));
  *tmdata = malloc(sizeof(complex double) * bs->n * bs->n * n_alloc);
  if (!*freqs || (!freqs_su != !*freqs_su) || !*tmdata)
      qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
          "malloc() failed.");
  size_t linebufsz = 256;
  char *linebuf = malloc(linebufsz);
  ssize_t readchars;
  double lastfreq_su = NAN;
  while((readchars = getline(&linebuf, &linebufsz, f)) != -1) {
    if (linebuf[0] == '#') continue;
    int Alpha, LAlpha, MAlpha, PAlpha, Beta, LBeta, MBeta, PBeta;
    double currentfreq_su, tr, ti;
    if (11 != sscanf(linebuf, "%lf %d %d %d %d %d %d %d %d %lf %lf",
          &currentfreq_su, &Alpha, &LAlpha, &MAlpha, &PAlpha,
          &Beta, &LBeta, &MBeta, &PBeta, &tr, &ti))
      abort(); // Malformed T-matrix file
    if (currentfreq_su != lastfreq_su) { // New frequency -> new T-matrix
      ++*n;
      lastfreq_su = currentfreq_su;
      if(*n > n_alloc) {
        n_alloc *= 2;
        *freqs = realloc(*freqs, n_alloc * sizeof(double));
        if (freqs_su) *freqs_su = realloc(*freqs_su, n_alloc * sizeof(double));
        *tmdata = realloc(*tmdata, sizeof(complex double) * bs->n * bs->n * n_alloc);
        if (!*freqs || (!freqs_su != !*freqs_su) || !*tmdata)
          qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
              "realloc() failed.");
      }
      if (freqs_su) (*freqs_su)[*n-1] = currentfreq_su;
      (*freqs)[*n-1] = qpms_SU2SI(currentfreq_su);

      for(size_t i = 0; i < bs->n * bs->n; ++i)
        (*tmdata)[(*n-1)*bs->n*bs->n + i] = NAN + I*NAN;
    }
    qpms_vswf_type_t TAlpha, TBeta;
    switch(PAlpha) {
      case 0: TAlpha = QPMS_VSWF_MAGNETIC; break;
      case 1: TAlpha = QPMS_VSWF_ELECTRIC; break;
      default: assert(0);
    }
    switch(PBeta) {
      case 0: TBeta = QPMS_VSWF_MAGNETIC; break;
      case 1: TBeta = QPMS_VSWF_ELECTRIC; break;
      default: assert(0);
    }
    qpms_uvswfi_t srcui = qpms_tmn2uvswfi(TAlpha, MAlpha, LAlpha),
                  destui = qpms_tmn2uvswfi(TBeta, MBeta, LBeta);
    ssize_t srci = qpms_vswf_set_spec_find_uvswfi(bs, srcui),
            desti = qpms_vswf_set_spec_find_uvswfi(bs, destui);
    if (srci == -1 || desti == -1)
      /* This element has not been requested in bs->ilist. */
      continue;
    else
      (*tmdata)[(*n-1)*bs->n*bs->n + desti*bs->n + srci] = tr + I*ti;
  }
  free(linebuf);
  // free some more memory
  n_alloc = *n;
  *freqs = realloc(*freqs, n_alloc * sizeof(double));
  if (freqs_su) *freqs_su = realloc(*freqs_su, n_alloc * sizeof(double));
  if (tmatrices_array) *tmatrices_array = realloc(*tmatrices_array, n_alloc * sizeof(qpms_tmatrix_t));
  *tmdata = realloc(*tmdata, sizeof(complex double) * bs->n * bs->n * n_alloc);
  if (!*freqs || (!freqs_su != !*freqs_su) || !*tmdata)
    qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
          "realloc() failed.");
  if (tmatrices_array) {
    *tmatrices_array = malloc(n_alloc * sizeof(qpms_tmatrix_t));
    if (!*tmatrices_array) 
      qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
          "malloc() failed.");
    for (size_t i = 0; i < *n; ++i) {
      (*tmatrices_array)[i].spec = bs;
      (*tmatrices_array)[i].m = (*tmdata) + i * bs->n * bs->n;
      (*tmatrices_array)[i].owns_m = false;
    }
  }
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_load_scuff_tmatrix(
    const char *path, ///< file path
    const qpms_vswf_set_spec_t * bs, ///< VSWF set spec
    size_t *const n, ///< Number of successfully loaded t-matrices
    double **const freqs, ///< Frequencies in SI units
    double ** const freqs_su, ///< Frequencies in SCUFF units (optional)
    qpms_tmatrix_t ** const tmatrices_array, ///< The resulting T-matrices (optional).
    complex double ** const tmdata
    ) {
  FILE *f = fopen(path, "r");
  if (!f)
    qpms_pr_error_at_line(__FILE__, __LINE__, __func__,
        "Could not open T-matrix file %s", path); 
  qpms_errno_t retval = 
    qpms_read_scuff_tmatrix(f, bs, n, freqs, freqs_su, tmatrices_array, tmdata);
  if(fclose(f)) qpms_pr_error_at_line(__FILE__, __LINE__, __func__,
        "Could not close the T-matrix file %s (well, that's weird, "
        "since it's read only).", path); 

  return retval;
}

