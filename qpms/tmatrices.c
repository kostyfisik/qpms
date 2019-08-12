#define _POSIX_C_SOURCE 200809L // for getline()
#define lapack_int int
#define lapack_complex_double complex double
#define lapack_complex_double_real(z) (creal(z))
#define lapack_complex_double_imag(z) (cimag(z))
#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <unistd.h>
#include "scatsystem.h"
#include "indexing.h"
#include "vswf.h"
#include "groups.h"
#include "symmetries.h"
#include <assert.h>
#include <unistd.h>
#include "vectors.h"
#include "quaternions.h"
#include <string.h>
#include "qpms_error.h"
#include "tmatrices.h"
#include "qpms_specfunc.h"
#include "normalisation.h"
#include <errno.h>
#include <gsl/gsl_integration.h>

// These are quite arbitrarily chosen constants for the quadrature in qpms_tmatrix_axialsym_fill()
#define TMATRIX_AXIALSYM_INTEGRAL_EPSREL (1e-5)
#define TMATRIX_AXIALSYM_INTEGRAL_EPSABS (1e-8)

#define HBAR (1.05457162825e-34)
#define ELECTRONVOLT (1.602176487e-19)
#define SCUFF_OMEGAUNIT (3e14)

#define SQ(x) ((x)*(x))
#define MAX(x,y) ((x) < (y) ? (y) : (x))
#define MIN(x,y) ((x) > (y) ? (y) : (x))

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
  QPMS_ENSURE_SUCCESS(qpms_tmatrix_interpolator_eval_fill(t, ip, freq));
  return t;
}

qpms_errno_t qpms_tmatrix_interpolator_eval_fill(qpms_tmatrix_t *t,
    const qpms_tmatrix_interpolator_t *ip, double freq) {
  QPMS_ENSURE(qpms_vswf_set_spec_isidentical(t->spec, ip->bspec), "Tried to fill a T-matrix with an incompatible interpolator!");
  const size_t n = ip->bspec->n;
  for (size_t i = 0; i < n*n; ++i){
    if (ip->splines_real[i]) t->m[i] = gsl_spline_eval(ip->splines_real[i], freq, NULL /*does this work?*/);
    if (ip->splines_imag[i]) t->m[i] += I* gsl_spline_eval(ip->splines_imag[i], freq, NULL /*does this work?*/);
  }
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_tmatrix_generator_interpolator(qpms_tmatrix_t *t,
    complex double omega, const void *interp)
{
  return qpms_tmatrix_interpolator_eval_fill(t, 
      (const qpms_tmatrix_interpolator_t *) interp, omega);
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
  if(bs->norm & (QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE 
        | QPMS_NORMALISATION_SPHARM_REAL))
    QPMS_NOT_IMPLEMENTED("Sorry, only standard complex-spherical harmonic based waves are supported right now");
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
      (*tmdata)[(*n-1)*bs->n*bs->n + desti*bs->n + srci] = (tr + I*ti)
        * qpms_normalisation_factor_uvswfi(bs->norm, srcui)
        / qpms_normalisation_factor_uvswfi(bs->norm, destui)
        * qpms_normalisation_factor_uvswfi(QPMS_NORMALISATION_CONVENTION_SCUFF, destui)
        / qpms_normalisation_factor_uvswfi(QPMS_NORMALISATION_CONVENTION_SCUFF, srcui);
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

bool qpms_load_scuff_tmatrix_crash_on_failure = true;

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
    if (qpms_load_scuff_tmatrix_crash_on_failure)
      qpms_pr_error_at_line(__FILE__, __LINE__, __func__,
          "Could not open T-matrix file %s", path); 
    else return errno;
  qpms_errno_t retval = 
    qpms_read_scuff_tmatrix(f, bs, n, freqs, freqs_su, tmatrices_array, tmdata);

  for (size_t i = 0; i < *n * bs->n * bs->n; ++i) 
    if(isnan(creal((*tmdata)[i])) || isnan(cimag((*tmdata)[i]))) {
      QPMS_WARN("Encountered NAN in a loaded T-matrix");
      retval |= QPMS_NAN_ENCOUNTERED;
      break;
    }

  if(fclose(f)) qpms_pr_error_at_line(__FILE__, __LINE__, __func__,
        "Could not close the T-matrix file %s (well, that's weird, "
        "since it's read only).", path); 

  return retval;
}

complex double *qpms_mie_coefficients_reflection(
		complex double *target, ///< Target array of length bspec->n. If NULL, a new one will be allocated.
		const qpms_vswf_set_spec_t *bspec, ///< Defines which of the coefficients are calculated.
		double a, ///< Radius of the sphere.
		complex double k_i, ///< Wave number of the internal material of the sphere.
		complex double k_e, ///< Wave number of the surrounding medium.
		complex double mu_i, ///< Relative permeability of the sphere material.
		complex double mu_e, ///< Relative permeability of the surrounding medium.
    qpms_bessel_t J_ext,
    qpms_bessel_t J_scat
		// TODO J_ext, J_scat?
		) {
  /*
   *  This implementation pretty much copies mie_coefficients() from qpms_p.py, so 
   *  any bugs there should affect this function as well and perhaps vice versa.
   *
   *  Compared to mie_coefficients() from qpms_p.py, this has opposite sign.
   */
  QPMS_ENSURE(J_ext != J_scat, "J_ext and J_scat must not be equal. Perhaps you want J_ext = 1 and J_scat = 3 anyways.");
  if (!target) QPMS_CRASHING_MALLOC(target, bspec->n * sizeof(complex double));
  qpms_l_t lMax = bspec->lMax;
  complex double x_i = k_i * a;
  complex double x_e = k_e * a;
  complex double m = k_i / k_e;
  complex double eta_inv_i = k_i / mu_i;
  complex double eta_inv_e = k_e / mu_e;

  complex double zi[lMax + 2];
  complex double ze[lMax + 2];
  complex double zs[lMax + 2];
  complex double RH[lMax + 1] /* MAGNETIC */, RV[lMax+1] /* ELECTRIC */;
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(QPMS_BESSEL_REGULAR, lMax+1, x_i, zi));
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(J_ext, lMax+1, x_e, ze));
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(J_scat, lMax+1, x_e, zs));
  for (qpms_l_t l = 0; l <= lMax; ++l) {
    // Bessel function derivatives as in DLMF 10.51.2
    complex double dzi = -zi[l+1] + l/x_i*zi[l];
    complex double dze = -ze[l+1] + l/x_e*ze[l];
    complex double dzs = -zs[l+1] + l/x_e*zs[l];
    complex double fi = zi[l]/x_i+dzi;
    complex double fs = zs[l]/x_e+dzs;
    complex double fe = ze[l]/x_e+dze;
    RV[l] = (-eta_inv_i * fe * zi[l] + eta_inv_e * ze[l] * fi)/(-eta_inv_e * fi * zs[l] + eta_inv_i * zi[l] * fs);
    RH[l] = (-eta_inv_e * fe * zi[l] + eta_inv_i * ze[l] * fi)/(-eta_inv_i * fi * zs[l] + eta_inv_e * zi[l] * fs);
  }

  for (size_t i = 0; i < bspec->n; ++i) {
    qpms_l_t l; qpms_m_t m; qpms_vswf_type_t t;
    QPMS_ENSURE_SUCCESS(qpms_uvswfi2tmn(bspec->ilist[i], &t, &m, &l));
    assert(l <= lMax);
    switch(t) {
      case QPMS_VSWF_ELECTRIC:
        target[i] = RV[l];
        break;
      case QPMS_VSWF_MAGNETIC:
        target[i] = RH[l];
        break;
      default: QPMS_WTF;
    }
  }
  return target;
}

/// Replaces the contents of an existing T-matrix with that of a spherical nanoparticle calculated using the Lorentz-mie theory.
qpms_errno_t qpms_tmatrix_spherical_fill(qpms_tmatrix_t *t, ///< T-matrix whose contents are to be replaced. Not NULL.
		double a, ///< Radius of the sphere.
		complex double k_i, ///< Wave number of the internal material of the sphere.
		complex double k_e, ///< Wave number of the surrounding medium.
		complex double mu_i, ///< Relative permeability of the sphere material.
		complex double mu_e ///< Relative permeability of the surrounding medium.
		) {
  qpms_l_t lMax = t->spec->lMax;
  complex double *miecoeffs = qpms_mie_coefficients_reflection(NULL, t->spec, a, k_i, k_e, mu_i, mu_e,
    QPMS_BESSEL_REGULAR, QPMS_HANKEL_PLUS);
  memset(t->m, 0, SQ(t->spec->n));
  for(size_t i = 0; i < t->spec->n; ++i) 
    t->m[t->spec->n*i + i] = miecoeffs[i];
  free(miecoeffs);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_tmatrix_generator_sphere(qpms_tmatrix_t *t,
    complex double omega, const void *params) 
{
  const qpms_tmatrix_generator_sphere_param_t *p = params;
  qpms_epsmu_t out = p->outside.function(omega, p->outside.params);
  qpms_epsmu_t in = p->inside.function(omega, p->inside.params);
  return qpms_tmatrix_spherical_fill(t, 
      p->radius,
      qpms_wavenumber(omega, in),
      qpms_wavenumber(omega, out),
      in.mu, out.mu);
}

/// Convenience function to calculate T-matrix of a non-magnetic spherical \
particle using the permittivity values, replacing existing T-matrix data.
qpms_errno_t qpms_tmatrix_spherical_mu0_fill(
		qpms_tmatrix_t *t, ///< T-matrix whose contents are to be replaced. Not NULL.
		double a, ///< Radius of the sphere.
		double omega, ///< Angular frequency.
		complex double epsilon_fg, ///< Permittivity of the sphere.
		complex double epsilon_bg ///< Permittivity of the background medium.
		) 
{
  complex double k_i = csqrt(epsilon_fg) * omega / SPEED_OF_LIGHT;
  complex double k_e = csqrt(epsilon_bg) * omega / SPEED_OF_LIGHT;
  return qpms_tmatrix_spherical_fill(t, a, k_i, k_e, 1, 1);
}

complex double *qpms_apply_tmatrix(
    complex double *f, const complex double *a, const qpms_tmatrix_t *T) 
{
  const size_t n = T->spec->n;
  if(!f)
    QPMS_CRASHING_CALLOC(f, n, sizeof(complex double));
  const complex double one = 1;
  const complex double zero = 0;
  cblas_zgemv(CblasRowMajor, CblasNoTrans, n, n, &one, T->m, n, a, 1, &zero, NULL, 1);
  return f;
}

qpms_arc_function_retval_t qpms_arc_sphere(double theta, const void *R) {
  qpms_arc_function_retval_t retval = {*(const double*)R, 0};
  return retval;
}

qpms_arc_function_retval_t qpms_arc_cylinder(double theta, const void *param) {
  const qpms_arc_cylinder_params_t *p = param;
  double thresh = atan(2 * p->R / p->h);
  QPMS_ENSURE(theta >= 0 && theta <= M_PI, 
      "theta = %g, but it must lie in interval [0, M_PI]", theta);
  qpms_arc_function_retval_t res;
  if (theta < thresh) {
    // Upper base
    res.r = 0.5 * p->h / cos(theta);
    res.beta = -theta;
  } else if (theta <= M_PI - thresh) {
    // Side
    res.r = p->R / cos(theta - M_PI_2);
    res.beta = -theta + M_PI_2;
  } else {
    // Lower base
    res.r = 0.5 * p->h / cos(theta - M_PI);
    res.beta = -theta + M_PI;
  }
  return res;
}


struct tmatrix_axialsym_integral_param_t {
  const qpms_vswf_set_spec_t *bspec;
  qpms_l_t l, l_in; 
  qpms_m_t m; // m_in = -m
  qpms_vswf_type_t t, t_in;
  qpms_arc_function_t f;
  complex double k_in, k, z_in, z;
  bool realpart; // Otherwise imaginary part
  qpms_bessel_t btype; // For Q QPMS_HANKEL_PLUS, for R QPMS_BESSEL_REGULAR
};

static double tmatrix_axialsym_integrand(double theta, void *param) {
  // Pretty inefficient; either real or imaginary part is thrown away etc.
  struct tmatrix_axialsym_integral_param_t *p = param;
  qpms_arc_function_retval_t rb = p->f.function(theta, p->f.params);
  csph_t kr = {rb.r * p->k,    theta, 0},
         kr_in = {rb.r * p->k_in, theta, 0};

  csphvec_t y_el = qpms_vswf_single_el_csph(p->m, p->l, kr, p->btype, p->bspec->norm);
  csphvec_t y_mg = qpms_vswf_single_mg_csph(p->m, p->l, kr, p->btype, p->bspec->norm);
  csphvec_t v_in_el = qpms_vswf_single_el_csph(-p->m, p->l_in, kr_in, QPMS_BESSEL_REGULAR, p->bspec->norm);
  csphvec_t v_in_mg = qpms_vswf_single_mg_csph(-p->m, p->l_in, kr_in, QPMS_BESSEL_REGULAR, p->bspec->norm);
  csphvec_t y1, y2, v_in1, v_in2;
  switch(p->t) {
    case QPMS_VSWF_ELECTRIC:  y1 = y_el; y2 = y_mg;  break;
    case QPMS_VSWF_MAGNETIC:  y1 = y_mg; y2 = y_el;  break;
    default: QPMS_WTF;
  }
  switch(p->t_in) {
    case QPMS_VSWF_ELECTRIC:  v_in1 = v_in_mg; v_in2 = v_in_el; break;
    case QPMS_VSWF_MAGNETIC:  v_in1 = v_in_el; v_in2 = v_in_mg; break;
    default: QPMS_WTF;
  }
  // Normal vector components
  double nrc = cos(rb.beta), nthetac = sin(rb.beta);
  // First triple product
  complex double tp1 = nrc * (y1.thetac * v_in1.phic - y1.phic * v_in1.thetac)
                 + nthetac * (y1.phic   * v_in1.rc   - y1.rc   * v_in1.phic);
  // Second thiple product
  complex double tp2 = nrc * (y2.thetac * v_in2.phic - y2.phic * v_in2.thetac)
                 + nthetac * (y2.phic   * v_in2.rc   - y2.rc   * v_in2.phic);
  double jac = SQ(rb.r) * sin(theta) / nrc; // Jacobian
  complex double res = p->z/p->z_in * tp1 + tp2;
  return p->realpart ? creal(res) : cimag(res);
}

qpms_errno_t qpms_tmatrix_axialsym_fill(
		qpms_tmatrix_t *t, complex double omega, qpms_epsmu_t outside,
		qpms_epsmu_t inside,qpms_arc_function_t shape, qpms_l_t lMaxQR)
{
  QPMS_UNTESTED; // TODO also check whether this is valid for all norms.
  const qpms_vswf_set_spec_t *bspec = t->spec;
  struct tmatrix_axialsym_integral_param_t p;
  p.k = qpms_wavenumber(omega, outside);
  p.k_in = qpms_wavenumber(omega, inside);
  p.z = qpms_waveimpedance(outside);
  p.z_in = qpms_waveimpedance(inside);
  p.f = shape;
  p.bspec = bspec;
  const gsl_function f = {tmatrix_axialsym_integrand, (void *) &p};

  if (lMaxQR < bspec->lMax) lMaxQR = bspec->lMax;
  qpms_vswf_set_spec_t *bspecQR = qpms_vswf_set_spec_from_lMax(lMaxQR, bspec->norm);
  size_t *reindex = qpms_vswf_set_reindex(bspec, bspecQR);

  // Q' and R' matrices
  complex double *Q, *R;
  QPMS_CRASHING_MALLOC(Q, sizeof(complex double) * SQ(bspecQR->n));
  QPMS_CRASHING_MALLOC(R, sizeof(complex double) * SQ(bspecQR->n));

  // Not sure what the size should be, but this should be more than enough.
  const size_t intlimit = 1024;
  const double epsrel = TMATRIX_AXIALSYM_INTEGRAL_EPSREL;
  const double epsabs = TMATRIX_AXIALSYM_INTEGRAL_EPSABS;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(intlimit);
  for(size_t i1 = 0; i1 < bspecQR->n; ++i1)
    for(size_t i2 = 0; i2 < bspecQR->n; ++i2) {
      // NOTE that the Q', R' matrices are !TRANSPOSED! here because of zgetrs
      const size_t iQR = i1 + i2 * bspecQR->n;
      double abserr; // gsl_integration_qag() needs this; thrown away for now
      qpms_m_t m1, m2;
      qpms_l_t l1, l2;
      qpms_vswf_type_t t1, t2;
      const qpms_uvswfi_t u1 = bspecQR->ilist[i1], u2 = bspecQR->ilist[i2];
      QPMS_ENSURE_SUCCESS(qpms_uvswfi2tmn(u1, &t1, &m1, &l1));
      QPMS_ENSURE_SUCCESS(qpms_uvswfi2tmn(u2, &t2, &m2, &l2));
      if (m1 + m2) {
        Q[iQR] = 0;
        R[iQR] = 0;
      } else {
        p.m = m1;
        p.l = l1; p.l_in = l2;
        p.t = t1; p.t_in = t2;

        double result; // We throw the quadrature error estimate away.
        // Re(R')
        p.btype = QPMS_BESSEL_REGULAR;
        p.realpart = true;
        QPMS_ENSURE_SUCCESS(gsl_integration_qag(&f, 0, M_PI, epsabs, epsrel,
              intlimit, 2, w, &result, &abserr));
        R[iQR] = result;
        
        // Im(R')
        p.realpart = false;
        QPMS_ENSURE_SUCCESS(gsl_integration_qag(&f, 0, M_PI, epsabs, epsrel,
              intlimit, 2, w, &result, &abserr));
        R[iQR] += I*result;

        // Re(Q')
        p.btype = QPMS_HANKEL_PLUS;
        p.realpart = true;
        QPMS_ENSURE_SUCCESS(gsl_integration_qag(&f, 0, M_PI, epsabs, epsrel,
              intlimit, 2, w, &result, &abserr));
        Q[iQR] = result;

        // Im(Q')
        p.realpart = false;
        QPMS_ENSURE_SUCCESS(gsl_integration_qag(&f, 0, M_PI, epsabs, epsrel,
              intlimit, 2, w, &result, &abserr));
        Q[iQR] += I*result;
      }
    }
  gsl_integration_workspace_free(w);

  // Compute T-matrix; maybe it would be better solved with some sparse matrix mechanism,
  // but fukkit.
  const size_t n = bspecQR->n;
  lapack_int *pivot;
  QPMS_CRASHING_MALLOC(pivot, sizeof(lapack_int) * n);
  QPMS_ENSURE_SUCCESS(LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, Q, n, pivot));
  // Solve Q'^T X = R'^T where X will be -T^T
  // Note that Q'^T, R'^T are already (transposed) in Q, R.
  const complex double minus1 = -1;
  QPMS_ENSURE_SUCCESS(LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', n, n /*nrhs*/,
        Q, n, pivot, R, n));
  // R now contains -T^T.
  for(size_t i1 = 0; i1 < bspec->n; ++i1) 
    for(size_t i2 = 0; i2 < bspec->n; ++i2) {
      if (reindex[i1] == ~(size_t) 0 && reindex[i2] == ~(size_t) 0) QPMS_WTF;
      const size_t it = i1 * bspec->n + i2;
      const size_t iQR = reindex[i1] + reindex[i2] * bspecQR->n;
      t->m[it] = -R[iQR];
    }

  free(pivot); free(R); free(Q); free(reindex);
  qpms_vswf_set_spec_free(bspecQR);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_tmatrix_generator_axialsym(qpms_tmatrix_t *t, complex double omega, const void *param) {
  const qpms_tmatrix_generator_axialsym_param_t *p = param;
  return qpms_tmatrix_axialsym_fill(t, omega, 
      p->outside.function(omega, p->outside.params),
      p->inside.function(omega, p->inside.params),
      p->shape,
      p->lMax_extend);
}

qpms_errno_t qpms_tmatrix_generator_constant(qpms_tmatrix_t *t,
    complex double omega, const void *tmorig) {
  QPMS_UNTESTED;
  // TODO more thorough check on bspec!
  const qpms_tmatrix_t *orig = tmorig;
  const size_t tocopy = SQ(MIN(t->spec->n, orig->spec->n));
  memcpy(t->m, orig->m, tocopy*sizeof(complex double));
  memset(t->m+tocopy, 0, (SQ(t->spec->n)-tocopy)*sizeof(complex double));
  return QPMS_SUCCESS;
}


