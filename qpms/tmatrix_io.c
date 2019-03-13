#define _POSIX_C_SOURCE 200809L // for getline()
#include "scatsystem.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "indexing.h" 
#include "vswf.h" // qpms_vswf_set_spec_find_uvswfi()
#include "qpms_error.h"


#define HBAR (1.05457162825e-34)
#define ELECTRONVOLT (1.602176487e-19)
#define SCUFF_OMEGAUNIT (3e14)

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


