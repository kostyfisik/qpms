// Functionality to speedup translation matrix computations in large finite arrays
// by caching Bessel function values etc.
#include <cblas.h>
#include "scatsystem.h"
#include "translations_inlines.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vectors.h"
#include "qpms_error.h"
#include "qpms_specfunc.h"

#define SQ(x) ((x)*(x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

typedef size_t ppid_t;
typedef size_t uoppid_t;


// Unordered exclusive pair ID count. 
static inline uoppid_t uopairarr_len(qpms_ss_pi_t pn_total) {
  return pn_total * (pn_total - 1) / 2;
}

// Unordered exclusive pair ID.
static inline uoppid_t uopairid(qpms_ss_pi_t pn_total, qpms_ss_pi_t p1, qpms_ss_pi_t p2) {
  qpms_ss_pi_t hi, lo;
  QPMS_ASSERT(p1 != p2);
  if (p1 < p2) {
    lo = p1;
    hi = p2;
  } else {
    lo = p2;
    hi = p1;
  }
  QPMS_ASSERT(lo >= 0);
  QPMS_ASSERT(hi < pn_total);
  return lo + hi * (hi - 1) / 2;
}

// Ordered pair ID.
static inline ppid_t pairid(qpms_ss_pi_t pn_total, qpms_ss_pi_t p1, qpms_ss_pi_t p2) {
  return pn_total * p1 + p2;
}

// Ordered pair ID count.
static inline ppid_t pairarr_len(qpms_ss_pi_t pn_total) {
  return SQ(pn_total);
}

typedef struct qpms_scatsys_translation_booster {
  double *r; // Unique distances array, indices are ppid_t
  qpms_l_t *lMax_r; // lMaxes associated with the r's (use same indices as with r)
  size_t r_count; // Number of different interparticle distances (length of r[])
  size_t *r_map; // maps pairs to the corresponding distances (index of uoppid_t type)
  // Offsets of the Bessel function values: bessel_offsets_r[i] has the cumulative sums
  // of 2*lMax_r[j]+2, j < i. bessel_offsets_r[r_count] is the total length of the 
  // Bessel function "cache"
  size_t *bessel_offsets_r; 
} booster_t;

/// Sorts an array and throws away duplicit elements.
/**
 * The unique elements (according to the compare function) 
 * are aligned to the left of the original array.
 *
 * \returns Number of kept unique array members
 */
static size_t sort_and_eliminate(void *base, size_t nmemb, size_t size,
    int (*compar)(const void *, const void *)) {
  if (nmemb == 0) return 0; // corner case
  qsort(base, nmemb, size, compar);
  size_t left = 0;
  for (size_t src = 1; src < nmemb; ++src) {
    assert(left <= src);
    if (compar((const char *)base + left*size, (const char *)base + src*size)) {
      left += 1;
      if (left < src) // Avoid copying to itself (when it's not needed and memcpy behaviour undef'd
        memcpy((char *)base + left*size, (char *)base + src*size, size);
    }
  }
  return left + 1;
}

static int cmp_double(const void *aa, const void *bb) {
  const double a = *(double*)aa;
  const double b = *(double*)aa;
  if (a < b) return -1;
  if (a == b) return 0;
  if (a > b) return 1;
  QPMS_WTF; // NaN or similar
} 

booster_t *qpms_scatsys_translation_booster_create(
    const qpms_scatsys_t *ss) {
  const qpms_ss_pi_t np = ss->p_count;
  booster_t *b;
  QPMS_CRASHING_MALLOC(b, sizeof(*b));
  QPMS_CRASHING_MALLOC(b->r, sizeof(double) * uopairarr_len(np));
  for(qpms_ss_pi_t i = 0; i < np; ++i)
    for(qpms_ss_pi_t j = 0; j < i; ++j) 
      b->r[uopairid(np, i, j)] = cart3_dist(ss->p[i].pos, ss->p[j].pos);
  b->r_count = sort_and_eliminate(b->r, uopairarr_len(np), sizeof(*b->r), cmp_double);

  QPMS_CRASHING_REALLOC(b->r, b->r_count * sizeof(*b->r));
  QPMS_CRASHING_CALLOC(b->lMax_r, b->r_count, sizeof(*b->lMax_r));
  for(qpms_ss_pi_t i = 0; i < np; ++i)
    for(qpms_ss_pi_t j = 0; j < i; ++j) { 
      const uoppid_t pid = uopairid(np, i, j);
      const double r = cart3_dist(ss->p[i].pos, ss->p[j].pos);
      double *rhit = bsearch(&r, b->r, b->r_count, sizeof(*b->r), cmp_double);
      QPMS_ASSERT(rhit != NULL);
      QPMS_ASSERT(rhit >= b->r);
      const size_t ri = b->r_map[pid] = rhit - b->r;
      b->lMax_r[ri] = MAX(b->lMax_r[ri], 
          MAX(qpms_ss_bspec_pi(ss, i)->lMax, qpms_ss_bspec_pi(ss, j)->lMax));
    }

  QPMS_CRASHING_MALLOC(b->bessel_offsets_r, (b->r_count + 1) * sizeof(*b->bessel_offsets_r));
  size_t bessel_offset = 0;
  for(size_t ri = 0; ri < b->r_count; ++ri) {
    b->bessel_offsets_r[ri] = bessel_offset;
    bessel_offset += 2 * b->lMax_r[ri] + 2;
  }
  b->bessel_offsets_r[b->r_count] = bessel_offset;

  return b;
}

int qpms_ss_create_translation_cache(qpms_scatsys_t *ss, qpms_ss_caching_mode_t m) {
QPMS_ASSERT(ss);
  if (ss->tbooster) {
    QPMS_WARN("Translation cache already created?");
    return 0;
  }
  switch(m) {
    case QPMS_SS_CACHE_NEVER:
      return 0;
    case QPMS_SS_CACHE_AUTO:
      QPMS_WARN("Translation operator cache heuristics not implemented, creating the cache");
    case QPMS_SS_CACHE_ALWAYS:
      ss->tbooster = qpms_scatsys_translation_booster_create(ss);
      if (ss->tbooster) return 0;
      else {
          QPMS_WARN("Failed to create tranlation operator cache");
          return -1;
      }
    default:
      QPMS_WTF;
  }
  QPMS_WTF;
}

static qpms_errno_t qpms_scatsys_translation_booster_eval_bessels(
    const booster_t *b, complex double *target, complex double k // includes ref. ind. effect
    ) {
  for(size_t ri = 0; ri < b->r_count; ++ri) {
    QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(QPMS_HANKEL_PLUS, // Is there a case for different J?
          2*b->lMax_r[ri]+1, k * b->r[ri], 
          target + b->bessel_offsets_r[ri]));
  }
  return QPMS_SUCCESS;
}

typedef struct qpms_scatsysw_translation_booster {
  // _Bool owned_by_ssw; // if False, this is not deallocated by parent ssw
  const booster_t *b;
  complex double *bessels;
} boosterw_t;

boosterw_t *qpms_scatsysw_translation_booster_create(
    const qpms_scatsys_at_omega_t *ssw) {
  QPMS_PARANOID_ASSERT(ssw->ss);
  const booster_t *const b = ssw->ss->tbooster;
  QPMS_ASSERT(b);
  boosterw_t *bw;
  QPMS_CRASHING_MALLOC(bw, sizeof(*bw));

  // Evaluate bessel functions
  QPMS_CRASHING_MALLOC(bw->bessels, b->bessel_offsets_r[b->r_count] * sizeof(*bw->bessels));
  for(size_t ri = 0; ri < b->r_count; ++ri) {
    QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(QPMS_HANKEL_PLUS, // is there a case for other J?
          b->bessel_offsets_r[ri+1]-b->bessel_offsets_r[ri]-1,
          b->r[ri] * ssw->wavenumber,
          bw->bessels + b->bessel_offsets_r[ri]));
  }

  bw->b = b;
  return bw;
}

void qpms_scatsysw_translation_booster_free(boosterw_t *bw) {
  if (bw) {
    free (bw->bessels);
  }
  free (bw);
}

complex double *qpms_scatsysw_build_modeproblem_matrix_full_boosted(
    /// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
    complex double *target,
    const qpms_scatsys_at_omega_t *ssw) 
{
  QPMS_ASSERT(ssw->translation_cache && ssw->ss->tbooster);
  const qpms_scatsys_t *const ss = ssw->ss;
  const booster_t *const b = ss->tbooster;
  const boosterw_t *const bw = ssw->translation_cache;
  const qpms_trans_calculator *const c = ss->c;

  const complex double k = ssw->wavenumber;
  const size_t full_len = ss->fecv_size;
  if(!target)
    QPMS_CRASHING_MALLOC(target, SQ(full_len) * sizeof(complex double));
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, SQ(ss->max_bspecn) * sizeof(complex double));
  memset(target, 0, SQ(full_len) * sizeof(complex double)); //unnecessary?
  double legendre_buf[gsl_sf_legendre_array_n(2*c->lMax + 1)]; //VLA, workspace for legendre arrays
  const complex double zero = 0, minusone = -1;
  { // Non-diagonal part; M[piR, piC] = -T[piR] S(piR<-piC)
    size_t fullvec_offsetR = 0;
    for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) {
      const qpms_vswf_set_spec_t *bspecR = ssw->tm[ss->p[piR].tmatrix_id]->spec;
      const cart3_t posR = ss->p[piR].pos;
      size_t fullvec_offsetC = 0;
      // dest particle T-matrix
      const complex double *tmmR = ssw->tm[ss->p[piR].tmatrix_id]->m;
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) {
        const qpms_vswf_set_spec_t *bspecC = ssw->tm[ss->p[piC].tmatrix_id]->spec;
        if(piC != piR) { // The diagonal will be dealt with later.
          uoppid_t pid = uopairid(ss->p_count, piC, piR);
          const cart3_t posC = ss->p[piC].pos;
          const sph_t dlj = cart2sph(cart3_substract(posR, posC));
          const size_t ri = b->r_map[pid];
          QPMS_PARANOID_ASSERT(dlj.r == b->r[ri]);
          const qpms_l_t pair_lMax = b->lMax_r[pid];
          const qpms_y_t pair_nelem = qpms_lMax2nelem(pair_lMax);
          { // this replaces qpms_trans_calculator_get_trans_array():
            // R is dest, C is src
            QPMS_PARANOID_ASSERT(c->normalisation == bspecC->norm && c->normalisation == bspecR->norm);
            QPMS_PARANOID_ASSERT(c->lMax >= bspecC->lMax && c->lMax >= bspecR->lMax);
            QPMS_PARANOID_ASSERT(bspecC->lMax_L < 0 && bspecR->lMax_L < 0);
            complex double A[pair_nelem][pair_nelem]; // VLA, TODO flatten
            complex double B[pair_nelem][pair_nelem]; // VLA, TODO flatten
            { // this replaces qpms_trans_calculator_get_AB_arrays() and ..._buf()
              const double costheta = cos(dlj.theta);
              QPMS_ENSURE_SUCCESS(gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,
                    2*c->lMax+1, costheta, -1, legendre_buf));
              const double * const legendres = legendre_buf;
              const complex double * const bessels = bw->bessels + b->bessel_offsets_r[ri];
              qpms_trans_calculator_get_AB_arrays_precalcbuf(c, pair_lMax, A[0], B[0],
                  /*deststride*/ pair_nelem, /*srcstride*/ 1, dlj.phi, bessels, legendres);
              qpms_trans_array_from_AB(tmp,// tmp is S(piR<-piC) 
                  bspecR, bspecC->n, bspecC, 1, A[0], B[0], pair_lMax);
            }
          }
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              bspecR->n /*m*/, bspecC->n /*n*/, bspecR->n /*k*/,
              &minusone/*alpha*/, tmmR/*a*/, bspecR->n/*lda*/,
              tmp/*b*/, bspecC->n/*ldb*/, &zero/*beta*/,
              target + fullvec_offsetR*full_len + fullvec_offsetC /*c*/,
              full_len /*ldc*/);
        }
        fullvec_offsetC += bspecC->n;
      }
      fullvec_offsetR += bspecR->n;
    }
  }
  // diagonal part M[pi,pi] = +1
  for (size_t i = 0; i < full_len; ++i) target[full_len * i + i] = +1;

  free(tmp);
  return target;
}
