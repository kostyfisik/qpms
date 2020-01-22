#include "scatsystem.h"
#include <assert.h>
#include "qpms_error.h"

#define SQ(x) ((x)*(x))

typedef size_t ppid_t;
typedef size_t uoppid_t;

// Unordered pair ID. TODO get rid of redundancies
static inline uoppid_t uopairid(qpms_ss_pi_t pn_total, qpms_ss_pi_t p1, qpms_ss_pi_t p2) {
  return pn_total * p1 + p2;
}

// Unordered pair ID count. TODO get rid of redundancies.
static inline uoppid_t uopairarr_len(qpms_ss_pi_t pn_total) {
  return SQ(pn_total);
}

// Ordered pair ID.
static inline uoppid_t pairid(qpms_ss_pi_t pn_total, qpms_ss_pi_t p1, qpms_ss_pi_t p2) {
  return pn_total * p1 + p2;
}

// Ordered pair ID count.
static inline uoppid_t pairarr_len(qpms_ss_pi_t pn_total) {
  return SQ(pn_total);
}

typedef struct qpms_scatsys_translation_booster {
    double *r; // Unique distances array, indices are ppid_t
    size_t r_count; // Number of different interparticle distances (length of r[])
    size_t *r_map; // maps pairs to the corresponding distances (index of uoppid_t type)
} qpms_scatsys_translation_booster_t;

struct uoppid_r_pair {
	double r;
	uoppid_t id;
}

/// Sorts an array and throws away duplicit elements.
/**
 * The unique elements (according to the compare function) 
 * are aligned to the left of the original array.
 *
 * \returns Number of kept unique array members
 */
static size_t sort_and_eliminate(void *base, size_t nmemb, size_t size,
		int (*compar)(const void *, const void *)) {
  if (nmemb = 0) return 0; // corner case
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

struct qpms_scatsys_translation_booster *qpms_scatsys_translation_booster_create(
        const qpms_scatsys_ss *ss) {
  const qpms_ss_pi_t np = ss->p_count;
  struct qpms_scatsys_translation_booster *b;
  QPMS_CRASHING_MALLOC(b, sizeof(struct qpms_scatsys_translation_booster));
  QPMS_CRASHING_MALLOC(b->r, sizeof(double) * uopairarr_len(np));
  for(qpms_ss_pi_t i = 0; i < np; ++i)
    for(qpms_ss_pi_t j = 0; j < np; ++i) // FIXME j < i when uopairid works as supposed
      b->r[uopairid(pn, i, j)] = cart3_dist(ss->p[i].pos, ss->p[j].pos);
  b->r_count = sort_and_eliminate(b->r, uopairrarr_len(np), sizeof(double));
  
  TODO;
}

