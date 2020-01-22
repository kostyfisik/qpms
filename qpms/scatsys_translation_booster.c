#include "scatsystem.h"

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

size_t sort_and_eliminate(void *base, size_t nmemb, size_t size,
		int (*compar)(const void *, const void *)) {
	qsort(base, nmemb, size, compar);
	_Bool eliminate;
	QPMS_CRASHING_CALLOC(eliminate, nmemb, sizeof(_Bool));
	size_t left = 0;
	
	TODO zde;

	free(eliminate);
	return left;
}


struct qpms_scatsys_translation_booster *qpms_scatsys_translation_booster_create(
        const qpms_scatsys_ss *ss) {
    const qpms_ss_pi_t np = ss->p_count;
    TODO;
}

