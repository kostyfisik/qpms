#ifndef QPMS_SCATSYS_PRIVATE_H
#define QPMS_SCATSYS_PRIVATE_H
/* 
 * This file includes some definitions shared between scatsystem.c
 * and scatsys_translation_booster.c that are not needed anywhere 
 * else.
 */

#include "scatsystem.h"

#include <pthread.h>

complex double *qpms_scatsysw_build_modeproblem_matrix_full_boosted(
		complex double *target, const qpms_scatsys_at_omega_t *ssw);
complex double *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_boosted(
		complex double *target_packed, const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri);
/// "private" destructor, called by qpms_scatsys_free()
void qpms_scatsys_translation_booster_free(struct qpms_scatsys_translation_booster *);
/// "private" constructor, use qpms_ss_create_translation_cache() instead.
struct qpms_scatsys_translation_booster *qpms_scatsys_translation_booster_create(
                const qpms_scatsys_t *ss);


struct qpms_scatsysw_translation_booster *
qpms_scatsysw_translation_booster_create(const qpms_scatsys_at_omega_t *ssw);

/// "private" destructor, called by qpms_scatsys_at_omega_free()
void qpms_scatsysw_translation_booster_free(struct qpms_scatsysw_translation_booster *);


struct qpms_scatsysw_build_modeproblem_matrix_irrep_packed_parallelR_thread_arg{
  const qpms_scatsys_at_omega_t *ssw;
  qpms_ss_pi_t *opistartR_ptr;
  pthread_mutex_t *opistartR_mutex;
  qpms_iri_t iri;
  complex double *target_packed;
};

void *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_parallelR_thread_boosted(void *arg);

#endif //QPMS_SCATSYS_PRIVATE_H
