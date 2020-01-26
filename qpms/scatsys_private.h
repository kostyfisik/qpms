#ifndef QPMS_SCATSYS_PRIVATE_H
#define QPMS_SCATSYS_PRIVATE_H

#include "scatsystem.h"

complex double *qpms_scatsysw_build_modeproblem_matrix_full_boosted(
		complex double *target, const qpms_scatsys_at_omega_t *ssw);

/// "private" destructor, called by qpms_scatsys_free()
void qpms_scatsys_translation_booster_free(struct qpms_scatsys_translation_booster *);
/// "private" constructor, use qpms_ss_create_translation_cache() instead.
struct qpms_scatsys_translation_booster *qpms_scatsys_translation_booster_create(
                const qpms_scatsys_t *ss);


struct qpms_scatsysw_translation_booster *
qpms_scatsysw_translation_booster_create(const qpms_scatsys_at_omega_t *ssw);

/// "private" destructor, called by qpms_scatsys_at_omega_free()
void qpms_scatsysw_translation_booster_free(struct qpms_scatsysw_translation_booster *);
#endif 
