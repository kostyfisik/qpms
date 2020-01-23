#ifndef TRANSLATIONS_INLINES_H
#define TRANSLATIONS_INLINES_H
#include "translations.h"
#include "indexing.h"

/// Rearranges the default-ordered "A,B" array elements into "bspec"-defined matrix.
// TODO DOC
static inline void qpms_trans_array_from_AB(
		complex double *t,
		const qpms_vswf_set_spec_t *const t_destspec,
		const size_t t_deststride,
		const qpms_vswf_set_spec_t *const t_srcspec,
		const size_t t_srcstride,
		const complex double *const A, const complex double *const B,
		/// A and B matrices' lMax.
		/** This also determines their size and stride: they are assumed to
		 * be square matrices of size `nelem * nelem` where
		 * `nelem = qpms_lMax2nelem(lMax_AB)`
		 */
		const qpms_l_t lMax_AB
		) {
	QPMS_PARANOID_ASSERT(lMax_AB >= t_srcspec->lMax && lMax_AB >= t_destspec->lMax);
	const qpms_y_t nelem_AB = qpms_lMax2nelem(lMax_AB);
	for (size_t desti = 0; desti < t_destspec->n; ++desti) {
		qpms_y_t desty; qpms_vswf_type_t destt;
		QPMS_ENSURE_SUCCESS_M(qpms_uvswfi2ty(t_destspec->ilist[desti], &destt, &desty),
				"Invalid u. vswf index %llx.", t_destspec->ilist[desti]);
		for (size_t srci = 0; srci < t_srcspec->n; ++srci){
			qpms_y_t srcy; qpms_vswf_type_t srct;
			QPMS_ENSURE_SUCCESS_M(qpms_uvswfi2ty(t_srcspec->ilist[srci], &srct, &srcy),
					"Invalid u. vswf index %llx.", t_srcspec->ilist[srci]);                                                                t[srci * t_srcstride + desti * t_deststride]
				= (srct == destt) ? A[desty*nelem_AB + srcy] : B[desty*nelem_AB + srcy];
		}
	}
}

int qpms_trans_calculator_get_AB_arrays_precalcbuf(const qpms_trans_calculator *c,
		qpms_y_t lMax, complex double *Adest, complex double *Bdest,
		size_t deststride, size_t srcstride, double kdlj_phi,
		const complex double *bessel_buf, const double *legendre_buf);
 
#endif
