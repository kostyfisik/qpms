#ifndef QPMS_INDEXING_H
#define QPMS_INDEXING_H

#include "qpms_types.h"
#include <math.h>

static inline qpms_y_t qpms_mn2y(qpms_m_t m, qpms_l_t n) {
        return n * (n + 1) + m - 1;
}

static inline qpms_lm_t qpms_y2n(qpms_y_t y) {
        //return (sqrt(5+y)-2)/2; // the cast will truncate the fractional part, which is what we want
        return sqrt(y+1);
}

static inline qpms_m_t qpms_yn2m(qpms_y_t y, qpms_l_t n) {
        return y-qpms_mn2y(0,n);
}

static inline void qpms_y2mn_p(qpms_y_t y, qpms_m_t *m, qpms_l_t *n){
        *m=qpms_yn2m(y,*n=qpms_y2n(y));
}

static inline qpms_y_t qpms_lMax2nelem(qpms_l_t lmax){
	return lmax * ((qpms_y_t)lmax + 2);
}

#endif //QPMS_INDEXING_H
