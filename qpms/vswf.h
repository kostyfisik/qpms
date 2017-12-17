#ifndef QPMS_VSWF_H
#define QPMS_VSWF_H
#include "qpms_types.h"

// array of pi, tau auxillary function (see [1,(37)])
// TODO put this elsewhere (no need 
typedef struct {
	//qpms_normalization_t norm;
	//qpms_l_t lMax;
	qpms_y_t nelem;
	double *pi, *tau;
} qpms_pitau_t;
void qpms_pitau_t_free(qpms_pitau_t);
void qpms_pitau_t_pfree(qpms_pitau_t*);

#endif // QPMS_VSWF_H
