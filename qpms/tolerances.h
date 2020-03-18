/*! \file tolerances.h */
#ifndef QPMS_TOLERANCES_H
#define QPMS_TOLERANCES_H

// TODO DOC

typedef struct qpms_tolerance_spec_t {
	double atol;
	double rtol;
} qpms_tolerance_spec_t;

/// A rather arbitrary default tolerance.
static const qpms_tolerance_spec_t QPMS_TOLERANCE_DEFAULT = {.atol = 1e-9, .rtol =  1e-8};

#endif // QPMS_TOLERANCES_H
