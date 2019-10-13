/** \file polynomials.h
 * \brief Basic operations with polynomials.
 *
 */
#ifndef QPMS_POLYNOMIALS_H
#include <gmp.h>

/// Polynomial with rational coeffs.
typedef struct qpq_t {
	int order;
	int offset;
	int capacity;
	mpq_t *coeffs;
} qpq_t;

/// Initiasise the coefficients array in qpq_t.
void qpq_init(qpq_t *x, int maxorder);

/// Deinitialise the coefficients array in qpq_t.
void qpq_clear(qpq_t *x);

/// Polynomial addition.
void qpq_add(qpq_t *sum, const qpq_t *addend1, const qpq_t *addend2);

/// Polynomial substraction.
void qpq_sub(qpq_t *difference, const qpq_t *minuend, const qpq_t *substrahend);

/// Polynomial multiplication.
void qpq_mul(qpq_t product, const qpq_t *multiplier, const qpq_t *multiplicand);


/// Polynomial with double coeffs.
typedef struct qpz_t {
	int order;
	int offset;
	int capacity;
	double *coeffs;
} qpz_t;

/// Initiasise the coefficients array in qpz_t.
void qpz_init(qpz_t *x, int maxorder);

/// Deinitialise the coefficients array in qpz_t.
void qpz_clear(qpz_t *x);

/// Polynomial addition.
void qpz_add(qpz_t *sum, const qpz_t *addend1, const qpz_t *addend2);

/// Polynomial substraction.
void qpz_sub(qpz_t *difference, const qpz_t *minuend, const qpz_t *substrahend);

/// Polynomial multiplication.
void qpz_mul(qpz_t product, const qpz_t *multiplier, const qpz_t *multiplicand);

/// Convert rational coefficient polynomial to double coefficient polynomial
void qpz_from_qpq(qpz_t *target, const qpq_t *src);



/// Table with pre-calculated Ferrers function coeffs.
typedef struct qpms_legendre_table_t {
	qpms_l_t lMax;
	// TODO
} qpms_legendre_table_t;

/// Constructor for qpms_legendre_table_t.
qpms_legendre_table_t *qpms_legendre_table_init(qpms_l_t lMax);

/// Destructor for qpms_legendre_table_t.
void qpms_legendre_table_free(qpms_legendre_table_t *);

/// Evaluates a Ferrers function.
double qpms_legendre_table_eval(const qpms_legendre_table_t *,
		qpms_l_t l, qpms_m_t m, double x);


// TODO pre-calculate also the products??

#endif //QPMS_POLYNOMIALS_H
