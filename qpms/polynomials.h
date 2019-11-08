/** \file polynomials.h
 * \brief Basic operations with polynomials.
 *
 */
#ifndef QPMS_POLYNOMIALS_H
#include <gmp.h>
#include "uthash.h"

/// Polynomial with rational coeffs.
// TODO more docs about initialisation etc.
typedef struct qpq_t {
	int order;
	int offset;
	int capacity;
	mpq_t *coeffs;
} qpq_t;

const static qpq_t QPQ_ZERO = {-1, 0, 0, NULL};

/// Initiasise the coefficients array in qpq_t.
/** Do not use on qpq_t that has already been initialised
 * (and not recently cleared),
 * otherwise you can get a memory leak.
 */
void qpq_init(qpq_t *p, int capacity);

/// Extend capacity of a qpq_t instance.
/** If the requested new_capacity is larger than the qpq_t's
 * capacity, the latter is extended to match new_capacity.
 * Otherwise, nothing happend (this function does _not_ trim
 * the capacity).
 */
void qpq_extend(qpq_t *p, int new_capacity);

/// Shrinks the capacity to the minimum that can store the current polynomial.
void qpq_shrink(qpq_t *p);

/// Canonicalises the coefficients and (re)sets the correct degree.
void qpq_canonicalise(qpq_t *p);

void qpq_set(qpq_t *copy, const qpq_t *orig);

void qpq_set_elem(qpq_t *p, int exponent, const mpq_t coeff);
void qpq_set_elem_si(qpq_t *p, int exponent, long numerator, unsigned long denominator);
void qpq_get_elem(mpq_t coeff, const qpq_t *p, int exponent);
/** \returns zero if the result fits into long / unsigned long; non-zero otherwise. */
int qpq_get_elem_si(long *numerator, unsigned long *denominator, const qpq_t *p, int exponent);

/// Deinitialise the coefficients array in qpq_t.
void qpq_clear(qpq_t *p);

/// Polynomial addition.
/** Supports operand and result pointer mixing. */
void qpq_add(qpq_t *sum, const qpq_t *addend1, const qpq_t *addend2);

/// Polynomial substraction.
/** Supports operand and result pointer mixing. */
void qpq_sub(qpq_t *difference, const qpq_t *minuend, const qpq_t *substrahend);

/// Polynomial multiplication.
/** Supports operand and result pointer mixing. */
void qpq_mul(qpq_t *product, const qpq_t *multiplier, const qpq_t *multiplicand);

/// Polynomial division with remainder.
/** Supports operand and result pointer mixing. */
void qpq_div(qpq_t *quotient, qpq_t *remainder, const qpq_t *dividend, const qpq_t *divisor);

/// Polynomial derivative.
/** Supports operand and result pointer mixing. */
void qpq_deriv(qpq_t *dPdx, const qpq_t *P);

/// Tests whether a polynomial is non-zero.
_Bool qpq_nonzero(const qpq_t *);



#if 0
/// Type representing a number of form \f$ a \sqrt{b}; a \in \ints, b \in \nats \f$.
typedef struct _qp_mpzs {
	mpz_t _1; ///< The integer factor \f$ a \f$.
	mpz_t _2; ///< The square-rooted factor \f$ b \f$. Always positive.
} mpzs_t[1];

void mpzs_init(mpzs_t x);
void mpzs_clear(mpzs_t x);
void mpzs_set(mpzs_t x, const mpzs_t y);
#endif


struct _qp_mpzs_hashed;

/// Sum of square roots of rational numbers.
/// Represented as \f$ \sum_s a_i \sqrt{b_i} / d \f$.
typedef struct _qp_mpqs {
	//int sz; ///< Used size of the numtree.
	struct _qp_mpzs_hashed *nt; ///< List of numerator components..
	mpq_t f; ///< Common rational factor. Always positive.
} mpqs_t[1];

void mpqs_init(mpqs_t x);
void mpqs_clear(mpqs_t x);
void mpqs_set(mpqs_t x, const mpqs_t y);
void mpqs_add(mpqs_t sum, const mpqs_t addend1, const mpqs_t addend2);
void mpqs_neg(mpqs_t negated_operand, const mpqs_t operand);
void mpqs_sub(mpqs_t difference, const mpqs_t minuend, const mpqs_t substrahend);
void mpqs_mul(mpqs_t product, const mpqs_t multiplier, const mpqs_t multiplicand);
//void mpqs_div_zs(mpqs_t quotient, const mpqs_t dividend, const mpzs_hh_t divisor);
void mpqs_clear_num(mpqs_t x); ///< Sets the numerator to zero, clearing the numerator tree.


/// A type representing a polynomial with rational coefficients times an optional factor \f$ \sqrt{1-x^2} \f$.
typedef struct qpq_legendroid_t {
	qpq_t p;
	_Bool f;
} qpq_legendroid_t;

void qpq_legendroid_init(qpq_legendroid_t *p);
void qpq_legendroid_clear(qpq_legendroid_t *p);

/// Polynomial multiplication.
void qpq_legendroid_mul(qpq_legendroid_t *product, const qpq_legendroid_t *multiplier, const qpq_legendroid_t *multiplicand);

/// Polynomial derivative.
void qpq_legendroid_deriv(qpq_legendroid_t *dP_dx, const qpq_legendroid_t *P);


/// Polynomial with double coeffs.
typedef struct qpz_t {
	int order;
	int offset;
	int capacity;
	double *coeffs;
} qpz_t;

/// Initiasise the coefficients array in qpz_t.
void qpz_init(qpz_t *p, int maxorder);

/// Deinitialise the coefficients array in qpz_t.
void qpz_clear(qpz_t *p);

/// Polynomial addition.
void qpz_add(qpz_t *sum, const qpz_t *addend1, const qpz_t *addend2);

/// Polynomial substraction.
void qpz_sub(qpz_t *difference, const qpz_t *minuend, const qpz_t *substrahend);

/// Polynomial multiplication.
void qpz_mul(qpz_t product, const qpz_t *multiplier, const qpz_t *multiplicand);

/// Convert rational coefficient polynomial to double coefficient polynomial
void qpz_from_qpq(qpz_t *target, const qpq_t *src);


#if 0 // This will go elsewhere
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
#endif

#endif //QPMS_POLYNOMIALS_H
