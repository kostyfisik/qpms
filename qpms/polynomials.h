/** \file polynomials.h
 * \brief Basic operations with polynomials.
 *
 */
#ifndef QPMS_POLYNOMIALS_H
#include <gmp.h>

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


/// Square root of a rational number.
typedef struct mpqs_struct_t {
	mpq_t _2; 
} mpqs_t[1];

static inline void mpqs_init(mpqs_t q) {mpq_init(q->_2);}
static inline void mpqs_clear(mpqs_t q) {mpq_clear(q->_2);}
static inline void mpqs_mul(mpqs_t product, const mpqs_t multiplier, const mpqs_t multiplicand) {
	mpq_mul(product->_2, multiplier->_2, multiplicand->_2);
}
static inline void mpqs_div(mpqs_t quotient, const mpqs_t dividend, const mpqs_t divisor) {
	mpq_div(quotient->_2, dividend->_2, divisor->_2);
}
static inline void mpqs_set(mpqs_t copy, const mpqs_t orig) { mpq_set(copy->_2, orig->_2); }
static inline void mpqs_set_mpq(mpqs_t qs, const mpq_t q) {
	mpq_mul(qs->_2, q, q);
}
static inline int mpqs_sgn(const mpqs_t q) { return mpq_sgn(q->_2); }
static inline void mpqs_canonicalize(mpqs_t q) { mpq_canonicalize(q->_2); }
static inline void mpqs_swap(mpqs_t a, mpqs_t b) { mpq_swap(a->_2, b->_2); }

/// Polynomial with rational square root coeffs.
// TODO more docs about initialisation etc.
typedef struct qpqs_t {
	int order;
	int offset;
	int capacity;
	mpqs_t *coeffs;
} qpqs_t;
const static qpqs_t QPQS_ZERO = {-1, 0, 0, NULL};

/// Initiasise the coefficients array in qpqs_t.
/** Do not use on qpqs_t that has already been initialised
 * (and not recently cleared),
 * otherwise you can get a memory leak.
 */
void qpqs_init(qpqs_t *p, int capacity);

/// Extend capacity of a qpqs_t instance.
/** If the requested new_capacity is larger than the qpqs_t's
 * capacity, the latter is extended to match new_capacity.
 * Otherwise, nothing happend (this function does _not_ trim
 * the capacity).
 */
void qpqs_extend(qpqs_t *p, int new_capacity);

/// Shrinks the capacity to the minimum that can store the current polynomial.
void qpqs_shrink(qpqs_t *p);

/// Canonicalises the coefficients and (re)sets the correct degree.
void qpqs_canonicalise(qpqs_t *p);

void qpqs_set(qpqs_t *copy, const qpqs_t *orig);

void qpqs_set_elem(qpqs_t *p, int exponent, const mpqs_t coeff);
void qpqs_set_elem_si(qpqs_t *p, int exponent, long numerator, unsigned long denominator);
void qpqs_get_elem(mpqs_t coeff, const qpqs_t *p, int exponent);
/** \returns zero if the result fits into long / unsigned long; non-zero otherwise. */
int qpqs_get_elem_si(long *numerator, unsigned long *denominator, const qpqs_t *p, int exponent);

/// Deinitialise the coefficients array in qpqs_t.
void qpqs_clear(qpqs_t *p);



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
