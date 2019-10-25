#include "polynomials.h"
#include <stdlib.h>
#include "qpms_error.h"
#include <stdbool.h>
#include <stdio.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))


static void qpq_dbgprint(const qpq_t *p) {
  for(int n = p->order; n >= p->offset; --n)
    if(mpq_sgn(p->coeffs[n - p->offset]))
      gmp_printf("%+Qdx**%d  ", p->coeffs[n - p->offset], n);
  gmp_printf("[%d, %d, %d]\n", p->capacity, p->order, p->offset);
}

// Auxillary function to set a mpq_t to 0/1
static inline void mpq_zero(mpq_t q) {
// Maybe not the best way to set it to zero.
// Alternatively, we could use just mpz_set_si(mpq_numref(sum->coeffs[i - minoffset]), 0);
  mpq_clear(q);
  mpq_init(q);
}

static inline void qpq_zero(qpq_t *q) {
  qpq_clear(q);
}

// qpq_t internal consistency check
static inline void qpq_cc(const qpq_t *p) {
  if (!p->coeffs) return;
  QPMS_ENSURE(p->capacity >= p->order - p->offset + 1, "qpq_t inconsistency detected; have you initialised it properly?");
}

_Bool qpq_nonzero(const qpq_t *p) {
  qpq_cc(p);
  if (p->capacity == 0) return false;
  
  for(int i = 0; i <= p->order - p->offset; ++i) 
    if (mpq_sgn(p->coeffs[i]))
      return true;
  return false;
}

void qpq_init(qpq_t *p, int capacity) {
  *p = QPQ_ZERO;
  if (capacity <= 0) 
    return;
  QPMS_CRASHING_MALLOC(p->coeffs, capacity * sizeof(mpq_t));
  for(int i = 0; i < capacity; ++i) 
    mpq_init(p->coeffs[i]);
  p->capacity = capacity;
}

void qpq_extend(qpq_t *p, int cap) {
  QPMS_ENSURE(p->capacity >= 0, "Got polynomial with negative capacity (%d). Is this a manually allocated one?", p->capacity);
  if (cap > 0 && cap > p->capacity) {
    QPMS_CRASHING_REALLOC(p->coeffs, sizeof(mpq_t) * cap);
    for(int i = p->capacity; i < cap; ++i) 
      mpq_init(p->coeffs[i]);
    p->capacity = cap;
  }
}

void qpq_shrink(qpq_t *p) {
  if (p->capacity > 0) {
      for(int n = p->capacity - 1; n > p->order - p->offset; --n)
        mpq_clear(p->coeffs[n]);
      p->capacity = p->order - p->offset + 1;
      if (p->capacity > 0) 
        QPMS_CRASHING_REALLOC(p->coeffs, p->capacity * sizeof(mpq_t));
  }
}

void qpq_canonicalise(qpq_t *p) {
  qpq_cc(p);
  // Lower the order if necessary (here one can get -1 if the polynomial is 0)
  while (p->order >= p->offset && !mpq_sgn(p->coeffs[p->order - p->offset])) --p->order;
  if (p->order < p->offset)
    qpq_zero(p);
  else { // remove the lowest-order coefficients which are in fact zero.
    int i = 0;
    while (i <= p->order - p->offset && !mpq_sgn(p->coeffs[i])) ++i;
    if (i > 0)
      for (int j = 0; j <= p->order - p->offset - i; ++j)
        mpq_swap(p->coeffs[j], p->coeffs[j+i]);
    p->offset += i;
    // canonicalise the fractions
    for (i = 0; i <= p->order - p->offset; ++i) 
      mpq_canonicalize(p->coeffs[i]);
  }
}

void qpq_set(qpq_t *p, const qpq_t *orig) {
  if(p != orig) {
    const int order = orig->order, offset = orig->offset;
    qpq_extend(p, order - offset + 1);
    for (int i = orig->offset; i <= order; ++i)
      mpq_set(p->coeffs[i - offset], orig->coeffs[i - offset]);
    p->offset = offset;
    p->order = order;
  }
}

void qpq_set_elem(qpq_t *p, const int exponent, const mpq_t coeff)  {
  QPMS_ENSURE(exponent >= 0, "Exponent must be non-negative, got %d", exponent);
  int offset = p->offset, order = p->order;
  if(mpq_sgn(coeff) == 0 && (exponent < offset || exponent > order)) 
    return; // exponent out of range, but zero needs not to be assigned explicitly
  if(exponent < p->offset) {
    qpq_extend(p, p->order - exponent + 1);
    offset = exponent;
    for(int i = order - offset; i >= p->offset - offset; --i) 
      mpq_swap(p->coeffs[i], p->coeffs[i - (p->offset - offset)]);
    for(int i = p->offset - offset - 1; i > 0; --i) 
      mpq_zero(p->coeffs[i]);
    p->offset = offset;
  }
  if(exponent > order) {
    qpq_extend(p, exponent - p->offset + 1);
    for(int i = p->order - p->offset + 1; i <= exponent - p->offset; ++i)
      mpq_zero(p->coeffs[i]);
    p->order = exponent;
  }
  mpq_set(p->coeffs[exponent - p->offset], coeff);
  return;
}

void qpq_set_elem_si(qpq_t *p, const int exponent, const long int num, const unsigned long int den) {
  mpq_t q;
  mpq_init(q);
  mpq_set_si(q, num, den);
  qpq_set_elem(p, exponent, q);
  mpq_clear(q);
}

void qpq_get_elem(mpq_t coeff, const qpq_t *p, const int exponent) {
  if (exponent < p->offset || exponent > p->order) 
    mpq_zero(coeff);
  else 
    mpq_set(coeff, p->coeffs[exponent-p->offset]);
}

int qpq_get_elem_si(long *num, unsigned long *den, const qpq_t *p, const int exponent) {
  mpq_t q;
  mpq_init(q);
  qpq_get_elem(q, p, exponent);
  int retval = 0;
  *num = mpz_get_si(mpq_numref(q));
  if (!mpz_fits_slong_p(mpq_numref(q)))
    retval += 1;
  *den = mpz_get_ui(mpq_denref(q));
  if (!mpz_fits_ulong_p(mpq_denref(q)))
    retval += 2;
  mpq_clear(q);
  return retval;
}

void qpq_clear(qpq_t *p) {
  if (p->capacity > 0) {
    for (int i = p->capacity - 1; i >= 0; --i)
      mpq_clear(p->coeffs[i]);
    free(p->coeffs);
  }
  *p = QPQ_ZERO;
}


// Auxillary function for lowering the offset
void qpq_lower_offset(qpq_t *p, int dec) {
  QPMS_ENSURE(dec >= 0, "Offset decrease must be positive, is %d!", dec);
  QPMS_ENSURE(dec <= p->offset, "Offset can't be pushed below 0, (offset=%d, decr.=%d!)", 
      p->offset,  dec);
  if(dec > 0) {
    qpq_extend(p, p->order - p->offset + dec + 1);
    for(int i = p->order; i >= p->offset; --i) 
      mpq_swap(p->coeffs[i+dec - p->offset], p->coeffs[i - p->offset]);
    p->offset -= dec;
    for(int i = dec - 1; i >= 0; --i)
      mpq_set_si(p->coeffs[i], 0, 1);
  }
}

void qpq_add(qpq_t *sum, const qpq_t *x, const qpq_t *y) {
  const int maxorder = MAX(x->order, y->order);
  const int minoffset = MIN(x->offset, y->offset);
  qpq_extend(sum, maxorder - minoffset + 1);
  /* if sum is actually some of the summands and that summand has higher
   * offset, we have to lower the offset.
   */
  if ((sum == x || sum == y) && sum->offset > minoffset) 
    qpq_lower_offset(sum, sum->offset - minoffset);

  for (int i = minoffset; i <= maxorder; ++i) {
    if (i - x->offset >= 0 && i <= x->order) {
      if (i - y->offset >= 0 && i <= y->order) 
        mpq_add(sum->coeffs[i - minoffset], 
        x->coeffs[i - x->offset], y->coeffs[i - y->offset]);
      else
        mpq_set(sum->coeffs[i - minoffset], x->coeffs[i - x->offset]);
    } else {
      if (i - y->offset >= 0 && i <= y->order)
        mpq_set(sum->coeffs[i - minoffset], y->coeffs[i - y->offset]);
      else {
        mpq_zero(sum->coeffs[i - minoffset]);
      }
    }
  }
  sum->offset = minoffset;
  sum->order = maxorder;
}
      
void qpq_sub(qpq_t *dif, const qpq_t *x, const qpq_t *y) {
  const int maxorder = MAX(x->order, y->order);
  const int minoffset = MIN(x->offset, y->offset);
  qpq_extend(dif, maxorder - minoffset + 1);
  /* if dif is actually some of the summands and that summand has higher
   * offset, we have to lower the offset.
   */
  if ((dif == x || dif == y) && dif->offset > minoffset) 
    qpq_lower_offset(dif, dif->offset - minoffset);

  for (int i = minoffset; i <= maxorder; ++i) {
    if (i - x->offset >= 0 && i <= x->order) {
      if (i - y->offset >= 0 && i <= y->order) 
        mpq_sub(dif->coeffs[i - minoffset], 
        x->coeffs[i - x->offset], y->coeffs[i - y->offset]);
      else
        mpq_set(dif->coeffs[i - minoffset], x->coeffs[i - x->offset]);
    } else {
      if (i - y->offset >= 0 && i <= y->order) {
        mpq_set(dif->coeffs[i - minoffset], y->coeffs[i - y->offset]);
        mpq_neg(dif->coeffs[i - minoffset], dif->coeffs[i - minoffset]);
      } else {
        mpq_zero(dif->coeffs[i - minoffset]);
      }
    }
  }
  dif->offset = minoffset;
  dif->order = maxorder;
}

void qpq_mul(qpq_t *p_orig, const qpq_t *x, const qpq_t *y) {

  const int maxorder = x->order + y->order;
  const int minoffset = x->offset + y->offset;
  
  qpq_t p[1];
  qpq_init(p, maxorder - minoffset + 1);
  mpq_t tmp; mpq_init(tmp);
  for (int xi = x->offset; xi <= x->order; ++xi)
    for (int yi = y->offset; yi <= y->order; ++yi) {
      mpq_mul(tmp, x->coeffs[xi - x->offset], y->coeffs[yi - y->offset]);
      mpq_add(p->coeffs[xi + yi - minoffset], p->coeffs[xi + yi - minoffset], tmp);
    }
  mpq_clear(tmp);
  p->order = maxorder;
  p->offset = minoffset;
  qpq_set(p_orig, p);
  qpq_clear(p);
} 

void qpq_div(qpq_t *q_orig, qpq_t *r_orig, const qpq_t *dend, const qpq_t *dor) {
  QPMS_ENSURE(q_orig != r_orig, 
      "Quotient and remainder must be _different_ instances of qpq_t.");


  // Split the divisor into "head" and "tail"
  qpq_t dor_tail[1]; qpq_init(dor_tail, dor->order - dor->offset); // divisor tail
  qpq_set(dor_tail, dor);
  qpq_canonicalise(dor_tail);
  QPMS_ENSURE(qpq_nonzero(dor_tail), "The divisor must be non-zero");

  const int dor_order = dor_tail->order;
  mpq_t dor_head; mpq_init(dor_head);
  mpq_set(dor_head, dor_tail->coeffs[dor_order - dor_tail->offset]);
  dor_tail->order--;

  qpq_t q[1];  qpq_init(q, dend->order - dor_order + 1);
  q->order = dend->order - dor_order;
  q->offset = 0;

  // Assign the dividend to r but with extended capacity (with zero offset)
  qpq_t r[1]; qpq_init(r, dend->order + 1);
  r->offset = 0;
  r->order = dend->order;
  for (int n = dend->offset; n <= dend->order; ++n) 
    mpq_set(r->coeffs[n], dend->coeffs[n - dend->offset]);

  qpq_t f[1]; qpq_init(f, 1);
  qpq_t ftail[1]; qpq_init(ftail, dor_tail->order - dor_tail->offset + 1);
  for(; r->order >= dor_order; --r->order) {
    // Compute the current order (r->order - dor_order) of q
    mpq_t * const hicoeff = &(q->coeffs[r->order - dor_order]);
    mpq_div(*hicoeff, r->coeffs[r->order], dor_head);
    mpq_canonicalize(*hicoeff);
    
    // Update the remainder
    f->offset = f->order = r->order - dor_order;
    mpq_set(f->coeffs[0], *hicoeff);
    qpq_mul(ftail, f, dor_tail);
    qpq_sub(r, r, ftail);
  }
  qpq_clear(ftail);
  qpq_clear(f);
  mpq_clear(dor_head);
  qpq_clear(dor_tail);

  qpq_canonicalise(r);
  qpq_canonicalise(q);
  qpq_set(q_orig, q);
  qpq_set(r_orig, r);
  qpq_clear(r); qpq_clear(f);
}


void qpq_deriv(qpq_t *dp, const qpq_t *p) {
  if (p->order <= 0) { // p is constant, dp is zero.
    qpq_clear(dp);
    return;
  }

  qpq_extend(dp, p->order - p->offset + (p->offset > 0));

  mpq_t qi;
  mpq_init(qi); // qi is now 0 / 1
  for(int i = p->offset + !p->offset; i <= p->order; ++i) {
    mpz_set_si(mpq_numref(qi), i); // qi is now i / 1
    mpq_mul(dp->coeffs[i-1], qi, p->coeffs[i]);
  }
  mpq_clear(qi);
  dp->order = p->order - 1;
  dp->offset = p->offset - 1 + !p->offset;
}


void qpq_legendroid_init(qpq_legendroid_t *p) {
  qpq_init(&p->p, 0);
  p->f = 0;
}

void qpq_legendroid_clear(qpq_legendroid_t *p) {
  qpq_clear(&p->p);
  p->f = 0;
}

void qpq_legendroid_mul(qpq_legendroid_t *p, const qpq_legendroid_t *a, const qpq_legendroid_t *b) {
  qpq_mul(&p->p, &a->p, &b->p);
  if (a->f && b->f) {
    // TODO make somehow a static representation of this constant polynomial
    qpq_t ff;
    qpq_init(&ff, 3);
    qpq_set_elem_si(&ff, 0, -1, 1);
    qpq_set_elem_si(&ff, 2, 1, 1);
    qpq_mul(&p->p, &p->p, &ff);
    qpq_clear(&ff);
  }
  p->f = !(a->f) != !(b->f);
}

