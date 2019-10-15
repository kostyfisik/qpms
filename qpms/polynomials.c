#include "polynomials.h"
#include <stdlib.h>
#include "qpms_error.h"
#include <stdbool.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))

// qpq_t internal consistency check
static inline void qpq_cc(const qpq_t *p) {
  if (!p->coeffs) return;
  QPMS_ENSURE(p->capacity >= p->order - p->offset + 1, "qpq_t inconsistency detected; have you initialised it properly?");
}

_Bool qpq_nonzero(const qpq_t *p) {
  qpq_cc(p);
  if (p->capacity <= 0) return false;
  
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
  if (cap > 0 && cap > p->capacity) {
    QPMS_CRASHING_REALLOC(p->coeffs, sizeof(mpq_t) * cap);
    for(int i = p->capacity; i < cap; ++i) 
      mpq_init(p->coeffs[i]);
    p->capacity = cap;
  }
}

void qpq_clear(qpq_t *p) {
  if (p->capacity > 0) {
    for (int i = p->capacity; i >= 0; --i)
      mpq_clear(p->coeffs[i]);
    free(p->coeffs);
  }
  *p = QPQ_ZERO;
}

void qpq_add(qpq_t *sum, const qpq_t *x, const qpq_t *y) {
  const int maxorder = MAX(x->order, y->order);
  const int minoffset = MIN(x->offset, y->offset);
  qpq_extend(sum, maxorder - minoffset + 1);
  for (int i = minoffset; i <= maxorder; ++i) {
    if (i - x->offset >= 0 && i <= x->order) {
      if (i - y->offset >= 0 && i <= y->order) 
        mpq_add(sum->coeffs[i - minoffset], 
        x->coeffs[i - x->offset], y->coeffs[i - y->offset]);
      else
        mpq_set(sum->coeffs[i - minoffset], x->coeffs[i - x->offset]);
    } else {
      if (i - y->offset >= 0 && i <= y->order)
        mpq_set(sum->coeffs[i - minoffset], y->coeffs[i - x->offset]);
      else {
        // Maybe not the best way to set it to zero.
        // Alternatively, we could use just mpz_set_si(mpq_numref(sum->coeffs[i - minoffset]), 0);
        mpq_clear(sum->coeffs[i - minoffset]);
        mpq_init(sum->coeffs[i - minoffset]);
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
  for (int i = minoffset; i <= maxorder; ++i) {
    if (i - x->offset >= 0 && i <= x->order) {
      if (i - y->offset >= 0 && i <= y->order) 
        mpq_sub(dif->coeffs[i - minoffset], 
        x->coeffs[i - x->offset], y->coeffs[i - y->offset]);
      else
        mpq_set(dif->coeffs[i - minoffset], x->coeffs[i - x->offset]);
    } else {
      if (i - y->offset >= 0 && i <= y->order) {
        mpq_set(dif->coeffs[i - minoffset], y->coeffs[i - x->offset]);
        mpq_neg(dif->coeffs[i - minoffset], dif->coeffs[i - minoffset]);
      } else {
        // Maybe not the best way to set it to zero.
        mpq_clear(dif->coeffs[i - minoffset]);
        mpq_init(dif->coeffs[i - minoffset]);
      }
    }
  }
  dif->offset = minoffset;
  dif->order = maxorder;
}

void qpq_mul(qpq_t *p, const qpq_t *x, const qpq_t *y) {
  const int maxorder = x->order + y->order;
  const int minoffset = x->offset + y->offset;
  // Easiest way to set p to a zero polynomial...
  qpq_clear(p);
  qpq_init(p, maxorder - minoffset + 1);
  for (int xi = x->offset; xi <= x->order; ++xi)
    for (int yi = y->offset; yi <= y->order; ++yi)
      mpq_mul(p->coeffs[xi + yi - minoffset],
          x->coeffs[xi - x->offset], y->coeffs[yi - y->offset]);
  p->order = maxorder;
  p->offset = minoffset;
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
    mpz_set_si(mpq_numref(qi), 1); // qi is now i / 1
    mpq_mul(dp->coeffs[i-1], qi, p->coeffs[i]);
  }
  mpq_clear(qi);
  dp->order = p->order - 1;
  dp->offset = p->offset - 1 + !p->offset;
}
