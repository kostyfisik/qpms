#include "polynomials.h"
#include <stdlib.h>
#include "qpms_error.h"
#include <stdbool.h>
#include <stdio.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))

#if 0
void mpzs_init(mpzs_t x) {
  mpz_init(x->_1);
  mpz_init(x->_2);
}

void mpzs_clear(mpzs_t x) {
  mpz_clear(x->_1);
  mpz_clear(x->_2);
}

void mpzs_set(mpzs_t x, const mpzs_t y) {
  mpz_set(x->_1, y->_1);
  mpz_set(x->_2, y->_2);
}

// Compares the square-rooted part of mpzs_t, so we can use it as a key in a search tree.
// Probably deprecated, since we now use hashes.
static int mpzs_cmp2(void *op1, void *op2) {
  mpz_cmpabs(((mpzs_t) op1)->_2, ((mpzs_t) op2)->_2);
}
#endif

/// Type representing a number of form \f$ a \sqrt{b}; a \in \ints, b \in \nats \f$.
/** This includes UT hash handle structure */
typedef struct _qp_mpzs_hashed {
        mpz_t _1; ///< The integer factor \f$ a \f$.
        mpz_t _2; ///< The square-rooted factor \f$ b \f$. Always positive.
        UT_hash_handle hh;
} mpzs_hh_t[1];

void mpzs_hh_init(mpzs_hh_t x);
void mpzs_hh_clear(mpzs_hh_t x);
void mpzs_hh_set(mpzs_hh_t x, const mpzs_hh_t y);
// void mpzs_hh_set_ui(mpzs_hh_t x, long integer_factor, unsigned long sqrt_part);


void mpzs_hh_init(mpzs_hh_t x) {
  mpz_init(x->_1);
  mpz_init(x->_2);
}

void mpzs_hh_clear(mpzs_hh_t x) {
  mpz_clear(x->_1);
  mpz_clear(x->_2);
}

void mpzs_hh_set(mpzs_hh_t x, const mpzs_hh_t y) {
  mpz_set(x->_1, y->_1);
  mpz_set(x->_2, y->_2);
}

//===== mpqs_t =====

void mpqs_init(mpqs_t x) {
  //x->sz = 0;
  x->nt = 0;
  mpq_init(x->f);
  mpq_set_ui(x->f, 1, 1);
}

void mpqs_clear_num(struct _qp_mpqs *x) {
  struct _qp_mpzs_hashed *current, *tmp;
  HASH_ITER(hh, x->nt, current, tmp) {
    HASH_DEL(x->nt, current);
    free(current);
    //x->sz--;
  }
}

void mpqs_clear(mpqs_t x) {
  mpqs_clear_num(x);
  mpq_clear(x->f);
  //x->sz = 0;
  x->nt = 0;
}

void mpqs_nt_append(mpqs_t x, const struct _qp_mpzs_hashed *numelem) {
  struct _qp_mpzs_hashed *n;
  QPMS_CRASHING_MALLOC(n, sizeof(mpzs_hh_t));
  mpzs_hh_init(n);
  mpzs_hh_set(n, numelem);
  HASH_ADD_KEYPTR(hh, x->nt, mpz_limbs_read(n->_2),
      mpz_size(n->_2) * sizeof(mp_limb_t), n);
}

void mpqs_set_z(mpqs_t x, const mpz_t num1, const mpz_t num2,
    const mpz_t den) {
  mpqs_clear_num(x);
  if (mpz_sgn(num1) == 0 || mpz_sgn(num2) == 0) return;
  mpz_abs(mpq_numref(x->f), num1);
  mpz_abs(mpq_denref(x->f), den);
  mpq_canonicalize(x->f);
  mpzs_hh_t tmp; mpzs_hh_init(tmp);
  mpz_set_si(tmp->_1, mpz_sgn(num1) * mpz_sgn(den));
  mpz_set(tmp->_2, num2);
  mpqs_nt_append(x, tmp);
  mpzs_hh_clear(tmp);
}

void mpqs_set_si(mpqs_t x, long num1, unsigned long num2,
    unsigned long den) {
  mpqs_clear_num(x);
  if(num1 == 0 || num2 == 0) return;
  mpz_set_ui(mpq_numref(x->f), labs(num1));
  mpz_set_ui(mpq_denref(x->f), den);
  mpq_canonicalize(x->f);
  mpzs_hh_t tmp; mpzs_hh_init(tmp);
  mpz_set_si(tmp->_1, (num1 < 0) ? -1 : 1);
  mpz_set_ui(tmp->_2, num2);
  mpqs_nt_append(x, tmp);
  mpzs_hh_clear(tmp);
}

void mpqs_nt_add(mpqs_t x, const struct _qp_mpzs_hashed *addend) {
  struct _qp_mpzs_hashed *s;
  HASH_FIND(hh, x->nt, mpz_limbs_read(addend->_2),
      mpz_size(addend->_2), s);
  if (!s) mpqs_nt_append(x, addend); // if not found
  else {
    mpz_add(s->_1, s->_1, addend->_1);
    if(!mpz_sgn(s->_1)) { // If zero, annihilate
      HASH_DEL(x->nt, s);
      mpzs_hh_clear(s);
      free(s);
    }
  }
}

void mpqs_init_set(mpqs_t dest, const mpqs_t src) {
  mpqs_init(dest);
  mpq_set(dest->f, src->f);
  struct _qp_mpzs_hashed *numitem, *tmp;
  HASH_ITER(hh, src->nt, numitem, tmp) {
    mpqs_nt_append(dest, numitem);
  }
}

void mpqs_set(mpqs_t x, const mpqs_t y) {
  mpqs_clear(x);
  mpqs_init_set(x, y);
}

void mpqs_neg_inplace(mpqs_t x) {
  for(struct _qp_mpzs_hashed *n = (x->nt)->hh.next; n != NULL; n = n->hh.next) {
    mpz_neg(n->_1, n->_1);
  }
}

void mpqs_neg(mpqs_t negated_operand, const mpqs_t operand) {
  mpqs_set(negated_operand, operand);
  mpqs_neg_inplace(negated_operand);
}

void mpqs_nt_gcd(mpz_t gcd, const mpqs_t x) {
  if(x->nt) {
    mpz_set(gcd, x->nt->_1);
    for(struct _qp_mpzs_hashed *n = (x->nt)->hh.next; n != NULL; n = n->hh.next) {
      mpz_gcd(gcd, gcd, n->_1);
    }
  } else 
    mpz_set_ui(gcd, 1);
}

void mpqs_canonicalise(mpqs_t x) {
  mpz_t gcd; mpz_init(gcd);
  mpqs_nt_gcd(gcd, x);
  mpz_mul(mpq_numref(x->f), mpq_numref(x->f), gcd);
  mpz_clear(gcd);
  mpq_canonicalize(x->f);
}

void mpqs_add(mpqs_t sum_final, const mpqs_t x, const mpqs_t y) {
  mpqs_t sum; mpqs_init(sum);
  mpqs_set(sum, x);

  // Denominators gcd
  mpz_t den_gcd; mpz_init(den_gcd);
  mpz_gcd(den_gcd, mpq_denref(x->f), mpq_denref(y->f));
  // Common denominator
  mpz_divexact(mpq_denref(sum->f), mpq_denref(sum->f), den_gcd);
  mpz_mul(mpq_denref(sum->f), mpq_denref(sum->f), mpq_denref(y->f));
  
  // Left operand numerator factor
  mpz_t tmp; mpz_init(tmp);
  mpz_divexact(tmp, mpq_denref(y->f), den_gcd);
  mpz_mul(tmp, tmp, mpq_numref(x->f));
  /* Distribute the factor to numerator elements; this approach might be
   * be suboptimal if x is a complicated expression
   * and y is quite simple. Maybe optimise later 
   */
  for(struct _qp_mpzs_hashed *n = sum->nt; n != NULL; n = n->hh.next)
    mpz_mul(n->_1, n->_1, tmp);
  mpz_set_ui(mpq_numref(sum->f), 1);
  
  // At this point, sum is equal to x but in a form prepared to
  // take summands from y.

  // Right operand numerator factor
  mpz_divexact(tmp, mpq_denref(x->f), den_gcd);
  mpz_mul(tmp, tmp, mpq_numref(y->f));

  mpzs_hh_t addend; mpzs_hh_init(addend);
  for(const struct _qp_mpzs_hashed *n = y->nt; n != NULL; n = n->hh.next) {
    mpz_mul(addend->_1, tmp, n->_1);
    mpz_set(addend->_2, n->_2);
    mpqs_nt_add(sum, addend);
  }
  mpzs_hh_clear(addend);
  mpz_clear(tmp);
  mpz_clear(den_gcd);
  
  mpqs_canonicalise(sum);
  mpqs_set(sum_final, sum);
  mpqs_clear(sum);
}

void mpqs_sub(mpqs_t dif, const mpqs_t x, const mpqs_t y) {
  // This is probably not optimal
  mpqs_t tmp; mpqs_init(tmp);
  mpqs_neg(tmp, y);
  mpqs_add(dif, x, tmp);
  mpqs_clear(tmp);
}

void mpqs_mul(mpqs_t product_final, const mpqs_t x, const mpqs_t y) {
  mpqs_t prod; mpqs_init(prod);

  mpz_mul(mpq_numref(prod->f), mpq_numref(x->f), mpq_numref(y->f));
  mpz_mul(mpq_denref(prod->f), mpq_denref(x->f), mpq_denref(y->f));

  mpz_t gcd; mpz_init(gcd); // common denominator of the sqrt to factor out
  mpz_t mx, my; mpz_init(mx); mpz_init(my);
  mpzs_hh_t addend; mpzs_hh_init(addend);

  for (const struct _qp_mpzs_hashed *nx = x->nt->hh.next; 
      nx != NULL; nx = nx->hh.next) {
    for (const struct _qp_mpzs_hashed *ny = x->nt->hh.next; 
        ny != NULL; ny = ny->hh.next) {
      mpz_gcd(gcd, nx->_2, ny->_2);
      mpz_divexact(mx, nx->_2, gcd);
      mpz_divexact(my, ny->_2, gcd);
      mpz_mul(addend->_2, mx, my);
      mpz_mul(addend->_1, ny->_1, nx->_1);
      mpz_mul(addend->_1, addend->_1, gcd);
      mpqs_nt_add(prod, addend);
    }
  }
  mpzs_hh_clear(addend);
  mpz_clear(mx); mpz_clear(my); mpz_clear(gcd);
  mpqs_canonicalise(prod);
  mpqs_set(product_final, prod);
  mpqs_clear(prod);
}


// Auxillary function to set a mpq_t to 0/1
static inline void mpq_zero(mpq_t q) {
// Maybe not the best way to set it to zero.
// Alternatively, we could use just mpz_set_si(mpq_numref(sum->coeffs[i - minoffset]), 0);
  mpq_clear(q);
  mpq_init(q);
}


// TODO to the template

static inline void qpq_zero(qpq_t *q) {
  qpq_clear(q);
}


#define TEMPLATE_QPQ
#include "polynomials.template"
#undef TEMPLATE_QPQ

#if 0
// Used in the template but perhaps not too useful to put in the header.
static inline void mpqs_set_si(mpqs_t q, int num, unsigned den) {mpq_set_si(q->_2, num, den) ;}

// TODO Put into the template
static inline void mpqs_zero(mpqs_t q) {
  mpqs_clear(q);
  mpqs_init(q);
}

static inline void qpqs_zero(qpqs_t *q) {
  qpqs_clear(q);
}

#define TEMPLATE_QPQS
#include "polynomials.template"
#undef TEMPLATE_QPQS


static void qpq_dbgprint(const qpq_t *p) {
  for(int n = p->order; n >= p->offset; --n)
    if(mpq_sgn(p->coeffs[n - p->offset]))
      gmp_printf("%+Qdx**%d  ", p->coeffs[n - p->offset], n);
  gmp_printf("[%d, %d, %d]\n", p->capacity, p->order, p->offset);
}
#endif

void qpq_set_elem_si(qpq_t *p, const int exponent, const long int num, const unsigned long int den) {
  mpq_t q;
  mpq_init(q);
  mpq_set_si(q, num, den);
  qpq_set_elem(p, exponent, q);
  mpq_clear(q);
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

