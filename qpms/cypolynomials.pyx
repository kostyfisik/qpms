from libc.limits cimport INT_MIN, INT_MAX, UINT_MAX, LONG_MIN, LONG_MAX, ULONG_MAX
from gmpy2 cimport * # requires gmpy2>=2.1.0 or later (pip has currently 2.0.something)

cdef extern from "gmp.h":
    # cdef struct __mpq_struct:
    #    pass
    # ctypedef __mpq_struct mpq_t[1]
    # cdef struct __mpz_struct:
    #     pass
    # ctypedef __mpz_struct mpz_t[1]
    double mpz_get_d(const mpz_t op)
    void mpq_init(mpq_t q)
    void mpq_clear(mpq_t q)
    int mpq_sgn(mpq_t q)

cdef extern from "polynomials.h":
    cdef struct qpq_t:
        int order
        int offset
        int capacity
        mpq_t *coeffs

    void qpq_init(qpq_t *p, int capacity)
    void qpq_extend(qpq_t *p, int new_capacity)
    void qpq_set(qpq_t *copy, const qpq_t *orig)
    void qpq_set_elem(qpq_t *p, int exponent, const mpq_t coeff)
    void qpq_set_elem_si(qpq_t *p, int exponent, long numerator, unsigned long denominator)
    void qpq_get_elem(mpq_t coeff, const qpq_t *p, int exponent)
    int qpq_get_elem_si(long *numerator, unsigned long *denominator, const qpq_t *p, int exponent)
    void qpq_clear(qpq_t *p)
    void qpq_add(qpq_t *sum, const qpq_t *addend1, const qpq_t *addend2)
    void qpq_sub(qpq_t *difference, const qpq_t *minuend, const qpq_t *substrahend)
    void qpq_mul(qpq_t *product, const qpq_t *multiplier, const qpq_t *multiplicand)
    void qpq_deriv(qpq_t *dPdx, const qpq_t *P)
    bint qpq_nonzero(const qpq_t *)


import_gmpy2() # needed to initialise the C-API

cdef class qpq:
    """ Polynomials with rational coefficients """
    cdef qpq_t p

    def __cinit__(self, *args, **kwargs):
        qpq_init(&self.p, 0)

    def __init__(self, *args, **kwargs):
        cdef int offset, order
        cdef dict thedict
        cdef mpq coeff
        if len(args) > 0 and isinstance(args[0], dict) and len(args[0]) > 0:
            thedict = args[0]
            keys = thedict.keys()
            for key in keys:
                if not isinstance(key, int) or key < 0 or key > INT_MAX:
                    raise TypeError("Coefficient dictionary keys must be non-negative integers.")
            offset = min(keys)
            order = max(keys)
            qpq_extend(&self.p, order - offset + 1)
            self.p.order = order
            self.p.offset = offset
            for key, val in thedict.items():
                if MPQ_Check(val):
                    coeff = val
                else:
                    coeff = mpq(val)
                qpq_set_elem(&self.p, key, MPQ(coeff))
        return
    
    def __dealloc__(self):
        qpq_clear(&self.p)

    def __getitem__(self, key):
        # Only one coefficient a time supported right now
        cdef mpq q
        cdef mpq_t cq
        if isinstance(key, int):
            if key >= 0 and key <= INT_MAX:
                """ Can we alternatively do it like:
                q = mpq() # or GMPy_MPQ_New(NULL) ?
                qpq_get_elem(MPQ(q), &self.p, key)
                return q
                ? """
                mpq_init(cq)
                qpq_get_elem(cq, &self.p, key)
                q = GMPy_MPQ_From_mpq(cq)
                mpq_clear(cq)
                return q
            else: raise IndexError("Only non-negative int exponents (indices) allowed.")
        else: raise TypeError("Only integer exponents (indices) allowed.")

    def __setitem__(self, key, value):
        # Only one coefficient a time supported right now
        if not MPQ_Check(value):
            value = mpq(value)
        if (isinstance(key, int)):
            if key >= 0 and key <= INT_MAX:
                qpq_set_elem(&self.p, key, MPQ(value))
            else: raise IndexError("Only non-negative int exponents (indices) allowed.")
        else: raise TypeError("Only integer exponents (indices) allowed.")

    def __add__(qpq self, qpq other):
        cdef qpq result = qpq()
        qpq_add(&result.p, &self.p, &other.p)
        return result
    def __sub__(qpq self, qpq other):
        cdef qpq result = qpq()
        qpq_sub(&result.p, &self.p, &other.p)
        return result
    def __mul__(qpq self, qpq other):
        cdef qpq result = qpq()
        qpq_mul(&result.p, &self.p, &other.p)
        return result
    def derivative(self):
        cdef qpq result = qpq()
        qpq_deriv(&result.p, &self.p)
        return result

    @property
    def order(self):
        return self.p.order
    @property
    def offset(self):
        return self.p.offset
    @property
    def capacity(self):
        return self.p.capacity

    def coeffdict(self):
        """ Returns a dictionary of all non-zero coefficients """
        cdef int exponent, i
        cdef dict result = dict()
        cdef mpq coeff
        if not qpq_nonzero(&self.p):
            return result
        for exponent in range(self.p.offset, self.p.order + 1):
            i = exponent - self.p.offset
            if mpq_sgn(self.p.coeffs[i]):
                result[i] = GMPy_MPQ_From_mpq(self.p.coeffs[i])
        return result

