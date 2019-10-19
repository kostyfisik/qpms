cdef extern from "<gmp.h>":
    cdef struct __mpq_struct:
        pass
    ctypedef __mpq_struct mpq_t[1]

    cdef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]

    double mpz_get_d(const mpz_t op)


cdef extern from "polynomials.h":
    cdef struct qpq_t:
        int order
        int offset
        int capacity
        mpq_t *coeffs

    void qpq_init(qpq_t *x, int capacity)
    void qpq_extend(qpq_t *x, int new_capacity)
    void qpq_set(qpq_t *copy, const qpq_t *orig)
    void qpq_clear(qpq_t 
