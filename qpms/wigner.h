/*! \file wigner.h
 * \brief Quaternions and Wigner matrices
 */
#ifndef QPMS_WIGNER_H
#define QPMS_WIGNER_H

#include "qpms_types.h"
#include "vectors.h"
#include "tiny_inlines.h"


/// Conversion from the 4*double to the 2*complex quaternion.
// TODO is this really correct? 
// I.e. do the axis from moble's text match this convention?
static inline qpms_quat_t qpms_quat_2c_from_4d (qpms_quat4d_t q) {
	qpms_quat_t q2c = {q.c1 + I * q.ck, q.cj + I * q.ci};
	return q2c;
}

/// Conversion from the 2*complex to the 4*double quaternion.
// TODO is this really correct? 
// I.e. do the axis from moble's text match this convention?
static inline qpms_quat4d_t qpms_quat_4d_from_2c (qpms_quat_t q) {
	qpms_quat4d_t q4d = {creal(q.a), cimag(q.b), creal(q.b), cimag(q.a)};
	return q4d;
}

/// Quaternion multiplication.
/**
 * \f[ (P Q)_a = P_a Q_a - \bar P_b Q_b, \f]
 * \f[ (P Q)_b = P_b Q_a + \bar P_a Q_b. \f]
 */
static inline qpms_quat_t qpms_quat_mult(qpms_quat_t p, qpms_quat_t q) {
	qpms_quat_t r;
	r.a = p.a * q.a - conj(p.b) * q.b;
	r.b = p.b * q.a + conj(p.a) * q.b;
	return r;
}

/// Quaternion addition.
static inline qpms_quat_t qpms_quat_add(qpms_quat_t p, qpms_quat_t q) {
	qpms_quat_t r;
	r.a = p.a+q.a;
	r.b = p.b+q.b;
	return r;
}

/// Exponential function of a quaternion \f$e^Q$\f.
static inline qpms_quat_t qpms_quat_exp(const qpms_quat_t q) {
	const qpms_quat4d_t q4 = qpms_quat_4d_from_2c(q);
	const double vn = sqrt(q4.ci*q4.ci + q4.cj*q4.cj + q4.ck *q4.ck);
	const double ea = exp(q4.c1);
	const double cv = vn ? (ea*sin(vn)/vn) : ea; // "vector" part common prefactor
	const qpms_quat4d_t r4 = {ea * cos(vn), cv*q4.ci, cv*q4.cj, cv*q4.ck};
	return qpms_quat_2c_from_4d(r4);
}

/// Quaternion scaling with a real number.
static inline qpms_quat_t qpms_quat_rscale(double s, qpms_quat_t q) {
	qpms_quat_t r = {s * q.a, s * q.b};
	return r;
}

// quaternion "basis"
static const qpms_quat_t qpms_quat_1 = {1, 0};
static const qpms_quat_t qpms_quat_i = {I, 0};
static const qpms_quat_t qpms_quat_j = {0, 1};
static const qpms_quat_t qpms_quat_k = {0, I};

/// Quaternion conjugation.
static inline qpms_quat_t qpms_quat_conj(const qpms_quat_t q) {
	qpms_quat_t r = {conj(q.a), -q.b};
	return r;
}

/// Quaternion norm.
static inline double qpms_quat_norm(const qpms_quat_t q) {
	return sqrt(creal(q.a * conj(q.a) + q.b * conj(q.b)));
}

/// Norm of the quaternion imaginary (vector) part.
static inline double qpms_quat_imnorm(const qpms_quat_t q) {
	const double z = cimag(q.a), x = cimag(q.b), y = creal(q.b);
	return sqrt(z*z + x*x + y*y);
}

/// Quaternion normalisation to unit norm.
static inline qpms_quat_t qpms_quat_normalise(qpms_quat_t q) {
	double n = qpms_quat_norm(q);
	return qpms_quat_rscale(1/n, q);
}

/// Logarithm of a quaternion.
static inline qpms_quat_t qpms_quat_log(const qpms_quat_t q) {
	const double n = qpms_quat_norm(q);
	const double imnorm = qpms_quat_imnorm(q);
	if (imnorm != 0.) {
		const double vc = acos(creal(q.a)/n) / imnorm;
		const qpms_quat_t r = {log(n) + cimag(q.a)*vc*I,
			q.b*vc};
		return r;
	}
	else {
		const qpms_quat_t r = {log(n), 0};
		return r;
	}
}

/// Quaternion power to a real exponent.
static inline qpms_quat_t qpms_quat_pow(const qpms_quat_t q, const double exponent) {
	const qpms_quat_t qe = qpms_quat_rscale(exponent, 
			qpms_quat_log(q));
	return qpms_quat_exp(qe);
}

/// Quaternion inversion.
/** \f[ q^{-1} = \frac{q*}{|q|}. \f] */
static inline qpms_quat_t qpms_quat_inv(const qpms_quat_t q) {
	const double norm = qpms_quat_norm(q);
	return qpms_quat_rscale(1./(norm*norm),
			qpms_quat_conj(q));
}
			
/// Make a pure imaginary quaternion from a 3d cartesian vector.
static inline qpms_quat_t qpms_quat_from_cart3(const cart3_t c) {
	const qpms_quat4d_t q4 = {0, c.x, c.y, c.z};
	return qpms_quat_2c_from_4d(q4);
}

/// Make a 3d cartesian vector from the imaginary part of a quaternion.
static inline cart3_t qpms_quat_to_cart3(const qpms_quat_t q) {
	const qpms_quat4d_t q4 = qpms_quat_4d_from_2c(q);
	const cart3_t c = {q4.ci, q4.cj, q4.ck};
	return c;
}

/// Rotate a 3-dimensional cartesian vector using the quaternion/versor representation.
static inline cart3_t qpms_quat_rot_cart3(qpms_quat_t q, const cart3_t v) {
	q = qpms_quat_normalise(q);
	//const qpms_quat_t qc = qpms_quat_normalise(qpms_quat_pow(q, -1)); // implementation of _pow wrong!
	const qpms_quat_t qc = qpms_quat_conj(q);
	const qpms_quat_t vv = qpms_quat_from_cart3(v);
	return qpms_quat_to_cart3(qpms_quat_mult(q,
				qpms_quat_mult(vv, qc)));
}

/// Wigner D matrix element from a rotator quaternion for integer \a l.
/**
 * The D matrix are calculated using formulae (3), (4), (6), (7) from
 * http://moble.github.io/spherical_functions/WignerDMatrices.html
 */
complex double qpms_wignerD_elem(qpms_quat_t q, qpms_l_t l,
	       qpms_m_t	mp, qpms_m_t m);

/// A VSWF representation element of the O(3) group.
/**
 * TODO more doc.
 */
complex double qpms_vswf_irot_elem_from_irot3(
		const qpms_irot3_t q, ///< The O(3) element in the quaternion representation.
		qpms_l_t l, qpms_m_t mp, qpms_m_t m,
		bool pseudo ///< Determines the sign of improper rotations. True for magnetic waves, false otherwise.
		);


static inline int qpms_irot3_checkdet(const qpms_irot3_t p) {
	if (p.det != 1 && p.det != -1) abort();
	return 0;
}

/// Improper rotation multiplication.
static inline qpms_irot3_t qpms_irot3_mult(const qpms_irot3_t p, const qpms_irot3_t q) {
#ifndef NDEBUG
	qpms_irot3_checkdet(p);
	qpms_irot3_checkdet(q);
#endif
	const qpms_irot3_t r = {qpms_quat_normalise(qpms_quat_mult(p.rot, q.rot)), p.det*q.det};
	return r;
}

/// Improper rotation power \f$ p^n \f$.
static inline qpms_irot3_t qpms_irot3_pow(const qpms_irot3_t p, int n) {
#ifndef NDEBUG
	qpms_irot3_checkdet(p);
#endif
	const qpms_irot3_t r = {qpms_quat_normalise(qpms_quat_pow(p.rot, n)),
		p.det == -1 ? min1pow(n) : 1};
	return r;
}

/// Apply an improper rotation onto a 3d cartesian vector.
static inline cart3_t qpms_irot3_apply_cart3(const qpms_irot3_t p, const cart3_t v) {
#ifndef NDEBUG
	qpms_irot3_checkdet(p);
#endif
	return cart3_scale(p.det, qpms_quat_rot_cart3(p.rot, v));
}

#endif //QPMS_WIGNER_H
