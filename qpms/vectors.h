#ifndef VECTORS_H
#define VECTORS_H
#include <math.h>
#ifndef M_PI_2
#define M_PI_2 (1.570796326794896619231321691639751442098584699687552910487)
#endif
#include "qpms_types.h"

//static inline double vectors_h_sq(double x) {return x*x;}

static inline double cart3norm(const cart3_t v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

static inline sph_t cart2sph(const cart3_t cart) {
	sph_t sph;
	sph.r = cart3norm(cart);
	sph.theta = sph.r ? acos(cart.z / sph.r) : M_PI_2;
	sph.phi = atan2(cart.y, cart.x);
	return sph;
}

static inline cart3_t sph2cart(const sph_t sph) {
	cart3_t cart;
	double sin_th = sin(sph.theta);
	cart.x = sph.r * sin_th * cos(sph.phi);
	cart.y = sph.r * sin_th * sin(sph.phi);
	cart.z = sph.r * cos(sph.theta);
	return cart;
}

static inline cart3_t cart3_add(const cart3_t a, const cart3_t b) {
	cart3_t res = {a.x+b.x, a.y+b.y, a.z+b.z};
	return res;
}

static inline cart3_t cart3_scale(const double c, const cart3_t v) {
	cart3_t res = {c * v.x, c * v.y, c * v.z};
	return res;
}

static inline ccart3_t ccart3_scale(const complex  double c, const ccart3_t v) {
	ccart3_t res = {c * v.x, c * v.y, c * v.z};
	return res;
}

static inline csphvec_t csphvec_add(const csphvec_t a, const csphvec_t b) {
	csphvec_t res = {a.rc + b.rc, a.thetac + b.thetac, a.phic + b.phic};
	return res;
}

static inline csphvec_t csphvec_scale(complex double c, const csphvec_t v) {
	csphvec_t res = {c * v.rc, c * v.thetac, c * v.phic};
	return res;
}

static inline complex double csphvec_dotnc(const csphvec_t a, const csphvec_t b) {
	//N.B. no complex conjugation done here
	return a.rc * b.rc + a.thetac * b.thetac + a.phic * b.phic;
}

static inline double cart3_dot(const cart3_t a, const cart3_t b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

// equivalent to sph_loccart2cart in qpms_p.py
static inline ccart3_t csphvec2ccart(const csphvec_t sphvec, const sph_t at) {
	const double st = sin(at.theta);
	const double ct = cos(at.theta);
	const double sf = sin(at.phi);
	const double cf = cos(at.phi);
	const double rx = st * cf;
	const double ry = st * sf;
	const double rz = ct;
	const double tx = ct * cf;
	const double ty = ct * sf;
	const double tz = -st;
	const double fx = -sf;
	const double fy = cf;
	const double fz = 0.;
	ccart3_t res;
	res.x = rx * sphvec.rc + tx * sphvec.thetac + fx * sphvec.phic;
	res.y = ry * sphvec.rc + ty * sphvec.thetac + fy * sphvec.phic;
	res.z = rz * sphvec.rc + tz * sphvec.thetac + fz * sphvec.phic;
	return res;
}

static inline csphvec_t ccart2csphvec(const ccart3_t cartvec, const sph_t at) {
	// this chunk is copy-pasted from csphvec2cart, so there should be a better way...
	const double st = sin(at.theta);
	const double ct = cos(at.theta);
	const double sf = sin(at.phi);
	const double cf = cos(at.phi);
	const double rx = st * cf;
	const double ry = st * sf;
	const double rz = ct;
	const double tx = ct * cf;
	const double ty = ct * sf;
	const double tz = -st;
	const double fx = -sf;
	const double fy = cf;
	const double fz = 0.;
	csphvec_t res;
	res.rc     = rx * cartvec.x + ry * cartvec.y + rz * cartvec.z;
	res.thetac = tx * cartvec.x + ty * cartvec.y + tz * cartvec.z;
	res.phic   = fx * cartvec.x + fy * cartvec.y + fz * cartvec.z;
	return res;
}

#endif //VECTORS_H
