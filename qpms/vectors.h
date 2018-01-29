#ifndef VECTORS_H
#define VECTORS_H
#include <math.h>
#define M_PI_2 (1.570796326794896619231321691639751442098584699687552910487)
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

// equivalent to sph_loccart2cart in qpms_p.py
static inline ccart3_t csphvec2cart(const csphvec_t sphvec, const sph_t at) {
	ccart3_t res = {0, 0, 0};
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
	res.x = rx * sphvec.rc + tx * sphvec.thetac + fx * sphvec.phic;
	res.y = ry * sphvec.rc + ty * sphvec.thetac + fy * sphvec.phic;
	res.z = rz * sphvec.rc + tz * sphvec.thetac + fz * sphvec.phic;
	return sphvec;
}



#endif //VECTORS_H
