#ifndef VECTORS_H
#define VECTORS_H
#include <math.h>
#define M_PI_2 (1.570796326794896619231321691639751442098584699687552910487)
#include "qpms_types.h"

//static inline double vectors_h_sq(double x) {return x*x;}

static inline double cart3norm(cart3_t v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

static inline sph_t cart2sph(cart3_t cart) {
	sph_t sph;
	sph.r = cart3norm(cart);
	sph.theta = sph.r ? acos(cart.z / sph.r) : M_PI_2;
	sph.phi = atan2(cart.y, cart.x);
	return sph;
}

static inline cart3_t sph2cart(sph_t sph) {
	cart3_t cart;
	double sin_th = sin(sph.theta);
	cart.x = sph.r * sin_th * cos(sph.phi);
	cart.y = sph.r * sin_th * sin(sph.phi);
	cart.z = sph.r * cos(sph.theta);
	return cart;
}


#endif //VECTORS_H
