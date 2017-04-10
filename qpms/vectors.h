#ifndef VECTORS_H
#define VECTORS_H
#include <math.h>

typedef struct {
	double x, y, z;
} cart3_t;

typedef struct {
	double x, y;
} cart2_t;

typedef struct {
	double r, theta, phi;
} sph_t;

typedef struct {
	double r, phi;
} pol_t;


//static inline double vectors_h_sq(double x) {return x*x;}

static inline cart3norm(cart3_t v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

static inline sph_t cart2sph(cart3_t cart) {
	sph_t sph;
	sph.r = cart3norm(cart);
	sph.theta = sph.r ? acos(cart.z / sph.r) : M_PI_2;
	sph.phi = atan2(cart.y, cart.x);
	return sph;
}

static inline cart_t sph2cart(sph_t sph) {
	cart_t cart;
	double sin_th = sin(sph.theta);
	cart.x = sph.r * sin_th * cos(sph.phi);
	cart.y = sph.r * sin_th * sin(sph.phi);
	cart.z = sph.r * cos(sph.theta);
	return cart;
}


#endif //VECTORS_H
