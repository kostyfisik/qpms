#ifndef VECTORS_H
#define VECTORS_H
#include <math.h>
#ifndef M_PI_2
#define M_PI_2 (1.570796326794896619231321691639751442098584699687552910487)
#endif
#include "qpms_types.h"

//static inline double vectors_h_sq(double x) {return x*x;}

static const cart2_t CART2_ZERO = {0, 0};
static const cart3_t CART3_ZERO = {0, 0, 0};

static inline cart2_t cart2_add(const cart2_t a, const cart2_t b) {
	cart2_t res = {a.x+b.x, a.y+b.y};
	return res;
}

static inline cart2_t cart2_substract(const cart2_t a, const cart2_t b) {
	cart2_t res = {a.x-b.x, a.y-b.y};
	return res;
}

static inline cart2_t cart2_scale(const double c, const cart2_t v) {
	cart2_t res = {c * v.x, c * v.y};
	return res;
}

static inline double cart2_dot(const cart2_t a, const cart2_t b) {
	return a.x * b.x + a.y * b.y;
}

static inline double cart2_normsq(const cart2_t a) {
	return cart2_dot(a, a);
}

static inline double cart2norm(const cart2_t v) {
	return hypot(v.x, v.y); //sqrt(v.x*v.x + v.y*v.y);
}

static inline pol_t cart2pol(const cart2_t cart) {
	pol_t pol;
	pol.r = cart2norm(cart);
	pol.phi = atan2(cart.y, cart.x);
	return pol;
}

static inline sph_t pol2sph_equator(const pol_t pol) {
	sph_t sph;
	sph.r = pol.r;
	sph.phi = pol.phi;
	sph.theta = M_PI_2;
	return sph;
}

static inline sph_t cart22sph(const cart2_t cart) {
	sph_t sph;
	sph.r = cart2norm(cart);
	sph.theta = M_PI_2;
	sph.phi = atan2(cart.y, cart.x);
	return sph;
}

static inline cart3_t cart22cart3xy(const cart2_t a) {
	cart3_t c;
	c.x = a.x;
	c.y = a.y;
	c.z = 0;
	return c;
}

static inline cart2_t cart3xy2cart2(const cart3_t a) {
	cart2_t c = {a.x, a.y};
	return c;
}

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
	double sin_th = 
#ifdef QPMS_VECTORS_NICE_TRANSFORMATIONS
	   (sph.theta == M_PI) ? 0 :
#endif
		sin(sph.theta);
	cart.x = sph.r * sin_th * cos(sph.phi);
	cart.y = sph.r * sin_th * sin(sph.phi);
	cart.z = sph.r * cos(sph.theta);
	return cart;
}

static inline cart2_t pol2cart(const pol_t pol) {
	cart2_t cart;
	cart.x = pol.r * cos(pol.phi);
	cart.y = pol.r * sin(pol.phi);
	return cart;
}

static inline cart3_t cart3_add(const cart3_t a, const cart3_t b) {
	cart3_t res = {a.x+b.x, a.y+b.y, a.z+b.z};
	return res;
}

static inline cart3_t cart3_substract(const cart3_t a, const cart3_t b) {
	cart3_t res = {a.x-b.x, a.y-b.y, a.z-b.z};
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

static inline ccart3_t ccart3_add(const ccart3_t a, const ccart3_t b) {
	ccart3_t res = {a.x+b.x, a.y+b.y, a.z+b.z};
	return res;
}

static inline ccart3_t ccart3_substract(const ccart3_t a, const ccart3_t b) {
	ccart3_t res = {a.x-b.x, a.y-b.y, a.z-b.z};
	return res;
}

static inline csphvec_t csphvec_add(const csphvec_t a, const csphvec_t b) {
	csphvec_t res = {a.rc + b.rc, a.thetac + b.thetac, a.phic + b.phic};
	return res;
}

static inline csphvec_t csphvec_substract(const csphvec_t a, const csphvec_t b) {
	csphvec_t res = {a.rc - b.rc, a.thetac - b.thetac, a.phic - b.phic};
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

void print_csphvec(csphvec_t);
void print_ccart3(ccart3_t);
void print_cart3(cart3_t);
void print_sph(sph_t);

// kahan sums for various types... TODO make generic code using macros

static inline void ccart3_kahaninit(ccart3_t *sum, ccart3_t *compensation) {
	sum->x = sum->y = sum->z = compensation->x = compensation->y = compensation->z = 0;
}
static inline void csphvec_kahaninit(csphvec_t *sum, csphvec_t *compensation) {
	sum->rc = sum->thetac = sum->phic = compensation->rc = compensation->thetac = compensation->phic = 0;
}

static inline void ccart3_kahanadd(ccart3_t *sum, ccart3_t *compensation, const ccart3_t input) {
	ccart3_t comped_input = ccart3_substract(input, *compensation);
	ccart3_t nsum = ccart3_add(*sum, comped_input);
	*compensation = ccart3_substract(ccart3_substract(nsum, *sum), comped_input);
	*sum = nsum;
}

static inline void csphvec_kahanadd(csphvec_t *sum, csphvec_t *compensation, const csphvec_t input) {
	csphvec_t comped_input = csphvec_substract(input, *compensation);
	csphvec_t nsum = csphvec_add(*sum, comped_input);
	*compensation = csphvec_substract(csphvec_substract(nsum, *sum), comped_input);
	*sum = nsum;
}

static inline double csphvec_norm(const csphvec_t a) {
	return  sqrt(creal(a.rc * conj(a.rc) + a.thetac * conj(a.thetac) + a.phic * conj(a.phic)));
}

static inline double csphvec_reldiff_abstol(const csphvec_t a, const csphvec_t b, double tolerance) {
	double anorm = csphvec_norm(a);
	double bnorm = csphvec_norm(b);
	if (anorm <= tolerance && bnorm <= tolerance) return 0;
	return csphvec_norm(csphvec_substract(a,b)) / (anorm + bnorm);
}

static inline double csphvec_reldiff(const csphvec_t a, const csphvec_t b) {
	return csphvec_reldiff_abstol(a, b, 0);
}

typedef double matrix3d[3][3];
typedef double matrix2d[2][2];
typedef complex double cmatrix3d[3][3];
typedef complex double cmatrix2d[2][2];
#endif //VECTORS_H
