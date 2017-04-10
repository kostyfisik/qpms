#include <math.h>
#include "gaunt.h"

static const double sqrtpi = 1.7724538509055160272981674833411451827975494561223871;
//static const double ln2 = 0.693147180559945309417232121458176568075500134360255254120;

// Associated Legendre polynomial at zero argument (DLMF 14.5.1)
double legendre0(int m, int n) {
	return pow(2,m) * sqrtpi / gamma(.5*n - .5*m + .5) / gamma(.5*n-.5*m);
}

// Derivative of associated Legendre polynomial at zero argument (DLMF 14.5.2)
double legendreD0(int m, int n) {
	return -2 * legendre0(m, n);
}


