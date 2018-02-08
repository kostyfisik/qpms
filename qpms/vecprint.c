#include <assert.h>
#include "vectors.h"
#include "complex.h"
#include <stdio.h>

void print_csphvec(csphvec_t v)
{
	printf("(%g+%gj)r̂ + (%g+%gj)θ̂ + (%g+%gj)φ̂",
			creal(v.rc), cimag(v.rc),
			creal(v.thetac), cimag(v.thetac),
			creal(v.phic), cimag(v.phic)
	      );
}


void print_cart3(cart3_t v)
{
	printf("%gx̂ + %gŷ + %gẑ", v.x, v.y, v.z);
}

void print_ccart3(ccart3_t v)
{
	printf("(%g+%gj)x̂ + (%g+%gj)ŷ + (%g+%gj)ẑ",
			creal(v.x), cimag(v.x),
			creal(v.y), cimag(v.y),
			creal(v.z), cimag(v.z)
	      );
}

void print_sph(sph_t r)
{
	printf("(r=%g, θ=%g, φ=%g)", r.r, r.theta, r.phi);
}
