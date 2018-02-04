#include <assert.h>
#include "vectors.h"
#include "complex.h"
#include <stdio.h>

void print_csphvec(csphvec_t v)
{
	printf("(%f+%fj)r̂ + (%f+%fj)θ̂ + (%f+%fj)φ̂",
			creal(v.rc), cimag(v.rc),
			creal(v.thetac), cimag(v.thetac),
			creal(v.phic), cimag(v.phic)
	      );
}


void print_cart3(cart3_t v)
{
	printf("%fx̂ + %fŷ + %fẑ", v.x, v.y, v.z);
}

void print_ccart3(ccart3_t v)
{
	printf("(%f+%fj)x̂ + (%f+%fj)ŷ + (%f+%fj)ẑ",
			creal(v.x), cimag(v.x),
			creal(v.y), cimag(v.y),
			creal(v.z), cimag(v.z)
	      );
}

void print_sph(sph_t r)
{
	printf("(r=%g, θ=%g, φ=%g)", r.r, r.theta, r.phi);
}
