#include "bessels.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static const double ln2 = 0.693147180559945309417;


// general; gives an array of size xxx with TODODESC
complex double * hankelcoefftable_init(size_t maxn) {
	complex double *hct = malloc((maxn+1)*(maxn+2)/2 * sizeof(complex double));
	for(size_t n = 0; n <= maxn; ++n) {
		complex double *hcs = hankelcoeffs_get(hct,n);
		for (size_t k = 0; k <= n; ++k) {
			double lcoeff = lgamma(n+k+1) - lgamma(n-k+1) - lgamma(k+1) - k*ln2;
			// for some reason, casting k-n to double does not work,so
			// cpow (I, k-n-1) cannot be used...
			complex double ifactor;
			switch ((n+1-k) % 4) {
				case 0:
					ifactor = 1;
					break;
				case 1:
					ifactor = -I;
					break;
				case 2:
					ifactor = -1;
					break;
				case 3:
					ifactor = I;
					break;
			}
			// the result should be integer, so round to remove inaccuracies
			hcs[k] = round(exp(lcoeff)) * ifactor;
		}
	}
	return hct;
}

void hankelparts_fill(complex double *lrt, complex double *srt, size_t maxn,
								size_t lrk_cutoff, complex double *hct,
								unsigned kappa, double c, double x) {
	if (lrt) memset(lrt, 0, (maxn+1)*sizeof(complex double));
	memset(srt, 0, (maxn+1)*sizeof(complex double));
	double regularisator = pow(1. - exp(-c * x), (double) kappa);
	double antiregularisator = 1. - regularisator;
	double xfrac = 1.; // x ** (-1-k)
	for (size_t k = 0; k <= maxn; ++k) {
	xfrac /= x;
	  for(size_t n = k; n <= maxn; ++n) 
	  	srt[n] += ((k<lrk_cutoff) ? antiregularisator : 1) 
		  				* xfrac * hankelcoeffs_get(hct,n)[k];
	  if (lrt && k < lrk_cutoff) for (size_t n = k; n <= maxn; ++n)
	  	lrt[n] += regularisator * xfrac * hankelcoeffs_get(hct,n)[k];
	}

	complex double expix = cexp(I * x);
	for(size_t n = 0; n <= maxn; ++n)
					srt[n] *= expix;
	if (lrt) for(size_t n = 0; n <= maxn; ++n)
					srt[n] *= expix;
}
