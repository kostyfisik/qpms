#include "bessels.h"
#include <stdlib.h>
#include <math.h>

static const double ln2 = 0.69314718055994531;

#include <stdio.h>

// general; gives an array of size 
complex double * hankelcoefftable_init(size_t maxn) {
	complex double *hct = malloc((maxn+1)*(maxn+2)/2 * sizeof(complex double));
	for(size_t n = 0; n <= maxn; ++n) {
		complex double *hcs = hankelcoeffs_get(hct,n);
		for (size_t k = 0; k <= n; ++k) {
			double lcoeff = lgamma(n+k+1) - lgamma(n-k+1) - lgamma(k+1) - k*ln2;
			printf("%f, %.16f\n", lcoeff, exp(lcoeff));
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

