#include "bessels.h"
#include <stdio.h>
#include <stdlib.h>
#if 0
int main() {
	size_t maxn = 5;
	complex double *hct = hankelcoefftable_init(maxn);
	for (size_t n = 0; n <= maxn; ++n) {
		printf("n = %zd\n", n);
		for(size_t k = 0; k<=n; ++k)
			printf("%p: %f + %fj,\n", hankelcoeffs_get(hct,n) + k, creal(hankelcoeffs_get(hct,n)[k]),
					cimag(hankelcoeffs_get(hct,n)[k]));
		printf("\n");
	}
	printf("%f+%fj\n",creal(cpow(I,(ptrdiff_t)-1)), cimag(cpow(I,(ptrdiff_t)-1)));
	return 0;
}
#endif

//#if 0
int main() {
	size_t maxn = 6;
	size_t lrk_cutoff = 2;
	size_t kappa = 4;
	double c = 0.1324;

	double xmin = 0.;
	double xstep = 0.001;
	double xmax = 20;

	complex double *hct = hankelcoefftable_init(maxn);
	complex double srhankel[maxn+1];
	complex double lrhankel[maxn+1];

	for(double x = xmin; x <= xmax; x += xstep) {
		hankelparts_fill(lrhankel, srhankel, maxn, lrk_cutoff, hct, kappa, c, x);
		printf("%f ", x);
		for(size_t n = 0; n <= maxn; ++n)
			printf("%.16e %.16e %.16e %.16e ", creal(lrhankel[n]), cimag(lrhankel[n]), creal(srhankel[n]), cimag(srhankel[n]));
		printf("\n");
	}
	
	free(hct);
	return 0;
}
//#endif

