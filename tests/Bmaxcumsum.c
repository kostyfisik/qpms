#include "../qpms/gaunt.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//const int lMax = 30;

#define MIN(x,y) (((x)>(y))?(y):(x))

static inline int Qmax(int m, int n, int mu, int nu) {
	return MIN(n, MIN(nu, (n+nu+1-abs(-m+mu))/2));
}

int main(int argc, char **argv) 
{
	int lMax;
	if (argc < 2)
		lMax = 30;
	else lMax = atoi(argv[1]);

	printf("// assuming lMax == %d:\n", lMax);
	size_t Bmaxcumsum = 0;
	printf("0x%zx,\n", Bmaxcumsum);
	for (int n = 0; n <= lMax; n++)
		for (int m = -n; m <= n; m++)
			for (int nu = 0; nu <= lMax; nu++)
				for (int mu = -nu; mu <= nu; mu++) {
					//FIXME toto je možná úplně blbě
					assert(gaunt_q_max(-m, n+1, mu, nu) == Qmax(m, n, mu, nu);
					Bmaxcumsum += gaunt_q_max(-m, n+1, mu, nu) + 1;
					printf("0x%zx,\n", Bmaxcumsum);
				}
	return 0;
}
