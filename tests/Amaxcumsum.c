#include "../qpms/gaunt.h"
#include <stdio.h>
#include <stdlib.h>

//const int lMax = 30;

int main(int argc, char **argv) 
{
	int lMax;
	if (argc < 2)
		lMax = 30;
	else lMax = atoi(argv[1]);

	printf("// assuming lMax == %d:\n", lMax);
	size_t Amaxcumsum = 0;
	printf("0x%zx,\n", Amaxcumsum);
	for (int n = 0; n <= lMax; n++)
		for (int m = -n; m <= n; m++)
			for (int nu = 0; nu <= lMax; nu++)
				for (int mu = -nu; mu <= nu; mu++) {
					Amaxcumsum += gaunt_q_max(-m, n, mu, nu) + 1;
					printf("0x%zx,\n", Amaxcumsum);
				}
	return 0;
}
