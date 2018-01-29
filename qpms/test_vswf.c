#include "vswf.h"
#include "indexing.h"
#include <stdio.h>
#include <gsl/gsl_math.h>
const double dtheta = 0.01 * M_PI;
const qpms_l_t lMax = 3;

int main() {
	qpms_y_t nelem = qpms_lMax2nelem(lMax);
	for (double theta = 0.; theta <= M_PI; theta += dtheta) {
		printf("%.5e %.5e ", theta, cos(theta));
		for(qpms_normalisation_t norm = 1; norm <= 3; ++norm) {//fujka :D
			qpms_pitau_t pt = qpms_pitau_get(theta, lMax, norm);
			for (qpms_y_t y = 0; y < nelem; ++y)
				printf("%.5e %.5e %.5e ", pt.leg[y], pt.pi[y], pt.tau[y]);
			qpms_pitau_free(pt);
		}
		printf("\n");
	}
}





