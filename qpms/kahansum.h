#ifndef KAHANSUM_H
#define KAHANSUM_H

static inline void kahaninit(double *sum, double *compensation) {
	sum = 0;
	compensation = 0;
}

static inline void kahaninc(double *sum, double *compensation, double input) {
	double compensated_input = input - *compensation;
	double nsum = *sum + compensated_input;
	*compensation = (nsum - *sum) - compensated_input;
	*sum = nsum;
}


static inline void ckahaninit(complex double *sum, complex double *compensation) {
	sum = 0;
	compensation = 0;
}

static inline void ckahaninc(complex double *sum, complex double *compensation, complex double input) {
	complex double compensated_input = input - *compensation;
	complex double nsum = *sum + compensated_input;
	*compensation = (nsum - *sum) - compensated_input;
	*sum = nsum;
}

#endif //KAHANSUM_H
