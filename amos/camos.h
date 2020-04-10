#ifndef CAMOS_H_
#define CAMOS_H_
#include "amos.h"

// TODO what about all the INTEGER_t and DOUBLE_PRECISION_t?

static inline int camos_zbesh(double zr, double zi, double fnu, int kode, int m,
          int n, double *cyr, double *cyi, int *nz) {
	int ierr;
	amos_zbesh(&zr, &zi, &fnu, &kode, &m, &n, cyr, cyi, nz, &ierr);
	return ierr;
}

static inline int camos_zbesj(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz) {
	int ierr;
	double cwrkr[n], cwrki[n];
	amos_zbesj(&zr, &zi, &fnu, &kode, &n, cyr, cyi, nz, &ierr);
	return ierr;
}

static inline int camos_zbesy(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz, double *cwrkr, double *cwrki) {
	int ierr;
	amos_zbesy(&zr, &zi, &fnu, &kode, &n, cyr, cyi, nz, cwrkr, cwrki, &ierr);
	return ierr;
}


#endif // CAMOS_H_
