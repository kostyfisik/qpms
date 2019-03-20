#ifndef AMOS_H
#define AMOS_H
#include "amos_mangling.h"

#define INTEGER_t int
#define DOUBLE_PRECISION_t double

void amos_zbesj(const DOUBLE_PRECISION_t *zr,
	      const DOUBLE_PRECISION_t *zi,
	      const DOUBLE_PRECISION_t *fnu,
	      const INTEGER_t *kode,
	      const INTEGER_t *n,
	      DOUBLE_PRECISION_t *cyr,
	      DOUBLE_PRECISION_t *cyi,
	      INTEGER_t *nz,
	      INTEGER_t *ierr);

void amos_zbesy(const DOUBLE_PRECISION_t *zr,
	      const DOUBLE_PRECISION_t *zi,
	      const DOUBLE_PRECISION_t *fnu,
	      const INTEGER_t *kode,
	      const INTEGER_t *n,
	      DOUBLE_PRECISION_t *cyr,
	      DOUBLE_PRECISION_t *cyi,
	      INTEGER_t *nz,
	      DOUBLE_PRECISION_t *cwrkr,
	      DOUBLE_PRECISION_t *cwrki,
	      INTEGER_t *ierr);

void amos_zbesh(const DOUBLE_PRECISION_t *zr,
	      const DOUBLE_PRECISION_t *zi,
	      const DOUBLE_PRECISION_t *fnu,
	      const INTEGER_t *kode,
	      const INTEGER_t *m,
	      const INTEGER_t *n,
	      DOUBLE_PRECISION_t *cyr,
	      DOUBLE_PRECISION_t *cyi,
	      INTEGER_t *nz,
	      INTEGER_t *ierr);


#endif
