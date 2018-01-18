#ifndef BESSELS_H
#define BESSELS_H

#include <stddef.h>
#include <complex.h>

complex double *hankelcoefftable_init(size_t maxn);

// general, gives the offset such that result[k] is 
// the coefficient corresponding to the e**(I * x) * x**(-k-1)
// term of the Hankel function; no boundary checks!
static inline complex double *
hankelcoeffs_get(complex double *hankelcoefftable, size_t n){
        return hankelcoefftable + n*(n+1)/2;
}


#endif //BESSELS_H
