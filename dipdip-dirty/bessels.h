#ifndef BESSELS_H
#define BESSELS_H

#include <stddef.h>
#include <complex.h>

complex double *hankelcoefftable_init(size_t maxn);



static inline complex double *
trindex_cd(complex double *arr, size_t n){
  return arr + n*(n+1)/2;
}

// general, gives the offset such that result[ql] is 
// the coefficient corresponding to the e**(I * x) * x**(-ql-1)
// term of the n-th Hankel function; no boundary checks!
static inline complex double *
hankelcoeffs_get(complex double *hankelcoefftable, size_t n){
  return trindex_cd(hankelcoefftable, n);
}

// general; target_longrange and target_shortrange are of size (maxn+1)
// if target_longrange is NULL, only the short-range part is calculated
void hankelparts_fill(complex double *target_longrange, complex double *target_shortrange,
		size_t maxn, size_t longrange_order_cutoff, // x**(-(order+1)-1) terms go completely to short-range part
		complex double *hankelcoefftable,
		unsigned kappa, double vc, double x); // x = k0 * r


// this declaration is general; however, the implementation 
// is so far only for kappa == ???, maxn == ??? TODO
void lrhankel_recpart_fill(complex  double *target_longrange_kspace /*Must be of size maxn*(maxn+1)/2*/,
	       size_t maxp, size_t longrange_k_cutoff /* terms e**(I x)/x**(k+1), k>= longrange_k_cutoff go
							 completely to the shortrange part
							 index with hankelcoeffs_get(target,p)l[delta_m] */, 
	       complex double  *hankelcoefftable,
	       unsigned kappa, double c, double k0, double k);

#endif //BESSELS_H
