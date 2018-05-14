#ifndef BESSELS_H
#define BESSELS_H
/* Short- and long-range parts of spherical Hankel functions
 * and (cylindrical) Hankel transforms of the long-range parts.
 * Currently, the implementation lies in bessels.c and 
 * lrhankel_recspace_dirty.c. The latter contains the implementation
 * of the Hankel transforms, but currenty only for a pretty limited
 * set of parameters. The general implementation is a BIG TODO here.
 */

#include <stddef.h>
#include <complex.h>

complex double *hankelcoefftable_init(size_t maxn);


// For navigating in the coefficients, maybe not for public use
static inline complex double *
trindex_cd(complex double const * const arr, size_t n){
  return (complex double *)(arr + n*(n+1)/2);
}

// general, gives the offset such that result[ql] is 
// the coefficient corresponding to the e**(I * x) * x**(-ql-1)
// term of the n-th Hankel function; no boundary checks!
static inline complex double *
hankelcoeffs_get(complex double const * const hankelcoefftable, size_t n){
  return trindex_cd(hankelcoefftable, n);
}

// general; target_longrange and target_shortrange are of size (maxn+1)
// if target_longrange is NULL, only the short-range part is calculated
void hankelparts_fill(complex double *target_longrange, complex double *target_shortrange,
		size_t maxn, size_t longrange_order_cutoff, /* terms e**(I x)/x**(k+1), 
							       k>= longrange_order_cutoff go 
							       completely to short-range part */
		complex double const * const hankelcoefftable,
		unsigned kappa, double vc, double x); // x = k0 * r



/* Hankel transforms of the long-range parts of the spherical Hankel functions */
// this declaration is general; however, the implementation 
// is so far only for kappa == 5, maxp == 5 TODO
void lrhankel_recpart_fill(complex  double *target_longrange_kspace /*Must be of size maxn*(maxn+1)/2*/,
	       size_t maxp /* Max. degree of transformed spherical Hankel function,
			      also the max. order of the Hankel transform */,
	       size_t longrange_order_cutoff /* terms e**(I x)/x**(k+1), k>= longrange_order_cutoff go
					        completely to the shortrange part
					        index with hankelcoeffs_get(target,p)l[delta_m] */, 
	       complex double const * const hankelcoefftable,
	       unsigned kappa, double c, double k0, double k);

#endif //BESSELS_H
