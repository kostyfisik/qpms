#include <cblas.h>
#include <string.h>
#include <math.h>
#include <qpms_error.h>

static inline size_t mu_index(size_t k, size_t j) {
  return k * (k - 1) / 2 + j;
}

/// Gram-Schmidt orthogonalisation.
/** Does not return the actual orthogonal basis (as it is not needed
 * for the LLL algorithm as such)  but rather only
 * the mu(i,j) coefficients and squared norms of the orthogonal vectors
 */
static void gram_schmidt(
    double *mu, ///< Array of \f[ \mu_{k,j} = \frac{\vect{v}_i \cdot\vect{v}_i^*}{|\vect v_j^*|^2}\f] of length mu_index(bsize, 0), 
    double *vstar_sqnorm, ///< Array of \f$ \abs{\vect v_i^*}^2 \f$ of length bsize.
    const double *v, ///< Vectors to orthogonalise, size [bsize][ndim], 
    const size_t bsize, ///< Size of the basis ( = dimensionality of the lattice)
    const size_t ndim ///< Dimensionality of the space.
    )
{
  double v_star[bsize][ndim];
  for (size_t i = 0; i < bsize; ++i) {
    memcpy(v_star[i], v+i*ndim, ndim*sizeof(double));
    double parallel_part[ndim /*???*/];
    memset(parallel_part, 0, sizeof(parallel_part));
    for (size_t j = 0; j < i; ++j) {
      double mu_numerator = cblas_ddot(ndim, v + i*ndim, 1, v_star[j], 1);
      mu[mu_index(i, j)] = mu_numerator / vstar_sqnorm[j];
      cblas_daxpy(ndim, mu[mu_index(i, j)], v_star[j], 1, parallel_part, 1);
    }
    cblas_daxpy(ndim, -1, parallel_part, 1, v_star[i], 1);
    vstar_sqnorm[i] = cblas_ddot(ndim, v_star[i], 1, v_star[i], 1);
  }
}

static inline double fsq(double x) { return x * x; };

// A naïve implementation of Lenstra-Lenstra-Lovász algorithm.
int qpms_reduce_lattice_basis(double *b, const size_t bsize, const size_t ndim, 
    double delta)
{
  QPMS_ENSURE(bsize <= ndim, 
      "The basis must have less elements (%zd) than the space dimension (%zd).",
      bsize, ndim);
  double mu[mu_index(bsize,0)];
  double bstar_sqnorm[bsize];
  gram_schmidt(mu, bstar_sqnorm, b, bsize, ndim);
  size_t k = 1;
  while (k < bsize) {
    // Step 1 of LLL, achieve mu(k, k-1) <= 0.5
    if (fabs(mu[mu_index(k, k-1)]) > 0.5) {
      double r = round(mu[mu_index(k, k-1)]);
      // "normalize" b(k), replacing it with b(k) - r b(k-1)
      cblas_daxpy(ndim, -r, b+(k-1)*ndim, 1, b+k*ndim, 1);
      // update mu to correspond to the new b(k)
      for(size_t j = 0; j < bsize; ++j)
        mu[mu_index(k, j)] -= r*mu[mu_index(k-1, j)];
    }
    // Step 2
    if (k > 0 && // Case 1
        bstar_sqnorm[k] < (delta - fsq(mu[mu_index(k, k-1)])) * bstar_sqnorm[k-1]) {
      // swap b(k) and b(k-1)
      cblas_dswap(ndim, &b[k*ndim], 1, &b[(k-1)*ndim], 1);
      double B_k = bstar_sqnorm[k];
      double mu_kkm1_old = mu[mu_index(k, k-1)];
      double C = B_k + fsq(mu_kkm1_old) * bstar_sqnorm[k-1];
      mu[mu_index(k, k-1)] *= bstar_sqnorm[k-1] / C;
      bstar_sqnorm[k] *= bstar_sqnorm[k-1] / C;
      bstar_sqnorm[k-1] = C;
      for(size_t j = k+1; j < bsize; ++j) {
        double m = mu[mu_index(j, k-1)];
        mu[mu_index(j, k-1)] = m*m + mu[mu_index(j, k)] * B_k / C;
        mu[mu_index(j, k)] = m - mu[mu_index(j, k)] * mu_kkm1_old;
      }
      for(size_t j = 0; j < k-1; ++j) {
        double m = mu[mu_index(k-1, j)];
        mu[mu_index(k-1, j)] = mu[mu_index(k, j)];
        mu[mu_index(k, j)] = m;
      }
      --k;
    } else { // Case 2
      size_t l = k;
      while(l > 0) {
        --l;
        if(fabs(mu[mu_index(k, l)] > 0.5)) {
          double r = round(mu[mu_index(k, l)]);
          cblas_daxpy(ndim, -r, b+l*ndim, 1, b+k*ndim, 1);
          for (size_t j = 0; j < l; ++j)
            mu[mu_index(k, j)] -= r * mu[mu_index(l, j)];
          mu[mu_index(k, l)] -= r;
          l = k;
        }
      }
      ++k;
    }
  }
  return 0;
}
