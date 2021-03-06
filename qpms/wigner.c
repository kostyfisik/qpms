#include "quaternions.h"
#include "tiny_inlines.h"
#include "kahansum.h"
#define WIGNER_ATOL (1e-15)

complex double qpms_wignerD_elem(const qpms_quat_t R,
    const qpms_l_t l, const qpms_m_t mp, const qpms_m_t m) {
  // TODO do some optimisations... The combinatoric coeffs could be precomputed.
  QPMS_ENSURE(abs(m) <= l || abs(mp) <= l, "Got invalid values l = %d, m = %d", (int)l, (int)m);
  const double mRa = cabs(R.a), mRb = cabs(R.b);
  if (mRa < WIGNER_ATOL) {
    if (m != -mp) return 0;
    else return min1pow(l+m) * cpow(R.b, 2*m);
  } else if (mRb < WIGNER_ATOL) {
    if (m != mp) return 0;
    else return cpow(R.a, 2*m);
  } else if (mRa >= mRb) {
    complex double prefac = exp(.5*(lgamma(l+m+1) + lgamma(l-m+1)
        - lgamma(l+mp+1) - lgamma(l-mp+1))) 
      * pow(mRa, 2*l-2*m)
      * cpow(R.a, mp+m)
      * cpow(R.b, -mp+m);
    double sum, sumc; kahaninit( &sum, &sumc);
    // FIXME this is probably quite stupid way to calculate the sum
    for(int rho=MAX(0, -m+mp); rho <= MIN(l-m, l+mp); ++rho)
      kahanadd(&sum, &sumc,
          exp( 
            lgamma(l+mp+1) - lgamma(rho+1) - lgamma(l+mp-rho+1)
            + lgamma(l-mp+1) - lgamma(l-rho-m+1) - lgamma(-mp+rho+m+1)
          )
          * pow(-mRb*mRb/(mRa*mRa), rho)
      );
    return prefac * sum;
  } else { // (mRa < mRb)
    complex double prefac = min1pow(l-m) * exp(.5*(lgamma(l+m+1) + lgamma(l-m+1)
        - lgamma(l+mp+1) - lgamma(l-mp+1))) 
      * pow(mRb, 2*l-2*m)
      * cpow(R.a, mp+m)
      * cpow(R.b, -mp+m);
    double sum, sumc; kahaninit( &sum, &sumc);
    // FIXME this is probably quite stupid way to calculate the sum
    for(int rho=MAX(0, -m-mp); rho <= MIN(l-mp, l-m); ++rho)
      kahanadd(&sum, &sumc,
          exp( 
            lgamma(l+mp+1) - lgamma(l-m-rho+1) - lgamma(mp+m+rho+1)
            + lgamma(l-mp+1) - lgamma(rho+1) - lgamma(l-mp-rho+1)
          )
          * pow(-mRa*mRa/(mRb*mRb), rho)
      );
    return prefac * sum;
  }
}

complex double qpms_vswf_irot_elem_from_irot3(const qpms_irot3_t q,
    const qpms_l_t l, const qpms_m_t mp, const qpms_m_t m, bool pseudo) {
#ifndef NDEBUG
  qpms_irot3_checkdet(q);
#endif
  complex double res = qpms_wignerD_elem(q.rot, l, mp, m);
  if (q.det == -1) {
    res *= min1pow(l);
    if (pseudo)
      res *= -1;
  }
  return res;
}
