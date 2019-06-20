#include <math.h>
#include <gsl/gsl_math.h>
#include "assert_cython_workaround.h"
#include "vswf.h"
#include "indexing.h"
#include "translations.h" // TODO move qpms_sph_bessel_fill elsewhere
#include "qpms_specfunc.h"
#include <stdlib.h>
#include <string.h>


qpms_vswf_set_spec_t *qpms_vswf_set_spec_init() {
  qpms_vswf_set_spec_t *s = calloc(1,sizeof(qpms_vswf_set_spec_t));
  if (s == NULL) return NULL; // TODO warn
  // The rest are zeroes automatically because of calloc:
  s->lMax_L = -1;
  return s;
}

#define MAX(x,y) (((x)<(y))?(y):(x))

qpms_errno_t qpms_vswf_set_spec_append(qpms_vswf_set_spec_t *s, const qpms_uvswfi_t u) {
  qpms_l_t l;
  qpms_m_t m;
  qpms_vswf_type_t t;
  if (qpms_uvswfi2tmn(u, &t, &m, &l)!=QPMS_SUCCESS) return QPMS_ERROR; // TODO WARN
  if (s->n + 1 > s->capacity) {
    size_t newcap = (s->capacity < 32) ? 32 : 2*s->capacity;
    qpms_uvswfi_t *newmem = realloc(s->ilist, newcap * sizeof(qpms_uvswfi_t));
    if (newmem == NULL) return QPMS_ENOMEM; // TODO WARN
    s->capacity = newcap;
    s->ilist = newmem;
  }
  s->ilist[s->n] = u;
  ++s->n;
  switch(t) {
    case QPMS_VSWF_ELECTRIC:
      s->lMax_N = MAX(s->lMax_N, l);
      break;
    case QPMS_VSWF_MAGNETIC:
      s->lMax_M = MAX(s->lMax_M, l);
      break;
    case QPMS_VSWF_LONGITUDINAL:
      s->lMax_L = MAX(s->lMax_L, l);
      break;
    default:
      abort();
  }
  s->lMax = MAX(s->lMax, l);
  return QPMS_SUCCESS;
}

bool qpms_vswf_set_spec_isidentical(const qpms_vswf_set_spec_t *a,
                    const qpms_vswf_set_spec_t *b) {
  if (a == b) return true;
  if (a->n != b->n) return false;
  for (size_t i = 0; i < a->n; ++i)
    if (a->ilist[i] != b->ilist[i])
      return false;
  return true;
}

qpms_vswf_set_spec_t *qpms_vswf_set_spec_copy(const qpms_vswf_set_spec_t *or){
  qpms_vswf_set_spec_t *c = malloc(sizeof(qpms_vswf_set_spec_t));
  if (!c) abort(); // return NULL
  *c = *or;
  c->ilist = malloc(sizeof(qpms_uvswfi_t) * c->n);
  memcpy(c->ilist, or->ilist, sizeof(qpms_uvswfi_t)*c->n);
  c->capacity = c->n;
  return c;
}

qpms_vswf_set_spec_t *qpms_vswf_set_spec_from_lMax(qpms_l_t lMax, 
    qpms_normalisation_t norm) {
  qpms_vswf_set_spec_t *c = malloc(sizeof(qpms_vswf_set_spec_t));
  if (!c) abort(); // return NULL
  c->n = c->capacity = 2 * qpms_lMax2nelem(lMax);
  c->ilist = malloc(sizeof(qpms_uvswfi_t) * c->capacity);
  size_t i = 0;
  for (int it = 0; it < 2; ++it)
    for (qpms_l_t n = 1; n <= lMax; ++n)
      for (qpms_m_t m = -n; m <= n; ++m) 
        c->ilist[i++] = 
          qpms_tmn2uvswfi(it ? QPMS_VSWF_MAGNETIC : QPMS_VSWF_ELECTRIC, m, n);
  c->norm = norm;
  c->lMax = c->lMax_M = c->lMax_N = lMax;
  c->lMax_L = -1;
  return c;
}

void qpms_vswf_set_spec_free(qpms_vswf_set_spec_t *s) {
  if(s) free(s->ilist);
  free(s);
}

csphvec_t qpms_vswf_single_el(qpms_m_t m, qpms_l_t l, sph_t kdlj,
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  lmcheck(l,m);
  csphvec_t N;
  complex double *bessel = malloc((l+1)*sizeof(complex double));
  if(qpms_sph_bessel_fill(btyp, l, kdlj.r, bessel)) abort();
  qpms_pitau_t pt = qpms_pitau_get(kdlj.theta, l, norm);
  complex double eimf = cexp(m * kdlj.phi * I);
  qpms_y_t y = qpms_mn2y(m,l);

  N.rc = l*(l+1) * pt.leg[y] * bessel[l] / kdlj.r * eimf;
  complex double besselfac = bessel[l-1] - l * bessel[l] / kdlj.r;
  N.thetac = pt.tau[y] * besselfac * eimf;
  N.phic = pt.pi[y] * besselfac * I * eimf;

  qpms_pitau_free(pt);
  free(bessel);
  return N;
}
csphvec_t qpms_vswf_single_mg(qpms_m_t m, qpms_l_t l, sph_t kdlj,
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  lmcheck(l,m);
  csphvec_t M;
  complex double *bessel = malloc((l+1)*sizeof(complex double));
  if(qpms_sph_bessel_fill(btyp, l, kdlj.r, bessel)) abort();
  qpms_pitau_t pt = qpms_pitau_get(kdlj.theta, l, norm);
  complex double eimf = cexp(m * kdlj.phi * I);
  qpms_y_t y = qpms_mn2y(m,l);

  M.rc = 0.;
  M.thetac = pt.pi[y] * bessel[l] * I * eimf;
  M.phic = -pt.tau[y] * bessel[l] * eimf;

  qpms_pitau_free(pt);
  free(bessel);
  return M;
}

qpms_vswfset_sph_t *qpms_vswfset_make(qpms_l_t lMax, sph_t kdlj, 
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  qpms_vswfset_sph_t *res = malloc(sizeof(qpms_vswfset_sph_t));
  res->lMax = lMax;
  qpms_y_t nelem = qpms_lMax2nelem(lMax);
  res->el = malloc(sizeof(csphvec_t)*nelem);
  res->mg = malloc(sizeof(csphvec_t)*nelem);
  if(QPMS_SUCCESS != qpms_vswf_fill(NULL, res->mg, res->el, lMax, kdlj, btyp, norm))
    abort(); // or return NULL? or rather assert?
  return res;
}

void qpms_vswfset_sph_pfree(qpms_vswfset_sph_t *w) {
  assert(NULL != w && NULL != w->el && NULL != w->mg);
  free(w->el);
  free(w->mg);
  free(w);
}

qpms_errno_t qpms_vswf_fill(csphvec_t *const longtarget, csphvec_t * const mgtarget, csphvec_t * const eltarget,
    qpms_l_t lMax, sph_t kr,
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  assert(lMax >= 1);
  complex double *bessel = malloc((lMax+1)*sizeof(complex double));
  if(qpms_sph_bessel_fill(btyp, lMax, kr.r, bessel)) abort();
  qpms_pitau_t pt = qpms_pitau_get(kr.theta, lMax, norm);
  complex double const *pbes = bessel + 1; // starting from l = 1
  double const *pleg = pt.leg;
  double const *ppi = pt.pi;
  double const *ptau = pt.tau;
  csphvec_t *plong = longtarget, *pmg = mgtarget, *pel = eltarget;
  for(qpms_l_t l = 1; l <= lMax; ++l) {
    complex double besfac;
    complex double besderfac;
    if (kr.r) {
      besfac = *pbes / kr.r;
    } else {
      besfac = (1 == l) ? 1/3. : 0;
    }
    besderfac = *(pbes-1) - l * besfac;
    for(qpms_m_t m = -l; m <= l; ++m) {
      complex double eimf = cexp(m * kr.phi * I);
      if (longtarget) {
        complex double longfac = sqrt(l*(l+1)) * eimf;
        plong->rc = // FATAL FIXME: I get wrong result here for plane wave re-expansion 
          // whenever kr.r > 0 (for waves with longitudinal component, ofcoz)
          /*(*(pbes-1) - (l+1)/kr.r* *pbes)*/
          (besderfac-besfac) 
          * (*pleg) * longfac;
        plong->thetac = *ptau * besfac * longfac;
        plong->phic = *ppi * I * besfac * longfac;
        ++plong;
      }
      if (eltarget) {
        pel->rc = l*(l+1) * (*pleg) * besfac * eimf;
        pel->thetac = *ptau * besderfac * eimf;
        pel->phic = *ppi * besderfac * I * eimf;
        ++pel;
      }
      if (mgtarget) {
        pmg->rc = 0.;
        pmg->thetac = *ppi * (*pbes) * I * eimf;
        pmg->phic = - *ptau * (*pbes) * eimf;
        ++pmg;
      }
      ++pleg; ++ppi; ++ptau;
    }
    ++pbes;
  }
  free(bessel);
  qpms_pitau_free(pt);
  return QPMS_SUCCESS;
}

// consistency check: this should give the same results as the above function (up to rounding errors)
qpms_errno_t qpms_vswf_fill_alternative(csphvec_t *const longtarget, csphvec_t * const mgtarget, csphvec_t * const eltarget,
    qpms_l_t lMax, sph_t kr,
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  assert(lMax >= 1);
  complex double *bessel = malloc((lMax+1)*sizeof(complex double));
  if(qpms_sph_bessel_fill(btyp, lMax, kr.r, bessel)) abort();
  complex double const *pbes = bessel + 1; // starting from l = 1

  qpms_y_t nelem = qpms_lMax2nelem(lMax);
  csphvec_t * const a1 = malloc(3*nelem*sizeof(csphvec_t)), * const a2 = a1 + nelem, * const a3 = a2 + nelem;
  if(qpms_vecspharm_fill(a1, a2, a3, lMax, kr, norm)) abort();
  const csphvec_t *p1 = a1; 
  const csphvec_t *p2 = a2;
  const csphvec_t *p3 = a3;

  csphvec_t *plong = longtarget, *pmg = mgtarget, *pel = eltarget;
  for(qpms_l_t l = 1; l <= lMax; ++l) {
    complex double besfac = *pbes / kr.r;
    complex double besderfac = *(pbes-1) - l * besfac;
    double sqrtlfac = sqrt(l*(l+1));
    for(qpms_m_t m = -l; m <= l; ++m) {
      complex double eimf = cexp(m * kr.phi * I); // FIXME unused variable?!!!
      if (longtarget) {
        *plong = csphvec_add(csphvec_scale(besderfac-besfac, *p3),
            csphvec_scale(sqrtlfac * besfac, *p2));
        ++plong;
      }
      if (eltarget) {
        *pel = csphvec_add(csphvec_scale(besderfac, *p2),
            csphvec_scale(sqrtlfac * besfac, *p3));
        ++pel;
      }
      if (mgtarget) {
        *pmg = csphvec_scale(*pbes, *p1);
        ++pmg;
      }
      ++p1; ++p2; ++p3;
    }
    ++pbes;
  }
  free(a1);
  free(bessel);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_vecspharm_fill(csphvec_t *const a1target, csphvec_t *const a2target, csphvec_t *const a3target,
    qpms_l_t lMax, sph_t dir, qpms_normalisation_t norm) {
  assert(lMax >= 1);
  qpms_pitau_t pt = qpms_pitau_get(dir.theta, lMax, norm);
  double const *pleg = pt.leg;
  double const *ppi = pt.pi;
  double const *ptau = pt.tau;
  csphvec_t *p1 = a1target, *p2 = a2target, *p3 = a3target;
  for (qpms_l_t l = 1; l <= lMax; ++l) {
    for(qpms_m_t m = -l; m <= l; ++m) {
      complex double eimf = cexp(m * dir.phi * I);
      if (a1target) {
        p1->rc = 0;
        p1->thetac = *ppi * I * eimf;
        p1->phic = -*ptau * eimf;
        ++p1;
      }
      if (a2target) {
        p2->rc = 0;
        p2->thetac = *ptau * eimf;
        p2->phic = *ppi * I * eimf;
        ++p2;
      }
      if (a3target) {
        p3->rc = sqrt(l*(l+1)) * (*pleg) * eimf;
        p3->thetac = 0;
        p3->phic = 0;
        ++p3;
      }
    }
    ++pleg; ++ppi; ++ptau;
  }
  qpms_pitau_free(pt);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_vecspharm_dual_fill(csphvec_t *const a1target, csphvec_t *const a2target, csphvec_t *const a3target,
    qpms_l_t lMax, sph_t dir, qpms_normalisation_t norm) {
  assert(lMax >= 1);
  qpms_pitau_t pt = qpms_pitau_get(dir.theta, lMax, norm);
  double const *pleg = pt.leg;
  double const *ppi = pt.pi;
  double const *ptau = pt.tau;
  csphvec_t *p1 = a1target, *p2 = a2target, *p3 = a3target;
  for(qpms_l_t l = 1; l <= lMax; ++l) {
    for(qpms_m_t m = -l; m <= l; ++m) {
      double normfac = 1./qpms_normalisation_t_factor_abssquare(norm, l, m); // factor w.r.t. Kristensson
      complex double eimf = cexp(m * dir.phi * I);
      if (a1target) {
        p1->rc = 0;
        p1->thetac = conj(*ppi * normfac * I * eimf);
        p1->phic = conj(-*ptau * normfac * eimf);
        ++p1;
      }
      if (a2target) {
        p2->rc = 0;
        p2->thetac = conj(*ptau * normfac * eimf);
        p2->phic = conj(*ppi * normfac * I * eimf);
        ++p2;
      }
      if (a3target) {
        p3->rc = conj(sqrt(l*(l+1)) * (*pleg) * normfac * eimf);
        p3->thetac = 0;
        p3->phic = 0;
        ++p3;
      }
      ++pleg; ++ppi; ++ptau;
    }
  }
  qpms_pitau_free(pt);
  return QPMS_SUCCESS;
}


static inline complex double ipowl(qpms_l_t l) {
  switch(l % 4) {
    case 0: return 1;
            break;
    case 1: return I;
            break;
    case 2: return -1;
            break;
    case 3: return -I;
            break;
    default: abort();
  }
  assert(0);
}

qpms_errno_t qpms_planewave2vswf_fill_sph(sph_t wavedir, csphvec_t amplitude,
    complex double *target_longcoeff, complex double *target_mgcoeff,
    complex double *target_elcoeff, qpms_l_t lMax, qpms_normalisation_t norm) {
  qpms_y_t nelem = qpms_lMax2nelem(lMax);
  csphvec_t * const dual_A1 = malloc(3*nelem*sizeof(csphvec_t)), *const dual_A2 = dual_A1 + nelem,
            * const dual_A3 = dual_A2 + nelem;
  if (QPMS_SUCCESS != qpms_vecspharm_dual_fill(dual_A1, dual_A2, dual_A3, lMax, wavedir, norm))
    abort();
  const csphvec_t *pA1 = dual_A1, *pA2 = dual_A2, *pA3 = dual_A3;
  complex double *plong = target_longcoeff, *pmg = target_mgcoeff, *pel = target_elcoeff;
  for (qpms_l_t l = 1; l <= lMax; ++l) {
    complex double prefac1 = 4 * M_PI * ipowl(l);
    complex double prefac23 = - 4 * M_PI * ipowl(l+1);
    for (qpms_m_t m = -l; m <= l; ++m) {
      *plong = prefac23 * csphvec_dotnc(*pA3, amplitude);
      *pmg = prefac1 * csphvec_dotnc(*pA1, amplitude);
      *pel = prefac23 * csphvec_dotnc(*pA2, amplitude);
      ++pA1; ++pA2; ++pA3; ++plong; ++pmg; ++pel;
    }

  }
  free(dual_A1);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_planewave2vswf_fill_cart(cart3_t wavedir_cart /*allow complex k?*/, ccart3_t amplitude_cart,
    complex double * const longcoeff, complex double * const mgcoeff,
    complex double * const elcoeff, qpms_l_t lMax, qpms_normalisation_t norm) 
{

  sph_t wavedir_sph = cart2sph(wavedir_cart);
  csphvec_t amplitude_sphvec = ccart2csphvec(amplitude_cart, wavedir_sph);
  return qpms_planewave2vswf_fill_sph(wavedir_sph, amplitude_sphvec,
      longcoeff, mgcoeff, elcoeff, lMax, norm);
}

csphvec_t qpms_eval_vswf(sph_t kr,
    complex double * const lc, complex double *const mc, complex double *const ec,
    qpms_l_t lMax, qpms_bessel_t btyp, qpms_normalisation_t norm)
{
  qpms_y_t nelem = qpms_lMax2nelem(lMax);
  csphvec_t lsum, msum, esum, lcomp, mcomp, ecomp;
  csphvec_kahaninit(&lsum, &lcomp);
  csphvec_kahaninit(&msum, &mcomp);
  csphvec_kahaninit(&esum, &ecomp);
  csphvec_t *lset = NULL, *mset = NULL, *eset = NULL;
  if(lc) lset = malloc(nelem * sizeof(csphvec_t));
  if(mc) mset = malloc(nelem * sizeof(csphvec_t));
  if(ec) eset = malloc(nelem * sizeof(csphvec_t));
  qpms_vswf_fill(lset, mset, eset, lMax, kr, btyp, norm);
  if(lc) for(qpms_y_t y = 0; y < nelem; ++y)
    csphvec_kahanadd(&lsum, &lcomp, csphvec_scale(lc[y], lset[y]));
  if(mc) for(qpms_y_t y = 0; y < nelem; ++y)
    csphvec_kahanadd(&msum, &mcomp, csphvec_scale(mc[y], mset[y]));
  if(ec) for(qpms_y_t y = 0; y < nelem; ++y)
    csphvec_kahanadd(&esum, &ecomp, csphvec_scale(ec[y], eset[y]));
  if(lc) free(lset);
  if(mc) free(mset);
  if(ec) free(eset);
  //return csphvec_add(esum, csphvec_add(msum, lsum));
  csphvec_kahanadd(&esum, &ecomp, msum);
  csphvec_kahanadd(&esum, &ecomp, lsum);
  return esum;
}


#if 0
csphvec_t qpms_eval_uvswf(const qpms_vswf_set_spec_t *setspec,
    const complex double *coeffs, sph_t evalpoint,
    qpms_bessel_t btyp) {
  const qpms_l_t lMax = b->lMax;
  double *M, *N, *L;

}
#endif 
