#include <math.h>
#include <gsl/gsl_math.h>
#include "assert_cython_workaround.h"
#include "vswf.h"
#include "indexing.h"
#include "translations.h" // TODO move qpms_sph_bessel_fill elsewhere
#include "qpms_specfunc.h"
#include <stdlib.h>
#include <string.h>
#include "qpms_error.h"
#include "normalisation.h"


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

struct bspec_reindex_pair {
  qpms_uvswfi_t ui;
  size_t i_orig;
};

static int cmp_bspec_reindex_pair(const void *aa, const void *bb) {
  const struct bspec_reindex_pair *a = aa, *b = bb;
  if (a->ui < b->ui) return -1;
  else if (a->ui == b->ui) return 0;
  else return 1;
}

size_t *qpms_vswf_set_reindex(const qpms_vswf_set_spec_t *small, const qpms_vswf_set_spec_t *big) {
  QPMS_UNTESTED;
  struct bspec_reindex_pair *small_pairs, *big_pairs;
  size_t *r;
  QPMS_CRASHING_MALLOC(small_pairs, sizeof(struct bspec_reindex_pair) * small->n);
  QPMS_CRASHING_MALLOC(big_pairs, sizeof(struct bspec_reindex_pair) * big->n);
  QPMS_CRASHING_MALLOC(r, sizeof(size_t) * small->n);
  for(size_t i = 0; i < small->n; ++i) {
    small_pairs[i].ui = small->ilist[i];
    small_pairs[i].i_orig = i;
  }
  for(size_t i = 0 ; i < big->n; ++i) {
    big_pairs[i].ui = big->ilist[i];
    big_pairs[i].i_orig = i;
  }
  qsort(small_pairs, small->n, sizeof(struct bspec_reindex_pair), cmp_bspec_reindex_pair);
  qsort(big_pairs, big->n, sizeof(struct bspec_reindex_pair), cmp_bspec_reindex_pair);

  size_t bi = 0;
  for(size_t si = 0; si < small->n; ++si) {
    while(big_pairs[bi].ui < small_pairs[si].ui) 
      ++bi;
    if(big_pairs[bi].ui == small_pairs[si].ui)
      r[small_pairs[si].i_orig] = big_pairs[si].i_orig;
    else
      r[small_pairs[si].i_orig] = ~(size_t)0;
  }

  free(small_pairs);
  free(big_pairs);
  return r;
}

csphvec_t qpms_vswf_single_el_csph(qpms_m_t m, qpms_l_t l, csph_t kdlj,
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  lmcheck(l,m);
  csphvec_t N;
  complex double *bessel;
  QPMS_CRASHING_MALLOC(bessel,(l+1)*sizeof(complex double));
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(btyp, l, kdlj.r, bessel));
  qpms_pitau_t pt = qpms_pitau_get(kdlj.theta, l, qpms_normalisation_t_csphase(norm));

  complex double eimf = qpms_spharm_azimuthal_part(norm, m, kdlj.phi);
  complex double d_eimf_dmf = qpms_spharm_azimuthal_part_derivative_div_m(norm, m, kdlj.phi);
  qpms_y_t y = qpms_mn2y(m,l);

  N.rc = l*(l+1) * pt.leg[y] * bessel[l] / kdlj.r * eimf;
  complex double besselfac = bessel[l-1] - l * bessel[l] / kdlj.r;
  N.thetac = pt.tau[y] * besselfac * eimf;
  N.phic = pt.pi[y] * besselfac * d_eimf_dmf;

  N = csphvec_scale(qpms_normalisation_factor_N_noCS(norm, l, m), N);

  qpms_pitau_free(pt);
  free(bessel);
  return N;
}

csphvec_t qpms_vswf_single_mg_csph(qpms_m_t m, qpms_l_t l, csph_t kdlj,
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  lmcheck(l,m);
  csphvec_t M;
  complex double *bessel;
  QPMS_CRASHING_MALLOC(bessel,(l+1)*sizeof(complex double));
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(btyp, l, kdlj.r, bessel));
  qpms_pitau_t pt = qpms_pitau_get(kdlj.theta, l, qpms_normalisation_t_csphase(norm));
  complex double eimf = qpms_spharm_azimuthal_part(norm, m, kdlj.phi);
  complex double d_eimf_dmf = qpms_spharm_azimuthal_part_derivative_div_m(norm, m, kdlj.phi);
  qpms_y_t y = qpms_mn2y(m,l);

  M.rc = 0.;
  M.thetac = pt.pi[y] * bessel[l] * d_eimf_dmf;
  M.phic = -pt.tau[y] * bessel[l] * eimf;

  M = csphvec_scale(qpms_normalisation_factor_M_noCS(norm, l, m), M);

  qpms_pitau_free(pt);
  free(bessel);
  return M;
}

csphvec_t qpms_vswf_single_el(qpms_m_t m, qpms_l_t l, sph_t kdlj,
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  return qpms_vswf_single_el_csph(m, l, sph2csph(kdlj), btyp, norm);
}

csphvec_t qpms_vswf_single_mg(qpms_m_t m, qpms_l_t l, sph_t kdlj,
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  return qpms_vswf_single_mg_csph(m, l, sph2csph(kdlj), btyp, norm);
}

qpms_vswfset_sph_t *qpms_vswfset_make(qpms_l_t lMax, sph_t kdlj, 
    qpms_bessel_t btyp, qpms_normalisation_t norm) {
  qpms_vswfset_sph_t *res = malloc(sizeof(qpms_vswfset_sph_t));
  res->lMax = lMax;
  qpms_y_t nelem = qpms_lMax2nelem(lMax);
  QPMS_CRASHING_MALLOC(res->el, sizeof(csphvec_t)*nelem);
  QPMS_CRASHING_MALLOC(res->mg, sizeof(csphvec_t)*nelem);
  QPMS_ENSURE_SUCCESS(qpms_vswf_fill(NULL, res->mg, res->el, lMax, kdlj, btyp, norm));
  return res;
}

void qpms_vswfset_sph_pfree(qpms_vswfset_sph_t *w) {
  assert(NULL != w && NULL != w->el && NULL != w->mg);
  free(w->el);
  free(w->mg);
  free(w);
}

csphvec_t qpms_vswf_L00(csph_t kr, qpms_bessel_t btyp,
    qpms_normalisation_t norm) { 
  QPMS_UNTESTED;
  // CHECKME Is it OK to ignore norm?? (Is L_0^0 the same for all conventions?)
  complex double bessel0;
  QPMS_ENSURE_SUCCESS(qpms_sph_bessel_fill(btyp, 0, kr.r, &bessel0));
  csphvec_t result = {0.25 * M_2_SQRTPI * bessel0, 0, 0};
  return result;
}

qpms_errno_t qpms_vswf_fill_csph(csphvec_t *const longtarget, 
    csphvec_t * const mgtarget, csphvec_t * const eltarget, qpms_l_t lMax,
    csph_t kr, qpms_bessel_t btyp, const qpms_normalisation_t norm) {
  assert(lMax >= 1);
  complex double *bessel = malloc((lMax+1)*sizeof(complex double));
  if(qpms_sph_bessel_fill(btyp, lMax, kr.r, bessel)) abort();
  qpms_pitau_t pt = qpms_pitau_get(kr.theta, lMax, qpms_normalisation_t_csphase(norm));
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
      complex double eimf = qpms_spharm_azimuthal_part(norm, m, kr.phi);
      complex double d_eimf_dmf = qpms_spharm_azimuthal_part_derivative_div_m(norm, m, kr.phi);
      if (longtarget) { QPMS_UNTESTED;
        double longfac = sqrt(l*(l+1));
        plong->rc = // FATAL FIXME: I get wrong result here for plane wave re-expansion 
          // whenever kr.r > 0 (for waves with longitudinal component, ofcoz)
          /*(*(pbes-1) - (l+1)/kr.r* *pbes)*/
          (besderfac-besfac) 
          * (*pleg) * longfac * eimf;
        plong->thetac = *ptau * besfac * longfac * eimf;
        plong->phic = *ppi * besfac * longfac * d_eimf_dmf;
        *plong = csphvec_scale(qpms_normalisation_factor_L_noCS(norm, l, m), *plong);
        ++plong;
      }
      if (eltarget) {
        pel->rc = l*(l+1) * (*pleg) * besfac * eimf;
        pel->thetac = *ptau * besderfac * eimf;
        pel->phic = *ppi * besderfac * d_eimf_dmf;
        *pel = csphvec_scale(qpms_normalisation_factor_N_noCS(norm, l, m), *pel);
        ++pel;
      }
      if (mgtarget) {
        pmg->rc = 0.;
        pmg->thetac = *ppi * (*pbes) * d_eimf_dmf;
        pmg->phic = - *ptau * (*pbes) * eimf;
        *pmg = csphvec_scale(qpms_normalisation_factor_M_noCS(norm, l, m), *pmg);
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

qpms_errno_t qpms_vswf_fill(csphvec_t *const longtarget, 
    csphvec_t * const mgtarget, csphvec_t * const eltarget, qpms_l_t lMax,
    sph_t kr, qpms_bessel_t btyp, qpms_normalisation_t norm) {
  csph_t krc = {kr.r, kr.theta, kr.phi};
  return qpms_vswf_fill_csph(longtarget, mgtarget, eltarget, lMax,
      krc, btyp, norm);
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
  csphvec_t *a;
  QPMS_CRASHING_MALLOC(a, 3*nelem*sizeof(csphvec_t))
  csphvec_t * const a1 = a, * const a2 = a1 + nelem, * const a3 = a2 + 2 * nelem;
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
      if (longtarget) {
        complex double L2Nfac = qpms_normalisation_factor_L_noCS(norm, l, m)
          / qpms_normalisation_factor_N_noCS(norm, l, m);
        *plong = csphvec_add(csphvec_scale(besderfac-besfac, *p3),
            csphvec_scale(sqrtlfac * besfac, *p2));
        *plong = csphvec_scale(L2Nfac, *plong);
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
  free(a);
  free(bessel);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_vecspharm_fill(csphvec_t *const a1target, csphvec_t *const a2target, csphvec_t *const a3target,
    qpms_l_t lMax, sph_t dir, qpms_normalisation_t norm) {
  assert(lMax >= 1);
  qpms_pitau_t pt = qpms_pitau_get(dir.theta, lMax, qpms_normalisation_t_csphase(norm));
  double const *pleg = pt.leg;
  double const *ppi = pt.pi;
  double const *ptau = pt.tau;
  csphvec_t *p1 = a1target, *p2 = a2target, *p3 = a3target;
  for (qpms_l_t l = 1; l <= lMax; ++l) {
    for(qpms_m_t m = -l; m <= l; ++m) {
      const complex double Mfac = qpms_normalisation_factor_M_noCS(norm, l, m);
      const complex double Nfac = qpms_normalisation_factor_N_noCS(norm, l, m);
      const complex double eimf = qpms_spharm_azimuthal_part(norm, m, dir.phi);
      const complex double deimf_dmf = qpms_spharm_azimuthal_part_derivative_div_m(norm, m, dir.phi);
      if (a1target) {
        p1->rc = 0;
        p1->thetac = *ppi * deimf_dmf * Mfac;
        p1->phic = -*ptau * eimf * Mfac;
        ++p1;
      }
      if (a2target) {
        p2->rc = 0;
        p2->thetac = *ptau * eimf * Nfac;
        p2->phic = *ppi * deimf_dmf * Nfac;
        ++p2;
      }
      if (a3target) {
        p3->rc = sqrt(l*(l+1)) * (*pleg) * eimf * Nfac;
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
#if 1
  return qpms_vecspharm_fill(a1target, a2target, a3target, lMax, dir, 
      qpms_normalisation_dual(norm));
#else
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
#endif
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
      if(target_longcoeff) *plong = prefac23 * csphvec_dotnc(*pA3, amplitude);
      if(target_mgcoeff) *pmg = prefac1 * csphvec_dotnc(*pA1, amplitude);
      if(target_elcoeff) *pel = prefac23 * csphvec_dotnc(*pA2, amplitude);
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

qpms_errno_t qpms_incfield_planewave(complex double *target, const qpms_vswf_set_spec_t *bspec,
    const cart3_t evalpoint, const void *args, bool add) {
  QPMS_UNTESTED;
  const qpms_incfield_planewave_params_t *p = args;

  const ccart3_t k_cart = p->use_cartesian ? p->k.cart : csph2ccart(p->k.sph);
  const complex double phase = ccart3_dotnc(k_cart, cart32ccart3(evalpoint));
  if(cimag(phase))
    QPMS_INCOMPLETE_IMPLEMENTATION("Complex-valued wave vector not implemented correctly; cf. docs.");
  const complex double phasefac = cexp(I*phase);

  // Throw away the imaginary component; TODO handle it correctly
  const sph_t k_sph = csph2sph(p->use_cartesian ? ccart2csph(p->k.cart) : p->k.sph);
  const csphvec_t E_sph = csphvec_scale(phasefac, 
      p->use_cartesian ? ccart2csphvec(p->E.cart, k_sph) : p->E.sph);

  complex double *lc = NULL, *mc = NULL, *nc = NULL;
  const qpms_y_t nelem = qpms_lMax2nelem(bspec->lMax);
  if (bspec->lMax_L > 0)  QPMS_CRASHING_MALLOC(lc, nelem * sizeof(complex double));
  if (bspec->lMax_M > 0)  QPMS_CRASHING_MALLOC(mc, nelem * sizeof(complex double));
  if (bspec->lMax_N > 0)  QPMS_CRASHING_MALLOC(nc, nelem * sizeof(complex double));

  qpms_errno_t retval = qpms_planewave2vswf_fill_sph(k_sph, E_sph, lc, mc, nc,
      bspec->lMax, bspec->norm);

  if (!add) memset(target, 0, bspec->n * sizeof(complex double));
  for (size_t i = 0; i < bspec->n; ++i) {
    const qpms_uvswfi_t ui = bspec->ilist[i];
    if (ui == QPMS_UI_L00) // for l = 0 the coefficient is zero due to symmetry (for real wave vector)
      target[i] = 0;
    else {
      qpms_vswf_type_t t; qpms_y_t y;
      QPMS_ENSURE_SUCCESS(qpms_uvswfi2ty_l(ui, &t, &y));
      switch(t) {
        case QPMS_VSWF_ELECTRIC:
          target[i] += nc[y];
          break;
        case QPMS_VSWF_MAGNETIC:
          target[i] += mc[y];
          break;
        case QPMS_VSWF_LONGITUDINAL:
          target[i] += lc[y];
          break;
        default:
          QPMS_WTF;
      }
    }
  }
  free(lc); free(mc); free(nc);
  return retval;
}
  

csphvec_t qpms_eval_vswf_csph(csph_t kr,
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
  qpms_vswf_fill_csph(lset, mset, eset, lMax, kr, btyp, norm);
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

csphvec_t qpms_eval_vswf(sph_t kr,
    complex double * const lc, complex double *const mc, complex double *const ec,
    qpms_l_t lMax, qpms_bessel_t btyp, qpms_normalisation_t norm) {
  csph_t krc = {kr.r, kr.theta, kr.phi};
  return qpms_eval_vswf_csph(krc, lc, mc, ec, lMax, btyp, norm);
}

qpms_errno_t qpms_uvswf_fill(csphvec_t *const target, const qpms_vswf_set_spec_t *bspec,
    csph_t kr, qpms_bessel_t btyp) {
  QPMS_UNTESTED;
  QPMS_ENSURE(target, "Target array pointer must not be NULL."); 
  csphvec_t *el = NULL, *mg = NULL, *lg = NULL;
  const qpms_y_t nelem = qpms_lMax2nelem(bspec->lMax);
  if (bspec->lMax_L > 0)
    QPMS_CRASHING_MALLOC(lg, nelem * sizeof(csphvec_t));
  if (bspec->lMax_M > 0)
    QPMS_CRASHING_MALLOC(mg, nelem * sizeof(csphvec_t));
  if (bspec->lMax_N > 0)
    QPMS_CRASHING_MALLOC(el, nelem * sizeof(csphvec_t));
  qpms_errno_t retval = qpms_vswf_fill_csph(lg, mg, el, bspec->lMax, kr, btyp, bspec->norm);
  for (size_t i = 0; i < bspec->n; ++i) {
    const qpms_uvswfi_t ui = bspec->ilist[i];
    if (ui == QPMS_UI_L00) // l = 0 longitudinal wave must be calculated separately.
      target[i] = qpms_vswf_L00(kr, btyp, bspec->norm);
    else {
      qpms_vswf_type_t t; qpms_y_t y;
      QPMS_ENSURE_SUCCESS(qpms_uvswfi2ty_l(ui, &t, &y));
      switch(t) {
        case QPMS_VSWF_ELECTRIC:
          target[i] = el[y];
          break;
        case QPMS_VSWF_MAGNETIC:
          target[i] = mg[y];
          break;
        case QPMS_VSWF_LONGITUDINAL:
          target[i] = lg[y];
          break;
        default:
          QPMS_WTF;
      }
    }
  }
  free(lg);
  free(mg);
  free(el);
  return retval;
}
  

csphvec_t qpms_eval_uvswf(const qpms_vswf_set_spec_t *bspec,
    const complex double *coeffs, const csph_t kr,
    const qpms_bessel_t btyp) {
  QPMS_UNTESTED;
  complex double *cM = NULL, *cN = NULL, *cL = NULL, cL00 = 0;
  if (bspec->lMax_L > 0)
    QPMS_CRASHING_CALLOC(cL, bspec->n, sizeof(complex double));
  if (bspec->lMax_M > 0)
    QPMS_CRASHING_CALLOC(cM, bspec->n, sizeof(complex double));
  if (bspec->lMax_N > 0)
    QPMS_CRASHING_CALLOC(cN, bspec->n, sizeof(complex double));
  for (size_t i = 0; i < bspec->n; ++i) {
    if (bspec->ilist[i] == 0) // L00, needs special care
      cL00 = coeffs[i];
    else {
      qpms_vswf_type_t t;
      qpms_y_t y;
      QPMS_ENSURE_SUCCESS(qpms_uvswfi2ty_l(bspec->ilist[i], &t, &y));
      switch(t) {
        case QPMS_VSWF_LONGITUDINAL:
          QPMS_ASSERT(cL);
          cL[y] = coeffs[i];
          break;
        case QPMS_VSWF_MAGNETIC:
          QPMS_ASSERT(cM);
          cM[y] = coeffs[i];
          break;
        case QPMS_VSWF_ELECTRIC:
          QPMS_ASSERT(cN);
          cN[y] = coeffs[i];
          break;
        default:
          QPMS_WTF;
      }
    }
  }
  csphvec_t result = qpms_eval_vswf_csph(kr, cL, cM, cN, bspec->lMax, btyp, bspec->norm);
  free(cM); free(cN); free(cL);
  if(cL00)
    result = csphvec_add(result,
       csphvec_scale(cL00, qpms_vswf_L00(kr, btyp, bspec->norm)));
  return result;
}

