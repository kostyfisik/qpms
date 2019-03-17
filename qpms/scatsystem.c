#include <stdlib.h>
#include <cblas.h>
#include <lapacke.h>
#include "scatsystem.h"
#include "indexing.h"
#include "vswf.h"
#include "groups.h"
#include "symmetries.h"
#include <gsl/gsl_spline.h>
#include <assert.h>
#include <unistd.h>
#include "vectors.h"
#include "wigner.h"
#include <string.h>
#include "qpms_error.h"
#include "translations.h"
#include <pthread.h>

#define SQ(x) ((x)*(x))
#define QPMS_SCATSYS_LEN_RTOL 1e-13
#define QPMS_SCATSYS_TMATRIX_ATOL 1e-14
#define QPMS_SCATSYS_TMATRIX_RTOL 1e-12
qpms_tmatrix_t *qpms_tmatrix_init(const qpms_vswf_set_spec_t *bspec) {
  qpms_tmatrix_t *t = malloc(sizeof(qpms_tmatrix_t));
  if (!t) abort();
  else {
    t->spec = bspec;
    size_t n = bspec->n;
    t->m = calloc(n*n, sizeof(complex double));
    if (!t->m) abort();
    t->owns_m = true;
  }
  return t;
}

qpms_tmatrix_t *qpms_tmatrix_copy(const qpms_tmatrix_t *T) {
  qpms_tmatrix_t *t = qpms_tmatrix_init(T->spec);
  size_t n = T->spec->n;
  for(size_t i = 0; i < n*n; ++i)
    t->m = T->m;
  return t;
}

void qpms_tmatrix_free(qpms_tmatrix_t *t){
  if(t && t->owns_m) free(t->m);
  free(t);
}

qpms_tmatrix_t *qpms_tmatrix_apply_symop_inplace(
                qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  //qpms_tmatrix_t *t = qpms_tmatrix_init(T->spec);
  const size_t n = T->spec->n;
  complex double tmp[n][n];
  // tmp = M T
  const complex double one = 1, zero = 0;
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasNoTrans,
      n, n, n, &one, M, n, T->m, n, &zero, tmp, n);
  // t->m = tmp M* = M T M*
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasConjTrans,
      n, n, n, &one, tmp, n, M, n, &zero, T->m, n);
  return T;
}

qpms_tmatrix_t *qpms_tmatrix_apply_symop(
                const qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  qpms_tmatrix_t *t = qpms_tmatrix_init(T->spec);
  const size_t n = T->spec->n;
  complex double tmp[n][n];
  // tmp = M T
  const complex double one = 1, zero = 0;
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasNoTrans,
      n, n, n, &one, M, n, T->m, n, &zero, tmp, n);
  // t->m = tmp M* = M T M*
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasConjTrans,
      n, n, n, &one, tmp, n, M, n, &zero, t->m, n);
  return t;
}

qpms_errno_t qpms_symmetrise_tmdata_irot3arr(
    complex double *tmdata, const size_t tmcount,
    const qpms_vswf_set_spec_t *bspec,
    const size_t n_symops, const qpms_irot3_t *symops) {
  const size_t n = bspec->n;
  qpms_tmatrix_t *tmcopy = qpms_tmatrix_init(bspec);
  complex double *symop_matrices = malloc(n*n*sizeof(complex double) * n_symops);
  if(!symop_matrices) qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
      "malloc() failed.");
  for (size_t i = 0; i < n_symops; ++i) 
    qpms_irot3_uvswfi_dense(symop_matrices + i*n*n, bspec, symops[i]);
  complex double tmp[n][n];
  const complex double one = 1, zero = 0;
  for (size_t tmi = 0; tmi < tmcount; ++tmi) {
    // Move the data in tmcopy; we will then write the sum directly into tmdata.
    memcpy(tmcopy->m, tmdata+n*n*tmi, n*n*sizeof(complex double));
    memset(tmdata+n*n*tmi, 0, n*n*sizeof(complex double));
    for (size_t i = 0; i < n_symops; ++i) {
      const complex double *const M = symop_matrices + i*n*n;
      // tmp = M T
      cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
          n, n, n, &one, M, n, tmcopy->m, n, &zero, tmp, n);
      // tmdata[...] += tmp M* = M T M*
      cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
          n, n, n, &one, tmp, n, M, n, &one, tmdata + tmi*n*n, n);
    }
    for (size_t ii = 0; ii < n*n; ++ii)
      tmdata[n*n*tmi+ii] /= n_symops;
  }
  free(symop_matrices);
  qpms_tmatrix_free(tmcopy);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_symmetrise_tmdata_finite_group(
    complex double *tmdata, const size_t tmcount,
    const qpms_vswf_set_spec_t *bspec,
    const qpms_finite_group_t *pointgroup) {
  if (!(pointgroup->rep3d)) qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
      "This function requires pointgroup->rep3d to be set correctly!");
  return qpms_symmetrise_tmdata_irot3arr(tmdata, tmcount, bspec,
      pointgroup->order, pointgroup->rep3d);
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_irot3arr_inplace(
    qpms_tmatrix_t *T,
    size_t n_symops,
    const qpms_irot3_t *symops
    ) {
  if(qpms_symmetrise_tmdata_irot3arr(T->m, 1,
        T->spec, n_symops, symops) != QPMS_SUCCESS)
    return NULL;
  else return T;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_finite_group_inplace(
    qpms_tmatrix_t *T,
    const qpms_finite_group_t *pointgroup
    ) {
  if(qpms_symmetrise_tmdata_finite_group(T->m, 1,
        T->spec, pointgroup) != QPMS_SUCCESS)
    return NULL;
  else return T;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_involution_inplace(
                qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  qpms_tmatrix_t *t = qpms_tmatrix_apply_symop(T, M);
  const size_t n = T->spec->n;
  for(size_t i = 0; i < n*n; ++i)
    T->m[i] = 0.5 * (t->m[i] + T->m[i]);
  qpms_tmatrix_free(t);
  return T;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_involution(
                const qpms_tmatrix_t *T, 
                const complex double *M 
                )
{
  qpms_tmatrix_t *t = qpms_tmatrix_init(T->spec);
  const size_t n = T->spec->n;
  complex double tmp[n][n];
  // tmp = M T
  const complex double one = 1, zero = 0;
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasNoTrans,
      n, n, n, &one, M, n, T->m, n, &zero, tmp, n);
  // t->m = tmp M* = M T M*
  cblas_zgemm(CblasRowMajor,
      CblasNoTrans,
      CblasConjTrans,
      n, n, n, &one, tmp, n, M, n, &zero, t->m, n);
  for(size_t i = 0; i < n*n; ++i)
    t->m[i] = 0.5 * (t->m[i] + T->m[i]);
  return t;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_inf(const qpms_tmatrix_t *T) {
  qpms_tmatrix_t *t = qpms_tmatrix_copy(T);
  return qpms_tmatrix_symmetrise_C_inf_inplace(t);
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_inf_inplace(qpms_tmatrix_t *T) {
  const size_t n = T->spec->n;
  for (size_t row = 0; row < n; row++) {
    qpms_m_t rm = qpms_uvswfi2m(T->spec->ilist[row]);
    for (size_t col = 0; col < n; col++) {
      qpms_m_t cm = qpms_uvswfi2m(T->spec->ilist[col]);
      if (rm == cm)
        ;// No-op // t->m[n*row + col] = T->m[n*row + col];
      else
        T->m[n*row + col] = 0;
    }
  }
  return T;
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_N(const qpms_tmatrix_t *T, int N) {
  qpms_tmatrix_t *t = qpms_tmatrix_copy(T);
  return qpms_tmatrix_symmetrise_C_N_inplace(t, N);
}

qpms_tmatrix_t *qpms_tmatrix_symmetrise_C_N_inplace(qpms_tmatrix_t *T, int N) {
  const size_t n = T->spec->n;
  for (size_t row = 0; row < n; row++) {
    qpms_m_t rm = qpms_uvswfi2m(T->spec->ilist[row]);
    for (size_t col = 0; col < n; col++) {
      qpms_m_t cm = qpms_uvswfi2m(T->spec->ilist[col]);
      if (((rm - cm) % N) == 0)
        ; // T->m[n*row + col] = T->m[n*row + col];
      else
        T->m[n*row + col] = 0;
    }
  }
  return T;
}

bool qpms_tmatrix_isclose(const qpms_tmatrix_t *A, const qpms_tmatrix_t *B,
    const double rtol, const double atol)
{
  if (!qpms_vswf_set_spec_isidentical(A->spec, B->spec))
    return false;
  if (A->m == B->m)
    return true;
  const size_t n = A->spec->n;
  for (size_t i = 0; i < n*n; ++i) {
    const double tol = atol + rtol * (cabs(B->m[i]));
    if ( cabs(B->m[i] - A->m[i]) > tol )
      return false;
  }
  return true;
}
    

qpms_tmatrix_interpolator_t *qpms_tmatrix_interpolator_create(const size_t incount,
    const double *freqs, const qpms_tmatrix_t *ta, const gsl_interp_type *iptype//, const bool copy_bspec
    ) {
  if (incount <= 0) return NULL;
  qpms_tmatrix_interpolator_t *ip = malloc(sizeof(qpms_tmatrix_interpolator_t));
  /*
  if (copy_bspec) {
    ip->bspec = qpms_vswf_set_spec_copy(ta[0].spec);
    ip->owns_bspec = true;
  }
  else {
  */
    ip->bspec = ta[0].spec;
  //  ip->owns_bspec = false;
  //}
  const size_t n = ip->bspec->n;

  // check if all matrices have the same bspec
  for (size_t i = 0; i < incount; ++i)
    if (!qpms_vswf_set_spec_isidentical(ip->bspec, ta[i].spec))
      abort();

  if (!(ip->splines_real = calloc(n*n,sizeof(gsl_spline *)))) abort();
  if (!(ip->splines_imag = calloc(n*n,sizeof(gsl_spline *)))) abort();
  for (size_t row = 0; row < n; ++row)
    for (size_t col = 0; col < n; ++col) {
      double y_real[incount], y_imag[incount];
      bool n0_real = false, n0_imag = false;
      for (size_t i = 0; i < incount; ++i) {
        complex double telem = ta[i].m[n * row + col];
        if ((y_real[i] = creal(telem))) n0_real = true;
        if ((y_imag[i] = cimag(telem))) n0_imag = true;
      }
      if (n0_real) {
        gsl_spline *s =
        ip->splines_real[n * row + col] = gsl_spline_alloc(iptype, incount);
        if (gsl_spline_init(s, freqs, y_real, incount) != 0 /*GSL_SUCCESS*/) abort();
      }
      else ip->splines_real[n * row + col] = NULL;
     if (n0_imag) {
        gsl_spline *s =
        ip->splines_imag[n * row + col] = gsl_spline_alloc(iptype, incount);
        if (gsl_spline_init(s, freqs, y_imag, incount) != 0 /*GSL_SUCCESS*/) abort();
      }
      else ip->splines_imag[n * row + col] = NULL;
    }
  return ip;
}

void qpms_tmatrix_interpolator_free(qpms_tmatrix_interpolator_t *ip) {
  if (ip) {
    const size_t n = ip->bspec->n;
    for (size_t i = 0; i < n*n; ++i) {
      if (ip->splines_real[i]) gsl_spline_free(ip->splines_real[i]);
      if (ip->splines_imag[i]) gsl_spline_free(ip->splines_imag[i]);
    }
    //if (ip->owns_bspec)
    //  qpms_vswf_set_spec_free(ip->bspec);
    free(ip);
  }
}

qpms_tmatrix_t *qpms_tmatrix_interpolator_eval(const qpms_tmatrix_interpolator_t *ip, double freq) {
  qpms_tmatrix_t *t = qpms_tmatrix_init(ip->bspec);
  const size_t n = ip->bspec->n;
  for (size_t i = 0; i < n*n; ++i){
    if (ip->splines_real[i]) t->m[i] = gsl_spline_eval(ip->splines_real[i], freq, NULL /*does this work?*/);
    if (ip->splines_imag[i]) t->m[i] += I* gsl_spline_eval(ip->splines_imag[i], freq, NULL /*does this work?*/);
  }
  return t;
}



// ------------ Stupid implementation of qpms_scatsys_apply_symmetry() -------------

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

// The following functions are just to make qpms_scatsys_apply_symmetry more readable.
// They are not to be used elsewhere, as they do not perform any array capacity checks etc.

/// Compare two orbit types in a scattering system.
static bool orbit_types_equal(const qpms_ss_orbit_type_t *a, const qpms_ss_orbit_type_t *b) {
  if (a->size != b->size) return false;
  if (memcmp(a->action, b->action, a->size*sizeof(qpms_ss_orbit_pi_t))) return false;
  if (memcmp(a->tmatrices, b->tmatrices, a->size*sizeof(qpms_ss_tmi_t))) return false;
  return true;
}

// Extend the action to all particles in orbit if only the action on the 0th 
// particle has been filled.
static void extend_orbit_action(qpms_scatsys_t *ss, qpms_ss_orbit_type_t *ot) {
  for(qpms_ss_orbit_pi_t src = 1; src < ot->size; ++src) {
    // find any element g that sends 0 to src:
    qpms_gmi_t g;
    for (g = 0; g < ss->sym->order; ++g)
      if (ot->action[g] == src) break;
    assert(g < ss->sym->order);
    // invg sends src to 0
    qpms_gmi_t invg = qpms_finite_group_inv(ss->sym, g);
    for (qpms_gmi_t f = 0; f < ss->sym->order; ++f)
      // if f sends 0 to dest, then f * invg sends src to dest
      ot->action[src * ss->sym->order + 
        qpms_finite_group_mul(ss->sym,f,invg)] = ot->action[f];
  }
}

//Add orbit type to a scattering system, updating the ss->otspace_end pointer accordingly
static void add_orbit_type(qpms_scatsys_t *ss, const qpms_ss_orbit_type_t *ot_current) {
  qpms_ss_orbit_type_t * const ot_new = & (ss->orbit_types[ss->orbit_type_count]);
  ot_new->size = ot_current->size;

  const qpms_vswf_set_spec_t *bspec = ss->tm[ot_current->tmatrices[0]]->spec;
  const size_t bspecn = bspec->n;
  ot_new->bspecn = bspecn;
  
  const size_t actionsiz = sizeof(ot_current->action[0]) * ot_current->size 
    * ss->sym->order;
  ot_new->action = (void *) (ss->otspace_end);
  memcpy(ot_new->action, ot_current->action, actionsiz);
  // N.B. we copied mostly garbage ^^^, most of it is initialized just now:
  extend_orbit_action(ss, ot_new);
#ifdef DUMP_ORBIT_ACTION
  fprintf(stderr, "Orbit action:\n");
  for (qpms_gmi_t gmi = 0; gmi < ss->sym->order; ++gmi)  {
    const qpms_quat4d_t q = qpms_quat_4d_from_2c(ss->sym->rep3d[gmi].rot);
    fprintf(stderr, "%+d[%g %g %g %g] ", (int)ss->sym->rep3d[gmi].det,
        q.c1, q.ci, q.cj, q.ck);
    fprintf(stderr, "%s\t", (ss->sym->permrep && ss->sym->permrep[gmi])?
      ss->sym->permrep[gmi] : "");
    for (qpms_ss_orbit_pi_t pi = 0; pi < ot_new->size; ++pi) 
      fprintf(stderr, "%d\t", (int) ot_new->action[gmi + pi*ss->sym->order]);
    fprintf(stderr, "\n");
  }
#endif
  ss->otspace_end += actionsiz;
  
  const size_t tmsiz = sizeof(ot_current->tmatrices[0]) * ot_current->size;
  ot_new->tmatrices = (void *) (ss->otspace_end);
  memcpy(ot_new->tmatrices, ot_current->tmatrices, tmsiz);
  ss->otspace_end += tmsiz;

  const size_t irbase_sizes_siz = sizeof(ot_new->irbase_sizes[0]) * ss->sym->nirreps;
  ot_new->irbase_sizes = (void *) (ss->otspace_end);
  ss->otspace_end += irbase_sizes_siz;
  ot_new->irbase_cumsizes = (void *) (ss->otspace_end);
  ss->otspace_end += irbase_sizes_siz;
  ot_new->irbase_offsets = (void *) (ss->otspace_end);
  ss->otspace_end += irbase_sizes_siz;
  
  const size_t irbases_siz = sizeof(ot_new->irbases[0]) * SQ(ot_new->size * bspecn);
  ot_new->irbases = (void *) (ss->otspace_end);
  ss->otspace_end += irbases_siz;

  size_t lastbs, bs_cumsum = 0;
  for(qpms_iri_t iri = 0; iri < ss->sym->nirreps; ++iri) {
    ot_new->irbase_offsets[iri] = bs_cumsum * bspecn * ot_new->size;
    qpms_orbit_irrep_basis(&lastbs, 
        ot_new->irbases + bs_cumsum*ot_new->size*bspecn,
        ot_new, bspec, ss->sym, iri);
    ot_new->irbase_sizes[iri] = lastbs;
    bs_cumsum += lastbs;
    ot_new->irbase_cumsizes[iri] = bs_cumsum;
  }
  if(bs_cumsum != ot_new->size * bspecn)
    qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
        "The cumulative size of the symmetry-adapted bases is wrong; "
        "expected %d = %d * %d, got %d.",
        ot_new->size * bspecn, ot_new->size, bspecn, bs_cumsum);
  ot_new->instance_count = 0;
  ss->orbit_type_count++;
}


// Almost 200 lines. This whole thing deserves a rewrite!
qpms_scatsys_t *qpms_scatsys_apply_symmetry(const qpms_scatsys_t *orig, const qpms_finite_group_t *sym) {
  // TODO check data sanity

  qpms_l_t lMax = 0; // the overall lMax of all base specs.
  qpms_normalisation_t normalisation = QPMS_NORMALISATION_UNDEF;

  // First, determine the rough radius of the array; it should be nonzero
  // in order to particle position equivalence work correctly
  double lenscale = 0;
  {
    double minx = +INFINITY, miny = +INFINITY, minz = +INFINITY;
    double maxx = -INFINITY, maxy = -INFINITY, maxz = -INFINITY;
    for (qpms_ss_pi_t i = 0; i < orig->p_count; ++i) {
      minx = MIN(minx,orig->p[i].pos.x);
      miny = MIN(miny,orig->p[i].pos.y);
      minz = MIN(minz,orig->p[i].pos.z);
      maxx = MAX(maxx,orig->p[i].pos.x);
      maxy = MAX(maxy,orig->p[i].pos.y);
      maxz = MAX(maxz,orig->p[i].pos.z);
    }
    lenscale = (fabs(maxx)+fabs(maxy)+fabs(maxz)+(maxx-minx)+(maxy-miny)+(maxz-minz)) / 3;
  }

  // Second, check that there are no duplicit positions in the input system.
  for (qpms_ss_pi_t i = 0; i < orig->p_count; ++i)
    for (qpms_ss_pi_t j = 0; j < i; ++j)
      assert(!cart3_isclose(orig->p[i].pos, orig->p[j].pos, 0, QPMS_SCATSYS_LEN_RTOL * lenscale));

  // Allocate T-matrix, particle and particle orbit info arrays
  qpms_scatsys_t *ss = malloc(sizeof(qpms_scatsys_t));
  ss->lenscale = lenscale;
  ss->sym = sym;

  ss->tm_capacity = sym->order * orig->tm_count;
  ss->tm = malloc(ss->tm_capacity * sizeof(qpms_tmatrix_t *));

  ss->p_capacity = sym->order * orig->p_count;
  ss->p = malloc(ss->p_capacity * sizeof(qpms_particle_tid_t));
  ss->p_orbitinfo = malloc(ss->p_capacity * sizeof(qpms_ss_particle_orbitinfo_t));
  for (qpms_ss_pi_t pi = 0; pi < ss->p_capacity; ++pi) {
    ss->p_orbitinfo[pi].t = QPMS_SS_P_ORBITINFO_UNDEF;
    ss->p_orbitinfo[pi].p = QPMS_SS_P_ORBITINFO_UNDEF;
  }

  // Copy T-matrices; checking for duplicities
  
  ss->max_bspecn = 0; // We'll need it later.for memory alloc estimates.

  qpms_ss_tmi_t tm_dupl_remap[ss->tm_capacity]; // Auxilliary array to label remapping the indices after ignoring t-matrix duplicities
  ss->tm_count = 0;
  for (qpms_ss_tmi_t i = 0; i < orig->tm_count; ++i) {
    qpms_ss_tmi_t j;
    for (j = 0; j < ss->tm_count; ++j) 
      if (qpms_tmatrix_isclose(orig->tm[i], ss->tm[j], QPMS_SCATSYS_TMATRIX_RTOL, QPMS_SCATSYS_TMATRIX_ATOL)) {
        break;
      }
    if (j == ss->tm_count) { // duplicity not found, copy the t-matrix
      ss->tm[j] = qpms_tmatrix_copy(orig->tm[i]);
      ss->max_bspecn = MAX(ss->tm[j]->spec->n, ss->max_bspecn);
      lMax = MAX(lMax, ss->tm[j]->spec->lMax);
      ++(ss->tm_count);
    } 
    tm_dupl_remap[i] = j;
    if (normalisation == QPMS_NORMALISATION_UNDEF)
      normalisation = ss->tm[i]->spec->norm;
    // We expect all bspec norms to be the same.
    else QPMS_ENSURE(normalisation == ss->tm[j]->spec->norm,
        "Normalisation convention must be the same for all T-matrices."
        " %d != %d\n", normalisation, ss->tm[j]->spec->norm);
  }

  // Copy particles, remapping the t-matrix indices
  for (qpms_ss_pi_t i = 0; i < orig->p_count; ++(i)) {
    ss->p[i] = orig->p[i];
    ss->p[i].tmatrix_id = tm_dupl_remap[ss->p[i].tmatrix_id];
  }
  ss->p_count = orig->p_count;

  // allocate t-matrix symmetry map
  ss->tm_sym_map = malloc(sizeof(qpms_ss_tmi_t) * sym->order * sym->order * ss->tm_count);

  // Extend the T-matrices list by the symmetry operations
  for (qpms_ss_tmi_t tmi = 0; tmi < ss->tm_count; ++tmi) 
    for (qpms_gmi_t gmi = 0; gmi < sym->order; ++gmi){
      const size_t d = ss->tm[tmi]->spec->n;
      complex double M[d][d]; // transformation matrix
      qpms_irot3_uvswfi_dense(M[0], ss->tm[tmi]->spec, sym->rep3d[gmi]);
      qpms_tmatrix_t *transformed = qpms_tmatrix_apply_symop(ss->tm[tmi], M[0]);
      qpms_ss_tmi_t tmj;
      for (tmj = 0; tmj < ss->tm_count; ++tmj)
        if (qpms_tmatrix_isclose(transformed, ss->tm[tmj], QPMS_SCATSYS_TMATRIX_RTOL, QPMS_SCATSYS_TMATRIX_ATOL))
          break;
      if (tmj < ss->tm_count) { // HIT, transformed T-matrix already exists
        qpms_tmatrix_free(transformed);
      } else { // MISS, save the matrix and increment the count
        ss->tm[ss->tm_count] = transformed;
        ++(ss->tm_count);
      }
      ss->tm_sym_map[gmi + tmi * sym->order] = tmj; // In any case, tmj now indexes the correct transformed matrix
    }
  // Possibly free some space using the new ss->tm_count instead of (old) ss->tm_count*sym->order
  ss->tm_sym_map = realloc(ss->tm_sym_map, sizeof(qpms_ss_tmi_t) * sym->order * ss->tm_count);
  // tm could be realloc'd as well, but those are just pointers, not likely many.
 

  // allocate particle symmetry map
  ss->p_sym_map = malloc(sizeof(qpms_ss_pi_t) * sym->order * sym->order * ss->p_count);
  // allocate orbit type array (TODO realloc later if too long)
  ss->orbit_type_count = 0;
  ss->orbit_types = calloc(ss->p_count, sizeof(qpms_ss_orbit_type_t));

  ss->otspace_end = ss->otspace = malloc( // reallocate later
      (sizeof(qpms_ss_orbit_pi_t) * sym->order * sym->order
       + sizeof(qpms_ss_tmi_t) * sym->order
       + 3*sizeof(size_t) * sym->nirreps
       + sizeof(complex double) * SQ(sym->order * ss->max_bspecn)) * ss->p_count
       );
  
  // Workspace for the orbit type determination
  qpms_ss_orbit_type_t ot_current;
  qpms_ss_orbit_pi_t ot_current_action[sym->order * sym->order];
  qpms_ss_tmi_t ot_current_tmatrices[sym->order];
  
  qpms_ss_pi_t current_orbit[sym->order];
  ot_current.action = ot_current_action;
  ot_current.tmatrices = ot_current_tmatrices;


  // Extend the particle list by the symmetry operations, check that particles mapped by symmetry ops on themselves
  // have the correct symmetry
  // TODO this could be sped up to O(npart * log(npart)); let's see later whether needed.
  for (qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    const bool new_orbit = (ss->p_orbitinfo[pi].t == QPMS_SS_P_ORBITINFO_UNDEF); // TODO build the orbit!!!
    if (new_orbit){
      current_orbit[0] = pi;
      ot_current.size = 1; 
      ot_current.tmatrices[0] = ss->p[pi].tmatrix_id;
#ifdef DUMP_PARTICLE_POSITIONS
      cart3_t pos = ss->p[pi].pos;
      fprintf(stderr, "An orbit [%.4g, %.4g, %.4g] => ", pos.x, pos.y, pos.z);
#endif
    }

    for (qpms_gmi_t gmi = 0; gmi < sym->order; ++gmi) {
      cart3_t newpoint = qpms_irot3_apply_cart3(sym->rep3d[gmi], ss->p[pi].pos);
      qpms_ss_tmi_t new_tmi = ss->tm_sym_map[gmi + ss->p[pi].tmatrix_id * sym->order]; // transformed t-matrix index
      qpms_ss_pi_t pj;
      for (pj = 0; pj < ss->p_count; ++pj)
        if (cart3_isclose(newpoint, ss->p[pj].pos, 0, ss->lenscale * QPMS_SCATSYS_LEN_RTOL)) {
          if (new_tmi != ss->p[pj].tmatrix_id)
            qpms_pr_error("The %d. particle with coords (%lg, %lg, %lg) "
                "is mapped to %d. another (or itself) with cords (%lg, %lg, %lg) "
                "without having the required symmetry", (int)pi,
                ss->p[pi].pos.x, ss->p[pi].pos.y, ss->p[pi].pos.z,
                (int)pj, ss->p[pj].pos.x, ss->p[pj].pos.y, ss->p[pj].pos.z);
          break;
        }
      if (pj < ss->p_count) { // HIT, the particle is transformed to an "existing" one.
        ;
      } else { // MISS, the symmetry transforms the particle to a new location, so add it.
        qpms_particle_tid_t newparticle = {newpoint, new_tmi};
        ss->p[ss->p_count] = newparticle;
        ++(ss->p_count);
#ifdef DUMP_PARTICLE_POSITIONS
        if(new_orbit)
          fprintf(stderr, "[%.4g, %.4g, %.4g] ", newpoint.x, newpoint.y, newpoint.z);
#endif
      }
      ss->p_sym_map[gmi + pi * sym->order] = pj;

      if (new_orbit) {
        // Now check whether the particle (result of the symmetry op) is already in the current orbit
        qpms_ss_orbit_pi_t opj;
        for (opj = 0; opj < ot_current.size; ++opj)
          if (current_orbit[opj] == pj) break; // HIT, pj already on current orbit
        if (opj == ot_current.size) { // MISS, pj is new on the orbit, extend the size and set the T-matrix id
          current_orbit[opj] = pj;
          ++ot_current.size;
          ot_current.tmatrices[opj] = ss->p[pj].tmatrix_id;
        }
        ot_current.action[gmi] = opj;
      }
    }
    if (new_orbit) { // Now compare if the new orbit corresponds to some of the existing types.
#ifdef DUMP_PARTICLE_POSITIONS
      fputc('\n', stderr);
#endif
      qpms_ss_oti_t oti;
      for(oti = 0; oti < ss->orbit_type_count; ++oti)
        if (orbit_types_equal(&ot_current, &(ss->orbit_types[oti]))) break; // HIT, orbit type already exists
      assert(0 == sym->order % ot_current.size);
      if (oti == ss->orbit_type_count)  // MISS, add the new orbit type
        add_orbit_type(ss, &ot_current);

      // Walk through the orbit again and set up the orbit info of the particles
      for (qpms_ss_orbit_pi_t opi = 0; opi < ot_current.size; ++opi) {
        const qpms_ss_pi_t pi_opi = current_orbit[opi];
        ss->p_orbitinfo[pi_opi].t = oti;
        ss->p_orbitinfo[pi_opi].p = opi;
        ss->p_orbitinfo[pi_opi].osn = ss->orbit_types[oti].instance_count;
      }
      ss->orbit_types[oti].instance_count++;
    }
  }
  // Possibly free some space using the new ss->p_count instead of (old) ss->p_count*sym->order
  ss->p_sym_map = realloc(ss->p_sym_map, sizeof(qpms_ss_pi_t) * sym->order * ss->p_count);
  ss->p = realloc(ss->p, sizeof(qpms_particle_tid_t) * ss->p_count); 
  ss->p_orbitinfo = realloc(ss->p_orbitinfo, sizeof(qpms_ss_particle_orbitinfo_t)*ss->p_count);
  ss->p_capacity = ss->p_count;

  {  // Reallocate the orbit type data space and update the pointers if needed.
    size_t otspace_sz = ss->otspace_end - ss->otspace;
    char *old_otspace = ss->otspace;
    ss->otspace = realloc(ss->otspace, otspace_sz);
    ptrdiff_t shift = ss->otspace - old_otspace;
    if(shift) {
      for (size_t oi = 0; oi < ss->orbit_type_count; ++oi) {
        ss->orbit_types[oi].action = (void *)(((char *) (ss->orbit_types[oi].action)) + shift);
        ss->orbit_types[oi].tmatrices = (void *)(((char *) (ss->orbit_types[oi].tmatrices)) + shift);
        ss->orbit_types[oi].irbase_sizes = (void *)(((char *) (ss->orbit_types[oi].irbase_sizes)) + shift);
        ss->orbit_types[oi].irbase_cumsizes = (void *)(((char *) (ss->orbit_types[oi].irbase_cumsizes)) + shift);
        ss->orbit_types[oi].irbase_offsets = (void *)(((char *) (ss->orbit_types[oi].irbase_offsets)) + shift);
        ss->orbit_types[oi].irbases = (void *)(((char *) (ss->orbit_types[oi].irbases)) + shift);
      }
      ss->otspace_end += shift;
    }
  }

  // Set ss->fecv_size and ss->fecv_pstarts
  ss->fecv_size = 0;
  ss->fecv_pstarts = malloc(ss->p_count * sizeof(size_t));
  for (qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    ss->fecv_pstarts[pi] = ss->fecv_size;
    ss->fecv_size += ss->tm[ss->p[pi].tmatrix_id]->spec->n; // That's a lot of dereferencing!
  }

  ss->saecv_sizes = malloc(sizeof(size_t) * sym->nirreps); if (!ss->saecv_sizes) abort();
  ss->saecv_ot_offsets = malloc(sizeof(size_t) * sym->nirreps * ss->orbit_type_count);
  if (!ss->saecv_ot_offsets) abort();
  for(qpms_iri_t iri = 0; iri < sym->nirreps; ++iri) {
    ss->saecv_sizes[iri] = 0;
    for(qpms_ss_oti_t oti = 0; oti < ss->orbit_type_count; ++oti) {
      ss->saecv_ot_offsets[iri * ss->orbit_type_count + oti] = ss->saecv_sizes[iri];
      const qpms_ss_orbit_type_t *ot = &(ss->orbit_types[oti]);
      ss->saecv_sizes[iri] += ot->instance_count * ot->irbase_sizes[iri];
    }
  }

  qpms_ss_pi_t p_ot_cumsum = 0;
  for (qpms_ss_oti_t oti = 0; oti < ss->orbit_type_count; ++oti) {
    qpms_ss_orbit_type_t *ot = ss->orbit_types + oti;
    ot->p_offset = p_ot_cumsum;
    p_ot_cumsum += ot->size * ot->instance_count;
  }

  // Set ss->p_by_orbit[]
  QPMS_CRASHING_MALLOC(ss->p_by_orbit, sizeof(qpms_ss_pi_t) * ss->p_count);
  for (qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    const qpms_ss_particle_orbitinfo_t *oi = ss->p_orbitinfo + pi;
    const qpms_ss_oti_t oti = oi->t;
    const qpms_ss_orbit_type_t *ot = ss->orbit_types + oti;
    ss->p_by_orbit[ot->p_offset + ot->size * oi->osn + oi->p] = pi;
  }
  
  ss->c = qpms_trans_calculator_init(lMax, normalisation);
  return ss;
}


void qpms_scatsys_free(qpms_scatsys_t *ss) {
  if(ss) {
    free(ss->tm);
    free(ss->p);
    free(ss->fecv_pstarts);
    free(ss->tm_sym_map);
    free(ss->p_sym_map);
    free(ss->otspace);
    free(ss->p_orbitinfo);
    free(ss->orbit_types);
    free(ss->saecv_sizes);
    free(ss->p_by_orbit);
    qpms_trans_calculator_free(ss->c);
  }
  free(ss);
}



// (copypasta from symmetries.c)
// TODO at some point, maybe support also other norms.
// (perhaps just use qpms_normalisation_t_factor() at the right places)
static inline void check_norm_compat(const qpms_vswf_set_spec_t *s)
{
  switch (qpms_normalisation_t_normonly(s->norm)) {
    case QPMS_NORMALISATION_POWER:
      break;
    case QPMS_NORMALISATION_SPHARM:
      break;
    default:
      abort(); // Only SPHARM and POWER norms are supported right now.
  }
}

complex double *qpms_orbit_action_matrix(complex double *target,
    const qpms_ss_orbit_type_t *ot, const qpms_vswf_set_spec_t *bspec,
    const qpms_finite_group_t *sym, const qpms_gmi_t g) {
  assert(sym); assert(g < sym->order);
  assert(sym->rep3d);
  assert(ot); assert(ot->size > 0);
  check_norm_compat(bspec);
  const size_t n = bspec->n;
  const qpms_gmi_t N = ot->size;
  if (target == NULL)
    target = malloc(n*n*N*N*sizeof(complex double));
  if (target == NULL) abort();
  memset(target, 0, n*n*N*N*sizeof(complex double));
  complex double tmp[n][n]; // this is the 'per-particle' action
  qpms_irot3_uvswfi_dense(tmp[0], bspec, sym->rep3d[g]); 
  for(qpms_gmi_t Col = 0; Col < ot->size; ++Col) {
    // Row is the 'destination' of the symmetry operation, Col is the 'source'
    const qpms_gmi_t Row = ot->action[sym->order * Col + g];
    for(size_t row = 0; row < bspec->n; ++row)
      for(size_t col = 0; col < bspec->n; ++col)
        target[n*n*N*Row + n*Col + n*N*row + col] = conj(tmp[row][col]); //CHECKCONJ
  }
#ifdef DUMP_ACTIONMATRIX
  fprintf(stderr,"%d: %s\n", 
      (int) g, (sym->permrep && sym->permrep[g])?
      sym->permrep[g] : "");
  for (size_t Row = 0; Row < ot->size; ++Row) {
    fprintf(stderr, "--------------------------\n");
    for (size_t row = 0; row < bspec->n; ++row) {
      for (size_t Col = 0; Col < ot->size; ++Col) {
        fprintf(stderr, "| ");
        for (size_t col = 0; col < bspec->n; ++col)
          fprintf(stderr, "%+2.3f%+2.3fj ", creal(target[n*n*N*Row + n*Col + n*N*row + col]),cimag(target[n*n*N*Row + n*Col + n*N*row + col]));
      }
      fprintf(stderr, "|\n");
    }
  }
  fprintf(stderr, "-------------------------------\n\n");
#endif

  return target;
}

complex double *qpms_orbit_irrep_projector_matrix(complex double *target,
    const qpms_ss_orbit_type_t *ot, const qpms_vswf_set_spec_t *bspec,
    const qpms_finite_group_t *sym, const qpms_iri_t iri) {
  assert(sym);
  assert(sym->rep3d);
  assert(ot); assert(ot->size > 0);
  assert(iri < sym->nirreps); assert(sym->irreps);
  check_norm_compat(bspec);
  const size_t n = bspec->n;
  const qpms_gmi_t N = ot->size;
  if (target == NULL)
    target = malloc(n*n*N*N*sizeof(complex double));
  if (target == NULL) abort();
  memset(target, 0, n*n*N*N*sizeof(complex double));
  // Workspace for the orbit group action matrices
  complex double *tmp = malloc(n*n*N*N*sizeof(complex double));
  const int d = sym->irreps[iri].dim;
  double prefac = d / (double) sym->order;
  for(int partner = 0; partner < d; ++partner) {
    for(qpms_gmi_t g = 0; g < sym->order; ++g) {
      // We use the diagonal elements of D_g
      complex double D_g_conj = sym->irreps[iri].m[g*d*d + partner*d + partner];
#ifdef DUMP_ACTIONMATRIX
      fprintf(stderr,"(factor %+g%+gj) ", creal(D_g_conj), cimag(D_g_conj));
#endif
      qpms_orbit_action_matrix(tmp, ot, bspec, sym, g);
      // TODO kahan sum?
      for(size_t i = 0; i < n*n*N*N; ++i)
        target[i] += prefac * D_g_conj * tmp[i];
    }
  }
  free(tmp);
#ifdef DUMP_PROJECTORMATRIX
  fprintf(stderr,"Projector %d (%s):\n", (int) iri, 
      sym->irreps[iri].name?sym->irreps[iri].name:"");
  for (size_t Row = 0; Row < ot->size; ++Row) {
    fprintf(stderr, "--------------------------\n");
    for (size_t row = 0; row < bspec->n; ++row) {
      for (size_t Col = 0; Col < ot->size; ++Col) {
        fprintf(stderr, "| ");
        for (size_t col = 0; col < bspec->n; ++col)
          fprintf(stderr, "%+2.3f%+2.3fj ", creal(target[n*n*N*Row + n*Col + n*N*row + col]),cimag(target[n*n*N*Row + n*Col + n*N*row + col]));
      }
      fprintf(stderr, "|\n");
    }
  }
  fprintf(stderr, "-------------------------------\n\n");
#endif
  return target;
}

#define SVD_ATOL 1e-8

complex double *qpms_orbit_irrep_basis(size_t *basis_size, 
    complex double *target,
    const qpms_ss_orbit_type_t *ot, const qpms_vswf_set_spec_t *bspec,
    const qpms_finite_group_t *sym, const qpms_iri_t iri) {
  assert(sym);
  assert(sym->rep3d);
  assert(ot); assert(ot->size > 0);
  assert(iri < sym->nirreps); assert(sym->irreps);
  check_norm_compat(bspec);
  const size_t n = bspec->n;
  const qpms_gmi_t N = ot->size;
  const bool newtarget = (target == NULL);
  if (newtarget)
    QPMS_CRASHING_MALLOC(target,n*n*N*N*sizeof(complex double));
  if (target == NULL) abort();
  memset(target, 0, n*n*N*N*sizeof(complex double));

  // Get the projector (also workspace for right sg. vect.)
  complex double *projector = qpms_orbit_irrep_projector_matrix(NULL,
    ot, bspec, sym, iri);
  if(!projector) abort();
  // Workspace for the right singular vectors.
  complex double *V_H = malloc(n*n*N*N*sizeof(complex double));
  if(!V_H) abort();
  // THIS SHOULD NOT BE NECESSARY
  complex double *U = malloc(n*n*N*N*sizeof(complex double));
  if(!U) abort();
  double *s = malloc(n*N*sizeof(double)); if(!s) abort();

  int info = LAPACKE_zgesdd(LAPACK_ROW_MAJOR,
      'A', // jobz; overwrite projector with left sg.vec. and write right into V_H
      n*N /* m */, n*N /* n */, projector /* a */, n*N /* lda */,
      s /* s */, U /* u */, n*N /* ldu, irrelev. */, V_H /* vt */,
      n*N /* ldvt */);
  if (info) qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
        "Something went wrong with the SVD.");


  size_t bs;
  for(bs = 0; bs < n*N; ++bs) {
#if 0
    qpms_pr_debug_at_flf(__FILE__, __LINE__, __func__,
        "%d. irrep, %zd. SV: %.16lf", (int) iri, bs, s[bs]);
#endif
    if(s[bs] > 1 + SVD_ATOL) qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
        "%zd. SV too large: %.16lf", bs, s[bs]);
    if(s[bs] > SVD_ATOL && fabs(1-s[bs]) > SVD_ATOL) 
      qpms_pr_error_at_flf(__FILE__, __LINE__, __func__,
        "%zd. SV in the 'wrong' interval: %.16lf", bs, s[bs]);
    if(s[bs] < SVD_ATOL) break;
  }

  memcpy(target, V_H, bs*N*n*sizeof(complex double));
  if(newtarget) QPMS_CRASHING_REALLOC(target, bs*N*n*sizeof(complex double));
  if(basis_size) *basis_size = bs;

  free(U);
  free(V_H);
  free(projector);
  return target;
}

complex double *qpms_scatsys_irrep_transform_matrix(complex double *U,
    const qpms_scatsys_t *ss, qpms_iri_t iri) {
  const size_t packedlen = ss->saecv_sizes[iri];
  const size_t full_len = ss->fecv_size;
  if (U == NULL)
    QPMS_CRASHING_MALLOC(U,full_len * packedlen * sizeof(complex double));
  memset(U, 0, full_len * packedlen * sizeof(complex double));

  size_t fullvec_offset = 0;

  for(qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    const qpms_ss_oti_t oti = ss->p_orbitinfo[pi].t;
    const qpms_ss_orbit_type_t *const ot = ss->orbit_types + oti;
    const qpms_ss_osn_t osn = ss->p_orbitinfo[pi].osn;
    const qpms_ss_orbit_pi_t opi = ss->p_orbitinfo[pi].p;
    // This is where the particle's orbit starts in the "packed" vector:
    const size_t packed_orbit_offset = ss->saecv_ot_offsets[iri*ss->orbit_type_count + oti]
      + osn * ot->irbase_sizes[iri];
    // Orbit coeff vector's full size:
    const size_t orbit_fullsize = ot->size * ot->bspecn;
    const size_t orbit_packedsize = ot->irbase_sizes[iri];
    // This is the orbit-level matrix projecting the whole orbit onto the irrep.
    const complex double *om = ot->irbases + ot->irbase_offsets[iri];

    for (size_t prow = 0; prow < orbit_packedsize; ++prow) 
      for (size_t pcol = 0; pcol < ot->bspecn; ++pcol) 
        U[full_len * (packed_orbit_offset + prow) + (fullvec_offset + pcol)]
          = om[orbit_fullsize * prow + (opi * ot->bspecn + pcol)];
    fullvec_offset += ot->bspecn;
  }

  return U;
}

complex double *qpms_scatsys_irrep_pack_matrix_stupid(complex double *target_packed,
		const complex double *orig_full, const qpms_scatsys_t *ss,
		qpms_iri_t iri){
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) // THIS IS A BIT PROBLEMATIC, TODO how to deal with empty irreps?
    return target_packed; 
  const size_t full_len = ss->fecv_size;
  if (target_packed == NULL)
    target_packed = malloc(SQ(packedlen)*sizeof(complex double));
  if (target_packed == NULL) abort();
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  // Workspace for the intermediate matrix
  complex double *tmp; 
  QPMS_CRASHING_MALLOC(tmp, full_len * packedlen * sizeof(complex double));

  complex double *U = qpms_scatsys_irrep_transform_matrix(NULL, ss, iri);

  const complex double one = 1, zero = 0;
  
  // tmp = F U*
  cblas_zgemm(
      CblasRowMajor, CblasNoTrans, CblasConjTrans,
      full_len /*M*/, packedlen /*N*/, full_len /*K*/,
      &one /*alpha*/, orig_full/*A*/, full_len/*ldA*/, 
      U /*B*/, full_len/*ldB*/, 
      &zero /*beta*/, tmp /*C*/, packedlen /*LDC*/);
  // target = U tmp
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
      packedlen /*M*/, packedlen /*N*/, full_len /*K*/,
      &one /*alpha*/, U/*A*/, full_len/*ldA*/,
      tmp /*B*/, packedlen /*ldB*/, &zero /*beta*/,
      target_packed /*C*/, packedlen /*ldC*/);
  
  free(tmp);
  free(U);
  return target_packed;
}
  
/// Transforms a big "packed" matrix into the full basis (trivial implementation).
complex double *qpms_scatsys_irrep_unpack_matrix_stupid(complex double *target_full,
		const complex double *orig_packed, const qpms_scatsys_t *ss,
		qpms_iri_t iri, bool add){
  const size_t packedlen = ss->saecv_sizes[iri];
  const size_t full_len = ss->fecv_size;
  if (target_full == NULL)
    target_full = malloc(SQ(full_len)*sizeof(complex double));
  if (target_full == NULL) abort();
  if(!add) memset(target_full, 0, SQ(full_len)*sizeof(complex double));

  if(!packedlen) return target_full; // Empty irrep, do nothing.

  // Workspace for the intermediate matrix result
  complex double *tmp; 
  QPMS_CRASHING_MALLOC(tmp, full_len * packedlen * sizeof(complex double));

  complex double *U = qpms_scatsys_irrep_transform_matrix(NULL, ss, iri);

  const complex double one = 1, zero = 0;
  
  // tmp = P U
  cblas_zgemm(
      CblasRowMajor, CblasNoTrans, CblasNoTrans,
      packedlen /*M*/, full_len /*N*/, packedlen /*K*/,
      &one /*alpha*/, orig_packed/*A*/, packedlen/*ldA*/, 
      U /*B*/, full_len/*ldB*/, 
      &zero /*beta*/, tmp /*C*/, full_len /*LDC*/);
  // target += U* tmp
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
      full_len /*M*/, full_len /*N*/, packedlen /*K*/,
      &one /*alpha*/, U/*A*/, full_len/*ldA*/,
      tmp /*B*/, full_len /*ldB*/, &one /*beta*/,
      target_full /*C*/, full_len /*ldC*/);
  free(tmp);
  free(U);
  return target_full;
}

complex double *qpms_scatsys_irrep_pack_matrix(complex double *target_packed,
		const complex double *orig_full, const qpms_scatsys_t *ss,
		qpms_iri_t iri){
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) // THIS IS A BIT PROBLEMATIC, TODO how to deal with empty irreps?
    return target_packed; 
  const size_t full_len = ss->fecv_size;
  if (target_packed == NULL)
    target_packed = malloc(SQ(packedlen)*sizeof(complex double));
  if (target_packed == NULL) abort();
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp = malloc(sizeof(complex double) * SQ(ss->max_bspecn)
      * ss->sym->order); if (!tmp) abort();

  const complex double one = 1, zero = 0;

  size_t fullvec_offsetR = 0;
  for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) { //Row loop
    const qpms_ss_oti_t otiR = ss->p_orbitinfo[piR].t;
    const qpms_ss_orbit_type_t *const otR = ss->orbit_types + otiR;
    const qpms_ss_osn_t osnR = ss->p_orbitinfo[piR].osn;
    const qpms_ss_orbit_pi_t opiR = ss->p_orbitinfo[piR].p;
    // This is where the particle's orbit starts in the "packed" vector:
    const size_t packed_orbit_offsetR =
      ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiR] 
      + osnR * otR->irbase_sizes[iri];
    // Orbit coeff vector's full size:
    const size_t orbit_fullsizeR = otR->size * otR->bspecn;
    const size_t particle_fullsizeR = otR->bspecn;
    const size_t orbit_packedsizeR = otR->irbase_sizes[iri];
    // This is the orbit-level matrix projecting the whole orbit onto the irrep.
    const complex double *omR = otR->irbases + otR->irbase_offsets[iri];

    size_t fullvec_offsetC = 0;
    if(orbit_packedsizeR) { // avoid zgemm crash on empty irrep
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) { //Column loop
        const qpms_ss_oti_t otiC = ss->p_orbitinfo[piC].t;
        const qpms_ss_orbit_type_t *const otC = ss->orbit_types + otiC;
        const qpms_ss_osn_t osnC = ss->p_orbitinfo[piC].osn;
        const qpms_ss_orbit_pi_t opiC = ss->p_orbitinfo[piC].p;
        // This is where the particle's orbit starts in the "packed" vector:
        const size_t packed_orbit_offsetC = 
          ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiC]
          + osnC * otC->irbase_sizes[iri];
        // Orbit coeff vector's full size:
        const size_t orbit_fullsizeC = otC->size * otC->bspecn;
        const size_t particle_fullsizeC = otC->bspecn;
        const size_t orbit_packedsizeC = otC->irbase_sizes[iri];
        // This is the orbit-level matrix projecting the whole orbit onto the irrep.
        const complex double *omC = otC->irbases + otC->irbase_offsets[iri];

        if(orbit_packedsizeC) { // avoid zgemm crash on empty irrep
          // tmp[oiR|piR,piC] = ∑_K M[piR,K] U*[K,piC]
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
              particle_fullsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeC /*K*/,
              &one /*alpha*/, orig_full + full_len*fullvec_offsetR + fullvec_offsetC/*A*/,
              full_len/*ldA*/, 
              omC + opiC*particle_fullsizeC /*B*/,
              orbit_fullsizeC/*ldB*/, &zero /*beta*/,
              tmp /*C*/, orbit_packedsizeC /*LDC*/);

          // target[oiR|piR,oiC|piC] += U[...] tmp[...]
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              orbit_packedsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeR /*K*/,
              &one /*alpha*/, omR + opiR*particle_fullsizeR/*A*/,
              orbit_fullsizeR/*ldA*/,
              tmp /*B*/, orbit_packedsizeC /*ldB*/, &one /*beta*/,
              target_packed + packedlen*packed_orbit_offsetR + packed_orbit_offsetC /*C*/,
              packedlen /*ldC*/);
        }
        fullvec_offsetC += otC->bspecn;
      }
    }
    fullvec_offsetR += otR->bspecn;
  }
  free(tmp);
  return target_packed;
}


/// Transforms a big "packed" matrix into the full basis.
/** TODO doc */
complex double *qpms_scatsys_irrep_unpack_matrix(complex double *target_full,
		const complex double *orig_packed, const qpms_scatsys_t *ss,
		qpms_iri_t iri, bool add){
  const size_t packedlen = ss->saecv_sizes[iri];
  const size_t full_len = ss->fecv_size;
  if (target_full == NULL)
    target_full = malloc(SQ(full_len)*sizeof(complex double));
  if (target_full == NULL) abort();
  if(!add) memset(target_full, 0, SQ(full_len)*sizeof(complex double));

  if(!packedlen) return target_full; // Empty irrep, do nothing.

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp = malloc(sizeof(complex double) * SQ(ss->max_bspecn)
      * ss->sym->order); if (!tmp) abort();

  const complex double one = 1, zero = 0;

  size_t fullvec_offsetR = 0;
  for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) { //Row loop
    const qpms_ss_oti_t otiR = ss->p_orbitinfo[piR].t;
    const qpms_ss_orbit_type_t *const otR = ss->orbit_types + otiR;
    const qpms_ss_osn_t osnR = ss->p_orbitinfo[piR].osn;
    const qpms_ss_orbit_pi_t opiR = ss->p_orbitinfo[piR].p;
    // This is where the particle's orbit starts in the "packed" vector:
    const size_t packed_orbit_offsetR =
      ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiR] 
      + osnR * otR->irbase_sizes[iri];
    // Orbit coeff vector's full size:
    const size_t orbit_fullsizeR = otR->size * otR->bspecn;
    const size_t particle_fullsizeR = otR->bspecn;
    const size_t orbit_packedsizeR = otR->irbase_sizes[iri];
    // This is the orbit-level matrix projecting the whole orbit onto the irrep.
    const complex double *omR = otR->irbases + otR->irbase_offsets[iri];

    size_t fullvec_offsetC = 0;
    if (orbit_packedsizeR) // avoid crash on empty irrep
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) { //Column loop
        const qpms_ss_oti_t otiC = ss->p_orbitinfo[piC].t;
        const qpms_ss_orbit_type_t *const otC = ss->orbit_types + otiC;
        const qpms_ss_osn_t osnC = ss->p_orbitinfo[piC].osn;
        const qpms_ss_orbit_pi_t opiC = ss->p_orbitinfo[piC].p;
        // This is where the particle's orbit starts in the "packed" vector:
        const size_t packed_orbit_offsetC = 
          ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiC]
          + osnC * otC->irbase_sizes[iri];
        // Orbit coeff vector's full size:
        const size_t orbit_fullsizeC = otC->size * otC->bspecn;
        const size_t particle_fullsizeC = otC->bspecn;
        const size_t orbit_packedsizeC = otC->irbase_sizes[iri];
        // This is the orbit-level matrix projecting the whole orbit onto the irrep.
        const complex double *omC = otC->irbases + otC->irbase_offsets[iri];
        if (orbit_packedsizeC) { // avoid crash on empty irrep
          // tmp = P U
          // tmp[oiR|piR,piC] = ∑_K M[piR,K] U[K,piC]
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            orbit_packedsizeR /*M*/, particle_fullsizeC /*N*/, orbit_packedsizeC /*K*/,
            &one /*alpha*/, orig_packed + packedlen*packed_orbit_offsetR + packed_orbit_offsetC/*A*/,
            packedlen/*ldA*/, 
            omC + opiC*particle_fullsizeC /*B*/,
            orbit_fullsizeC/*ldB*/, &zero /*beta*/,
            tmp /*C*/, particle_fullsizeC /*LDC*/);

          // target[oiR|piR,oiC|piC] += U*[...] tmp[...]
          cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
            particle_fullsizeR /*M*/, particle_fullsizeC /*N*/, orbit_packedsizeR /*K*/,
            &one /*alpha*/, omR + opiR*particle_fullsizeR/*A*/,
            orbit_fullsizeR/*ldA*/,
            tmp /*B*/, particle_fullsizeC /*ldB*/, &one /*beta*/,
            target_full + full_len*fullvec_offsetR + fullvec_offsetC /*C*/,
            full_len /*ldC*/);
        }
      fullvec_offsetC += otC->bspecn;
    }
    fullvec_offsetR += otR->bspecn;
  }

  free(tmp);
  return target_full;
}

/// Projects a "big" vector onto an irrep.
/** TODO doc */
complex double *qpms_scatsys_irrep_pack_vector(complex double *target_packed,
		const complex double *orig_full, const qpms_scatsys_t *ss,
		const qpms_iri_t iri) {
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) return target_packed; // Empty irrep
  if (target_packed == NULL)
    target_packed = malloc(packedlen*sizeof(complex double));
  if (target_packed == NULL) abort();
  memset(target_packed, 0, packedlen*sizeof(complex double));

  const complex double one = 1;

  size_t fullvec_offset = 0;
  for(qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    const qpms_ss_oti_t oti = ss->p_orbitinfo[pi].t;
    const qpms_ss_orbit_type_t *const ot = ss->orbit_types + oti;
    const qpms_ss_osn_t osn = ss->p_orbitinfo[pi].osn;
    const qpms_ss_orbit_pi_t opi = ss->p_orbitinfo[pi].p;
    // This is where the particle's orbit starts in the "packed" vector:
    const size_t packed_orbit_offset = ss->saecv_ot_offsets[iri*ss->orbit_type_count + oti] 
      + osn * ot->irbase_sizes[iri];
    // Orbit coeff vector's full size:
    const size_t orbit_fullsize = ot->size * ot->bspecn;
    const size_t particle_fullsize = ot->bspecn;
    const size_t orbit_packedsize = ot->irbase_sizes[iri];
    // This is the orbit-level matrix projecting the whole orbit onto the irrep.
    const complex double *om = ot->irbases + ot->irbase_offsets[iri];
    if (orbit_packedsize) // avoid crash on empty irrep
      cblas_zgemv(CblasRowMajor/*order*/, CblasNoTrans/*transA*/,
        orbit_packedsize/*M*/, particle_fullsize/*N*/, &one/*alpha*/, 
        om + opi*particle_fullsize/*A*/, orbit_fullsize/*lda*/,
        orig_full+fullvec_offset/*X*/, 1/*incX*/,
        &one/*beta*/, target_packed+packed_orbit_offset/*Y*/, 1/*incY*/);
    fullvec_offset += ot->bspecn;
  }
  return target_packed;
}

/// Transforms a big "packed" vector into the full basis.
/** TODO doc */
complex double *qpms_scatsys_irrep_unpack_vector(complex double *target_full,
		const complex double *orig_packed, const qpms_scatsys_t *ss,
		const qpms_iri_t iri, bool add) {
  const size_t full_len = ss->fecv_size;
  if (target_full == NULL)
    target_full = malloc(full_len*sizeof(complex double));
  if (target_full == NULL) abort();
  if (!add) memset(target_full, 0, full_len*sizeof(complex double));

  const complex double one = 1;
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) return target_full; // Completely empty irrep
 
  size_t fullvec_offset = 0;
  for(qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    const qpms_ss_oti_t oti = ss->p_orbitinfo[pi].t;
    const qpms_ss_orbit_type_t *const ot = ss->orbit_types + oti;
    const qpms_ss_osn_t osn = ss->p_orbitinfo[pi].osn;
    const qpms_ss_orbit_pi_t opi = ss->p_orbitinfo[pi].p;
    // This is where the particle's orbit starts in the "packed" vector:
    const size_t packed_orbit_offset = ss->saecv_ot_offsets[iri*ss->orbit_type_count + oti] 
      + osn * ot->irbase_sizes[iri];
    // Orbit coeff vector's full size:
    const size_t orbit_fullsize = ot->size * ot->bspecn;
    const size_t particle_fullsize = ot->bspecn;
    const size_t orbit_packedsize = ot->irbase_sizes[iri];
    // This is the orbit-level matrix projecting the whole orbit onto the irrep.
    const complex double *om = ot->irbases + ot->irbase_offsets[iri];

    if (orbit_packedsize) // empty irrep, avoid zgemv crashing.
      cblas_zgemv(CblasRowMajor/*order*/, CblasConjTrans/*transA*/,
        orbit_packedsize/*M*/, particle_fullsize/*N*/, &one/*alpha*/, om + opi*particle_fullsize/*A*/,
        orbit_fullsize/*lda*/, orig_packed+packed_orbit_offset /*X*/, 1/*incX*/, &one/*beta*/, 
        target_full+fullvec_offset/*Y*/, 1/*incY*/);

    fullvec_offset += ot->bspecn;
  }
  return target_full; 
}

complex double *qpms_scatsys_build_translation_matrix_full(
    /// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
    complex double *target,
    const qpms_scatsys_t *ss,
    double k ///< Wave number to use in the translation matrix.
    )
{
  const size_t full_len = ss->fecv_size;
  if(!target)
    QPMS_CRASHING_MALLOC(target, SQ(full_len) * sizeof(complex double));
  memset(target, 0, SQ(full_len) * sizeof(complex double)); //unnecessary?
  { // Non-diagonal part; M[piR, piC] = T[piR] S(piR<-piC)
    size_t fullvec_offsetR = 0;
    for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) {
      const qpms_vswf_set_spec_t *bspecR = ss->tm[ss->p[piR].tmatrix_id]->spec;
      const cart3_t posR = ss->p[piR].pos;
      size_t fullvec_offsetC = 0;
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) {
        const qpms_vswf_set_spec_t *bspecC = ss->tm[ss->p[piC].tmatrix_id]->spec;
        if(piC != piR)  { // The diagonal will be dealt with later.
          const cart3_t posC = ss->p[piC].pos;
          QPMS_ENSURE_SUCCESS(qpms_trans_calculator_get_trans_array_lc3p(ss->c,
                target + fullvec_offsetR*full_len + fullvec_offsetC,
                bspecR, full_len, bspecC, 1,
                k, posR, posC));      
        }
        fullvec_offsetC += bspecC->n;
      }
      assert(fullvec_offsetC = full_len);
      fullvec_offsetR += bspecR->n;
    }
    assert(fullvec_offsetR == full_len);
  }
  
  return target;
}


complex double *qpms_scatsys_build_modeproblem_matrix_full(
    /// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
    complex double *target,
    const qpms_scatsys_t *ss,
    double k ///< Wave number to use in the translation matrix.
    )
{
  const size_t full_len = ss->fecv_size;
  if(!target)
    QPMS_CRASHING_MALLOC(target, SQ(full_len) * sizeof(complex double));
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, SQ(ss->max_bspecn) * sizeof(complex double));
  memset(target, 0, SQ(full_len) * sizeof(complex double)); //unnecessary?
  complex double one = 1, zero = 0;
  { // Non-diagonal part; M[piR, piC] = T[piR] S(piR<-piC)
    size_t fullvec_offsetR = 0;
    for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) {
      const qpms_vswf_set_spec_t *bspecR = ss->tm[ss->p[piR].tmatrix_id]->spec;
      const cart3_t posR = ss->p[piR].pos;
      size_t fullvec_offsetC = 0;
      // dest particle T-matrix
      const complex double *tmmR = ss->tm[ss->p[piR].tmatrix_id]->m;
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) {
        const qpms_vswf_set_spec_t *bspecC = ss->tm[ss->p[piC].tmatrix_id]->spec;
        if(piC != piR) { // The diagonal will be dealt with later.
          const cart3_t posC = ss->p[piC].pos;
          QPMS_ENSURE_SUCCESS(qpms_trans_calculator_get_trans_array_lc3p(ss->c,
                tmp, // tmp is S(piR<-piC)
                bspecR, bspecC->n, bspecC, 1,
                k, posR, posC));      
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              bspecR->n /*m*/, bspecC->n /*n*/, bspecR->n /*k*/, 
              &one/*alpha*/, tmmR/*a*/, bspecR->n/*lda*/,
              tmp/*b*/, bspecC->n/*ldb*/, &zero/*beta*/,
              target + fullvec_offsetR*full_len + fullvec_offsetC /*c*/,
              full_len /*ldc*/);
        }
        fullvec_offsetC += bspecC->n;
      }
      fullvec_offsetR += bspecR->n;
    }
  }
  // diagonal part M[pi,pi] = -1
  for (size_t i = 0; i < full_len; ++i) target[full_len * i + i] = -1;
  
  free(tmp);
  return target;
}

complex double *qpms_scatsys_build_modeproblem_matrix_irrep_packed(
    /// Target memory with capacity for ss->saecv_sizes[iri]**2 elements. If NULL, new will be allocated.
    complex double *target_packed,
    const qpms_scatsys_t *ss, qpms_iri_t iri,
    double k ///< Wave number to use in the translation matrix.
    )
{
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) // THIS IS A BIT PROBLEMATIC, TODO how to deal with empty irreps?
    return target_packed; 
  if (target_packed == NULL)
    target_packed = malloc(SQ(packedlen)*sizeof(complex double));
  if (target_packed == NULL) abort();
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  // some of the following workspaces are probably redundant; TODO optimize later.

  // workspaces for the uncompressed particle<-particle tranlation matrix block
  // and the result of multiplying with a T-matrix
  complex double *Sblock, *TSblock;
  QPMS_CRASHING_MALLOC(Sblock, sizeof(complex double)*SQ(ss->max_bspecn));
  QPMS_CRASHING_MALLOC(TSblock, sizeof(complex double)*SQ(ss->max_bspecn));

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, sizeof(complex double) * SQ(ss->max_bspecn) * ss->sym->order);

  const complex double one = 1, zero = 0;

  for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) { //Row loop
    const qpms_ss_oti_t otiR = ss->p_orbitinfo[piR].t;
    const qpms_ss_orbit_type_t *const otR = ss->orbit_types + otiR;
    const qpms_ss_osn_t osnR = ss->p_orbitinfo[piR].osn;
    const qpms_ss_orbit_pi_t opiR = ss->p_orbitinfo[piR].p;
    // This is where the particle's orbit starts in the "packed" vector:
    const size_t packed_orbit_offsetR =
      ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiR] 
      + osnR * otR->irbase_sizes[iri];
    const qpms_vswf_set_spec_t *bspecR = ss->tm[ss->p[piR].tmatrix_id]->spec;
    // Orbit coeff vector's full size:
    const size_t orbit_fullsizeR = otR->size * otR->bspecn;
    const size_t particle_fullsizeR = otR->bspecn; // == bspecR->n
    const size_t orbit_packedsizeR = otR->irbase_sizes[iri];
    // This is the orbit-level matrix projecting the whole orbit onto the irrep.
    const complex double *omR = otR->irbases + otR->irbase_offsets[iri];
    const cart3_t posR = ss->p[piR].pos;
    if(orbit_packedsizeR) { // avoid zgemm crash on empty irrep
      // dest particle T-matrix
      const complex double *tmmR = ss->tm[ss->p[piR].tmatrix_id]->m;
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) { //Column loop
        const qpms_ss_oti_t otiC = ss->p_orbitinfo[piC].t;
        const qpms_ss_orbit_type_t *const otC = ss->orbit_types + otiC;
        const qpms_ss_osn_t osnC = ss->p_orbitinfo[piC].osn;
        const qpms_ss_orbit_pi_t opiC = ss->p_orbitinfo[piC].p;
        // This is where the particle's orbit starts in the "packed" vector:
        const size_t packed_orbit_offsetC = 
          ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiC]
          + osnC * otC->irbase_sizes[iri];
        const qpms_vswf_set_spec_t *bspecC = ss->tm[ss->p[piC].tmatrix_id]->spec;
        // Orbit coeff vector's full size:
        const size_t orbit_fullsizeC = otC->size * otC->bspecn;
        const size_t particle_fullsizeC = otC->bspecn; // == bspecC->n
        const size_t orbit_packedsizeC = otC->irbase_sizes[iri];
        // This is the orbit-level matrix projecting the whole orbit onto the irrep.
        const complex double *omC = otC->irbases + otC->irbase_offsets[iri];

        if(orbit_packedsizeC) { // avoid zgemm crash on empty irrep
          if(piC != piR) { // non-diagonal, calculate TS
            const cart3_t posC = ss->p[piC].pos;
            QPMS_ENSURE_SUCCESS(qpms_trans_calculator_get_trans_array_lc3p(ss->c,
                  Sblock, // Sblock is S(piR->piC)
                  bspecR, bspecC->n, bspecC, 1,
                  k, posR, posC));

						cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
								bspecR->n /*m*/, bspecC->n /*n*/, bspecR->n /*k*/,
								&one/*alpha*/, tmmR/*a*/, bspecR->n/*lda*/,
								Sblock/*b*/, bspecC->n/*ldb*/, &zero/*beta*/,
								TSblock /*c*/, bspecC->n /*ldc*/);
          } else { // diagonal, fill with diagonal -1
            for (size_t row = 0; row < bspecR->n; ++row)
              for (size_t col = 0; col < bspecC->n; ++col)
                TSblock[row * bspecC->n + col] = (row == col)? -1 : 0;
          }

          // tmp[oiR|piR,piC] = ∑_K M[piR,K] U*[K,piC]
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
              particle_fullsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeC /*K*/,
              &one /*alpha*/, TSblock/*A*/, particle_fullsizeC/*ldA*/, 
              omC + opiC*particle_fullsizeC /*B*/,
              orbit_fullsizeC/*ldB*/, &zero /*beta*/,
              tmp /*C*/, orbit_packedsizeC /*LDC*/);

          // target[oiR|piR,oiC|piC] += U[...] tmp[...]
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              orbit_packedsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeR /*K*/,
              &one /*alpha*/, omR + opiR*particle_fullsizeR/*A*/, orbit_fullsizeR/*ldA*/,
              tmp /*B*/, orbit_packedsizeC /*ldB*/, &one /*beta*/,
              target_packed + packedlen*packed_orbit_offsetR + packed_orbit_offsetC /*C*/,
              packedlen /*ldC*/);
        }
      }
    }
  }
  free(tmp);
  free(Sblock);
  free(TSblock);
  return target_packed;
}

complex double *qpms_scatsys_build_modeproblem_matrix_irrep_packed_orbitorderR(
    /// Target memory with capacity for ss->saecv_sizes[iri]**2 elements. If NULL, new will be allocated.
    complex double *target_packed,
    const qpms_scatsys_t *ss, qpms_iri_t iri,
    double k ///< Wave number to use in the translation matrix.
    )
{
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) // THIS IS A BIT PROBLEMATIC, TODO how to deal with empty irreps?
    return target_packed; 
  if (target_packed == NULL)
    target_packed = malloc(SQ(packedlen)*sizeof(complex double));
  if (target_packed == NULL) abort();
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  // some of the following workspaces are probably redundant; TODO optimize later.

  // workspaces for the uncompressed particle<-particle tranlation matrix block
  // and the result of multiplying with a T-matrix
  complex double *Sblock, *TSblock;
  QPMS_CRASHING_MALLOC(Sblock, sizeof(complex double)*SQ(ss->max_bspecn));
  QPMS_CRASHING_MALLOC(TSblock, sizeof(complex double)*SQ(ss->max_bspecn));

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, sizeof(complex double) * SQ(ss->max_bspecn) * ss->sym->order);

  const complex double one = 1, zero = 0;

  for(qpms_ss_pi_t opistartR = 0; opistartR < ss->p_count; 
      opistartR += ss->orbit_types[ss->p_orbitinfo[opistartR].t].size //orbit_p_countR; might write a while() instead
     ) {
    const qpms_ss_pi_t orbitstartpiR = ss->p_by_orbit[opistartR];
    const qpms_ss_oti_t otiR = ss->p_orbitinfo[orbitstartpiR].t;
    const qpms_ss_osn_t osnR = ss->p_orbitinfo[orbitstartpiR].osn;
    const qpms_ss_orbit_type_t *const otR = ss->orbit_types + otiR;
    const qpms_ss_orbit_pi_t orbit_p_countR = otR->size;
    const size_t orbit_packedsizeR = otR->irbase_sizes[iri];

    if(orbit_packedsizeR) { // avoid zgemm crash on empty irrep
      const size_t particle_fullsizeR = otR->bspecn; // == bspecR->n
      const qpms_vswf_set_spec_t *bspecR = ss->tm[ss->p[orbitstartpiR].tmatrix_id]->spec;
      // This is the orbit-level matrix projecting the whole orbit onto the irrep.
      const complex double *omR = otR->irbases + otR->irbase_offsets[iri];
      // Orbit coeff vector's full size:
      const size_t orbit_fullsizeR = otR->size * otR->bspecn;
      // This is where the orbit starts in the "packed" vector:
      const size_t packed_orbit_offsetR =
        ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiR] 
        + osnR * otR->irbase_sizes[iri];
      for(qpms_ss_orbit_pi_t opiR = 0; opiR < orbit_p_countR; ++opiR) {
        qpms_ss_pi_t piR = ss->p_by_orbit[opistartR + opiR];
        assert(opiR == ss->p_orbitinfo[piR].p);
        const qpms_ss_oti_t otiR = ss->p_orbitinfo[piR].t;
        assert(ss->p_orbitinfo[piR].osn == osnR);
        const cart3_t posR = ss->p[piR].pos;
        // dest particle T-matrix
        const complex double *tmmR = ss->tm[ss->p[piR].tmatrix_id]->m;
        for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) { //Column loop
          const qpms_ss_oti_t otiC = ss->p_orbitinfo[piC].t;
          const qpms_ss_orbit_type_t *const otC = ss->orbit_types + otiC;
          const qpms_ss_osn_t osnC = ss->p_orbitinfo[piC].osn;
          const qpms_ss_orbit_pi_t opiC = ss->p_orbitinfo[piC].p;
          // This is where the particle's orbit starts in the "packed" vector:
          const size_t packed_orbit_offsetC = 
            ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiC]
            + osnC * otC->irbase_sizes[iri];
          const qpms_vswf_set_spec_t *bspecC = ss->tm[ss->p[piC].tmatrix_id]->spec;
          // Orbit coeff vector's full size:
          const size_t orbit_fullsizeC = otC->size * otC->bspecn;
          const size_t particle_fullsizeC = otC->bspecn; // == bspecC->n
          const size_t orbit_packedsizeC = otC->irbase_sizes[iri];
          // This is the orbit-level matrix projecting the whole orbit onto the irrep.
          const complex double *omC = otC->irbases + otC->irbase_offsets[iri];

          if(orbit_packedsizeC) { // avoid zgemm crash on empty irrep
            if(piC != piR) { // non-diagonal, calculate TS
              const cart3_t posC = ss->p[piC].pos;
              QPMS_ENSURE_SUCCESS(qpms_trans_calculator_get_trans_array_lc3p(ss->c,
                    Sblock, // Sblock is S(piR->piC)
                    bspecR, bspecC->n, bspecC, 1,
                    k, posR, posC));

              cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  bspecR->n /*m*/, bspecC->n /*n*/, bspecR->n /*k*/,
                  &one/*alpha*/, tmmR/*a*/, bspecR->n/*lda*/,
                  Sblock/*b*/, bspecC->n/*ldb*/, &zero/*beta*/,
                  TSblock /*c*/, bspecC->n /*ldc*/);
            } else { // diagonal, fill with diagonal -1
              for (size_t row = 0; row < bspecR->n; ++row)
                for (size_t col = 0; col < bspecC->n; ++col)
                  TSblock[row * bspecC->n + col] = (row == col)? -1 : 0;
            }

            // tmp[oiR|piR,piC] = ∑_K M[piR,K] U*[K,piC]
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                particle_fullsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeC /*K*/,
                &one /*alpha*/, TSblock/*A*/, particle_fullsizeC/*ldA*/, 
                omC + opiC*particle_fullsizeC /*B*/,
                orbit_fullsizeC/*ldB*/, &zero /*beta*/,
                tmp /*C*/, orbit_packedsizeC /*LDC*/);

            // target[oiR|piR,oiC|piC] += U[...] tmp[...]
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                orbit_packedsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeR /*K*/,
                &one /*alpha*/, omR + opiR*particle_fullsizeR/*A*/, orbit_fullsizeR/*ldA*/,
                tmp /*B*/, orbit_packedsizeC /*ldB*/, &one /*beta*/,
                target_packed + packedlen*packed_orbit_offsetR + packed_orbit_offsetC /*C*/,
                packedlen /*ldC*/);
          }
        }
      }
    }
  }
  free(tmp);
  free(Sblock);
  free(TSblock);
  return target_packed;
}
