#include <stdlib.h>
#include <cblas.h>
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
  const size_t actionsiz = sizeof(ot_current->action[0]) * ot_current->size 
    * ss->sym->order;
  ot_new->action = (void *) (ss->otspace_end);
  memcpy(ot_new->action, ot_current->action, actionsiz);
  // N.B. we copied mostly garbage ^^^, most of it is initialized just now:
  extend_orbit_action(ss, ot_new);
  ss->otspace_end += actionsiz;
  const size_t tmsiz = sizeof(ot_current->tmatrices[0]) * ot_current->size;
  ot_new->tmatrices = (void *) (ss->otspace_end);
  memcpy(ot_new->tmatrices, ot_current->tmatrices, tmsiz);
  ss->otspace_end += tmsiz;
}


// Almost 200 lines. This whole thing deserves a rewrite!
qpms_scatsys_t *qpms_scatsys_apply_symmetry(const qpms_scatsys_t *orig, const qpms_finite_group_t *sym) {
  // TODO check data sanity

  // First, determine the rough radius of the array
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
    lenscale = ((maxx-minx)+(maxy-miny)+(maxz-minz)) / 3;
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

  qpms_ss_tmi_t tm_dupl_remap[ss->tm_capacity]; // Auxilliary array to label remapping the indices after ignoring t-matrix duplicities
  ss->tm_count = 0;
  for (qpms_ss_tmi_t i = 0; i < orig->tm_count; ++i) {
    qpms_ss_tmi_t j;
    for (j = 0; j < ss->tm_count; ++j) 
      if (qpms_tmatrix_isclose(orig->tm[i], ss->tm[j], QPMS_SCATSYS_TMATRIX_RTOL, QPMS_SCATSYS_TMATRIX_ATOL)) {
        break;
      }
    if (j == ss->tm_count) { // duplicity not found, copy the t-matrix
      ss->tm[i] = qpms_tmatrix_copy(orig->tm[i]);
      ++(ss->tm_count);
    } 
    tm_dupl_remap[i] = j;
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
       + sizeof(qpms_ss_tmi_t) * sym->order) * ss->p_count);
  
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
    }

    for (qpms_gmi_t gmi = 0; gmi < sym->order; ++gmi) {
      cart3_t newpoint = qpms_irot3_apply_cart3(sym->rep3d[gmi], ss->p[pi].pos);
      qpms_ss_tmi_t new_tmi = ss->tm_sym_map[gmi + ss->p[pi].tmatrix_id * sym->order]; // transformed t-matrix index
      qpms_ss_pi_t pj;
      for (pj = 0; pj < ss->p_count; ++pj)
        if (cart3_isclose(newpoint, ss->p[pj].pos, 0, ss->lenscale * QPMS_SCATSYS_LEN_RTOL)) {
          if (new_tmi != ss->p[pj].tmatrix_id)
            abort(); // The particle is mapped to another (or itself) without having the required symmetry.  TODO give some error message.
          break;
        }
      if (pj < ss->p_count) { // HIT, the particle is transformed to an "existing" one.
        ;
      } else { // MISS, the symmetry transforms the particle to a new location, so add it.
        qpms_particle_tid_t newparticle = {newpoint, new_tmi};
        ss->p[ss->p_count] = newparticle;
        ++(ss->p_count);
      }
      ss->p_sym_map[gmi + pi * sym->order] = pi;
      
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
      qpms_ss_oti_t oti;
      for(oti = 0; oti < ss->orbit_type_count; ++oti)
        if (orbit_types_equal(&ot_current, &(ss->orbit_types[oti]))) break; // HIT, orbit type already exists
      if (oti == ss->orbit_type_count)  // MISS, add the new orbit type
        add_orbit_type(ss, &ot_current);
      
      // Walk through the orbit again and set up the orbit info of the particles
      for (qpms_ss_orbit_pi_t opi = 0; opi < ot_current.size; ++opi) {
        const qpms_ss_pi_t pi_opi = current_orbit[opi];
        ss->p_orbitinfo[pi_opi].t = oti;
        ss->p_orbitinfo[pi_opi].p = opi;
      }
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
    if(shift)
      for (size_t oi = 0; oi < ss->orbit_type_count; ++oi) {
        ss->orbit_types[oi].action = (void *)(((char *) (ss->orbit_types[oi].action)) + shift);
        ss->orbit_types[oi].tmatrices = (void *)(((char *) (ss->orbit_types[oi].tmatrices)) + shift);
      }
  }

  // Set ss->fecv_size and ss->fecv_pstarts
  ss->fecv_size = 0;
  ss->fecv_pstarts = malloc(ss->p_count * sizeof(size_t));
  for (qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    ss->fecv_pstarts[pi] = ss->fecv_size;
    ss->fecv_size += ss->tm[ss->p[pi].tmatrix_id]->spec->n; // That's a lot of dereferencing!
  }

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
        target[n*n*N*Row + n*Col + n*N*row + col] = tmp[row][col]; //CHECKCONJ
  }
  return target;
}

