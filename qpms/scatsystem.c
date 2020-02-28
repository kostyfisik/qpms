#include <stdlib.h>
#define lapack_int int
#define lapack_complex_double complex double
#define lapack_complex_double_real(z) (creal(z))
#define lapack_complex_double_imag(z) (cimag(z))
#include <lapacke.h>
#include <cblas.h>
#include <lapacke.h>
#include "scatsystem.h"
#include "indexing.h"
#include "vswf.h"
#include "groups.h"
#include "symmetries.h"
#include <assert.h>
#include <unistd.h>
#include "vectors.h"
#include "quaternions.h"
#include <string.h>
#include "qpms_error.h"
#include "translations.h"
#include "tmatrices.h"
#include <pthread.h>
#include "kahansum.h"
#include "tolerances.h"
#include "beyn.h"

#ifdef QPMS_SCATSYSTEM_USE_OWN_BLAS
#include "qpmsblas.h"
#define SERIAL_ZGEMM qpms_zgemm
#else
#define SERIAL_ZGEMM cblas_zgemm
#endif

#define SQ(x) ((x)*(x))
#define QPMS_SCATSYS_LEN_RTOL 1e-13
#define QPMS_SCATSYS_TMATRIX_ATOL 1e-12
#define QPMS_SCATSYS_TMATRIX_RTOL 1e-12

long qpms_scatsystem_nthreads_default = 4;
long qpms_scatsystem_nthreads_override = 0;

void qpms_scatsystem_set_nthreads(long n) {
  qpms_scatsystem_nthreads_override = n;
}

static inline void qpms_ss_ensure_periodic(const qpms_scatsys_t *ss) {
  QPMS_ENSURE(ss->lattice_dimension > 0, "This method is applicable only to periodic systems.");
}

static inline void qpms_ss_ensure_periodic_a(const qpms_scatsys_t *ss, const char *s) {
  QPMS_ENSURE(ss->lattice_dimension > 0, "This method is applicable only to periodic systems. Use %s instead.", s);
}

static inline void qpms_ss_ensure_nonperiodic(const qpms_scatsys_t *ss) {
  QPMS_ENSURE(ss->lattice_dimension == 0, "This method is applicable only to nonperiodic systems.");
}

static inline void qpms_ss_ensure_nonperiodic_a(const qpms_scatsys_t *ss, const char *s) {
  QPMS_ENSURE(ss->lattice_dimension == 0, "This method is applicable only to nonperiodic systems. Use %s instead.", s);
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

  const qpms_vswf_set_spec_t *bspec = qpms_ss_bspec_tmi(ss, ot_current->tmatrices[0]);
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
  QPMS_ENSURE(bs_cumsum == ot_new->size * bspecn,
        "The cumulative size of the symmetry-adapted bases is wrong; "
        "expected %d = %d * %d, got %d.",
        ot_new->size * bspecn, ot_new->size, bspecn, bs_cumsum);
  ot_new->instance_count = 0;
  ss->orbit_type_count++;
}

// Standardises the lattice base vectors, and fills the other contents of ss->per[0].
// LPTODO split the functionality in smaller functions, these might be useful elsewhere.
static void process_lattice_bases(qpms_scatsys_t *ss, const qpms_tolerance_spec_t *tol) {
  switch(ss->lattice_dimension) {
    case 1: {
              double normsq;
              normsq = cart3_normsq(ss->per.lattice_basis[0]);
              ss->per.unitcell_volume = sqrt(normsq);
              ss->per.reciprocal_basis[0] = cart3_divscale(ss->per.lattice_basis[0], normsq);
              ss->per.eta = 2. / M_2_SQRTPI / ss->per.unitcell_volume;
      } break;
    case 2: {
              // (I) Create orthonormal basis of the plane in which the basis vectors lie
              cart3_t onbasis[2];
              double norm0 = cart3norm(ss->per.lattice_basis[0]);
              // 0: Just normalised 0. basis vector
              onbasis[0] = cart3_divscale(ss->per.lattice_basis[0], norm0); 
              // 1: Gram-Schmidt complement of the other basis vector
              const double b0norm_dot_b1 = cart3_dot(onbasis[0], ss->per.lattice_basis[1]);
              onbasis[1] = cart3_substract(ss->per.lattice_basis[1],
                  cart3_scale(b0norm_dot_b1, onbasis[0]));
              onbasis[1] = cart3_divscale(onbasis[1], cart3norm(onbasis[1]));
              // (II) Express the lattice basis in terms of the new 2d plane vector basis onbasis
              cart2_t b2d[2] = {{.x = norm0, .y = 0}, 
                {.x = b0norm_dot_b1, .y = cart3_dot(onbasis[1], ss->per.lattice_basis[1])}};
              // (III) Reduce lattice basis and get reciprocal bases in the 2D plane
              l2d_reduceBasis(b2d[0], b2d[1], &b2d[0], &b2d[1]);
              ss->per.unitcell_volume = l2d_unitcell_area(b2d[0], b2d[1]);
              ss->per.eta = 2. / M_2_SQRTPI / sqrt(ss->per.unitcell_volume);
              cart2_t r2d[2];
              QPMS_ENSURE_SUCCESS(l2d_reciprocalBasis1(b2d[0], b2d[1], &r2d[0], &r2d[1]));
              // (IV) Rotate everything back to the original 3D space.
              for(int i = 0; i < 2; ++i) {
                ss->per.lattice_basis[i] = cart3_add(
                    cart3_scale(b2d[i].x, onbasis[0]),
                    cart3_scale(b2d[i].y, onbasis[1]));
                ss->per.reciprocal_basis[i] = cart3_add(
                    cart3_scale(r2d[i].x, onbasis[0]),
                    cart3_scale(r2d[i].y, onbasis[1]));
              }
      } break;
    case 3: {
              l3d_reduceBasis(ss->per.lattice_basis, ss->per.lattice_basis); 
              ss->per.unitcell_volume = fabs(l3d_reciprocalBasis1(ss->per.lattice_basis,
                  ss->per.reciprocal_basis));
              ss->per.eta = 2. / M_2_SQRTPI / cbrt(ss->per.unitcell_volume);
              // TODO check unitcell_volume for sanity w.r.t tolerance
      } break;
    default:
      QPMS_WTF;
  }
}


// Almost 200 lines. This whole thing deserves a rewrite!
qpms_scatsys_at_omega_t *qpms_scatsys_apply_symmetry(const qpms_scatsys_t *orig,
    const qpms_finite_group_t *sym, complex double omega, 
    const qpms_tolerance_spec_t *tol) {
  if (sym == NULL) {
    QPMS_WARN("No symmetry group pointer provided, assuming trivial symmetry"
        " (this is currently the only option for periodic systems).");
    sym = &QPMS_FINITE_GROUP_TRIVIAL;
  }

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
  qpms_scatsys_t *ss;
  QPMS_CRASHING_MALLOC(ss, sizeof(*ss));
  ss->lattice_dimension = orig->lattice_dimension;
  // TODO basic periodic lattices related stuff here.
  if (ss->lattice_dimension) {
    ss->per = orig->per;
    process_lattice_bases(ss, tol);
  }
  for(int i = 0; i < ss->lattice_dimension; ++i)  // extend lenscale by basis vectors
    lenscale = MAX(lenscale, cart3norm(ss->per.lattice_basis[i]));
  ss->lenscale = lenscale;
  ss->sym = sym;
  ss->medium = orig->medium;

  // Copy the qpms_tmatrix_fuction_t from orig
  ss->tmg_count = orig->tmg_count;
  QPMS_CRASHING_MALLOC(ss->tmg, ss->tmg_count * sizeof(*(ss->tmg)));
  memcpy(ss->tmg, orig->tmg, ss->tmg_count * sizeof(*(ss->tmg)));

  ss->tm_capacity = sym->order * orig->tm_count;
  QPMS_CRASHING_MALLOC(ss->tm, ss->tm_capacity * sizeof(*(ss->tm)));

  ss->p_capacity = sym->order * orig->p_count;
  QPMS_CRASHING_MALLOC(ss->p, ss->p_capacity * sizeof(qpms_particle_tid_t));
  QPMS_CRASHING_MALLOC(ss->p_orbitinfo, ss->p_capacity * sizeof(qpms_ss_particle_orbitinfo_t));
  for (qpms_ss_pi_t pi = 0; pi < ss->p_capacity; ++pi) {
    ss->p_orbitinfo[pi].t = QPMS_SS_P_ORBITINFO_UNDEF;
    ss->p_orbitinfo[pi].p = QPMS_SS_P_ORBITINFO_UNDEF;
  }

  // Evaluate the original T-matrices at omega
  qpms_tmatrix_t **tm_orig_omega; 
  QPMS_CRASHING_MALLOC(tm_orig_omega, orig->tmg_count * sizeof(*tm_orig_omega));
  for(qpms_ss_tmgi_t i = 0; i < orig->tmg_count; ++i) 
    tm_orig_omega[i] = qpms_tmatrix_init_from_function(orig->tmg[i], omega);

  // Evaluate the medium and derived T-matrices at omega.
  qpms_scatsys_at_omega_t *ssw;
  QPMS_CRASHING_MALLOC(ssw, sizeof(*ssw)); // returned
  ssw->ss = ss;
  ssw->omega = omega;
  ssw->medium = qpms_epsmu_generator_eval(ss->medium, omega);
  ssw->wavenumber = qpms_wavenumber(omega, ssw->medium);
  // we will be using ss->tm_capacity also for ssw->tm
  QPMS_CRASHING_MALLOC(ssw->tm, ss->tm_capacity * sizeof(*(ssw->tm))); // returned

  // Evaluate T-matrices at omega; checking for duplicities
  
  ss->max_bspecn = 0; // We'll need it later.for memory alloc estimates.

  qpms_ss_tmi_t tm_dupl_remap[ss->tm_capacity]; // Auxilliary array to label remapping the indices after ignoring t-matrix duplicities; VLA!
  ss->tm_count = 0;
  for (qpms_ss_tmi_t i = 0; i < orig->tm_count; ++i) {
    qpms_tmatrix_t *ti = qpms_tmatrix_apply_operation(&(orig->tm[i].op), tm_orig_omega[orig->tm[i].tmgi]);
    qpms_ss_tmi_t j;
    for (j = 0; j < ss->tm_count; ++j) 
      if (qpms_tmatrix_isclose(ti, ssw->tm[j], tol->rtol, tol->atol)) {
        break;
      }
    if (j == ss->tm_count) { // duplicity not found, save both the "abstract" and "at omega" T-matrices
      qpms_tmatrix_operation_copy(&ss->tm[j].op, &orig->tm[i].op);
      ss->tm[j].tmgi = orig->tm[i].tmgi; // T-matrix functions are preserved.
      ssw->tm[j] = ti;
      ss->max_bspecn = MAX(ssw->tm[j]->spec->n, ss->max_bspecn);
      lMax = MAX(lMax, ssw->tm[j]->spec->lMax);
      ++(ss->tm_count);
    } 
    else qpms_tmatrix_free(ti);
    tm_dupl_remap[i] = j;
    if (normalisation == QPMS_NORMALISATION_UNDEF)
      normalisation = ssw->tm[i]->spec->norm;
    // We expect all bspec norms to be the same.
    else QPMS_ENSURE(normalisation == ssw->tm[j]->spec->norm,
        "Normalisation convention must be the same for all T-matrices."
        " %d != %d\n", normalisation, ssw->tm[j]->spec->norm);
  }

  // Free the original T-matrices at omega
  for(qpms_ss_tmgi_t i = 0; i < orig->tmg_count; ++i)
    qpms_tmatrix_free(tm_orig_omega[i]);
  free(tm_orig_omega);

  // Copy particles, remapping the t-matrix indices
  for (qpms_ss_pi_t i = 0; i < orig->p_count; ++(i)) {
    ss->p[i] = orig->p[i];
    ss->p[i].tmatrix_id = tm_dupl_remap[ss->p[i].tmatrix_id];
  }
  ss->p_count = orig->p_count;

  // allocate t-matrix symmetry map
  QPMS_CRASHING_MALLOC(ss->tm_sym_map, sizeof(qpms_ss_tmi_t) * sym->order * sym->order * ss->tm_count);

  // Extend the T-matrices list by the symmetry operations
  for (qpms_ss_tmi_t tmi = 0; tmi < ss->tm_count; ++tmi) 
    for (qpms_gmi_t gmi = 0; gmi < sym->order; ++gmi){
      const size_t d = ssw->tm[tmi]->spec->n;
      complex double *m;
      QPMS_CRASHING_MALLOC(m, d*d*sizeof(complex double)); // ownership passed to ss->tm[ss->tm_count].op
      qpms_irot3_uvswfi_dense(m, ssw->tm[tmi]->spec, sym->rep3d[gmi]);
      qpms_tmatrix_t *transformed = qpms_tmatrix_apply_symop(ssw->tm[tmi], m);
      qpms_ss_tmi_t tmj;
      for (tmj = 0; tmj < ss->tm_count; ++tmj)
        if (qpms_tmatrix_isclose(transformed, ssw->tm[tmj], tol->rtol, tol->atol))
          break;
      if (tmj < ss->tm_count) { // HIT, transformed T-matrix already exists
        //TODO some "rounding error cleanup" symmetrisation could be performed here?
        qpms_tmatrix_free(transformed);
        free(m);
      } else { // MISS, save the matrix (also the "abstract" one)
        ssw->tm[ss->tm_count] = transformed;
	ss->tm[ss->tm_count].tmgi = ss->tm[tmi].tmgi;
        qpms_tmatrix_operation_compose_chain_init(&(ss->tm[ss->tm_count].op), 2, 1); // FIXME MEMLEAK
        struct qpms_tmatrix_operation_compose_chain * const o = &(ss->tm[ss->tm_count].op.op.compose_chain);
        o->ops[0] = & ss->tm[tmi].op; // Let's just borrow this
        o->ops_owned[0] = false;
        o->opmem[0].typ = QPMS_TMATRIX_OPERATION_LRMATRIX;
        o->opmem[0].op.lrmatrix.m = m;
        o->opmem[0].op.lrmatrix.owns_m = true;
        o->opmem[0].op.lrmatrix.m_size = d * d;
        o->ops[1] = o->opmem;
        o->ops_owned[1] = true;
        ++(ss->tm_count);
      }
      ss->tm_sym_map[gmi + tmi * sym->order] = tmj; // In any case, tmj now indexes the correct transformed matrix
    }
  // Possibly free some space using the new ss->tm_count instead of (old) ss->tm_count*sym->order
  QPMS_CRASHING_REALLOC(ss->tm_sym_map, sizeof(qpms_ss_tmi_t) * sym->order * ss->tm_count);
  // tm could be realloc'd as well, but those are just pointers, not likely many.
 

  // allocate particle symmetry map
  QPMS_CRASHING_MALLOC(ss->p_sym_map, sizeof(qpms_ss_pi_t) * sym->order * sym->order * ss->p_count);
  // allocate orbit type array (TODO realloc later if too long)
  ss->orbit_type_count = 0;
  QPMS_CRASHING_CALLOC(ss->orbit_types, ss->p_count, sizeof(qpms_ss_orbit_type_t));

  QPMS_CRASHING_MALLOC(ss->otspace, // reallocate later
      (sizeof(qpms_ss_orbit_pi_t) * sym->order * sym->order
       + sizeof(qpms_ss_tmi_t) * sym->order
       + 3*sizeof(size_t) * sym->nirreps
       + sizeof(complex double) * SQ(sym->order * ss->max_bspecn)) * ss->p_count
       );
  ss->otspace_end = ss->otspace;
  
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
  QPMS_CRASHING_REALLOC(ss->p_sym_map, sizeof(qpms_ss_pi_t) * sym->order * ss->p_count);
  QPMS_CRASHING_REALLOC(ss->p, sizeof(qpms_particle_tid_t) * ss->p_count); 
  QPMS_CRASHING_REALLOC(ss->p_orbitinfo, sizeof(qpms_ss_particle_orbitinfo_t)*ss->p_count);
  ss->p_capacity = ss->p_count;

  {  // Reallocate the orbit type data space and update the pointers if needed.
    size_t otspace_sz = ss->otspace_end - ss->otspace;
    char *old_otspace = ss->otspace;
    QPMS_CRASHING_REALLOC(ss->otspace, otspace_sz);
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
  QPMS_CRASHING_MALLOC(ss->fecv_pstarts, ss->p_count * sizeof(size_t));
  for (qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    ss->fecv_pstarts[pi] = ss->fecv_size;
    ss->fecv_size += ssw->tm[ss->p[pi].tmatrix_id]->spec->n; // That's a lot of dereferencing!
  }

  QPMS_CRASHING_MALLOC(ss->saecv_sizes, sizeof(size_t) * sym->nirreps);
  QPMS_CRASHING_MALLOC(ss->saecv_ot_offsets, sizeof(size_t) * sym->nirreps * ss->orbit_type_count);
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
  return ssw;
}


void qpms_scatsys_free(qpms_scatsys_t *ss) {
  if(ss) {
    QPMS_ASSERT(ss->tm);
    for(qpms_ss_tmi_t tmi = 0; tmi < ss->tm_count; ++tmi) 
        qpms_tmatrix_operation_clear(&ss->tm[tmi].op); 
    free(ss->tm);
    free(ss->tmg);
    free(ss->p);
    free(ss->fecv_pstarts);
    free(ss->tm_sym_map);
    free(ss->p_sym_map);
    free(ss->otspace);
    free(ss->p_orbitinfo);
    free(ss->orbit_types);
    free(ss->saecv_ot_offsets);
    free(ss->saecv_sizes);
    free(ss->p_by_orbit);
    qpms_trans_calculator_free(ss->c);
  }
  free(ss);
}

void qpms_scatsys_at_omega_refill(qpms_scatsys_at_omega_t *ssw, 
    complex double omega) {
  const qpms_scatsys_t * const ss = ssw->ss;
  ssw->omega = omega;
  ssw->medium = qpms_epsmu_generator_eval(ss->medium, omega);
  ssw->wavenumber = qpms_wavenumber(omega, ssw->medium);
  qpms_tmatrix_t **tmatrices_preop;
  QPMS_CRASHING_CALLOC(tmatrices_preop, ss->tmg_count, sizeof(*tmatrices_preop));
  for (qpms_ss_tmgi_t tmgi = 0; tmgi < ss->tmg_count; ++tmgi) 
    tmatrices_preop[tmgi] = qpms_tmatrix_init_from_function(ss->tmg[tmgi], omega);
  for (qpms_ss_tmi_t tmi = 0; tmi < ss->tm_count; ++tmi) 
    qpms_tmatrix_apply_operation_replace(ssw->tm[tmi], &ss->tm[tmi].op,
        tmatrices_preop[ss->tm[tmi].tmgi]);
  for (qpms_ss_tmgi_t tmgi = 0; tmgi < ss->tmg_count; ++tmgi)
    qpms_tmatrix_free(tmatrices_preop[tmgi]);
  free(tmatrices_preop);
}

qpms_scatsys_at_omega_t *qpms_scatsys_at_omega(const qpms_scatsys_t *ss,
    complex double omega) {
  // TODO
  qpms_scatsys_at_omega_t *ssw;
  QPMS_CRASHING_MALLOC(ssw, sizeof(*ssw));
  ssw->omega = omega;
  ssw->ss = ss;
  ssw->medium = qpms_epsmu_generator_eval(ss->medium, omega);
  ssw->wavenumber = qpms_wavenumber(omega, ssw->medium);
  QPMS_CRASHING_CALLOC(ssw->tm, ss->tm_count, sizeof(*ssw->tm));
  qpms_tmatrix_t **tmatrices_preop;
  QPMS_CRASHING_CALLOC(tmatrices_preop, ss->tmg_count, sizeof(*tmatrices_preop));
  for (qpms_ss_tmgi_t tmgi = 0; tmgi < ss->tmg_count; ++tmgi) 
    tmatrices_preop[tmgi] = qpms_tmatrix_init_from_function(ss->tmg[tmgi], omega);
  for (qpms_ss_tmi_t tmi = 0; tmi < ss->tm_count; ++tmi) {
    ssw->tm[tmi] = qpms_tmatrix_apply_operation(&ss->tm[tmi].op,
        tmatrices_preop[ss->tm[tmi].tmgi]); //<- main difference to .._refill()
    QPMS_ENSURE(ssw->tm[tmi], 
        "Got NULL pointer from qpms_tmatrix_apply_operation");
  }
  for (qpms_ss_tmgi_t tmgi = 0; tmgi < ss->tmg_count; ++tmgi)
    qpms_tmatrix_free(tmatrices_preop[tmgi]);
  free(tmatrices_preop);
  return ssw;
}

void qpms_scatsys_at_omega_free(qpms_scatsys_at_omega_t *ssw) {
  if (ssw) {
    if(ssw->tm) 
     for(qpms_ss_tmi_t i = 0; i < ssw->ss->tm_count; ++i) 
      qpms_tmatrix_free(ssw->tm[i]); 
    free(ssw->tm);
  }
  free(ssw);
}

// (copypasta from symmetries.c)
// TODO at some point, maybe support also other norms.
// (perhaps just use qpms_normalisation_t_factor() at the right places)
static inline void check_norm_compat(const qpms_vswf_set_spec_t *s)
{
  switch (s->norm & QPMS_NORMALISATION_NORM_BITS) {
    case QPMS_NORMALISATION_NORM_POWER:
      break;
    case QPMS_NORMALISATION_NORM_SPHARM:
      break;
    default:
      QPMS_WTF; // Only SPHARM and POWER norms are supported right now.
  }
}

complex double *qpms_orbit_action_matrix(complex double *target,
    const qpms_ss_orbit_type_t *ot, const qpms_vswf_set_spec_t *bspec,
    const qpms_finite_group_t *sym, const qpms_gmi_t g) {
  assert(sym); assert(g < sym->order);
  assert(sym->rep3d);
  assert(ot); assert(ot->size > 0);
  // check_norm_compat(bspec); not needed here, the qpms_irot3_uvswfi_dense should complain if necessary
  const size_t n = bspec->n;
  const qpms_gmi_t N = ot->size;
  if (target == NULL) 
    QPMS_CRASHING_MALLOC(target, n*n*N*N*sizeof(complex double));
  memset(target, 0, n*n*N*N*sizeof(complex double));
  complex double tmp[n][n]; // this is the 'per-particle' action
  qpms_irot3_uvswfi_dense(tmp[0], bspec, sym->rep3d[g]); 
  for(qpms_ss_orbit_pi_t Col = 0; Col < ot->size; ++Col) {
    // Row is the 'destination' of the symmetry operation, Col is the 'source'
    const qpms_ss_orbit_pi_t Row = ot->action[sym->order * Col + g];
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
  // check_norm_compat(bspec); // probably not needed here, let the called functions complain if necessary, but CHEKME
  const size_t n = bspec->n;
  const qpms_gmi_t N = ot->size;
  if (target == NULL)
    QPMS_CRASHING_MALLOC(target, n*n*N*N*sizeof(complex double));
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
  check_norm_compat(bspec); // Here I'm not sure; CHECKME
  const size_t n = bspec->n;
  const qpms_gmi_t N = ot->size;
  const bool newtarget = (target == NULL);
  if (newtarget)
    QPMS_CRASHING_MALLOC(target,n*n*N*N*sizeof(complex double));
  memset(target, 0, n*n*N*N*sizeof(complex double));

  // Get the projector (also workspace for right sg. vect.)
  complex double *projector = qpms_orbit_irrep_projector_matrix(NULL,
    ot, bspec, sym, iri);
  if(!projector) abort();
  // Workspace for the right singular vectors.
  complex double *V_H; QPMS_CRASHING_MALLOC(V_H, n*n*N*N*sizeof(complex double));
  // THIS SHOULD NOT BE NECESSARY
  complex double *U; QPMS_CRASHING_MALLOC(U, n*n*N*N*sizeof(complex double));
  double *s; QPMS_CRASHING_MALLOC(s, n*N*sizeof(double)); 

  QPMS_ENSURE_SUCCESS(LAPACKE_zgesdd(LAPACK_ROW_MAJOR,
      'A', // jobz; overwrite projector with left sg.vec. and write right into V_H
      n*N /* m */, n*N /* n */, projector /* a */, n*N /* lda */,
      s /* s */, U /* u */, n*N /* ldu, irrelev. */, V_H /* vt */,
      n*N /* ldvt */));

  size_t bs;
  for(bs = 0; bs < n*N; ++bs) {
    QPMS_ENSURE(s[bs] <= 1 + SVD_ATOL, "%zd. SV too large: %.16lf", bs, s[bs]);
    QPMS_ENSURE(!(s[bs] > SVD_ATOL && fabs(1-s[bs]) > SVD_ATOL),
        "%zd. SV in the 'wrong' interval: %.16lf", bs, s[bs]);
    if(s[bs] < SVD_ATOL) break;
  }

  memcpy(target, V_H, bs*N*n*sizeof(complex double));
  if(newtarget) QPMS_CRASHING_REALLOC(target, bs*N*n*sizeof(complex double));
  if(basis_size) *basis_size = bs;

  free(s);
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
    QPMS_CRASHING_MALLOC(target_packed, SQ(packedlen)*sizeof(complex double));
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
    QPMS_CRASHING_MALLOC(target_full, SQ(full_len)*sizeof(complex double));
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
    QPMS_CRASHING_MALLOC(target_packed, SQ(packedlen)*sizeof(complex double));
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, sizeof(complex double) * SQ(ss->max_bspecn) * ss->sym->order);

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
    QPMS_CRASHING_MALLOC(target_full, SQ(full_len)*sizeof(complex double));
  if(!add) memset(target_full, 0, SQ(full_len)*sizeof(complex double));

  if(!packedlen) return target_full; // Empty irrep, do nothing.

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, sizeof(complex double) * SQ(ss->max_bspecn) * ss->sym->order);

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
    QPMS_CRASHING_MALLOC(target_packed, packedlen*sizeof(complex double));
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
    QPMS_CRASHING_MALLOC(target_full, full_len*sizeof(complex double));
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
    complex double k ///< Wave number to use in the translation matrix.
    )
{
  return qpms_scatsys_build_translation_matrix_e_full(
      target, ss, k, QPMS_HANKEL_PLUS);
}

complex double *qpms_scatsyswk_build_translation_matrix_full(
    /// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
    complex double *target,
    const qpms_scatsys_at_omega_k_t *sswk
    )
{
  const qpms_scatsys_at_omega_t *ssw = sswk->ssw;
  const complex double wavenumber = ssw->wavenumber;
  const qpms_scatsys_t *ss = ssw->ss;
  qpms_ss_ensure_periodic(ss);
  const cart3_t k_cart3 = cart3_from_double_array(sswk->k);
  return qpms_scatsys_periodic_build_translation_matrix_full(target, ss, wavenumber, &k_cart3);
}

complex double *qpms_scatsys_build_translation_matrix_e_full(
    /// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
    complex double *target,
    const qpms_scatsys_t *ss,
    complex double k, ///< Wave number to use in the translation matrix.
    qpms_bessel_t J ///< Bessel function type.
    )
{
  qpms_ss_ensure_nonperiodic(ss);
  const size_t full_len = ss->fecv_size;
  if(!target)
    QPMS_CRASHING_MALLOC(target, SQ(full_len) * sizeof(complex double));
  memset(target, 0, SQ(full_len) * sizeof(complex double)); //unnecessary?
  { // Non-diagonal part; M[piR, piC] = T[piR] S(piR<-piC)
    size_t fullvec_offsetR = 0;
    for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) {
      const qpms_vswf_set_spec_t *bspecR = qpms_ss_bspec_pi(ss, piR); 
      const cart3_t posR = ss->p[piR].pos;
      size_t fullvec_offsetC = 0;
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) {
        const qpms_vswf_set_spec_t *bspecC = qpms_ss_bspec_pi(ss, piC);
        if(piC != piR)  { // The diagonal will be dealt with later.
          const cart3_t posC = ss->p[piC].pos;
          QPMS_ENSURE_SUCCESS(qpms_trans_calculator_get_trans_array_lc3p(ss->c,
                target + fullvec_offsetR*full_len + fullvec_offsetC,
                bspecR, full_len, bspecC, 1,
                k, posR, posC, J));      
        }
        fullvec_offsetC += bspecC->n;
      }
      assert(fullvec_offsetC == full_len);
      fullvec_offsetR += bspecR->n;
    }
    assert(fullvec_offsetR == full_len);
  }
  
  return target;
}

static inline int qpms_ss_ppair_W32xy(const qpms_scatsys_t *ss,
    qpms_ss_pi_t pdest, qpms_ss_pi_t psrc, complex double wavenumber, const cart2_t kvector,
    complex double *target, const ptrdiff_t deststride, const ptrdiff_t srcstride,
    qpms_ewald_part parts) {
  const qpms_vswf_set_spec_t *srcspec = qpms_ss_bspec_pi(ss, psrc);
  const qpms_vswf_set_spec_t *destspec = qpms_ss_bspec_pi(ss, pdest);

  // This might be a bit arbitrary, roughly "copied" from Unitcell constructor. TODO review.
  const double maxR = sqrt(ss->per.unitcell_volume) * 32;
  const double maxK = 1024 * 2 * M_PI / maxR;

  return qpms_trans_calculator_get_trans_array_e32_e(ss->c,
      target, NULL /*err*/, destspec, deststride, srcspec, srcstride,
      ss->per.eta, wavenumber, 
      cart3xy2cart2(ss->per.lattice_basis[0]), cart3xy2cart2(ss->per.lattice_basis[1]),
      kvector, 
      cart2_substract(cart3xy2cart2(ss->p[pdest].pos), cart3xy2cart2(ss->p[psrc].pos)),
      maxR, maxK, parts);
}

static inline int qpms_ss_ppair_W(const qpms_scatsys_t *ss,
    qpms_ss_pi_t pdest, qpms_ss_pi_t psrc, complex double wavenumber, const double wavevector[],
    complex double *target, const ptrdiff_t deststride, const ptrdiff_t srcstride,
    qpms_ewald_part parts) {
  if(ss->lattice_dimension == 2 && // Currently, we can only the xy-plane
      !ss->per.lattice_basis[0].z && !ss->per.lattice_basis[1].z &&
      !wavevector[2]) 
    return qpms_ss_ppair_W32xy(ss, pdest, psrc, wavenumber, cart2_from_double_array(wavevector), 
              target + deststride * ss->fecv_pstarts[pdest] + srcstride * ss->fecv_pstarts[psrc],
              deststride, srcstride, parts);
  else 
    QPMS_NOT_IMPLEMENTED("Only 2D xy-lattices currently supported");
}

complex double *qpms_scatsys_periodic_build_translation_matrix_full(
    complex double *target, const qpms_scatsys_t *ss,
    complex double wavenumber, const cart3_t *wavevector) {
  QPMS_UNTESTED;
  qpms_ss_ensure_periodic(ss);
  const size_t full_len = ss->fecv_size;
  if(!target)
    QPMS_CRASHING_MALLOC(target, SQ(full_len) * sizeof(complex double));
  const ptrdiff_t deststride = ss->fecv_size, srcstride = 1;
  // We have some limitations in the current implementation
  if(ss->lattice_dimension == 2 && // Currently, we can only the xy-plane
      !ss->per.lattice_basis[0].z && !ss->per.lattice_basis[1].z &&
      !wavevector->z) {
    for (qpms_ss_pi_t pd = 0; pd < ss->p_count; ++pd) 
      for (qpms_ss_pi_t ps = 0; ps < ss->p_count; ++ps) {
        QPMS_ENSURE_SUCCESS(qpms_ss_ppair_W32xy(ss, pd, ps, wavenumber, cart3xy2cart2(*wavevector), 
              target + deststride * ss->fecv_pstarts[pd] + srcstride * ss->fecv_pstarts[ps],
              deststride, srcstride, QPMS_EWALD_FULL));
      }
  } else 
    QPMS_NOT_IMPLEMENTED("Only 2D xy-lattices currently supported");
  return target;
}

// Common implementation of qpms_scatsysw[k]_build_modeproblem_matrix_full
static inline complex double *qpms_scatsysw_scatsyswk_build_modeproblem_matrix_full(
    /// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
    complex double *target,
    const qpms_scatsys_at_omega_t *ssw,
    const double k[] // NULL if non-periodic
    )
{
  const complex double wavenumber = ssw->wavenumber;
  const qpms_scatsys_t *ss = ssw->ss;
  qpms_ss_ensure_nonperiodic(ss);
  const size_t full_len = ss->fecv_size;
  if(!target)
    QPMS_CRASHING_MALLOC(target, SQ(full_len) * sizeof(complex double));
  complex double *tmp; // translation matrix, S or W
  QPMS_CRASHING_MALLOC(tmp, SQ(ss->max_bspecn) * sizeof(complex double));
  memset(target, 0, SQ(full_len) * sizeof(complex double));
  const complex double zero = 0, minusone = -1;
  { // Non-diagonal part; M[piR, piC] = -T[piR] S(piR<-piC)
    size_t fullvec_offsetR = 0;
    for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) {
      const qpms_vswf_set_spec_t *bspecR = ssw->tm[ss->p[piR].tmatrix_id]->spec;
      const cart3_t posR = ss->p[piR].pos;
      size_t fullvec_offsetC = 0;
      // dest particle T-matrix
      const complex double *tmmR = ssw->tm[ss->p[piR].tmatrix_id]->m;
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) {
        const qpms_vswf_set_spec_t *bspecC = ssw->tm[ss->p[piC].tmatrix_id]->spec;
        if (k == NULL) { // non-periodic case
          if(piC != piR) { // No "self-interaction" in non-periodic case
            const cart3_t posC = ss->p[piC].pos;
            QPMS_ENSURE_SUCCESS(qpms_trans_calculator_get_trans_array_lc3p(ss->c,
                  tmp, // tmp is S(piR<-piC)
                  bspecR, bspecC->n, bspecC, 1,
                  wavenumber, posR, posC, QPMS_HANKEL_PLUS));      
          }
        } else { // periodic case 
          QPMS_ENSURE_SUCCESS(qpms_ss_ppair_W(ss, piR, piC, wavenumber, k,
                tmp /*target*/, bspecC->n /*deststride*/, 1 /*srcstride*/,
                QPMS_EWALD_FULL));
        }

        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            bspecR->n /*m*/, bspecC->n /*n*/, bspecR->n /*k*/, 
            &minusone/*alpha*/, tmmR/*a*/, bspecR->n/*lda*/,
            tmp/*b*/, bspecC->n/*ldb*/, &zero/*beta*/,
            target + fullvec_offsetR*full_len + fullvec_offsetC /*c*/,
            full_len /*ldc*/);
        
        fullvec_offsetC += bspecC->n;
      }
      fullvec_offsetR += bspecR->n;
    }
  }

  // Add the identity, diagonal part M[pi,pi] += 1
  for (size_t i = 0; i < full_len; ++i) target[full_len * i + i] += 1;

  free(tmp);
  return target;
}

complex double *qpms_scatsysw_build_modeproblem_matrix_full(
    complex double *target, const qpms_scatsys_at_omega_t *ssw) {
  qpms_ss_ensure_nonperiodic_a(ssw->ss, "qpms_scatsyswk_build_modeproblem_matrix_full()");
  return qpms_scatsysw_scatsyswk_build_modeproblem_matrix_full(
      target, ssw, NULL);
}

complex double *qpms_scatsyswk_build_modeproblem_matrix_full(
    complex double *target, const qpms_scatsys_at_omega_k_t *sswk) 
{
  qpms_ss_ensure_periodic_a(sswk->ssw->ss, "qpms_scatsysw_build_modeproblem_matrix_full()");
  return qpms_scatsysw_scatsyswk_build_modeproblem_matrix_full(target, sswk->ssw, sswk->k);
}


// Serial reference implementation.
complex double *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_serial(
    /// Target memory with capacity for ss->saecv_sizes[iri]**2 elements. If NULL, new will be allocated.
    complex double *target_packed,
    const qpms_scatsys_at_omega_t *ssw,
    qpms_iri_t iri
    )
{
  const qpms_scatsys_t *ss = ssw->ss;
  qpms_ss_ensure_nonperiodic(ss);
  const complex double k = ssw->wavenumber;
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) // THIS IS A BIT PROBLEMATIC, TODO how to deal with empty irreps?
    return target_packed; 
  if (target_packed == NULL)
    QPMS_CRASHING_MALLOC(target_packed, SQ(packedlen)*sizeof(complex double));
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  // some of the following workspaces are probably redundant; TODO optimize later.

  // workspaces for the uncompressed particle<-particle tranlation matrix block
  // and the result of multiplying with a T-matrix (times -1)
  complex double *Sblock, *TSblock;
  QPMS_CRASHING_MALLOC(Sblock, sizeof(complex double)*SQ(ss->max_bspecn));
  QPMS_CRASHING_MALLOC(TSblock, sizeof(complex double)*SQ(ss->max_bspecn));

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, sizeof(complex double) * SQ(ss->max_bspecn) * ss->sym->order);

  const complex double one = 1, zero = 0, minusone = -1;

  for(qpms_ss_pi_t piR = 0; piR < ss->p_count; ++piR) { //Row loop
    const qpms_ss_oti_t otiR = ss->p_orbitinfo[piR].t;
    const qpms_ss_orbit_type_t *const otR = ss->orbit_types + otiR;
    const qpms_ss_osn_t osnR = ss->p_orbitinfo[piR].osn;
    const qpms_ss_orbit_pi_t opiR = ss->p_orbitinfo[piR].p;
    // This is where the particle's orbit starts in the "packed" vector:
    const size_t packed_orbit_offsetR =
      ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiR] 
      + osnR * otR->irbase_sizes[iri];
    const qpms_vswf_set_spec_t *bspecR = ssw->tm[ss->p[piR].tmatrix_id]->spec;
    // Orbit coeff vector's full size:
    const size_t orbit_fullsizeR = otR->size * otR->bspecn;
    const size_t particle_fullsizeR = otR->bspecn; // == bspecR->n
    const size_t orbit_packedsizeR = otR->irbase_sizes[iri];
    // This is the orbit-level matrix projecting the whole orbit onto the irrep.
    const complex double *omR = otR->irbases + otR->irbase_offsets[iri];
    const cart3_t posR = ss->p[piR].pos;
    if(orbit_packedsizeR) { // avoid zgemm crash on empty irrep
      // dest particle T-matrix
      const complex double *tmmR = ssw->tm[ss->p[piR].tmatrix_id]->m;
      for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) { //Column loop
        const qpms_ss_oti_t otiC = ss->p_orbitinfo[piC].t;
        const qpms_ss_orbit_type_t *const otC = ss->orbit_types + otiC;
        const qpms_ss_osn_t osnC = ss->p_orbitinfo[piC].osn;
        const qpms_ss_orbit_pi_t opiC = ss->p_orbitinfo[piC].p;
        // This is where the particle's orbit starts in the "packed" vector:
        const size_t packed_orbit_offsetC = 
          ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiC]
          + osnC * otC->irbase_sizes[iri];
        const qpms_vswf_set_spec_t *bspecC = ssw->tm[ss->p[piC].tmatrix_id]->spec;
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
                  k, posR, posC, QPMS_HANKEL_PLUS));

						cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
								bspecR->n /*m*/, bspecC->n /*n*/, bspecR->n /*k*/,
								&minusone/*alpha*/, tmmR/*a*/, bspecR->n/*lda*/,
								Sblock/*b*/, bspecC->n/*ldb*/, &zero/*beta*/,
								TSblock /*c*/, bspecC->n /*ldc*/);
          } else { // diagonal, fill with diagonal +1
            for (size_t row = 0; row < bspecR->n; ++row)
              for (size_t col = 0; col < bspecC->n; ++col)
                TSblock[row * bspecC->n + col] = (row == col)? +1 : 0;
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

complex double *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_orbitorderR(
    /// Target memory with capacity for ss->saecv_sizes[iri]**2 elements. If NULL, new will be allocated.
    complex double *target_packed,
    const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri
    )
{
  const qpms_scatsys_t *ss = ssw->ss;
  qpms_ss_ensure_nonperiodic(ss);
  const complex double k = ssw->wavenumber;
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) // THIS IS A BIT PROBLEMATIC, TODO how to deal with empty irreps?
    return target_packed; 
  if (target_packed == NULL)
    QPMS_CRASHING_MALLOC(target_packed, SQ(packedlen)*sizeof(complex double));
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  // some of the following workspaces are probably redundant; TODO optimize later.

  // workspaces for the uncompressed particle<-particle tranlation matrix block
  // and the result of multiplying with a T-matrix (times -1)
  complex double *Sblock, *TSblock;
  QPMS_CRASHING_MALLOC(Sblock, sizeof(complex double)*SQ(ss->max_bspecn));
  QPMS_CRASHING_MALLOC(TSblock, sizeof(complex double)*SQ(ss->max_bspecn));

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, sizeof(complex double) * SQ(ss->max_bspecn) * ss->sym->order);

  const complex double one = 1, zero = 0, minusone = -1;

  for(qpms_ss_pi_t opistartR = 0; opistartR < ss->p_count; 
      opistartR += ss->orbit_types[ss->p_orbitinfo[ss->p_by_orbit[opistartR]].t].size //orbit_p_countR; might write a while() instead
     ) {
    const qpms_ss_pi_t orbitstartpiR = ss->p_by_orbit[opistartR];
    const qpms_ss_oti_t otiR = ss->p_orbitinfo[orbitstartpiR].t;
    const qpms_ss_osn_t osnR = ss->p_orbitinfo[orbitstartpiR].osn;
    const qpms_ss_orbit_type_t *const otR = ss->orbit_types + otiR;
    const qpms_ss_orbit_pi_t orbit_p_countR = otR->size;
    const size_t orbit_packedsizeR = otR->irbase_sizes[iri];

    if(orbit_packedsizeR) { // avoid zgemm crash on empty irrep
      const size_t particle_fullsizeR = otR->bspecn; // == bspecR->n
      const qpms_vswf_set_spec_t *bspecR = ssw->tm[ss->p[orbitstartpiR].tmatrix_id]->spec;
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
        assert(otiR == ss->p_orbitinfo[piR].t);
        assert(ss->p_orbitinfo[piR].osn == osnR);
        const cart3_t posR = ss->p[piR].pos;
        // dest particle T-matrix
        const complex double *tmmR = ssw->tm[ss->p[piR].tmatrix_id]->m;
        for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) { //Column loop
          const qpms_ss_oti_t otiC = ss->p_orbitinfo[piC].t;
          const qpms_ss_orbit_type_t *const otC = ss->orbit_types + otiC;
          const qpms_ss_osn_t osnC = ss->p_orbitinfo[piC].osn;
          const qpms_ss_orbit_pi_t opiC = ss->p_orbitinfo[piC].p;
          // This is where the particle's orbit starts in the "packed" vector:
          const size_t packed_orbit_offsetC = 
            ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiC]
            + osnC * otC->irbase_sizes[iri];
          const qpms_vswf_set_spec_t *bspecC = ssw->tm[ss->p[piC].tmatrix_id]->spec;
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
                    k, posR, posC, QPMS_HANKEL_PLUS));

              cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  bspecR->n /*m*/, bspecC->n /*n*/, bspecR->n /*k*/,
                  &minusone/*alpha*/, tmmR/*a*/, bspecR->n/*lda*/,
                  Sblock/*b*/, bspecC->n/*ldb*/, &zero/*beta*/,
                  TSblock /*c*/, bspecC->n /*ldc*/);
            } else { // diagonal, fill with diagonal +1
              for (size_t row = 0; row < bspecR->n; ++row)
                for (size_t col = 0; col < bspecC->n; ++col)
                  TSblock[row * bspecC->n + col] = (row == col)? +1 : 0;
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

struct qpms_scatsysw_build_modeproblem_matrix_irrep_packed_parallelR_thread_arg{
  const qpms_scatsys_at_omega_t *ssw;
  qpms_ss_pi_t *opistartR_ptr;
  pthread_mutex_t *opistartR_mutex;
  qpms_iri_t iri;
  complex double *target_packed;
};

static void *qpms_scatsysw_build_modeproblem_matrix_irrep_packed_parallelR_thread(void *arg)
{
  const struct qpms_scatsysw_build_modeproblem_matrix_irrep_packed_parallelR_thread_arg 
    *a = arg;
  const qpms_scatsys_at_omega_t *ssw = a->ssw;
  const complex double k = ssw->wavenumber;
  const qpms_scatsys_t *ss = ssw->ss;
  const qpms_iri_t iri = a->iri;
  const size_t packedlen = ss->saecv_sizes[iri];

  // some of the following workspaces are probably redundant; TODO optimize later.

  // workspaces for the uncompressed particle<-particle tranlation matrix block
  // and the result of multiplying with a T-matrix (times -1)
  complex double *Sblock, *TSblock;
  QPMS_CRASHING_MALLOC(Sblock, sizeof(complex double)*SQ(ss->max_bspecn));
  QPMS_CRASHING_MALLOC(TSblock, sizeof(complex double)*SQ(ss->max_bspecn));

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, sizeof(complex double) * SQ(ss->max_bspecn) * ss->sym->order);

  const complex double one = 1, zero = 0, minusone = -1;

  while(1) {
    // In the beginning, pick a target (row) orbit for this thread
    QPMS_ENSURE_SUCCESS(pthread_mutex_lock(a->opistartR_mutex));
    if(*(a->opistartR_ptr) >= ss->p_count) {// Everything is already done, end
      QPMS_ENSURE_SUCCESS(pthread_mutex_unlock(a->opistartR_mutex));
      break; 
    }
    const qpms_ss_pi_t opistartR = *(a->opistartR_ptr);
    // Now increment it for another thread:
    *(a->opistartR_ptr) += ss->orbit_types[ss->p_orbitinfo[ss->p_by_orbit[opistartR]].t].size;
    QPMS_ENSURE_SUCCESS(pthread_mutex_unlock(a->opistartR_mutex));
    
    // Orbit picked (defined by opistartR), do the work.
    const qpms_ss_pi_t orbitstartpiR = ss->p_by_orbit[opistartR];
    const qpms_ss_oti_t otiR = ss->p_orbitinfo[orbitstartpiR].t;
    const qpms_ss_osn_t osnR = ss->p_orbitinfo[orbitstartpiR].osn;
    const qpms_ss_orbit_type_t *const otR = ss->orbit_types + otiR;
    const qpms_ss_orbit_pi_t orbit_p_countR = otR->size;
    const size_t orbit_packedsizeR = otR->irbase_sizes[iri];

    if(orbit_packedsizeR) { // avoid zgemm crash on empty irrep
      const size_t particle_fullsizeR = otR->bspecn; // == bspecR->n
      const qpms_vswf_set_spec_t *bspecR = ssw->tm[ss->p[orbitstartpiR].tmatrix_id]->spec;
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
        assert(otiR == ss->p_orbitinfo[piR].t);
        assert(ss->p_orbitinfo[piR].osn == osnR);
        const cart3_t posR = ss->p[piR].pos;
        // dest particle T-matrix
        const complex double *tmmR = ssw->tm[ss->p[piR].tmatrix_id]->m;
        for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) { //Column loop
          const qpms_ss_oti_t otiC = ss->p_orbitinfo[piC].t;
          const qpms_ss_orbit_type_t *const otC = ss->orbit_types + otiC;
          const qpms_ss_osn_t osnC = ss->p_orbitinfo[piC].osn;
          const qpms_ss_orbit_pi_t opiC = ss->p_orbitinfo[piC].p;
          // This is where the particle's orbit starts in the "packed" vector:
          const size_t packed_orbit_offsetC = 
            ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiC]
            + osnC * otC->irbase_sizes[iri];
          const qpms_vswf_set_spec_t *bspecC = ssw->tm[ss->p[piC].tmatrix_id]->spec;
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
                    k, posR, posC, QPMS_HANKEL_PLUS));

              SERIAL_ZGEMM(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  bspecR->n /*m*/, bspecC->n /*n*/, bspecR->n /*k*/,
                  &minusone/*alpha*/, tmmR/*a*/, bspecR->n/*lda*/,
                  Sblock/*b*/, bspecC->n/*ldb*/, &zero/*beta*/,
                  TSblock /*c*/, bspecC->n /*ldc*/);
            } else { // diagonal, fill with diagonal +1
              for (size_t row = 0; row < bspecR->n; ++row)
                for (size_t col = 0; col < bspecC->n; ++col)
                  TSblock[row * bspecC->n + col] = (row == col)? +1 : 0;
            }

            // tmp[oiR|piR,piC] = ∑_K M[piR,K] U*[K,piC]
            SERIAL_ZGEMM(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                particle_fullsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeC /*K*/,
                &one /*alpha*/, TSblock/*A*/, particle_fullsizeC/*ldA*/, 
                omC + opiC*particle_fullsizeC /*B*/,
                orbit_fullsizeC/*ldB*/, &zero /*beta*/,
                tmp /*C*/, orbit_packedsizeC /*LDC*/);

            // target[oiR|piR,oiC|piC] += U[...] tmp[...]
            SERIAL_ZGEMM(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                orbit_packedsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeR /*K*/,
                &one /*alpha*/, omR + opiR*particle_fullsizeR/*A*/, orbit_fullsizeR/*ldA*/,
                tmp /*B*/, orbit_packedsizeC /*ldB*/, &one /*beta*/,
                a->target_packed + packedlen*packed_orbit_offsetR + packed_orbit_offsetC /*C*/,
                packedlen /*ldC*/);
          }
        }
      }
    }



  }
  free(tmp);
  free(Sblock);
  free(TSblock);
  return NULL;
}

// this differs from the ...build_modeproblem_matrix... only by the `J`
// maybe I should use this one there as well to save lines... TODO
struct qpms_scatsys_build_translation_matrix_e_irrep_packed_parallelR_thread_arg{
  const qpms_scatsys_t *ss;
  qpms_ss_pi_t *opistartR_ptr;
  pthread_mutex_t *opistartR_mutex;
  qpms_iri_t iri;
  complex double *target_packed;
  complex double k;
  qpms_bessel_t J;
};

static void *qpms_scatsys_build_translation_matrix_e_irrep_packed_parallelR_thread(void *arg)
{
  const struct qpms_scatsys_build_translation_matrix_e_irrep_packed_parallelR_thread_arg 
    *a = arg;
  const qpms_scatsys_t *ss = a->ss;
  const qpms_iri_t iri = a->iri;
  const size_t packedlen = ss->saecv_sizes[iri];
  const qpms_bessel_t J = a->J;

  // some of the following workspaces are probably redundant; TODO optimize later.

  // workspace for the uncompressed particle<-particle tranlation matrix block
  complex double *Sblock;
  QPMS_CRASHING_MALLOC(Sblock, sizeof(complex double)*SQ(ss->max_bspecn));

  // Workspace for the intermediate particle-orbit matrix result
  complex double *tmp;
  QPMS_CRASHING_MALLOC(tmp, sizeof(complex double) * SQ(ss->max_bspecn) * ss->sym->order);

  const complex double one = 1, zero = 0;

  while(1) {
    // In the beginning, pick a target (row) orbit for this thread
    QPMS_ENSURE_SUCCESS(pthread_mutex_lock(a->opistartR_mutex));
    if(*(a->opistartR_ptr) >= ss->p_count) {// Everything is already done, end
      QPMS_ENSURE_SUCCESS(pthread_mutex_unlock(a->opistartR_mutex));
      break; 
    }
    const qpms_ss_pi_t opistartR = *(a->opistartR_ptr);
    // Now increment it for another thread:
    *(a->opistartR_ptr) += ss->orbit_types[ss->p_orbitinfo[ss->p_by_orbit[opistartR]].t].size;
    QPMS_ENSURE_SUCCESS(pthread_mutex_unlock(a->opistartR_mutex));
    
    // Orbit picked (defined by opistartR), do the work.
    const qpms_ss_pi_t orbitstartpiR = ss->p_by_orbit[opistartR];
    const qpms_ss_oti_t otiR = ss->p_orbitinfo[orbitstartpiR].t;
    const qpms_ss_osn_t osnR = ss->p_orbitinfo[orbitstartpiR].osn;
    const qpms_ss_orbit_type_t *const otR = ss->orbit_types + otiR;
    const qpms_ss_orbit_pi_t orbit_p_countR = otR->size;
    const size_t orbit_packedsizeR = otR->irbase_sizes[iri];

    if(orbit_packedsizeR) { // avoid zgemm crash on empty irrep
      const size_t particle_fullsizeR = otR->bspecn; // == bspecR->n
      const qpms_vswf_set_spec_t *bspecR = qpms_ss_bspec_pi(ss, orbitstartpiR);
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
        assert(otiR == ss->p_orbitinfo[piR].t);
        assert(ss->p_orbitinfo[piR].osn == osnR);
        const cart3_t posR = ss->p[piR].pos;
        for(qpms_ss_pi_t piC = 0; piC < ss->p_count; ++piC) { //Column loop
          const qpms_ss_oti_t otiC = ss->p_orbitinfo[piC].t;
          const qpms_ss_orbit_type_t *const otC = ss->orbit_types + otiC;
          const qpms_ss_osn_t osnC = ss->p_orbitinfo[piC].osn;
          const qpms_ss_orbit_pi_t opiC = ss->p_orbitinfo[piC].p;
          // This is where the particle's orbit starts in the "packed" vector:
          const size_t packed_orbit_offsetC = 
            ss->saecv_ot_offsets[iri*ss->orbit_type_count + otiC]
            + osnC * otC->irbase_sizes[iri];
          const qpms_vswf_set_spec_t *bspecC = qpms_ss_bspec_pi(ss, piC);
          // Orbit coeff vector's full size:
          const size_t orbit_fullsizeC = otC->size * otC->bspecn;
          const size_t particle_fullsizeC = otC->bspecn; // == bspecC->n
          const size_t orbit_packedsizeC = otC->irbase_sizes[iri];
          // This is the orbit-level matrix projecting the whole orbit onto the irrep.
          const complex double *omC = otC->irbases + otC->irbase_offsets[iri];

          if(orbit_packedsizeC) { // avoid zgemm crash on empty irrep
					// THIS IS THE MAIN PART DIFFERENT FROM ...modeproblem...() TODO unify 
          // somehow to save lines
            if(piC != piR) { // non-diagonal, calculate S
              const cart3_t posC = ss->p[piC].pos;
              QPMS_ENSURE_SUCCESS(qpms_trans_calculator_get_trans_array_lc3p(ss->c,
                    Sblock, // Sblock is S(piR->piC)
                    bspecR, bspecC->n, bspecC, 1,
                    a->k, posR, posC, J));
            } else { // diagonal, fill with zeros; TODO does this make sense?
              // would unit matrix be better? or unit only for QPMS_BESSEL_REGULAR?
              for (size_t row = 0; row < bspecR->n; ++row)
                for (size_t col = 0; col < bspecC->n; ++col)
                  Sblock[row * bspecC->n + col] = 0; //(row == col)? 1 : 0;
            }

            // tmp[oiR|piR,piC] = ∑_K M[piR,K] U*[K,piC]
            SERIAL_ZGEMM(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                particle_fullsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeC /*K*/,
                &one /*alpha*/, Sblock/*A*/, particle_fullsizeC/*ldA*/, 
                omC + opiC*particle_fullsizeC /*B*/,
                orbit_fullsizeC/*ldB*/, &zero /*beta*/,
                tmp /*C*/, orbit_packedsizeC /*LDC*/);

            // target[oiR|piR,oiC|piC] += U[...] tmp[...]
            SERIAL_ZGEMM(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                orbit_packedsizeR /*M*/, orbit_packedsizeC /*N*/, particle_fullsizeR /*K*/,
                &one /*alpha*/, omR + opiR*particle_fullsizeR/*A*/, orbit_fullsizeR/*ldA*/,
                tmp /*B*/, orbit_packedsizeC /*ldB*/, &one /*beta*/,
                a->target_packed + packedlen*packed_orbit_offsetR + packed_orbit_offsetC /*C*/,
                packedlen /*ldC*/);
          }
        }
      }
    }



  }
  free(tmp);
  free(Sblock);
  return NULL;
}

// Almost the same as ...build_modeproblem_matrix_...parallelR
// --> TODO write this in a more generic way to save LoC.
complex double *qpms_scatsys_build_translation_matrix_e_irrep_packed(
		/// Target memory with capacity for ss->fecv_size**2 elements. If NULL, new will be allocated.
		complex double *target_packed,
		const qpms_scatsys_t *ss,
		qpms_iri_t iri,
    const complex double k,
		qpms_bessel_t J
		)
{
  qpms_ss_ensure_nonperiodic(ss);
  QPMS_UNTESTED;
  const size_t packedlen = ss->saecv_sizes[iri];
  if (!packedlen) // THIS IS A BIT PROBLEMATIC, TODO how to deal with empty irreps?
    return target_packed; 
  if (target_packed == NULL)
    QPMS_CRASHING_MALLOC(target_packed, SQ(packedlen)*sizeof(complex double));
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  qpms_ss_pi_t opistartR = 0;
  pthread_mutex_t opistartR_mutex;
  QPMS_ENSURE_SUCCESS(pthread_mutex_init(&opistartR_mutex, NULL));
  const struct qpms_scatsys_build_translation_matrix_e_irrep_packed_parallelR_thread_arg
    arg = {ss, &opistartR, &opistartR_mutex, iri, target_packed, k, J};

  // FIXME THIS IS NOT PORTABLE:
  long nthreads;
  if (qpms_scatsystem_nthreads_override > 0) {
    nthreads = qpms_scatsystem_nthreads_override;
    QPMS_DEBUG(QPMS_DBGMSG_THREADS, "Using overriding value of %ld thread(s).", 
        nthreads);
  } else {
    nthreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (nthreads < 1) {
      QPMS_DEBUG(QPMS_DBGMSG_THREADS, "_SC_NPROCESSORS_ONLN returned %ld, using %ld thread(s) instead.",
         nthreads, qpms_scatsystem_nthreads_default);
      nthreads = qpms_scatsystem_nthreads_default; 
    } else {
      QPMS_DEBUG(QPMS_DBGMSG_THREADS, "_SC_NRPOCESSORS_ONLN returned %ld.", nthreads);
    }
  }
  pthread_t thread_ids[nthreads];
  for(long thi = 0; thi < nthreads; ++thi)
    QPMS_ENSURE_SUCCESS(pthread_create(thread_ids + thi, NULL,
      qpms_scatsys_build_translation_matrix_e_irrep_packed_parallelR_thread,
      (void *) &arg));
  for(long thi = 0; thi < nthreads; ++thi) {
    void *retval;
    QPMS_ENSURE_SUCCESS(pthread_join(thread_ids[thi], &retval));
  }

  QPMS_ENSURE_SUCCESS(pthread_mutex_destroy(&opistartR_mutex));
  return target_packed;
}


// Parallel implementation, now default
complex double *qpms_scatsysw_build_modeproblem_matrix_irrep_packed(
    /// Target memory with capacity for ss->saecv_sizes[iri]**2 elements. If NULL, new will be allocated.
    complex double *target_packed,
    const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri
    )
{
  qpms_ss_ensure_nonperiodic(ssw->ss);
  const size_t packedlen = ssw->ss->saecv_sizes[iri];
  if (!packedlen) // THIS IS A BIT PROBLEMATIC, TODO how to deal with empty irreps?
    return target_packed; 
  if (target_packed == NULL)
    QPMS_CRASHING_MALLOC(target_packed,SQ(packedlen)*sizeof(complex double));
  memset(target_packed, 0, SQ(packedlen)*sizeof(complex double));

  qpms_ss_pi_t opistartR = 0;
  pthread_mutex_t opistartR_mutex;
  QPMS_ENSURE_SUCCESS(pthread_mutex_init(&opistartR_mutex, NULL));
  const struct qpms_scatsysw_build_modeproblem_matrix_irrep_packed_parallelR_thread_arg
    arg = {ssw, &opistartR, &opistartR_mutex, iri, target_packed};

  // FIXME THIS IS NOT PORTABLE:
  long nthreads;
  if (qpms_scatsystem_nthreads_override > 0) {
    nthreads = qpms_scatsystem_nthreads_override;
    QPMS_DEBUG(QPMS_DBGMSG_THREADS, "Using overriding value of %ld thread(s).", 
        nthreads);
  } else {
    nthreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (nthreads < 1) {
      QPMS_DEBUG(QPMS_DBGMSG_THREADS, "_SC_NPROCESSORS_ONLN returned %ld, using %ld thread(s) instead.",
         nthreads, qpms_scatsystem_nthreads_default);
      nthreads = qpms_scatsystem_nthreads_default; 
    } else {
      QPMS_DEBUG(QPMS_DBGMSG_THREADS, "_SC_NRPOCESSORS_ONLN returned %ld.", nthreads);
    }
  }
  pthread_t thread_ids[nthreads];
  for(long thi = 0; thi < nthreads; ++thi)
    QPMS_ENSURE_SUCCESS(pthread_create(thread_ids + thi, NULL,
      qpms_scatsysw_build_modeproblem_matrix_irrep_packed_parallelR_thread,
      (void *) &arg));
  for(long thi = 0; thi < nthreads; ++thi) {
    void *retval;
    QPMS_ENSURE_SUCCESS(pthread_join(thread_ids[thi], &retval));
  }

  QPMS_ENSURE_SUCCESS(pthread_mutex_destroy(&opistartR_mutex));
  return target_packed;
}


complex double *qpms_scatsys_incident_field_vector_full(
    complex double *target_full, const qpms_scatsys_t *ss,
    qpms_incfield_t f,	const void *args, bool add ) {
  QPMS_UNTESTED;
  if (!target_full) QPMS_CRASHING_CALLOC(target_full, ss->fecv_size,
      sizeof(complex double));
  for(qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    complex double *ptarget = target_full + ss->fecv_pstarts[pi];
    const qpms_vswf_set_spec_t *bspec = qpms_ss_bspec_pi(ss, pi);
    const cart3_t pos = ss->p[pi].pos;
    QPMS_ENSURE_SUCCESS(f(ptarget, bspec, pos, args, add));
  }
  return target_full;
}


#if 0
complex double *qpms_scatsys_incident_field_vector_irrep_packed(
    complex double *target_full, const qpms_scatsys_t *ss,
    const qpms_iri_t iri, qpms_incfield_t f,
    const void *args, bool add) {
  TODO;
}
#endif


complex double *qpms_scatsysw_apply_Tmatrices_full(
		complex double *target_full, const complex double *inc_full, 
		const qpms_scatsys_at_omega_t *ssw) {
  QPMS_UNTESTED;
  const qpms_scatsys_t *ss = ssw->ss;
  if (!target_full) QPMS_CRASHING_CALLOC(target_full, ss->fecv_size,
      sizeof(complex double));
  for(qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    complex double *ptarget = target_full + ss->fecv_pstarts[pi];
    const complex double *psrc = inc_full + ss->fecv_pstarts[pi];
    // TODO check whether T-matrix is non-virtual after virtual t-matrices are implemented.
    const qpms_tmatrix_t *T = ssw->tm[ss->p[pi].tmatrix_id];
    qpms_apply_tmatrix(ptarget, psrc, T);
  }
  return target_full;
}



ccart3_t qpms_scatsys_eval_E(const qpms_scatsys_t *ss, 
    const complex double *cvf, const cart3_t where,
    const complex double k) {
  QPMS_UNTESTED;
  ccart3_t res = {0,0,0};
  ccart3_t res_kc = {0,0,0}; // kahan sum compensation

  for (qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    const qpms_vswf_set_spec_t *bspec = qpms_ss_bspec_pi(ss, pi);
    const cart3_t particle_pos = ss->p[pi].pos;
    const complex double *particle_cv = cvf + ss->fecv_pstarts[pi];

    const csph_t kr = sph_cscale(k, cart2sph(
          cart3_substract(where, particle_pos)));
    const csphvec_t E_sph = qpms_eval_uvswf(bspec, particle_cv, kr, 
        QPMS_HANKEL_PLUS);
    const ccart3_t E_cart = csphvec2ccart_csph(E_sph, kr);
    ckahanadd(&(res.x), &(res_kc.x), E_cart.x);
    ckahanadd(&(res.y), &(res_kc.y), E_cart.y);
    ckahanadd(&(res.z), &(res_kc.z), E_cart.z);
  }
  return res;
}

#if 0
ccart3_t qpms_scatsys_eval_E_irrep(const qpms_scatsys_t *ss,
   qpms_iri_t iri, const complex double *cvr, cart3_t where) {
  TODO;
}
#endif

void qpms_ss_LU_free(qpms_ss_LU lu) {
  free(lu.a);
  free(lu.ipiv);
}

qpms_ss_LU qpms_scatsysw_modeproblem_matrix_full_factorise(complex double *mpmatrix_full,
    int *target_piv, const qpms_scatsys_at_omega_t *ssw, const qpms_scatsys_at_omega_k_t *sswk) {
  if (sswk) {
    QPMS_ASSERT(sswk->ssw == ssw || !ssw);
    ssw = sswk->ssw;
    QPMS_ASSERT(ssw->ss->lattice_dimension > 0);
  } else {
    QPMS_ASSERT(ssw->ss->lattice_dimension == 0);
  }
  const qpms_scatsys_t *ss = ssw->ss;
  QPMS_ENSURE(mpmatrix_full, "A non-NULL pointer to the pre-calculated mode matrix is required");
  if (!target_piv) QPMS_CRASHING_MALLOC(target_piv, ss->fecv_size * sizeof(int));
  QPMS_ENSURE_SUCCESS(LAPACKE_zgetrf(LAPACK_ROW_MAJOR, ss->fecv_size, ss->fecv_size,
        mpmatrix_full, ss->fecv_size, target_piv));
  qpms_ss_LU lu;
  lu.a = mpmatrix_full;
  lu.ipiv = target_piv;
  lu.ssw = ssw;
  lu.sswk = sswk;
  lu.full = true;
  lu.iri = -1;
  return lu;
}

qpms_ss_LU qpms_scatsysw_modeproblem_matrix_irrep_packed_factorise(complex double *mpmatrix_packed,
    int *target_piv, const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri) {
  QPMS_ENSURE(mpmatrix_packed, "A non-NULL pointer to the pre-calculated mode matrix is required");
  qpms_ss_ensure_nonperiodic(ssw->ss);
  size_t n = ssw->ss->saecv_sizes[iri];
  if (!target_piv) QPMS_CRASHING_MALLOC(target_piv, n * sizeof(int));
  QPMS_ENSURE_SUCCESS(LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n,
        mpmatrix_packed, n, target_piv));
  qpms_ss_LU lu;
  lu.a = mpmatrix_packed;
  lu.ipiv = target_piv;
  lu.ssw = ssw;
  lu.full = false;
  lu.iri = iri;
  return lu;
}

qpms_ss_LU qpms_scatsysw_build_modeproblem_matrix_full_LU(
    complex double *target, int *target_piv,
    const qpms_scatsys_at_omega_t *ssw){
  qpms_ss_ensure_nonperiodic_a(ssw->ss, "qpms_scatsyswk_build_modeproblem_matrix_full_LU()"); 
  target = qpms_scatsysw_build_modeproblem_matrix_full(target, ssw);
  return qpms_scatsysw_modeproblem_matrix_full_factorise(target, target_piv, ssw, NULL);
}

qpms_ss_LU qpms_scatsyswk_build_modeproblem_matrix_full_LU(
    complex double *target, int *target_piv,
    const qpms_scatsys_at_omega_k_t *sswk){
  target = qpms_scatsyswk_build_modeproblem_matrix_full(target, sswk);
  return qpms_scatsysw_modeproblem_matrix_full_factorise(target, target_piv, sswk->ssw, sswk);
}

qpms_ss_LU qpms_scatsysw_build_modeproblem_matrix_irrep_packed_LU(
    complex double *target, int *target_piv,
    const qpms_scatsys_at_omega_t *ssw, qpms_iri_t iri){
  target = qpms_scatsysw_build_modeproblem_matrix_irrep_packed(target, ssw, iri);
  return qpms_scatsysw_modeproblem_matrix_irrep_packed_factorise(target, target_piv, ssw, iri);
}

complex double *qpms_scatsys_scatter_solve(
    complex double *f, const complex double *a_inc, qpms_ss_LU lu) {
  const size_t n = lu.full ? lu.ssw->ss->fecv_size : lu.ssw->ss->saecv_sizes[lu.iri];
  if (!f) QPMS_CRASHING_MALLOC(f, n * sizeof(complex double));
  memcpy(f, a_inc, n*sizeof(complex double)); // It will be rewritten by zgetrs
  QPMS_ENSURE_SUCCESS(LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N' /*trans*/,  n /*n*/, 1 /*nrhs number of right hand sides*/,
        lu.a /*a*/, n /*lda*/, lu.ipiv /*ipiv*/, f/*b*/, 1 /*ldb; CHECKME*/));
  return f;
}

struct qpms_scatsys_finite_eval_Beyn_ImTS_param {
  const qpms_scatsys_t *ss;
  qpms_iri_t iri;
};

/// Wrapper for Beyn algorithm (non-periodic system)
static int qpms_scatsys_finite_eval_Beyn_ImTS(complex double *target,
    size_t m, complex double omega, void *params) {
  const struct qpms_scatsys_finite_eval_Beyn_ImTS_param *p = params;
  qpms_scatsys_at_omega_t *ssw = qpms_scatsys_at_omega(p->ss, omega);
  QPMS_ENSURE(ssw != NULL, "qpms_scatsys_at_omega() returned NULL");
  if (p->iri == QPMS_NO_IRREP) {
    QPMS_ASSERT(m == p->ss->fecv_size);
    QPMS_ENSURE(NULL != qpms_scatsysw_build_modeproblem_matrix_full(
          target, ssw),
      "qpms_scatsysw_build_modeproblem_matrix_full() returned NULL");
  } else {
    QPMS_ASSERT(m == p->ss->saecv_sizes[p->iri]);
    QPMS_ENSURE(NULL != qpms_scatsysw_build_modeproblem_matrix_irrep_packed(
          target, ssw, p->iri),
      "qpms_scatsysw_build_modeproblem_matrix_irrep_packed() returned NULL");
  }
  qpms_scatsys_at_omega_free(ssw);
  return QPMS_SUCCESS;
}

beyn_result_t *qpms_scatsys_finite_find_eigenmodes(
    const qpms_scatsys_t * const ss, const qpms_iri_t iri,
    complex double omega_centre, double omega_rr, double omega_ri, 
    size_t contour_npoints, 
    double rank_tol, size_t rank_sel_min, double res_tol) {
  qpms_ss_ensure_nonperiodic_a(ss, "qpms_scatsys_periodic_find_eigenmodes()");
  size_t n; // matrix dimension
  if (qpms_iri_is_valid(ss->sym, iri)) {
    n = ss->saecv_sizes[iri];
  } else if (iri == QPMS_NO_IRREP) {
    n = ss->fecv_size;
  } else QPMS_WTF;

  beyn_contour_t *contour = beyn_contour_ellipse(omega_centre,
      omega_rr, omega_ri, contour_npoints);
  
  struct qpms_scatsys_finite_eval_Beyn_ImTS_param p = {ss, iri};
  beyn_result_t *result = beyn_solve(n, n /* possibly make smaller? */,
      qpms_scatsys_finite_eval_Beyn_ImTS, NULL, (void *) &p,
      contour, rank_tol, rank_sel_min, res_tol);
  QPMS_ENSURE(result != NULL, "beyn_solve() returned NULL");

  free(contour);
  return result;
}

struct qpms_scatsys_periodic_eval_Beyn_ImTW_param {
  const qpms_scatsys_t *ss;
  const double *k; ///< Wavevector in cartesian coordinates.
};

/// Wrapper for Beyn algorithm (periodic system)
static int qpms_scatsys_periodic_eval_Beyn_ImTW(complex double *target,
    size_t m, complex double omega, void *params){
  const struct qpms_scatsys_periodic_eval_Beyn_ImTW_param *p = params;
  qpms_scatsys_at_omega_t *ssw = qpms_scatsys_at_omega(p->ss, omega);
  QPMS_ENSURE(ssw != NULL, "qpms_scatsys_at_omega() returned NULL");
  qpms_scatsys_at_omega_k_t sswk = {
    .ssw = ssw,
    .k = {p->k[0], p->k[1], p->k[2]}
  };
  QPMS_ASSERT(m == p->ss->fecv_size);
  QPMS_ENSURE(NULL != 
      qpms_scatsyswk_build_modeproblem_matrix_full(target, &sswk),
      "qpms_scatsyswk_build_modeproblem_matrix_full() returned NULL");
  qpms_scatsys_at_omega_free(ssw);
  return QPMS_SUCCESS;
}

beyn_result_t *qpms_scatsys_periodic_find_eigenmodes(
    const qpms_scatsys_t * const ss, const double k[3], 
    complex double omega_centre, double omega_rr, double omega_ri, 
    size_t contour_npoints, 
    double rank_tol, size_t rank_sel_min, double res_tol) {
  qpms_ss_ensure_nonperiodic_a(ss, "qpms_scatsys_finite_find_eigenmodes()");
  size_t n = ss->fecv_size; // matrix dimension

  beyn_contour_t *contour = beyn_contour_ellipse(omega_centre,
      omega_rr, omega_ri, contour_npoints);

  struct qpms_scatsys_periodic_eval_Beyn_ImTW_param p = {ss, k};
  beyn_result_t *result = beyn_solve(n, n, 
      qpms_scatsys_periodic_eval_Beyn_ImTW, NULL, (void *) &p,
      contour, rank_tol, rank_sel_min, res_tol);
  QPMS_ENSURE(result != NULL, "beyn_solve() returned NULL");

  free(contour);
  return result;
}

