//  c99 -g -I.. ss_syms_packs.c staticgroups.c ../qpms/scatsystem.c ../qpms/vswf.c ../qpms/error.c  ../qpms/translations.c ../qpms/symmetries.c ../qpms/legendre.c ../qpms/gaunt.c  ../qpms/wigner.c -lm -lgsl -lblas -llapacke
typedef int qpms_gmi_t;// There is something wrong in the includes, apparently.
#include <qpms/qpms_types.h>
#include <qpms/scatsystem.h>
#include <stdlib.h>
#include <qpms/vswf.h>
#include <qpms/indexing.h>
#include <stdio.h>
#include "staticgroups.h"

const qpms_finite_group_t *D3h = &QPMS_FINITE_GROUP_D3h;
const qpms_finite_group_t *C4v = &QPMS_FINITE_GROUP_C4v;
const qpms_finite_group_t *TRIVG = &QPMS_FINITE_GROUP_trivial_g;
const qpms_finite_group_t *C2v = &QPMS_FINITE_GROUP_C2v;
const qpms_finite_group_t *D2h = &QPMS_FINITE_GROUP_D2h;
const qpms_finite_group_t *D4h = &QPMS_FINITE_GROUP_D4h;

double uniform_random(double min, double max) {
  double random_value = min + (max-min)*(double)rand()/RAND_MAX;
  return random_value;
}

int main()
{
  srand(666);
#if 0
  qpms_vswf_set_spec_t 
    *b1 = qpms_vswf_set_spec_from_lMax(1,QPMS_NORMALISATION_POWER_CS),
    *b2 = qpms_vswf_set_spec_from_lMax(2,QPMS_NORMALISATION_POWER_CS);
#else
  // Only electric waves
  qpms_vswf_set_spec_t *b1 = qpms_vswf_set_spec_init(), 
                       *b2 = qpms_vswf_set_spec_init();
  b1->norm = b2-> norm = QPMS_NORMALISATION_POWER_CS;
  for(qpms_l_t l = 1; l <= 1; ++l)
    for (qpms_m_t m = -l; m <= l; ++m)
      qpms_vswf_set_spec_append(b1, qpms_tmn2uvswfi(QPMS_VSWF_ELECTRIC, m, l));
  for(qpms_l_t l = 1; l <= 1; ++l)
    for (qpms_m_t m = -l; m <= l; ++m)
      qpms_vswf_set_spec_append(b2, qpms_tmn2uvswfi(QPMS_VSWF_ELECTRIC, m, l));
#endif
  qpms_tmatrix_t *t1 = qpms_tmatrix_init(b1);
  qpms_tmatrix_t *t2 = qpms_tmatrix_init(b2);

#if 0
  // Random diagonal T-matrices
  for(size_t i = 0; i < b1->n; ++i)
    t1->m[i + i*b1->n] = uniform_random(-1,1) + I*uniform_random(-1,1);
  for(size_t i = 0; i < b2->n; ++i)
    t2->m[i + i*b2->n] = uniform_random(-1,1) + I*uniform_random(-1,1);
#else
  for(size_t i = 0; i < b1->n; ++i)
    t1->m[i + i*b1->n] = 1;
  for(size_t i = 0; i < b2->n; ++i)
    t2->m[i + i*b2->n] = 1;
#endif

  const cart3_t pp1 = {0, 0, 1}, pp2 = {0,0, 2}, pp3 = {0,0 , 0};
  qpms_tmatrix_t * tmlist[] = {t1, t2};
  qpms_particle_tid_t plist[] = {{pp1, 0}, {pp2, 0}, {pp3, 1}};

  qpms_scatsys_t protoss;
  protoss.tm = tmlist;
  protoss.tm_count=2;
  protoss.p = plist;
  protoss.p_count=3;

  qpms_scatsys_t *ss = qpms_scatsys_apply_symmetry(&protoss, D3h);

  printf("p_count: %d, tm_count: %d, nirreps: %d, orbit_type_count: %d\n",
      (int)ss->p_count, (int)ss->tm_count, (int)ss->sym->nirreps,
      (int)ss->orbit_type_count);

  const double k = 1.7;

  complex double *S_full = qpms_scatsys_build_translation_matrix_full(
      NULL, ss, k); 
  complex double *S_packed[ss->sym->nirreps];
  for (qpms_iri_t iri = 0; iri < ss->sym->nirreps; ++iri)
    S_packed[iri] = qpms_scatsys_irrep_pack_matrix(NULL,
        S_full, ss, iri);

  complex double *S_recfull = qpms_scatsys_irrep_unpack_matrix(NULL,
      S_packed[0], ss, 0, false);
  for (qpms_iri_t iri = 1; iri < ss->sym->nirreps; ++iri)
    qpms_scatsys_irrep_unpack_matrix(S_recfull, S_packed[iri],
        ss, iri, true);

  double maxerr = 0;
  for (size_t i = 0; i < ss->fecv_size; ++i) {
    double err = cabs(S_full[i] - S_recfull[i]);
    maxerr = (err > maxerr) ? err : maxerr;
  }

  printf("maxerr: %lg\n", maxerr);
  
  for (qpms_iri_t iri = 0; iri < ss->sym->nirreps; ++iri) free(S_packed[iri]);
  free(S_full);
  qpms_scatsys_free(ss);
  qpms_tmatrix_free(t1);
  qpms_tmatrix_free(t2);
  qpms_vswf_set_spec_free(b1);
  qpms_vswf_set_spec_free(b2);
  return 0;
}
