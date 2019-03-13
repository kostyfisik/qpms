// c99 -g -DNLINE -DDAGRUP=C4v -DDUMP_PARTICLE_POSITIONS -DDUMP_ORBIT_ACTION -DDUMP_PROJECTORMATRIX -DDUMP_ACTIONMATRIX -I.. sss3.c staticgroups.c ../qpms/scatsystem.c ../qpms/vswf.c ../qpms/error.c  ../qpms/translations.c ../qpms/symmetries.c ../qpms/legendre.c ../qpms/gaunt.c  ../qpms/wigner.c -lm -lgsl  -llapacke ~/repo/CBLAS/lib/cblas_LINUX.a ~/repo/BLAS-3.8.0/blas_LINUX.a
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
const qpms_finite_group_t *C2 = &QPMS_FINITE_GROUP_C2;
const qpms_finite_group_t *C4 = &QPMS_FINITE_GROUP_C4;
const qpms_finite_group_t *D2h = &QPMS_FINITE_GROUP_D2h;
const qpms_finite_group_t *D4h = &QPMS_FINITE_GROUP_D4h;
const qpms_finite_group_t *x_and_z_flip = &QPMS_FINITE_GROUP_x_and_z_flip;
const qpms_finite_group_t *y_and_z_flip = &QPMS_FINITE_GROUP_y_and_z_flip;

#ifndef DAGRUP
#define DAGRUP D4h
#endif

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
    for (qpms_m_t m = -0l; m <= l; m += 2)
      qpms_vswf_set_spec_append(b1, qpms_tmn2uvswfi(QPMS_VSWF_ELECTRIC, m, l));
  for(qpms_l_t l = 1; l <= 1; ++l)
    for (qpms_m_t m = -0l; m <= l; m += 2)
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

#ifdef YLINE
  const cart3_t pp1 = {0, 1.1, 0};
  const cart3_t pp2 = {0, 1.4, 0};
#elif defined XLINE
  const cart3_t pp1 = {1.1, 0, 0};
  const cart3_t pp2 = {1.4, 0, 0};
#elif defined ZLINE
  const cart3_t pp1 = {0, 0, 1.1};
  const cart3_t pp2 = {0, 0, 1.4};
#else
  const cart3_t pp1 = {1.1, 1, 0};
  const cart3_t pp2 = {0, 1.4, 0};
#endif
  const cart3_t pp3 = {0, 0, 1};
  qpms_tmatrix_t * tmlist[] = {t1, t2};
  qpms_particle_tid_t plist[] = {{pp1,1}, {pp2, 0},  {pp3, 1},
  };

  qpms_scatsys_t protoss;
  protoss.tm = tmlist;
  protoss.tm_count=sizeof(tmlist)/sizeof(qpms_tmatrix_t *);
  protoss.p = plist;
  protoss.p_count=sizeof(plist)/sizeof(qpms_particle_tid_t);

  qpms_scatsys_t *ss = qpms_scatsys_apply_symmetry(&protoss, DAGRUP);

  printf("p_count: %d, tm_count: %d, nirreps: %d, orbit_type_count: %d\n",
      (int)ss->p_count, (int)ss->tm_count, (int)ss->sym->nirreps,
      (int)ss->orbit_type_count);

  fputs("Orbit projection matrices:\n", stderr);
  for (qpms_ss_oti_t oti = 0; oti < ss->orbit_type_count; ++oti) {
    fprintf(stderr, "Orbit type %d:\n", (int)oti);
    const qpms_ss_orbit_type_t *ot = &(ss->orbit_types[oti]);
    size_t row = 0;
    for (qpms_iri_t iri = 0; iri < ss->sym->nirreps; ++iri) {
      assert(row*ot->size*ot->bspecn == ot->irbase_offsets[iri]);
      fprintf(stderr, "---------- IR %d (%s) -------------\n", (int)iri, ss->sym->irreps[iri].name);
      for (size_t irbi = 0; irbi < ot->irbase_sizes[iri]; ++irbi) {
        for(qpms_ss_orbit_pi_t opi = 0; opi < ot->size; ++opi){
          fputs("| ", stderr);
          for(size_t i = 0; i < ot->bspecn; ++i) {
            const complex double elem = ot->irbases[row * ot->size * ot->bspecn
            + opi * ot->bspecn + i];
            fprintf(stderr, "%+.3f%+.3fj ", creal(elem), cimag(elem));
          }
        }
        fputs("|\n", stderr);
        ++row;
      }
    }
    fputs("------------------------\n\n",stderr);
  }

  const double k = 1.7;

  complex double *S_full = qpms_scatsys_build_translation_matrix_full(
      NULL, ss, k); 
  {
    const size_t full_len = ss->fecv_size;
    size_t fullvec_offset_dest = 0;
    for (qpms_ss_pi_t pdest = 0; pdest < ss->p_count; pdest++) {
      size_t fullvec_offset_src = 0;
      const size_t bspecn_dest = ss->tm[ss->p[pdest].tmatrix_id]->spec->n;
      for (qpms_ss_pi_t psrc = 0; psrc < ss->p_count; psrc++) {
        const size_t bspecn_src = ss->tm[ss->p[psrc].tmatrix_id]->spec->n;
        fprintf(stderr, "Translation matrix element %d<-%d; (%g %g %g)<-(%g %g %g):\n",
            (int)pdest, (int)psrc, ss->p[pdest].pos.x, ss->p[pdest].pos.y, ss->p[pdest].pos.z,
            ss->p[psrc].pos.x, ss->p[psrc].pos.y, ss->p[psrc].pos.z);

        for(size_t row = 0; row < bspecn_dest; ++row) {
          for(size_t col = 0; col < bspecn_src; ++col)
            fprintf(stderr, "%+2.3f%+2.3fj ", creal(S_full[full_len * (fullvec_offset_dest+row) + fullvec_offset_src+col]),
                cimag(S_full[full_len * (fullvec_offset_dest+row) + fullvec_offset_src+col]));
          fputc('\n', stderr);
        }
        fullvec_offset_src += bspecn_src;
      }
      fullvec_offset_dest += bspecn_dest;
    }
  }
  {
    fputs("\n\n", stderr);
    const size_t full_len = ss->fecv_size;
    for (size_t row = 0 ; row < full_len; ++row) {
      for (size_t col = 0 ; col < full_len; ++col)
        fprintf(stderr, "%+2.3f%+2.3fj ", creal(S_full[full_len * row + col]), cimag(S_full[full_len * row + col]));
      fputc('\n', stderr);
    }
  }

  complex double *S_packed[ss->sym->nirreps];
  for (qpms_iri_t iri = 0; iri < ss->sym->nirreps; ++iri) {
    S_packed[iri] = qpms_scatsys_irrep_pack_matrix(NULL,
        S_full, ss, iri);
    fprintf(stderr, "--- Packed matrix for irrep %d (%s):\n", (int) iri, ss->sym->irreps[iri].name);
    for (size_t row = 0; row < ss->saecv_sizes[iri]; ++row) {
      for (size_t col = 0; col < ss->saecv_sizes[iri]; ++col) {
        complex double elem = S_packed[iri][row * ss->saecv_sizes[iri] + col];
        fprintf(stderr, "%+.3f+%.3fj ", creal(elem), cimag(elem));
      }
      fputc('\n', stderr);
    }
  }
  {
    complex double *S_partrecfull = qpms_scatsys_irrep_unpack_matrix(NULL,
        S_packed[0], ss, 0, false);
    for (qpms_iri_t iri = 0; iri < ss->sym->nirreps; ++iri) {
      qpms_scatsys_irrep_unpack_matrix(S_partrecfull, S_packed[iri],
          ss, iri, false);
      fprintf(stderr, "\nPartial reconstruction %d (%s):\n", (int)iri, ss->sym->irreps[iri].name);
      const size_t full_len = ss->fecv_size;
      for (size_t row = 0 ; row < full_len; ++row) {
        for (size_t col = 0 ; col < full_len; ++col)
          fprintf(stderr, "%+2.3f%+2.3fj ", creal(S_partrecfull[full_len * row + col]), cimag(S_partrecfull[full_len * row + col]));
        fputc('\n', stderr);
      }
    }

    double maxerr = 0;
    for (size_t i = 0; i < ss->fecv_size; ++i) {
      double err = cabs(S_full[i] - S_partrecfull[i]);
      maxerr = (err > maxerr) ? err : maxerr;
    }
    free(S_partrecfull);
  }


  complex double *S_recfull = qpms_scatsys_irrep_unpack_matrix(NULL,
      S_packed[0], ss, 0, false);
  for (qpms_iri_t iri = 1; iri < ss->sym->nirreps; ++iri)
    qpms_scatsys_irrep_unpack_matrix(S_recfull, S_packed[iri],
        ss, iri, true);
  {
    fputs("\n\n", stderr);
    const size_t full_len = ss->fecv_size;
    for (size_t row = 0 ; row < full_len; ++row) {
      for (size_t col = 0 ; col < full_len; ++col)
        fprintf(stderr, "%+2.3f%+2.3fj ", creal(S_recfull[full_len * row + col]), cimag(S_recfull[full_len * row + col]));
      fputc('\n', stderr);
    }
  }

  double maxerr = 0;
  for (size_t i = 0; i < ss->fecv_size; ++i) {
    double err = cabs(S_full[i] - S_recfull[i]);
    maxerr = (err > maxerr) ? err : maxerr;
  }


  printf("maxerr: %lg\n", maxerr);
  
  fprintf(stderr, "pi\tpos\toti\tosn\tp\n");
  for(qpms_ss_pi_t pi = 0; pi < ss->p_count; ++pi) {
    cart3_t pos = ss->p[pi].pos;
    qpms_ss_oti_t oti = ss->p_orbitinfo[pi].t;
    qpms_ss_osn_t osn = ss->p_orbitinfo[pi].osn;
    qpms_ss_orbit_pi_t p = ss->p_orbitinfo[pi].p;
    fprintf(stderr, "%d\t(%.3g,%.3g,%.3g)\t%d\t%d\t%d\n",
        (int)pi, pos.x, pos.y, pos.z, (int)oti, (int)osn, (int)p);
  }

  for (qpms_iri_t iri = 0; iri < ss->sym->nirreps; ++iri) free(S_packed[iri]);
  free(S_full);
  qpms_scatsys_free(ss);
  qpms_tmatrix_free(t1);
  qpms_tmatrix_free(t2);
  qpms_vswf_set_spec_free(b1);
  qpms_vswf_set_spec_free(b2);
  return 0;
}
