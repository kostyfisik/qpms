#define STATIC_ASSERT(COND,MSG) typedef char static_assertion_##MSG[(COND)?1:-1]

#define cg2s(x) gsl_complex_tostd(x)
#define cs2g(x) gsl_complex_fromstd(x)

#include <complex.h>
#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "qpms_error.h"
#include <string.h>
#include <cblas.h>

#include "beyn.h"
#define SQ(x) ((x)*(x))

// matrix access
#define MAT(mat_, n_rows_, n_cols_, i_row_, i_col_) (mat_[(n_cols_) * (i_row_) + (i_col_)])

typedef struct BeynSolver
{
  int M;   // dimension of matrices
  int L;   // number of columns of VHat matrix

  complex double *eigenvalues, *eigenvalue_errors; // [L]
  complex double *eigenvectors; // [L][M] (!!!)
  complex double *A0, *A1, *A0_coarse, *A1_coarse, *MInvVHat; // [M][L]
  complex double *VHat; // [M][L]
  double *Sigma, *residuals; // [L]
} BeynSolver;

// constructor, destructor
BeynSolver *BeynSolver_create(int M, int L);
void BeynSolver_free(BeynSolver *solver);

// reset the random matrix VHat used in Beyn's algorithm
void BeynSolver_srandom(BeynSolver *solver, unsigned int RandSeed);

// Uniformly random number from interval [a, b].
static double randU(double a, double b) { return a + (b-a) * random() * (1. / RAND_MAX); }

// Random number from normal distribution (via Box-Muller transform, which is enough for our purposes).
static double randN(double Sigma, double Mu) {
  double u1 = randU(0, 1);
  double u2 = randU(0, 1);
  return Mu + Sigma*sqrt(-2*log(u1))*cos(2.*M_PI*u2);
}

static complex double zrandN(double sigma, double mu){
  return randN(sigma, mu) + I*randN(sigma, mu);
}

static inline double dsq(double a) { return a * a; }

static _Bool beyn_contour_ellipse_inside_test(struct beyn_contour_t *c, complex double z) {
  double rRe = c->z_dz[c->n][0];
  double rIm = c->z_dz[c->n][1];
  complex double zn = z - c->centre;
  return dsq(creal(zn)/rRe) + dsq(cimag(zn)/rIm) < 1;
}

beyn_contour_t *beyn_contour_ellipse(complex double centre, double rRe, double rIm, size_t n) 
{
  beyn_contour_t *c;
  QPMS_CRASHING_MALLOC(c, sizeof(beyn_contour_t) + (n+1)*sizeof(c->z_dz[0]));
  c->centre = centre;
  c->n = n;
  for(size_t i = 0; i < n; ++i) {
    double t = i*2*M_PI/n;
    double st = sin(t), ct = cos(t);
    c->z_dz[i][0] = centre + ct*rRe + I*st*rIm;
    c->z_dz[i][1] = (-rRe*st + I*rIm*ct) * (2*M_PI / n);
  }
  // We hide the half-axes metadata after the discretisation points.
  c->z_dz[n][0] = rRe;
  c->z_dz[n][1] = rIm;
  c->inside_test = beyn_contour_ellipse_inside_test;
  return c;
}


// Sets correct sign to zero for a given branch cut orientation
static inline complex double 
align_zero(complex double z, beyn_contour_halfellipse_orientation or) 
{
  // Maybe redundant, TODO check the standard.
  const double positive_zero = copysign(0., +1.);
  const double negative_zero = copysign(0., -1.);
  switch(or) { // ensure correct zero signs; CHECKME!!!
    case BEYN_CONTOUR_HALFELLIPSE_RE_PLUS:
      if(creal(z) == 0 && signbit(creal(z)))
        z = positive_zero + I * cimag(z);
      break;
    case BEYN_CONTOUR_HALFELLIPSE_RE_MINUS:
      if(creal(z) == 0 && !signbit(creal(z)))
        z = negative_zero + I * cimag(z);
      break;
    case BEYN_CONTOUR_HALFELLIPSE_IM_PLUS:
      if(cimag(z) == 0 && signbit(cimag(z)))
        z = creal(z) + I * positive_zero;
      break;
    case BEYN_CONTOUR_HALFELLIPSE_IM_MINUS:
      if(cimag(z) == 0 && !signbit(cimag(z)))
        z = creal(z) + I * negative_zero;
      break;
    default: QPMS_WTF;
  }
  return z;
}


beyn_contour_t *beyn_contour_halfellipse(complex double centre, double rRe,
    double rIm, size_t n, beyn_contour_halfellipse_orientation or)
{
  beyn_contour_t *c;
  QPMS_CRASHING_MALLOC(c, sizeof(beyn_contour_t) + (n+1)*sizeof(c->z_dz[0])
      + sizeof(beyn_contour_halfellipse_orientation));
  c->centre = centre;
  c->n = n;
  const size_t nline = n/2;
  const size_t narc = n - nline;
  complex double faktor;
  double l = rRe, h = rIm;
  switch(or) {
    case BEYN_CONTOUR_HALFELLIPSE_RE_PLUS:
      faktor = -I;
      l = rIm; h = rRe;
      break;
    case BEYN_CONTOUR_HALFELLIPSE_RE_MINUS:
      faktor = I;
      l = rIm; h = rRe;
      break;
    case BEYN_CONTOUR_HALFELLIPSE_IM_PLUS:
      faktor = 1;
      break;
    case BEYN_CONTOUR_HALFELLIPSE_IM_MINUS:
      faktor = -1;
      break;
    default: QPMS_WTF;
  }

  for(size_t i = 0; i < narc; ++i) {
    double t = (i+0.5)*M_PI/narc;
    double st = sin(t), ct = cos(t);
    c->z_dz[i][0] = centre + faktor*(ct*l + I*st*h);
    c->z_dz[i][1] = faktor * (-l*st + I*h*ct) * (M_PI / narc);
  }
  for(size_t i = 0; i < nline; ++i) {
    double t = 0.5 * (1 - (double) nline) + i;
    c->z_dz[narc + i][0] = align_zero(centre + faktor * t * 2 * l / nline, or);
    c->z_dz[narc + i][1] = faktor * 2 * l / nline;
  }
  // We hide the half-axes metadata after the discretisation points.
  c->z_dz[n][0] = rRe;
  c->z_dz[n][1] = rIm;
  // ugly...
  *((beyn_contour_halfellipse_orientation *) &c->z_dz[n+1][0]) = or;
  c->inside_test = NULL; // TODO beyn_contour_halfellipse_inside_test;
  return c;
}

beyn_contour_t *beyn_contour_kidney(complex double centre, double rRe,
    double rIm, const double rounding, const size_t n, beyn_contour_halfellipse_orientation or)
{
  beyn_contour_t *c;
  QPMS_ENSURE(rounding >= 0 && rounding < .5, "rounding must lie in the interval [0, 0.5)");
  QPMS_CRASHING_MALLOC(c, sizeof(beyn_contour_t) + (n+1)*sizeof(c->z_dz[0])
      + sizeof(beyn_contour_halfellipse_orientation));
  c->centre = centre;
  c->n = n;
  complex double faktor;
  double l = rRe, h = rIm;
  switch(or) {
    case BEYN_CONTOUR_HALFELLIPSE_RE_PLUS:
      faktor = -I;
      l = rIm; h = rRe;
      break;
    case BEYN_CONTOUR_HALFELLIPSE_RE_MINUS:
      faktor = I;
      l = rIm; h = rRe;
      break;
    case BEYN_CONTOUR_HALFELLIPSE_IM_PLUS:
      faktor = 1;
      break;
    case BEYN_CONTOUR_HALFELLIPSE_IM_MINUS:
      faktor = -1;
      break;
    default: QPMS_WTF;
  }

  // Small circle centre coordinates.
  const double y = rounding * h; // distance from the cut / straight line
  const double x = sqrt(SQ(h - y) - SQ(y)); 

  const double alpha = asin(y/(h-y));
  const double ar = l/h; // aspect ratio

  // Parameter range (equal to the contour length if ar == 1)
  const double tmax = 2 * (x + (M_PI_2 + alpha) * y + (M_PI_2 - alpha) * h);
  const double dt = tmax / n;

  size_t i = 0;
  double t;
  // Straight line, first part
  double t_lo = 0, t_hi = x;
  for(; t = i * dt, t <= t_hi; ++i) { 
    c->z_dz[i][0] = align_zero(centre + (t - t_lo) * ar * faktor, or);
    c->z_dz[i][1] = dt * ar * faktor;
  }
  // First small arc
  t_lo = t_hi; t_hi = t_lo + (M_PI_2 + alpha) * y;
  for(; t = i * dt, t < t_hi; ++i) {
    double phi = (t - t_lo) / y - M_PI_2;
    c->z_dz[i][0] = centre + (ar * (x + y * cos(phi)) + y * (1 + sin(phi)) * I) * faktor;
    c->z_dz[i][1] = dt * (- ar * sin(phi) + cos(phi) * I) * faktor;
  }
  // Big arc
  t_lo = t_hi; t_hi = t_lo + (M_PI - 2 * alpha) * h;
  for(; t = i * dt, t < t_hi; ++i) {
    double phi = (t - t_lo) / h + alpha;
    c->z_dz[i][0] = centre + (ar * (h * cos(phi)) + h * sin(phi) * I) * faktor;
    c->z_dz[i][1] = dt * (- ar * sin(phi) + cos(phi) * I) * faktor;
  }
  // Second small arc
  t_lo = t_hi; t_hi = t_lo + (M_PI_2 + alpha) * y;
  for(; t = i * dt, t < t_hi; ++i) {
    double phi = (t - t_lo) / y + M_PI - alpha;
    c->z_dz[i][0] = centre + (ar * (- x + y * cos(phi)) + y * (1 + sin(phi)) * I) * faktor;
    c->z_dz[i][1] = dt * (- ar * sin(phi) + cos(phi) * I) * faktor;
  }
  // Straight line, second part
  t_lo = t_hi; t_hi = tmax;
  for(; t = i * dt, i < n; ++i) {
    c->z_dz[i][0] = align_zero(centre + (t - t_lo - x) * ar * faktor, or);
    c->z_dz[i][1] = dt * ar * faktor;
  }
  
#if 0 // TODO later
  // We hide the half-axes metadata after the discretisation points.
  c->z_dz[n][0] = rRe;
  c->z_dz[n][1] = rIm;
  // ugly...
  *((beyn_contour_halfellipse_orientation *) &c->z_dz[n+1][0]) = or;
#endif
  c->inside_test = NULL; // TODO beyn_contour_halfellipse_inside_test;
  return c;
}

void beyn_result_free(beyn_result_t *r) {
  if(r) {
    free(r->eigval);
    free(r->eigval_err);
    free(r->residuals);
    free(r->eigvec);
    free(r->ranktest_SV);
    free(r);
  }
}

BeynSolver *BeynSolver_create(int M, int L)
{
  BeynSolver *solver;
  QPMS_CRASHING_CALLOC(solver, 1, sizeof(*solver));

  solver->M = M;
  solver->L = L;
  QPMS_ENSURE(L <= M, "We expect L <= M, but we got L = %d, M = %d ", L, M);

  // storage for eigenvalues and eigenvectors
  //solver->eigenvalues  = gsl_vector_complex_calloc(L);
  QPMS_CRASHING_CALLOC(solver->eigenvalues, L, sizeof(complex double));
  //solver->eigenvalue_errors     = gsl_vector_complex_calloc(L);
  QPMS_CRASHING_CALLOC(solver->eigenvalue_errors, L, sizeof(complex double));
  //solver->residuals    = gsl_vector_calloc(L);
  QPMS_CRASHING_CALLOC(solver->residuals, L, sizeof(double));
  //solver->eigenvectors = gsl_matrix_complex_calloc(L, M);
  QPMS_CRASHING_CALLOC(solver->eigenvectors, L * M, sizeof(complex double));

  // storage for singular values, random VHat matrix, etc. used in algorithm
  //solver->A0           = gsl_matrix_complex_calloc(M,L);
  QPMS_CRASHING_CALLOC(solver->A0, M * L, sizeof(complex double));
  //solver->A1           = gsl_matrix_complex_calloc(M,L);
  QPMS_CRASHING_CALLOC(solver->A1, M * L, sizeof(complex double));
  //solver->A0_coarse     = gsl_matrix_complex_calloc(M,L);
  QPMS_CRASHING_CALLOC(solver->A0_coarse, M * L, sizeof(complex double));
  //solver->A1_coarse     = gsl_matrix_complex_calloc(M,L);
  QPMS_CRASHING_CALLOC(solver->A1_coarse, M * L, sizeof(complex double));
  //solver->MInvVHat     = gsl_matrix_complex_calloc(M,L);
  QPMS_CRASHING_CALLOC(solver->MInvVHat, M * L, sizeof(complex double));
  //solver->VHat         = gsl_matrix_complex_calloc(M,L);
  QPMS_CRASHING_CALLOC(solver->VHat, M * L, sizeof(complex double));
  //solver->Sigma        = gsl_vector_calloc(L);
  QPMS_CRASHING_CALLOC(solver->Sigma, L, sizeof(double));
  // Beyn Step 1: Generate random matrix VHat
  BeynSolver_srandom(solver,(unsigned)time(NULL));

  return solver;

}

void BeynSolver_free(BeynSolver *solver)
{
  free(solver->eigenvalues);
  free(solver->eigenvalue_errors);  
  free(solver->eigenvectors);

  free(solver->A0);
  free(solver->A1);
  free(solver->A0_coarse);
  free(solver->A1_coarse);
  free(solver->MInvVHat);
  free(solver->Sigma);
  free(solver->residuals);
  free(solver->VHat);

  free(solver);
}

void BeynSolver_free_partial(BeynSolver *solver) 
{
  free(solver->A0);
  free(solver->A1);
  free(solver->A0_coarse);
  free(solver->A1_coarse);
  free(solver->MInvVHat);
  free(solver->VHat);

  free(solver);
}

void BeynSolver_srandom(BeynSolver *solver, unsigned int RandSeed)
{
  if (RandSeed==0) 
    RandSeed=time(0);
  srandom(RandSeed);
  for(size_t i = 0; i < solver->M * solver->L; ++i)
    solver->VHat[i] = zrandN(1, 0);
}


/*
 * linear-algebra manipulations on the A0 and A1 matrices 
 * (obtained via numerical quadrature) to extract eigenvalues 
 * and eigenvectors                                
 */

static int beyn_process_matrices(BeynSolver *solver, beyn_function_M_t M_function, 
    void *Params,
    complex double *A0, complex double *A1, double complex z0,
    complex double *eigenvalues, complex double *eigenvectors, const double rank_tol, size_t rank_sel_min, const double res_tol)
{
  const size_t m = solver->M;
  const size_t l = solver->L;
  double *Sigma = solver->Sigma;

  int verbose = 1; // TODO

  // Beyn Step 3: Compute SVD: A0 = V0_full * Sigma * W0T_full
  if(verbose) printf(" Beyn: computing SVD...\n");
  //gsl_matrix_complex* V0_full = gsl_matrix_complex_alloc(m,l);
  //gsl_matrix_complex_memcpy(V0_full,A0);
  complex double *V0_full; 
  QPMS_CRASHING_MALLOCPY(V0_full, A0, m * l * sizeof(complex double));
  //gsl_matrix_complex* W0T_full = gsl_matrix_complex_alloc(l,l);
  complex double *W0T_full; QPMS_CRASHING_MALLOC(W0T_full, l * l * sizeof(complex double));
  //QPMS_ENSURE(Sigma->stride == 1, "Sigma vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  //QPMS_ENSURE(V0_full->size1 >= V0_full->size2, "m = %zd, l = %zd, what the hell?");
  QPMS_ENSURE_SUCCESS(LAPACKE_zgesdd(LAPACK_ROW_MAJOR, // A = U*Σ*conjg(V') 
        'O' /*jobz, 'O' overwrites a with U and saves conjg(V') in vt if m >= n, i.e. if M >= L, which holds*/,
        m, // V0_full->size1 /* m, number of rows */,
        l, // V0_full->size2 /* n, number of columns */,
        V0_full, //(lapack_complex_double *)(V0_full->data) /*a*/,
        l, //V0_full->tda /*lda*/,
        Sigma, //Sigma->data /*s*/, 
        NULL /*u; not used*/, 
        m /*ldu; should not be really used but lapacke requires it for some obscure reason*/,
        W0T_full, //(lapack_complex_double *)W0T_full->data /*vt*/,
        l //W0T_full->tda /*ldvt*/
        ));


  // Beyn Step 4: Rank test for Sigma
  // compute effective rank K (number of eigenvalue candidates)
  int K=0;
  for (int k=0; k < l/* k<Sigma->size*/ /* this is l, actually */; k++) { 
    if (verbose) printf("Beyn: SV(%d)=%e\n",k, Sigma[k] /*gsl_vector_get(Sigma, k)*/);
    if (k < rank_sel_min || Sigma[k] /*gsl_vector_get(Sigma, k)*/ > rank_tol)
      K++;
  }
  if (verbose)printf(" Beyn: %d/%zd relevant singular values\n",K,l);
  if (K==0) { 
    QPMS_WARN("no singular values found in Beyn eigensolver\n");
    return 0;
  }

  // Beyn step 5:  B = V0' * A1 * W0 * Sigma^-1
  // set V0, W0T = matrices of first K right, left singular vectors
  //gsl_matrix_complex *V0 = gsl_matrix_complex_alloc(m,K);
  //gsl_matrix_complex *W0T= gsl_matrix_complex_alloc(K,l);
  complex double *V0, *W0T;
  QPMS_CRASHING_MALLOC(V0, m * K * sizeof(complex double));
  QPMS_CRASHING_MALLOC(W0T, K * l * sizeof(complex double));

  // TODO this is stupid, some parts could be handled simply by realloc.
  for (int k = 0; k < K; ++k) {
    //gsl_vector_complex_view tmp;
    //tmp = gsl_matrix_complex_column(V0_full, k);
    //gsl_matrix_complex_set_col(V0, k, &(tmp.vector));
    for(int i = 0; i < m; ++i)
      MAT(V0, m, K, i, k) = MAT(V0_full, m, l, i, k);
    //tmp = gsl_matrix_complex_row(W0T_full, k);
    //gsl_matrix_complex_set_row(W0T, k, &(tmp.vector));
    for(int j = 0; j < K; ++j)
      MAT(W0T, K, l, k, j) = MAT(W0T_full, m, l, k, j);
  }
  //gsl_matrix_complex_free(V0_full);
  free(V0_full);
  //gsl_matrix_complex_free(W0T_full);
  free(W0T_full);

  //gsl_matrix_complex *TM2 = gsl_matrix_complex_calloc(K,l);
  //gsl_matrix_complex *B = gsl_matrix_complex_calloc(K,K);
  complex double *TM2, *B;
  QPMS_CRASHING_CALLOC(TM2, K * l, sizeof(complex double));
  QPMS_CRASHING_CALLOC(B, K * K, sizeof(complex double));

  if(verbose) printf(" Multiplying V0*A1->TM...\n");

  //const gsl_complex one = gsl_complex_rect(1,0);
  //const gsl_complex zero = gsl_complex_rect(0,0);
  //gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, one,
  //    V0, A1, zero, TM2);
  // dims:  V0[m,K], A1[m,l], TM2[K,l]
  const complex double one = 1, zero = 0;
  cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, K, l, m,
      &one, V0, K, A1, l, &zero, TM2, l);


  if(verbose) printf(" Multiplying TM*W0T->B...\n");

  //gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, one,
  //    TM2, W0T, zero, B);
  // DIMS: TM2[K,l], W0T[K,l], B[K,K]
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, K, K, l,
        &one, TM2, l, W0T, l, &zero, B, K);

  //gsl_matrix_complex_free(W0T);
  //gsl_matrix_complex_free(TM2);
  free(W0T);
  free(TM2);

  if(verbose) printf(" Scaling B <- B*Sigma^{-1}\n");
  //gsl_vector_complex *tmp = gsl_vector_complex_calloc(K);
  for(int i = 0; i < K; i++) {
    //gsl_matrix_complex_get_col(tmp, B, i);
    //gsl_vector_complex_scale(tmp, gsl_complex_rect(1.0/gsl_vector_get(Sigma,i), 0));
    //gsl_matrix_complex_set_col(B,i,tmp);
    for(int j = 0; j < K; j++)
      MAT(B, K, K, j, i) /= Sigma[i]; 
  }
  //gsl_vector_complex_free(tmp);

  //for(int m=0; m<K; m++)                  //  B <- B * Sigma^{-1}
  
  // Beyn step 6.
  // Eigenvalue decomposition B -> S*Lambda*S'
  /* According to Beyn's paper (Algorithm 1), one should check conditioning
   * of the eigenvalues; if they are ill-conditioned, one should perform
   * a procedure involving Schur decomposition and reordering.
   *
   * Beyn refers to MATLAB routines eig, condeig, schur, ordschur.
   * I am not sure about the equivalents in LAPACK, TODO check zgeevx, zgeesx.
   */
  if(verbose) printf(" Eigensolving (%i,%i)\n",K,K);

  //gsl_vector_complex *Lambda = gsl_vector_complex_alloc(K); // eigenvalues
  //gsl_matrix_complex *S = gsl_matrix_complex_alloc(K,K); // eigenvectors
  complex double *Lambda /* eigenvalues */ , *S /* eigenvector */;
  QPMS_CRASHING_MALLOC(Lambda, K * sizeof(complex double));
  QPMS_CRASHING_MALLOC(S, K * K * sizeof(complex double));

  //QPMS_ENSURE(Sigma->stride == 1, "Sigma vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  //QPMS_ENSURE(Lambda->stride == 1, "Lambda vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  // dims: B[K,K], S[K,K], Lambda[K]
  QPMS_ENSURE_SUCCESS(LAPACKE_zgeev(
        LAPACK_ROW_MAJOR,
        'N' /* jobvl; don't compute left eigenvectors */,
        'V' /* jobvr; do compute right eigenvectors */,
        K   /* n */,
        B, //(lapack_complex_double *)(B->data) /* a */,
        K, //B->tda  /* lda */,
        Lambda, //(lapack_complex_double *) Lambda->data /* w */,
        NULL /* vl */,
        m /* ldvl, not used but for some reason needed */,
        S, //(lapack_complex_double *)(S->data)/* vr */,
        K //S->tda/* ldvr */
        ));

  //gsl_matrix_complex_free(B);
  free(B);

  // V0S <- V0*S
  printf("Multiplying V0*S...\n");
  //gsl_matrix_complex *V0S = gsl_matrix_complex_alloc(m, K);
  //QPMS_ENSURE_SUCCESS(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
  //      one, V0, S, zero, V0S));
  complex double *V0S;
  QPMS_CRASHING_MALLOC(V0S, m * K * sizeof(complex double));
  // dims: V0[m,K], S[K,K], V0S[m,K]
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, K, K,
        &one, V0, K, S, K, &zero, V0S, K);

  //gsl_matrix_complex_free(V0);
  free(V0);

  // FIXME!!! make clear relation between KRetained and K in the results!
  // (If they differ, there are possibly some spurious eigenvalues.)
  int KRetained = 0;
  //gsl_matrix_complex *Mmat = gsl_matrix_complex_alloc(m, m);
  //gsl_vector_complex *MVk = gsl_vector_complex_alloc(m);
  complex double *Mmat, *MVk;
  QPMS_CRASHING_MALLOC(Mmat, m * m * sizeof(complex double));
  QPMS_CRASHING_MALLOC(MVk, m * sizeof(complex double));
  for (int k = 0; k < K; ++k) {
    //const gsl_complex zgsl = gsl_complex_add(gsl_complex_rect(creal(z0), cimag(z0)), gsl_vector_complex_get(Lambda, k));
    //const complex double z = GSL_REAL(zgsl) + I*GSL_IMAG(zgsl);
    const complex double z = z0 + Lambda[k];
    //gsl_vector_complex_const_view Vk = gsl_matrix_complex_const_column(V0S, k);

    double residual = 0;
    if(res_tol > 0) {
      QPMS_ENSURE_SUCCESS(M_function(Mmat, m,  z, Params));
      //QPMS_ENSURE_SUCCESS(gsl_blas_zgemv(CblasNoTrans, one, Mmat, &(Vk.vector), zero, MVk));
      // Vk[i] == V0S[i, k];  dims: Mmat[m,m], Vk[m] (V0S[m, K]), MVk[m]
      cblas_zgemv(CblasRowMajor, CblasNoTrans, m, m,
            &one, Mmat, m, &MAT(V0S, m, K, 0, k), m /* stride of Vk in V0S */,
            &zero, MVk, 1);

      //residual = gsl_blas_dznrm2(MVk);
      residual = cblas_dznrm2(m, MVk, 1);
      if (verbose) printf("Beyn: Residual(%i)=%e\n",k,residual);
    }
    if (res_tol > 0 && residual > res_tol) continue;

    //gsl_vector_complex_set(eigenvalues, KRetained, zgsl);
    eigenvalues[KRetained] = z;
    if(eigenvectors) {
      //gsl_matrix_complex_set_row(eigenvectors, KRetained, &(Vk.vector));
      for(int j = 0; j < m; ++j) 
        MAT(eigenvectors, l, m, KRetained, j) = MAT(V0S, m, K, j, k);
      //gsl_vector_set(solver->residuals, KRetained, residual);

    }
    ++KRetained;
  }

  free(V0S);
  free(Mmat);
  free(MVk);
  free(S);
  free(Lambda);

  return KRetained;
}


beyn_result_t *beyn_solve(const size_t m, const size_t l,
    beyn_function_M_t M_function, beyn_function_M_inv_Vhat_t M_inv_Vhat_function,
    void *params, const beyn_contour_t *contour,
    double rank_tol, size_t rank_sel_min, double res_tol)
{  
  BeynSolver *solver = BeynSolver_create(m, l);

  complex double *A0           = solver->A0;
  complex double *A1           = solver->A1;
  complex double *A0_coarse     = solver->A0_coarse;
  complex double *A1_coarse     = solver->A1_coarse;
  complex double *MInvVHat     = solver->MInvVHat;
  complex double *VHat         = solver->VHat;

  /***************************************************************/
  /* evaluate contour integrals by numerical quadrature to get   */
  /* A0 and A1 matrices                                          */
  /***************************************************************/

  // TODO zeroing probably redundant (Used calloc...)
  memset(A0, 0, m * l * sizeof(complex double));
  memset(A1, 0, m * l * sizeof(complex double));
  memset(A0_coarse, 0, m * l * sizeof(complex double));
  memset(A1_coarse, 0, m * l * sizeof(complex double));
  const size_t N = contour->n;
  if(N & 1) QPMS_WARN("Contour discretisation point number is odd (%zd),"
     " the error estimates might be a bit off.", N);


  // Beyn Step 2: Computa A0, A1 
  const complex double z0 = contour->centre;
  for(int n=0; n<N; n++) { 
    const complex double z = contour->z_dz[n][0];
    const complex double dz = contour->z_dz[n][1];

    memcpy(MInvVHat, VHat, m * l * sizeof(complex double));

    if(M_inv_Vhat_function) {
      QPMS_ENSURE_SUCCESS(M_inv_Vhat_function(MInvVHat, m, l, VHat, z, params));
    } else {
      lapack_int *pivot;
      //gsl_matrix_complex *Mmat = gsl_matrix_complex_alloc(m,m);
      complex double *Mmat;
      QPMS_CRASHING_MALLOC(Mmat, m * m * sizeof(complex double));
      QPMS_ENSURE_SUCCESS(M_function(Mmat, m, z, params));
      QPMS_CRASHING_MALLOC(pivot, sizeof(lapack_int) * m);
#if 0
      QPMS_ENSURE_SUCCESS(LAPACKE_zgetrf(LAPACK_ROW_MAJOR,
            m /*m*/, m /*n*/,(lapack_complex_double *) Mmat->data /*a*/, Mmat->tda /*lda*/, pivot /*ipiv*/));
      QPMS_ENSURE(MInvVHat->tda == l, "wut?");
      QPMS_ENSURE_SUCCESS(LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N' /*trans*/,
            m /*n*/, l/*nrhs*/, (lapack_complex_double *)(Mmat->data) /*a*/, Mmat->tda /*lda*/, pivot/*ipiv*/, 
            (lapack_complex_double *)(MInvVHat->data) /*b*/, MInvVHat->tda/*ldb*/));
#endif
      QPMS_ENSURE_SUCCESS(LAPACKE_zgetrf(LAPACK_ROW_MAJOR,
            m /*m*/, m /*n*/, Mmat /*a*/, m /*lda*/, pivot /*ipiv*/));
      QPMS_ENSURE_SUCCESS(LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N' /*trans*/,
            m /*n*/, l/*nrhs*/, Mmat /*a*/, m /*lda*/, pivot/*ipiv*/, 
            MInvVHat /*b*/, l /*ldb*/));

      free(pivot);
      free(Mmat);
    }

    // gsl_matrix_complex_scale(MInvVHat, cs2g(dz));
    for(size_t i = 0; i < m * l; ++i)
      MInvVHat[i] *= dz;
    //gsl_matrix_complex_add(A0, MInvVHat);
    for(size_t i = 0; i < m * l; ++i)
      A0[i] += MInvVHat[i];
    if((n%2)==0) {
      //gsl_matrix_complex_add(A0_coarse, MInvVHat);
      //gsl_matrix_complex_add(A0_coarse, MInvVHat);
      for(size_t i = 0; i < m * l; ++i)
        A0_coarse[i] += 2 * MInvVHat[i];
    }

    // A_1 scaling as in Beyn's Remark 3.2(b) for numerical stability.
    //gsl_matrix_complex_scale(MInvVHat, cs2g(z - z0)); 
    for(size_t i = 0; i < m * l; ++i)
      MInvVHat[i] *= (z - z0);
    //gsl_matrix_complex_add(A1, MInvVHat);
    for(size_t i = 0; i < m * l; ++i)
      A1[i] += MInvVHat[i];
    if((n%2)==0) {
      for(size_t i = 0; i < m * l; ++i)
      //gsl_matrix_complex_add(A1_coarse, MInvVHat);
      //gsl_matrix_complex_add(A1_coarse, MInvVHat);
        A1_coarse[i] += 2 * MInvVHat[i];
    }
  }

  complex double *eigenvalues  = solver->eigenvalues;
  complex double *eigenvalue_errors     = solver->eigenvalue_errors;
  complex double *eigenvectors = solver->eigenvectors;

  // Repeat Steps 3–6 with rougher contour approximation to get an error estimate.
  int K_coarse = beyn_process_matrices(solver, M_function, params, A0_coarse, A1_coarse, z0, eigenvalue_errors, /*eigenvectors_coarse*/ NULL, rank_tol, rank_sel_min, res_tol);
  // Reid did also fabs on the complex and real parts ^^^.

  // Beyn Steps 3–6
  int K = beyn_process_matrices(solver, M_function, params, A0, A1, z0, eigenvalues, eigenvectors, rank_tol, rank_sel_min, res_tol);

  //gsl_blas_zaxpy(gsl_complex_rect(-1,0), eigenvalues, eigenvalue_errors);
  const complex double minusone = -1.;
  //TODO maybe change the sizes to correspend to retained eigval count K, not l
  cblas_zaxpy(l, &minusone, eigenvalues, 1, eigenvalue_errors, 1);

  beyn_result_t *result;
  QPMS_CRASHING_MALLOC(result, sizeof(*result));
  result->eigval = solver->eigenvalues;
  result->eigval_err = solver->eigenvalue_errors;
  result->residuals = solver->residuals;
  result->eigvec = solver->eigenvectors;
  result->ranktest_SV = solver->Sigma;
  result->neig = K;

  BeynSolver_free_partial(solver);

  return result;
}

#if 0
// Wrapper of pure C array M-matrix function to GSL.

struct beyn_function_M_carr2gsl_param {
  beyn_function_M_t M_function;
  beyn_function_M_inv_Vhat_t M_inv_Vhat_function;
  void *param;
};

static int beyn_function_M_carr2gsl(gsl_matrix_complex *target_M, complex double z, void *params) 
{
  struct beyn_function_M_carr2gsl_param *p = params;
  // These could rather be asserts.
  QPMS_ENSURE(target_M->size2 == target_M->tda, "Target GSL matrix is not a C-contiguous array. This is a bug, please report!");
  QPMS_ENSURE(target_M->size1 == target_M->size2, "Target is not a square matrix. This is a bug, please report!");
  return p->M_function((complex double *) target_M->data, target_M->size1, z, p->param);
}

static int beyn_function_M_inv_Vhat_carr2gsl(gsl_matrix_complex *target,
    const gsl_matrix_complex *Vhat, complex double z, void *params) 
{
  QPMS_UNTESTED;
  struct beyn_function_M_carr2gsl_param *p = params;
  // These could rather be asserts.
  QPMS_ENSURE(target->size2 == target->tda, "Target GSL matrix is not a C-contiguous array. This is a bug, please report!");
  QPMS_ENSURE(Vhat->size2 == Vhat->tda, "Source GSL matrix is not a C-contiguous array. This is a bug, please report!");
  // TODO here I could also check whether Vhat and target have compatible sizes.
  return p->M_inv_Vhat_function((complex double *) target->data, target->size1, target->size2, 
      (complex double *) Vhat->data, z, p->param);
}

beyn_result_t *beyn_solve(size_t m, size_t l, beyn_function_M_t M, beyn_function_M_inv_Vhat_t M_inv_Vhat,
    void *params, const beyn_contour_t *contour, double rank_tol, size_t rank_sel_min, double res_tol) {
  struct beyn_function_M_carr2gsl_param p = {M, M_inv_Vhat, params};
  return beyn_result_from_beyn_result_gsl(
    beyn_solve_gsl(m, l, beyn_function_M_carr2gsl,
      (p.M_inv_Vhat_function) ? beyn_function_M_inv_Vhat_carr2gsl : NULL, 
      (void *) &p, contour, rank_tol, rank_sel_min, res_tol)
    );
}

beyn_result_t *beyn_result_from_beyn_result_gsl(beyn_result_gsl_t *g) {
  struct beyn_result_t *result;
  QPMS_CRASHING_MALLOC(result, sizeof(beyn_result_t));
  result->gsl = g;
  result->neig = result->gsl->neig;
  result->vlen = result->gsl->eigvec->size2;
  result->eigval = (complex double *) result->gsl->eigval->data;
  result->eigval_err = (complex double *) result->gsl->eigval_err->data;
  result->residuals = result->gsl->residuals->data;
  result->eigvec = (complex double *) result->gsl->eigvec->data;
  result->ranktest_SV = result->gsl->ranktest_SV->data;
  return result;
}
  
void beyn_result_free(beyn_result_t *result) {
  if(result)
    beyn_result_gsl_free(result->gsl);
  free(result);
}
#endif

