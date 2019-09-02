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

// Maybe GSL works?
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include "beyn.h"

STATIC_ASSERT((sizeof(lapack_complex_double) == sizeof(gsl_complex)), lapack_and_gsl_complex_must_be_consistent);


typedef struct BeynSolver
{
  int M;   // dimension of matrices
  int L;   // number of columns of VHat matrix

  gsl_vector_complex *eigenvalues, *eigenvalue_errors;
  gsl_matrix_complex *eigenvectors;
  gsl_matrix_complex *A0, *A1, *A0_coarse, *A1_coarse, *MInvVHat;
  gsl_matrix_complex *VHat;
  gsl_vector *Sigma, *residuals;
  double complex *Workspace;
} BeynSolver;

// constructor, destructor
BeynSolver *BeynSolver_create(int M, int L);
void BeynSolver_free(BeynSolver *solver);

// reset the random matrix VHat used in the Beyn algorithm
//
void BeynSolver_srandom(BeynSolver *solver, unsigned int RandSeed);

// for both of the following routines,
// the return value is the number of eigenvalues found,
// and the eigenvalues and eigenvectors are stored in the
// Lambda and eigenvectors fields of the BeynSolver structure

// Beyn method for circular contour of radius R,
// centered at z0, using N quadrature points
//int BeynSolve(BeynSolver *solver,
//              BeynFunction UserFunction, void *Params,
//              double complex z0, double R, int N);

// Beyn method for elliptical contour of horizontal, vertical
// radii Rx, Ry, centered at z0, using N quadrature points
int BeynSolve(BeynSolver *solver,
    beyn_function_M_gsl_t M_function, beyn_function_M_inv_Vhat_gsl_t M_inv_Vhat_function, void *params,
    double complex z0, double Rx, double Ry, int N);




static double randU(double A, double B) { return A + (B-A) * random() * (1. / RAND_MAX); }

static double randN(double Sigma, double Mu) {
  double u1 = randU(0, 1);
  double u2 = randU(0, 1);
  return Mu + Sigma*sqrt(-2*log(u1))*cos(2.*M_PI*u2);
}

static complex double zrandN(double sigma, double mu){
  return randN(sigma, mu) + I*randN(sigma, mu);
}


beyn_contour_t *beyn_contour_ellipse(complex double centre, double rRe, double rIm, size_t n) 
{
  beyn_contour_t *c;
  QPMS_CRASHING_MALLOC(c, sizeof(beyn_contour_t) + n*sizeof(c->z_dz[0]));
  c->centre = centre;
  c->n = n;
  for(size_t i = 0; i < n; ++i) {
    double t = i*2*M_PI/n;
    double st = sin(t), ct = cos(t);
    c->z_dz[i][0] = centre + ct*rRe + I*st*rIm;
    c->z_dz[i][1] = (I*rRe*st + rIm*ct) / n;
  }
  return c;
}

void beyn_result_gsl_free(beyn_result_gsl_t *r) {
  if(r) {
    gsl_vector_complex_free(r->eigval);
    gsl_vector_complex_free(r->eigval_err);
    gsl_vector_free(r->residuals);
    gsl_matrix_complex_free(r->eigvec);
    free(r);
  }
}

BeynSolver *BeynSolver_create(int M, int L)
{
  BeynSolver *solver= (BeynSolver *)malloc(sizeof(*solver));

  solver->M = M;
  solver->L = L;
  QPMS_ENSURE(L <= M, "We expect L <= M, but we got L = %d, M = %d ", L, M);

  // storage for eigenvalues and eigenvectors
  solver->eigenvalues  = gsl_vector_complex_calloc(L);
  solver->eigenvalue_errors     = gsl_vector_complex_calloc(L);
  solver->residuals    = gsl_vector_calloc(L);
  solver->eigenvectors = gsl_matrix_complex_calloc(M, L);

  // storage for singular values, random VHat matrix, etc. used in algorithm
  solver->A0           = gsl_matrix_complex_calloc(M,L);
  solver->A1           = gsl_matrix_complex_calloc(M,L);
  solver->A0_coarse     = gsl_matrix_complex_calloc(M,L);
  solver->A1_coarse     = gsl_matrix_complex_calloc(M,L);
  solver->MInvVHat     = gsl_matrix_complex_calloc(M,L);
  solver->VHat         = gsl_matrix_complex_calloc(M,L);
  solver->Sigma        = gsl_vector_calloc(L);
  BeynSolver_srandom(solver,(unsigned)time(NULL));

  return solver;

}

void BeynSolver_free(BeynSolver *solver)
{
  gsl_vector_complex_free(solver->eigenvalues);
  gsl_vector_complex_free(solver->eigenvalue_errors);  
  gsl_matrix_complex_free(solver->eigenvectors);

  gsl_matrix_complex_free(solver->A0);
  gsl_matrix_complex_free(solver->A1);
  gsl_matrix_complex_free(solver->A0_coarse);
  gsl_matrix_complex_free(solver->A1_coarse);
  gsl_matrix_complex_free(solver->MInvVHat);
  gsl_vector_free(solver->Sigma);
  gsl_vector_free(solver->residuals);
  gsl_matrix_complex_free(solver->VHat);

  free(solver);
}

void BeynSolver_free_partial(BeynSolver *solver) 
{
  gsl_matrix_complex_free(solver->A0);
  gsl_matrix_complex_free(solver->A1);
  gsl_matrix_complex_free(solver->A0_coarse);
  gsl_matrix_complex_free(solver->A1_coarse);
  gsl_matrix_complex_free(solver->MInvVHat);
  gsl_vector_free(solver->Sigma);
  gsl_matrix_complex_free(solver->VHat);

  free(solver);
}

void BeynSolver_srandom(BeynSolver *solver, unsigned int RandSeed)
{
  if (RandSeed==0) 
    RandSeed=time(0);
  srandom(RandSeed);
  gsl_matrix_complex *VHat=solver->VHat;
  for(int nr=0; nr<VHat->size1; nr++)
    for(int nc=0; nc<VHat->size2; nc++)
      gsl_matrix_complex_set(VHat,nr,nc,cs2g(zrandN(1, 0)));

}


/***************************************************************/
/* perform linear-algebra manipulations on the A0 and A1       */
/* matrices (obtained via numerical quadrature) to extract     */
/* eigenvalues and eigenvectors                                */
/***************************************************************/

static int beyn_process_matrices(BeynSolver *solver, beyn_function_M_gsl_t M_function, 
    void *Params,
    gsl_matrix_complex *A0, gsl_matrix_complex *A1, double complex z0,
    gsl_vector_complex *eigenvalues, gsl_matrix_complex *eigenvectors, const double RankTol, const double ResTol)
{
  int M = solver->M;
  int L = solver->L;
  gsl_vector *Sigma = solver->Sigma;


  int verbose = 1; // TODO

  // A0 -> V0_full * Sigma * W0T_full' 
  printf(" Beyn: computing SVD...\n");
  gsl_matrix_complex* V0_full = gsl_matrix_complex_alloc(M,L);
  gsl_matrix_complex_memcpy(V0_full,A0);

  gsl_matrix_complex* W0T_full = gsl_matrix_complex_alloc(L,L);
  //A0->SVD(Sigma, &V0_full, &W0T_full);

  QPMS_ENSURE(Sigma->stride == 1, "Sigma vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  QPMS_ENSURE(V0_full->size1 >= V0_full->size2, "m = %zd, l = %zd, what the hell?");
  QPMS_ENSURE_SUCCESS(LAPACKE_zgesdd(LAPACK_ROW_MAJOR, // A = U*Σ*conjg(V') 
        'O' /*jobz, 'O' overwrites a with U and saves conjg(V') in vt if m >= n, i.e. if M >= L, which holds*/,
        V0_full->size1 /* m, number of rows */,
        V0_full->size2 /* n, number of columns */,
        (lapack_complex_double *)(V0_full->data) /*a*/,
        V0_full->tda /*lda*/,
        Sigma->data /*s*/, 
        NULL /*u; not used*/, 
        M /*ldu; should not be really used but lapacke requires it for some obscure reason*/,
        (lapack_complex_double *)W0T_full->data /*vt*/,
        W0T_full->tda /*ldvt*/
        ));


  // compute effective rank K (number of eigenvalue candidates)
  int K=0;
  for(int k=0; k<Sigma->size /* this is L, actually */; k++)
  { if (verbose) printf("Beyn: SV(%i)=%e\n",k,gsl_vector_get(Sigma, k));
    if (gsl_vector_get(Sigma, k) > RankTol )
      K++;
  }
  printf(" Beyn: %i/%i relevant singular values\n",K,L);
  if (K==0)
  { printf("no singular values found in Beyn eigensolver\n");
    return 0;
  }

  // set V0, W0T = matrices of first K right, left singular vectors
  gsl_matrix_complex *V0 = gsl_matrix_complex_alloc(M,K);
  gsl_matrix_complex *W0T= gsl_matrix_complex_alloc(K,L);

  for (int k = 0; k < K; ++k) {
    gsl_vector_complex_view tmp;
    tmp = gsl_matrix_complex_column(V0_full, k);
    gsl_matrix_complex_set_col(V0, k, &(tmp.vector));
    tmp = gsl_matrix_complex_row(W0T_full, k);
    gsl_matrix_complex_set_row(W0T, k, &(tmp.vector));
  }

  gsl_matrix_complex_free(V0_full);
  gsl_matrix_complex_free(W0T_full);

  // B <- V0' * A1 * W0 * Sigma^-1
  gsl_matrix_complex *TM2 = gsl_matrix_complex_calloc(K,L);
  gsl_matrix_complex *B = gsl_matrix_complex_calloc(K,K);

  printf(" Multiplying V0*A1->TM...\n");
  //V0.Multiply(A1, &TM2, "--transA C");   // TM2 <- V0' * A1
  const gsl_complex one = gsl_complex_rect(1,0);
  const gsl_complex zero = gsl_complex_rect(0,0);
  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, one,
      V0, A1, zero, TM2);

  printf(" Multiplying TM*W0T->B...\n");
  //TM2.Multiply(&W0T, &B, "--transB C");   //  B <- TM2 * W0

  gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, one,
      TM2, W0T, zero, B);

  gsl_matrix_complex_free(W0T);
  gsl_matrix_complex_free(TM2);

  printf(" Scaling B <- B*Sigma^{-1}\n");
  gsl_vector_complex *tmp = gsl_vector_complex_calloc(K);
  for(int i = 0; i < K; i++){
    gsl_matrix_complex_get_col(tmp, B, i);
    gsl_vector_complex_scale(tmp, gsl_complex_rect(1.0/gsl_vector_get(Sigma,i), 0));
    gsl_matrix_complex_set_col(B,i,tmp);
  }
  gsl_vector_complex_free(tmp);

  //for(int m=0; m<K; m++)                  //  B <- B * Sigma^{-1}

  // for(int n=0; n<K; n++)
  //  B.ScaleEntry(m,n,1.0/Sigma->GetEntry(n));

  // B -> S*Lambda*S'
  printf(" Eigensolving (%i,%i)\n",K,K);

  gsl_vector_complex *Lambda = gsl_vector_complex_alloc(K); // eigenvalues
  gsl_matrix_complex *S = gsl_matrix_complex_alloc(K,K); // eigenvectors

  QPMS_ENSURE(Sigma->stride == 1, "Sigma vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  QPMS_ENSURE(Lambda->stride == 1, "Lambda vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  QPMS_ENSURE_SUCCESS(LAPACKE_zgeev(
        LAPACK_ROW_MAJOR,
        'N' /* jobvl; don't compute left eigenvectors */,
        'V' /* jobvr; do compute right eigenvectors */,
        K   /* n */,
        (lapack_complex_double *)(B->data) /* a */,
        B->tda  /* lda */,
        (lapack_complex_double *) Lambda->data /* w */,
        NULL /* vl */,
        M /* ldvl, not used by for some reason needed */,
        (lapack_complex_double *)(S->data)/* vr */,
        S->tda/* ldvr */
        ));

  gsl_matrix_complex_free(B);

  // V0S <- V0*S
  printf("Multiplying V0*S...\n");
  gsl_matrix_complex *V0S = gsl_matrix_complex_alloc(M, K);
  QPMS_ENSURE_SUCCESS(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
        one, V0, S, zero, V0S));

  gsl_matrix_complex_free(V0);

  int KRetained = 0;
  gsl_matrix_complex *Mmat = gsl_matrix_complex_alloc(M,M);
  gsl_vector_complex *MVk = gsl_vector_complex_alloc(M);
  for (int k = 0; k < K; ++k) {
    const gsl_complex zgsl = gsl_complex_add(gsl_complex_rect(creal(z0), cimag(z0)), gsl_vector_complex_get(Lambda, k));
    const complex double z = GSL_REAL(zgsl) + I*GSL_IMAG(zgsl);
    gsl_vector_complex_const_view Vk = gsl_matrix_complex_const_column(V0S, k);

    double residual = 0;
    if(ResTol > 0) {
      QPMS_ENSURE_SUCCESS(M_function(Mmat, z, Params));
      QPMS_ENSURE_SUCCESS(gsl_blas_zgemv(CblasNoTrans, one, Mmat, &(Vk.vector), zero, MVk));
      residual = gsl_blas_dznrm2(MVk);
      if (verbose) printf("Beyn: Residual(%i)=%e\n",k,residual);
    }
    if (ResTol > 0 && residual > ResTol) continue;

    gsl_vector_complex_set(eigenvalues, KRetained, zgsl);
    if(eigenvectors) {
      gsl_matrix_complex_set_col(eigenvectors, KRetained, &(Vk.vector));
      gsl_vector_set(solver->residuals, KRetained, residual);
    }
    ++KRetained;
  }

  gsl_matrix_complex_free(V0S);
  gsl_matrix_complex_free(Mmat);
  gsl_vector_complex_free(MVk);
  gsl_matrix_complex_free(S);
  gsl_vector_complex_free(Lambda);

  return KRetained;
}


int beyn_solve_gsl(beyn_result_gsl_t **result, size_t m, size_t l,
    beyn_function_M_gsl_t M_function, beyn_function_M_inv_Vhat_gsl_t M_inv_Vhat_function,
    void *params, const beyn_contour_t *contour,
    double rank_tol, double res_tol)
{  
  BeynSolver *solver = BeynSolver_create(m, l);

  const int M = solver->M;
  const int L = solver->L;
  gsl_matrix_complex *A0           = solver->A0;
  gsl_matrix_complex *A1           = solver->A1;
  gsl_matrix_complex *A0_coarse     = solver->A0_coarse;
  gsl_matrix_complex *A1_coarse     = solver->A1_coarse;
  gsl_matrix_complex *MInvVHat     = solver->MInvVHat;
  gsl_matrix_complex *VHat         = solver->VHat;

  /***************************************************************/
  /* evaluate contour integrals by numerical quadrature to get   */
  /* A0 and A1 matrices                                          */
  /***************************************************************/

  gsl_matrix_complex_set_zero(A0);
  gsl_matrix_complex_set_zero(A1);
  gsl_matrix_complex_set_zero(A0_coarse);
  gsl_matrix_complex_set_zero(A1_coarse);
  const size_t N = contour->n;
  if(N & 1) QPMS_WARN("Contour discretisation point number is odd (%zd),"
     " the error estimates might be a bit off.", N);

  const complex double z0 = contour->centre;
  for(int n=0; n<N; n++) { 
    const complex double z = contour->z_dz[n][0];
    const complex double dz = contour->z_dz[n][1];

    gsl_matrix_complex_memcpy(MInvVHat, VHat);

    // Tän pitäis kutsua eval_WT
    // Output ilmeisesti tallentuun neljänteen parametriin

    if(M_inv_Vhat_function) {
      QPMS_ENSURE_SUCCESS(M_inv_Vhat_function(MInvVHat, VHat, z, params));
    } else {
      lapack_int *pivot;
      gsl_matrix_complex *Mmat = gsl_matrix_complex_alloc(M,M);
      QPMS_ENSURE_SUCCESS(M_function(Mmat, z, params));
      QPMS_CRASHING_MALLOC(pivot, sizeof(lapack_int) * M);
      QPMS_ENSURE_SUCCESS(LAPACKE_zgetrf(LAPACK_ROW_MAJOR,
            M /*m*/, M /*n*/,(lapack_complex_double *) Mmat->data /*a*/, Mmat->tda /*lda*/, pivot /*ipiv*/));
      QPMS_ENSURE(MInvVHat->tda == L, "wut?");
      QPMS_ENSURE_SUCCESS(LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N' /*trans*/,
            M /*n*/, L/*nrhs*/, (lapack_complex_double *)(Mmat->data) /*a*/, Mmat->tda /*lda*/, pivot/*ipiv*/, 
            (lapack_complex_double *)(MInvVHat->data) /*b*/, MInvVHat->tda/*ldb*/));

      free(pivot);
      gsl_matrix_complex_free(Mmat);
    }

    gsl_matrix_complex_scale(MInvVHat, cs2g(dz));
    gsl_matrix_complex_add(A0, MInvVHat);
    if((n%2)==0) {
      gsl_matrix_complex_add(A0_coarse, MInvVHat);
      gsl_matrix_complex_add(A0_coarse, MInvVHat);
    }

    gsl_matrix_complex_scale(MInvVHat, cs2g(z - z0));
    gsl_matrix_complex_add(A1, MInvVHat);
    if((n%2)==0) {
      gsl_matrix_complex_add(A1_coarse, MInvVHat);
      gsl_matrix_complex_add(A1_coarse, MInvVHat);
    }
  }

  gsl_vector_complex *eigenvalues  = solver->eigenvalues;
  gsl_vector_complex *eigenvalue_errors     = solver->eigenvalue_errors;
  gsl_matrix_complex *eigenvectors = solver->eigenvectors;

  int K = beyn_process_matrices(solver, M_function, params, A0, A1, z0, eigenvalues, eigenvectors, rank_tol, res_tol);
  int K_coarse = beyn_process_matrices(solver, M_function, params, A0_coarse, A1_coarse, z0, eigenvalue_errors, eigenvectors, rank_tol, res_tol);

  gsl_blas_zaxpy(gsl_complex_rect(-1,0), eigenvalues, eigenvalue_errors);
  // TODO Original did also fabs on the complex and real parts ^^^.

  QPMS_CRASHING_MALLOC(*result, sizeof(beyn_result_gsl_t));
  (*result)->eigval = solver->eigenvalues;
  (*result)->eigval_err = solver->eigenvalue_errors;
  (*result)->residuals = solver->residuals;
  (*result)->eigvec = solver->eigenvectors;

  BeynSolver_free_partial(solver);

  return K;
}

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

int beyn_solve(beyn_result_t **result, size_t m, size_t l, beyn_function_M_t M, beyn_function_M_inv_Vhat_t M_inv_Vhat,
    void *params, const beyn_contour_t *contour, double rank_tol, double res_tol) {
  struct beyn_function_M_carr2gsl_param p = {M, M_inv_Vhat, params};
  QPMS_CRASHING_MALLOC(*result, sizeof(beyn_result_t));
  int retval = beyn_solve_gsl(&((*result)->gsl), m, l, beyn_function_M_carr2gsl,
      (p.M_inv_Vhat_function) ? beyn_function_M_inv_Vhat_carr2gsl : NULL, 
      (void *) &p, contour, rank_tol, res_tol);
  (*result)->neig = (*result)->gsl->neig;
  (*result)->eigval = (complex double *) (*result)->gsl->eigval->data;
  (*result)->eigval_err = (complex double *) (*result)->gsl->eigval_err->data;
  (*result)->residuals = (*result)->gsl->residuals->data;
  (*result)->eigvec = (complex double *) (*result)->gsl->eigvec->data;
  return retval;
}

void beyn_result_free(beyn_result_t *result) {
  if(result)
    beyn_result_gsl_free(result->gsl);
  free(result);
}

