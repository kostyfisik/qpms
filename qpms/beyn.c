/* 
 * This file was originally part of SCUFF-EM by M. T. Homer Reid.
 * Modified by Kristian Arjas and Marek Nečada to work without libhmat and libhrutil.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#define STATIC_ASSERT(COND,MSG) typedef char static_assertion_##MSG[(COND)?1:-1]

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

#include <beyn.h>

STATIC_ASSERT((sizeof(lapack_complex_double) == sizeof(gsl_complex)), lapack_and_gsl_complex_must_be_consistent);

// Uniformly random number between -2 and 2
gsl_complex zrandN(){
    double a = (rand()*4.0/RAND_MAX) - 2;
    double b = (rand()*4.0/RAND_MAX) - 2;
    return gsl_complex_rect(a, b);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
BeynSolver *CreateBeynSolver(int M, int L)
{
  BeynSolver *Solver= (BeynSolver *)malloc(sizeof(*Solver));

  Solver->M = M;
  Solver->L = L;
  QPMS_ENSURE(L <= M, "We expect L <= M, but we got L = %d, M = %d ", L, M);

#if 0
  int MLMax = (M>L) ? M : L;
#endif
  int MLMin = (M<L) ? M : L;

  // storage for eigenvalues and eigenvectors
  Solver->Eigenvalues  = gsl_vector_complex_calloc(L);
  Solver->EVErrors     = gsl_vector_complex_calloc(L);
  Solver->Residuals    = gsl_vector_calloc(L);
  Solver->Eigenvectors = gsl_matrix_complex_calloc(M, L);

  // storage for singular values, random VHat matrix, etc. used in algorithm
  Solver->A0           = gsl_matrix_complex_calloc(M,L);
  Solver->A1           = gsl_matrix_complex_calloc(M,L);
  Solver->A0Coarse     = gsl_matrix_complex_calloc(M,L);
  Solver->A1Coarse     = gsl_matrix_complex_calloc(M,L);
  Solver->MInvVHat     = gsl_matrix_complex_calloc(M,L);
  Solver->VHat         = gsl_matrix_complex_calloc(M,L);
  Solver->Sigma        = gsl_vector_calloc(MLMin);
  ReRandomize(Solver,(unsigned)time(NULL));

#if 0
  // internal workspace: need storage for 2 MxL matrices
  // plus 3 LxL matrices
  #define MLBUFFERS 2
  #define LLBUFFERS 3
  size_t ML = MLMax*L, LL = L*L;
#endif
 
  return Solver;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DestroyBeynSolver(BeynSolver *Solver)
{
  gsl_vector_complex_free(Solver->Eigenvalues);
  gsl_vector_complex_free(Solver->EVErrors);  
  gsl_matrix_complex_free(Solver->Eigenvectors);

  gsl_matrix_complex_free(Solver->A0);
  gsl_matrix_complex_free(Solver->A1);
  gsl_matrix_complex_free(Solver->A0Coarse);
  gsl_matrix_complex_free(Solver->A1Coarse);
  gsl_matrix_complex_free(Solver->MInvVHat);
  gsl_vector_free(Solver->Sigma);
  gsl_vector_free(Solver->Residuals);
  gsl_matrix_complex_free(Solver->VHat);

  free(Solver);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

void ReRandomize(BeynSolver *Solver, unsigned int RandSeed)
{
  if (RandSeed==0) 
   RandSeed=time(0);
  srandom(RandSeed);
  gsl_matrix_complex *VHat=Solver->VHat;
  for(int nr=0; nr<VHat->size1; nr++)
   for(int nc=0; nc<VHat->size2; nc++)
    gsl_matrix_complex_set(VHat,nr,nc,zrandN());

}


/***************************************************************/
/* perform linear-algebra manipulations on the A0 and A1       */
/* matrices (obtained via numerical quadrature) to extract     */
/* eigenvalues and eigenvectors                                */
/***************************************************************/

int ProcessAMatrices(BeynSolver *Solver, beyn_function_M_t M_function, 
    void *Params,
    gsl_matrix_complex *A0, gsl_matrix_complex *A1, double complex z0,
    gsl_vector_complex *Eigenvalues, gsl_matrix_complex *Eigenvectors)
{
  int M = Solver->M;
  int L = Solver->L;
  gsl_vector *Sigma = Solver->Sigma;


  int Verbose = 0;//CheckEnv("SCUFF_BEYN_VERBOSE");
  double RankTol=1.0e-4; //CheckEnv("SCUFF_BEYN_RANK_TOL",&RankTol);
  double ResTol=0.0;    // CheckEnv("SCUFF_BEYN_RES_TOL",&ResTol);
 
  // A0 -> V0Full * Sigma * W0TFull' 
  printf(" Beyn: computing SVD...\n");
  gsl_matrix_complex* V0Full = gsl_matrix_complex_calloc(M,L);
  gsl_matrix_complex_memcpy(V0Full,A0);

  gsl_matrix_complex* W0TFull = gsl_matrix_complex_calloc(L,L);
  //A0->SVD(Sigma, &V0Full, &W0TFull);
 
  QPMS_ENSURE(Sigma->stride == 1, "Sigma vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  // FIXME not supported by GSL; use LAPACKE_zgesdd 
  //gsl_linalg_complex_SV_decomp(V0Full, W0TFull, Sigma, work);
  QPMS_ENSURE_SUCCESS(LAPACKE_zgesdd(LAPACK_ROW_MAJOR, // A = U*Σ*conjg(V') 
      'O' /*jobz, 'O' overwrites a with U and saves conjg(V') in vt if m >= n, i.e. if M >= L, which holds*/,
      V0Full->size1 /* m, number of rows */,
      V0Full->size2 /* n, number of columns */,
      (lapack_complex_double *)(V0Full->data) /*a*/,
      V0Full->tda /*lda*/,
      Sigma->data /*s*/, 
      NULL /*u; not used*/, 
      0 /*ldu; not used*/,
      (lapack_complex_double *)W0TFull->data /*vt*/,
      W0TFull->tda /*ldvt*/
      ));

  
  // compute effective rank K (number of eigenvalue candidates)
  int K=0;
  for(int k=0; k<Sigma->size /* this is L, actually */; k++)
   { if (Verbose) printf("Beyn: SV(%i)=%e",k,gsl_vector_get(Sigma, k));
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
  
  gsl_vector_complex *TempM = gsl_vector_complex_calloc(M);
  gsl_vector_complex *TempL = gsl_vector_complex_calloc(L);
  for(int k=0; k<K; k++){
      // It should be rows and cols like this, right..?
      gsl_matrix_complex_get_row(TempM, V0Full,k);
      gsl_matrix_complex_set_row(V0, k, TempM);
      gsl_matrix_complex_get_col(TempL,W0TFull,k);
      gsl_matrix_complex_set_col(W0T,k, TempL);
   }

  gsl_matrix_complex_free(V0Full);
  gsl_matrix_complex_free(W0TFull);
  gsl_vector_complex_free(TempM);
  gsl_vector_complex_free(TempL);
  

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
  
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, one,
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

  gsl_vector_complex *Lambda = gsl_vector_complex_alloc(K); // Eigenvalues
  gsl_matrix_complex *S = gsl_matrix_complex_alloc(K,K); // Eigenvectors

  // FIXME general complex eigensystems not supported by GSL (use LAPACKE_zgee?)
  //gsl_eigen_genv_workspace * W = gsl_eigen_genv_alloc(K);
  //gsl_eigen_genv(B, Eye, alph, beta, S,W);
  //gsl_eigen_genv_free(W);
  QPMS_ENSURE(Sigma->stride == 1, "Sigma vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  QPMS_ENSURE(Lambda->stride == 1, "Lambda vector stride must be 1 for LAPACKE_zgesdd, got %zd.", Sigma->stride);
  QPMS_ENSURE_SUCCESS(LAPACKE_zgeev(
        LAPACK_ROW_MAJOR,
        'N' /* jobvl; don't compute left eigenvectors */,
        'V' /* jobvr; do compute right eigenvectors */,
        K   /* n */,
        (lapack_complex_double *)(S->data) /* a */,
        S->tda  /* lda */,
        (lapack_complex_double *) Lambda->data /* w */,
        NULL /* vl */,
        0 /* ldvl */,
        (lapack_complex_double *)(B->data)/* vr */,
        B->tda/* ldvr */
        ));


 
  //B.NSEig(&Lambda, &S); 

  // V0S <- V0*S
  printf("Multiplying V0*S...\n");
  gsl_matrix_complex *V0S = gsl_matrix_complex_alloc(M, K);
  QPMS_ENSURE_SUCCESS(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
        one, V0, S, zero, V0S));



  int KRetained = 0;
  gsl_matrix_complex *Mmat = gsl_matrix_complex_alloc(M,L);
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
      if (Verbose) printf("Beyn: Residual(%i)=%e\n",k,residual);
    }
    if (ResTol > 0 && residual > ResTol) continue;

    gsl_vector_complex_set(Eigenvalues, KRetained, zgsl);
    if(Eigenvectors) {
      gsl_matrix_complex_set_col(Eigenvectors, KRetained, &(Vk.vector));
      gsl_vector_set(Solver->Residuals, KRetained, residual);
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

/***************************************************************/
/***************************************************************/
/***************************************************************/

int BeynSolve(BeynSolver *Solver, beyn_function_M_t M_function,
   beyn_function_M_inv_Vhat_t M_inv_Vhat_function, void *Params,
              double complex z0, double Rx, double Ry, int N)
{  
  /***************************************************************/
  /* force N to be even so we can simultaneously evaluate        */
  /* the integral with N/2 quadrature points                     */
  /***************************************************************/
    
  if ( (N%2)==1 ) N++;

  /*if (Rx==Ry)
   printf("Applying Beyn method with z0=%s,R=%e,N=%i...\n",z2s(z0),Rx,N);
  else
   printf("Applying Beyn method with z0=%s,Rx=%e,Ry=%e,N=%i...\n",z2s(z0),Rx,Ry,N);
*/
  //int M = Solver->M;
  //int L = Solver->L;
  gsl_matrix_complex *A0           = Solver->A0;
  gsl_matrix_complex *A1           = Solver->A1;
  gsl_matrix_complex *A0Coarse     = Solver->A0Coarse;
  gsl_matrix_complex *A1Coarse     = Solver->A1Coarse;
  gsl_matrix_complex *MInvVHat     = Solver->MInvVHat;
  gsl_matrix_complex *VHat         = Solver->VHat;

  /***************************************************************/
  /* evaluate contour integrals by numerical quadrature to get   */
  /* A0 and A1 matrices                                          */
  /***************************************************************/
  
  gsl_matrix_complex_set_zero(A0);
  gsl_matrix_complex_set_zero(A1);
  gsl_matrix_complex_set_zero(A0Coarse);
  gsl_matrix_complex_set_zero(A1Coarse);
  double DeltaTheta = 2.0*M_PI / ((double)N);
  printf(" Evaluating contour integral (%i points)...\n",N);
  for(int n=0; n<N; n++)
   { 
     double Theta = ((double)n)*DeltaTheta;
     double CT    = cos(Theta), ST=sin(Theta);
     gsl_complex z1   = gsl_complex_rect(Rx*CT, Ry*ST);
     gsl_complex dz   = gsl_complex_rect(Ry*CT/((double)N),(Rx*ST/((double)N)));
     
     double complex zz = Rx*CT + Ry*ST*I;

     //MInvVHat->Copy(VHat);
     gsl_matrix_complex_memcpy(MInvVHat, VHat);

     // Tän pitäis kutsua eval_WT
     // Output ilmeisesti tallentuun neljänteen parametriin
     
     if(M_inv_Vhat_function) {
       QPMS_ENSURE_SUCCESS(M_inv_Vhat_function(MInvVHat, VHat, z0+zz, Params));
     } else {
       QPMS_NOT_IMPLEMENTED("TODO");
     }
     //UserFunc(z0+zz, Params, VHat, MInvVHat);

     gsl_matrix_complex_scale(MInvVHat, dz);
     gsl_matrix_complex_add(A0, MInvVHat);
     if((n%2)==0) {
         gsl_matrix_complex_add(A0Coarse, MInvVHat);
         gsl_matrix_complex_add(A0Coarse, MInvVHat);
     }
     
     gsl_matrix_complex_scale(MInvVHat, z1);
     gsl_matrix_complex_add(A1, MInvVHat);
     if((n%2)==0) {
         gsl_matrix_complex_add(A1Coarse, MInvVHat);
         gsl_matrix_complex_add(A1Coarse, MInvVHat);
     }
   }

  gsl_vector_complex *Eigenvalues  = Solver->Eigenvalues;
  //gsl_vector_complex *EVErrors     = Solver->EVErrors;
  gsl_matrix_complex *Eigenvectors = Solver->Eigenvectors;
  
  int K = ProcessAMatrices(Solver, M_function, Params, A0, A1, z0, Eigenvalues, Eigenvectors);
  //int KCoarse = ProcessAMatrices(Solver, UserFunc, Params, A0Coarse, A1Coarse, z0, EVErrors, Eigenvectors);
//  Log("{K,KCoarse}={%i,%i}",K,KCoarse);
    /*
  for(int k=0; k<EVErrors->N && k<Eigenvalues->N; k++)
   { EVErrors->ZV[k] -= Eigenvalues->ZV[k];
     EVErrors->ZV[k] = cdouble( fabs(real(EVErrors->ZV[k])),
                                fabs(imag(EVErrors->ZV[k]))
                              );
   }

*/
  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
/*
int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunction, void *Params,
              cdouble z0, double R, int N)
{ return BeynSolve(Solver, UserFunction, Params, z0, R, R, N); }
*/
