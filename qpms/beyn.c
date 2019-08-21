/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
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

/*
 * libBeyn.cc      -- implementation of Beyn's method for
 *                 -- nonlinear eigenvalue problems
 *
 * Homer Reid      -- 6/2016
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <complex.h>

// Maybe GSL works?
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include <beyn.h>

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

  int MLMax = (M>L) ? M : L;
  int MLMin = (M<L) ? M : L;

  // storage for eigenvalues and eigenvectors
  Solver->Eigenvalues  = gsl_vector_complex_calloc(L);
  Solver->EVErrors     = gsl_vector_complex_calloc(L);
  Solver->Residuals    = gsl_vector_complex_calloc(L);
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

  // internal workspace: need storage for 2 MxL matrices
  // plus 3 LxL matrices
  #define MLBUFFERS 2
  #define LLBUFFERS 3
  size_t ML = MLMax*L, LL = L*L;
  size_t WorkspaceSize = (MLBUFFERS*ML + LLBUFFERS*LL)*sizeof(double complex);
 
  Solver->Workspace = (double complex*)malloc( WorkspaceSize );

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
  gsl_matrix_complex_free(Solver->VHat);

  free(Solver->Workspace);

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

int ProcessAMatrices(BeynSolver *Solver, BeynFunction UserFunc, 
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
  gsl_vector_complex *work = gsl_vector_complex_alloc(M);
 
  // FIXME not supported by GSL; use LAPACKE_zgesdd 
  gsl_linalg_complex_SV_decomp(V0Full, W0TFull, Sigma, work);

  
  // compute effective rank K (number of eigenvalue candidates)
  int K=0;
  for(int k=0; k<Sigma->size; k++)
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
  gsl_vector_complex_free(work);
  gsl_vector_complex_free(TempM);
  gsl_vector_complex_free(TempL);
  

  // B <- V0' * A1 * W0 * Sigma^-1
  gsl_matrix_complex *TM2 = gsl_matrix_complex_calloc(K,L);
  gsl_matrix_complex *B = gsl_matrix_complex_calloc(K,K);

  printf(" Multiplying V0*A1->TM...\n");
  //V0.Multiply(A1, &TM2, "--transA C");   // TM2 <- V0' * A1
  gsl_complex one = gsl_complex_rect(1,0);
  gsl_complex zero = gsl_complex_rect(0,0);
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

  gsl_vector_complex *Lambda = gsl_vector_complex_calloc(K); // Eigenvalues
  gsl_matrix_complex *S = gsl_matrix_complex_calloc(K,K); // Eigenvectors
  gsl_matrix_complex *Eye = gsl_matrix_complex_alloc(K,K);

  gsl_vector_complex *alph = gsl_vector_complex_calloc(K);
  gsl_vector_complex *beta = gsl_vector_complex_calloc(K);

  gsl_matrix_complex_set_identity(Eye);

  // FIXME general complex eigensystems not supported by GSL (use LAPACKE_zgee?)
  gsl_eigen_genv_workspace * W = gsl_eigen_genv_alloc(K);
  gsl_eigen_genv(B, Eye, alph, beta, S,W);
  gsl_eigen_genv_free(W);

  gsl_complex tmpa;
  gsl_complex tmpb;
  for(int i = 0; i < K; i++){
      tmpb = gsl_vector_complex_get(beta,i);
      tmpa = gsl_vector_complex_get(alph,i);
      if(gsl_complex_abs(tmpb)){
          gsl_vector_complex_set(Lambda, i, gsl_complex_div(tmpa,tmpb));
          printf("Eigenvalue %e + %e i found\n",GSL_REAL(gsl_complex_div(tmpa,tmpb)), GSL_IMAG(gsl_complex_div(tmpa,tmpb)));
      } else{
          printf("Beta %d is zero \n",i);
      }
      if(!gsl_complex_abs(tmpa)){
          printf("Alpha %d is zero \n",i);
      }
  }

  gsl_vector_complex_free(alph);
  gsl_vector_complex_free(beta);
  gsl_matrix_complex_free(Eye);
 
  //B.NSEig(&Lambda, &S);

  
  // V0S <- V0*S
  printf("Multiplying V0*S...\n");

  gsl_vector_complex *V = gsl_vector_complex_alloc(K);
  gsl_vector_complex *s = gsl_vector_complex_alloc(K);

  printf("Evaluating retained values \n");
  int KRetained=0;
  gsl_vector_complex * om = gsl_vector_complex_alloc(1);
  for(int k=0; k<K; k++)
   { 
     if(gsl_complex_abs(gsl_vector_complex_get(Lambda,k))){
         gsl_complex tmp_c = gsl_vector_complex_get(Lambda, k);
         double complex  z = z0 + GSL_REAL(tmp_c) + GSL_IMAG(tmp_c)*I;
         //gsl_matrix_get_col(V, V0S, k);
         gsl_matrix_complex_get_col(s, S, k);
         gsl_blas_zgemv(CblasNoTrans, one, V0, s, zero, V);

 
     double Residual=0.0;
     if (ResTol>0.0)
      { /*gsl_matrix_complex Vk(M,1,V);
        gsl_matrix_complex MVk(M,1,MLBuffers[0]);
        UserFunc(z, Params, &Vk, &MVk);
        Residual=VecNorm(MVk.ZM, M);
        */
        gsl_vector_complex_set(om,1,tmp_c);
        Residual = min_sing(om,Params); // in unitcell.c
        if (1) printf("Beyn: Residual(%i)=%e\n",k,Residual);
      }
     if (ResTol>0.0 && Residual>ResTol) continue;
    
    //Eigenvalues->SetEntry(KRetained, z);
    gsl_vector_complex_set(Eigenvalues, KRetained, tmp_c);
    gsl_matrix_complex_set_col(Eigenvectors, KRetained, V);
    /*if (Eigenvectors) 
     { 
         //Eigenvectors->SetEntries(":", KRetained, V);
       //Solver->Residuals->SetEntry(KRetained,Residual);
     }
     */
    KRetained++;
    }
   }

  printf("%d eigenvalues found \n",KRetained);

  gsl_matrix_complex_free(S);
  gsl_matrix_complex_free(V0);
  gsl_vector_complex_free(Lambda);
  gsl_vector_complex_free(om);
  return KRetained;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunc, void *Params,
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
     
     UserFunc(z0+zz, Params, VHat, MInvVHat);

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
  
  int K = ProcessAMatrices(Solver, UserFunc, Params, A0, A1, z0, Eigenvalues, Eigenvectors);
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
