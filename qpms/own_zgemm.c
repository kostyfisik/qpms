/* IMPORTANT! This code is partially taken from GSL, so everything must be GPL'd
 * or this has to be rewritten (or removed; the only reason to use this are problems
 * with OpenBLAS) when distributed.
 */

#include "qpmsblas.h"
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

void
cblas_xerbla (int p, const char *rout, const char *form, ...)
{
  va_list ap;

  va_start (ap, form);

  if (p)
    {
      fprintf (stderr, "Parameter %d to routine %s was incorrect\n", p, rout);
    }

  vfprintf (stderr, form, ap);
  va_end (ap);

  abort ();
}


#define BASE double

#define INDEX QPMS_BLAS_INDEX_T
#define OFFSET(N, incX) ((incX) > 0 ?  0 : ((N) - 1) * (-(incX)))
#define BLAS_ERROR(x)  cblas_xerbla(0, __FILE__, x); 

#define MAX(x,y) (((x) < (y)) ? (y) : (x))

#define CONJUGATE(x) ((x) == CblasConjTrans)
#define TRANSPOSE(x) ((x) == CblasTrans || (x) == CblasConjTrans)
#define UPPER(x) ((x) == CblasUpper)
#define LOWER(x) ((x) == CblasLower)

/* Handling of packed complex types... */

#define REAL(a,i) (((BASE *) a)[2*(i)])
#define IMAG(a,i) (((BASE *) a)[2*(i)+1])

#define REAL0(a) (((BASE *)a)[0])
#define IMAG0(a) (((BASE *)a)[1])

#define CONST_REAL(a,i) (((const BASE *) a)[2*(i)])
#define CONST_IMAG(a,i) (((const BASE *) a)[2*(i)+1])

#define CONST_REAL0(a) (((const BASE *)a)[0])
#define CONST_IMAG0(a) (((const BASE *)a)[1])


#define GB(KU,KL,lda,i,j) ((KU+1+(i-j))*lda + j)

#define TRCOUNT(N,i) ((((i)+1)*(2*(N)-(i)))/2)

/* #define TBUP(N,i,j) */
/* #define TBLO(N,i,j) */

#define TPUP(N,i,j) (TRCOUNT(N,(i)-1)+(j)-(i))
#define TPLO(N,i,j) (((i)*((i)+1))/2 + (j))


/* check if CBLAS_ORDER is correct */
#define CHECK_ORDER(pos,posIfError,order) \
if(((order)!=CblasRowMajor)&&((order)!=CblasColMajor)) \
     pos = posIfError;

/* check if CBLAS_TRANSPOSE is correct */
#define CHECK_TRANSPOSE(pos,posIfError,Trans) \
if(((Trans)!=CblasNoTrans)&&((Trans)!=CblasTrans)&&((Trans)!=CblasConjTrans)) \
    pos = posIfError;

/* check if a dimension argument is correct */
#define CHECK_DIM(pos,posIfError,dim) \
if((dim)<0) \
    pos = posIfError;

/* cblas_xgemm() */
#define CBLAS_ERROR_GEMM(pos,Order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc) \
{ \
    CBLAS_TRANSPOSE __transF=CblasNoTrans,__transG=CblasNoTrans; \
    if((Order)==CblasRowMajor) { \
        __transF = ((TransA)!=CblasConjTrans) ? (TransA) : CblasTrans; \
        __transG = ((TransB)!=CblasConjTrans) ? (TransB) : CblasTrans; \
    } else { \
        __transF = ((TransB)!=CblasConjTrans) ? (TransB) : CblasTrans; \
        __transG = ((TransA)!=CblasConjTrans) ? (TransA) : CblasTrans; \
    } \
    CHECK_ORDER(pos,1,Order); \
    CHECK_TRANSPOSE(pos,2,TransA); \
    CHECK_TRANSPOSE(pos,3,TransB); \
    CHECK_DIM(pos,4,M); \
    CHECK_DIM(pos,5,N); \
    CHECK_DIM(pos,6,K); \
    if((Order)==CblasRowMajor) { \
        if(__transF==CblasNoTrans) { \
            if((lda)<MAX(1,(K))) { \
                (pos) = 9; \
            } \
        } else { \
            if((lda)<MAX(1,(M))) { \
                (pos) = 9; \
            } \
        } \
        if(__transG==CblasNoTrans) { \
            if((ldb)<MAX(1,(N))) { \
                (pos) = 11; \
            } \
        } else { \
            if((ldb)<MAX(1,(K))) { \
                (pos) = 11; \
            } \
        } \
        if((ldc)<MAX(1,(N))) { \
            (pos) = 14; \
        } \
    } else if((Order)==CblasColMajor) { \
        if(__transF==CblasNoTrans) { \
            if((ldb)<MAX(1,(K))) { \
                (pos) = 11; \
            } \
        } else { \
            if((ldb)<MAX(1,(N))) { \
                (pos) = 11; \
            } \
        } \
        if(__transG==CblasNoTrans) { \
            if((lda)<MAX(1,(M))) { \
                (pos) = 9; \
            } \
        } else { \
            if((lda)<MAX(1,(K))) { \
                (pos) = 9; \
            } \
        } \
        if((ldc)<MAX(1,(M))) { \
            (pos) = 14; \
        } \
    } \
}



#define CHECK_ARGS_X(FUNCTION,VAR,ARGS) do { int VAR = 0 ;      \
    CBLAS_ERROR_##FUNCTION ARGS ; \
    if (VAR) cblas_xerbla(pos,__FILE__,""); } while (0)

#define CHECK_ARGS14(FUNCTION,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14))

void qpms_zgemm(CBLAS_LAYOUT Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, 
    const INDEX M, const INDEX N, const INDEX K, 
    const _Complex double *alpha, const _Complex double *A, const INDEX lda,
    const _Complex double *B, const INDEX ldb,
    const _Complex double *beta, _Complex double *C, const INDEX ldc)
{
  INDEX i, j, k;
  INDEX n1, n2;
  INDEX ldf, ldg;
  int conjF, conjG, TransF, TransG;
  const BASE *F, *G;

  CHECK_ARGS14(GEMM,Order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);

  {
    const BASE alpha_real = CONST_REAL0(alpha);
    const BASE alpha_imag = CONST_IMAG0(alpha);

    const BASE beta_real = CONST_REAL0(beta);
    const BASE beta_imag = CONST_IMAG0(beta);

    if ((alpha_real == 0.0 && alpha_imag == 0.0)
        && (beta_real == 1.0 && beta_imag == 0.0))
      return;

    if (Order == CblasRowMajor) {
      n1 = M;
      n2 = N;
      F = (const BASE *)A;
      ldf = lda;
      conjF = (TransA == CblasConjTrans) ? -1 : 1;
      TransF = (TransA == CblasNoTrans) ? CblasNoTrans : CblasTrans;
      G = (const BASE *)B;
      ldg = ldb;
      conjG = (TransB == CblasConjTrans) ? -1 : 1;
      TransG = (TransB == CblasNoTrans) ? CblasNoTrans : CblasTrans;
    } else {
      n1 = N;
      n2 = M;
      F = (const BASE *)B;
      ldf = ldb;
      conjF = (TransB == CblasConjTrans) ? -1 : 1;
      TransF = (TransB == CblasNoTrans) ? CblasNoTrans : CblasTrans;
      G = (const BASE *)A;
      ldg = lda;
      conjG = (TransA == CblasConjTrans) ? -1 : 1;
      TransG = (TransA == CblasNoTrans) ? CblasNoTrans : CblasTrans;
    }

    /* form  y := beta*y */
    if (beta_real == 0.0 && beta_imag == 0.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          REAL(C, ldc * i + j) = 0.0;
          IMAG(C, ldc * i + j) = 0.0;
        }
      }
    } else if (!(beta_real == 1.0 && beta_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const BASE Cij_real = REAL(C, ldc * i + j);
          const BASE Cij_imag = IMAG(C, ldc * i + j);
          REAL(C, ldc * i + j) = beta_real * Cij_real - beta_imag * Cij_imag;
          IMAG(C, ldc * i + j) = beta_real * Cij_imag + beta_imag * Cij_real;
        }
      }
    }

    if (alpha_real == 0.0 && alpha_imag == 0.0)
      return;

    if (TransF == CblasNoTrans && TransG == CblasNoTrans) {

      /* form  C := alpha*A*B + C */

      for (k = 0; k < K; k++) {
        for (i = 0; i < n1; i++) {
          const BASE Fik_real = CONST_REAL(F, ldf * i + k);
          const BASE Fik_imag = conjF * CONST_IMAG(F, ldf * i + k);
          const BASE temp_real = alpha_real * Fik_real - alpha_imag * Fik_imag;
          const BASE temp_imag = alpha_real * Fik_imag + alpha_imag * Fik_real;
          if (!(temp_real == 0.0 && temp_imag == 0.0)) {
            for (j = 0; j < n2; j++) {
              const BASE Gkj_real = CONST_REAL(G, ldg * k + j);
              const BASE Gkj_imag = conjG * CONST_IMAG(G, ldg * k + j);
              REAL(C, ldc * i + j) += temp_real * Gkj_real - temp_imag * Gkj_imag;
              IMAG(C, ldc * i + j) += temp_real * Gkj_imag + temp_imag * Gkj_real;
            }
          }
        }
      }

    } else if (TransF == CblasNoTrans && TransG == CblasTrans) {

      /* form  C := alpha*A*B' + C */

      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          BASE temp_real = 0.0;
          BASE temp_imag = 0.0;
          for (k = 0; k < K; k++) {
            const BASE Fik_real = CONST_REAL(F, ldf * i + k);
            const BASE Fik_imag = conjF * CONST_IMAG(F, ldf * i + k);
            const BASE Gjk_real = CONST_REAL(G, ldg * j + k);
            const BASE Gjk_imag = conjG * CONST_IMAG(G, ldg * j + k);
            temp_real += Fik_real * Gjk_real - Fik_imag * Gjk_imag;
            temp_imag += Fik_real * Gjk_imag + Fik_imag * Gjk_real;
          }
          REAL(C, ldc * i + j) += alpha_real * temp_real - alpha_imag * temp_imag;
          IMAG(C, ldc * i + j) += alpha_real * temp_imag + alpha_imag * temp_real;
        }
      }

    } else if (TransF == CblasTrans && TransG == CblasNoTrans) {

      for (k = 0; k < K; k++) {
        for (i = 0; i < n1; i++) {
          const BASE Fki_real = CONST_REAL(F, ldf * k + i);
          const BASE Fki_imag = conjF * CONST_IMAG(F, ldf * k + i);
          const BASE temp_real = alpha_real * Fki_real - alpha_imag * Fki_imag;
          const BASE temp_imag = alpha_real * Fki_imag + alpha_imag * Fki_real;
          if (!(temp_real == 0.0 && temp_imag == 0.0)) {
            for (j = 0; j < n2; j++) {
              const BASE Gkj_real = CONST_REAL(G, ldg * k + j);
              const BASE Gkj_imag = conjG * CONST_IMAG(G, ldg * k + j);
              REAL(C, ldc * i + j) += temp_real * Gkj_real - temp_imag * Gkj_imag;
              IMAG(C, ldc * i + j) += temp_real * Gkj_imag + temp_imag * Gkj_real;
            }
          }
        }
      }

    } else if (TransF == CblasTrans && TransG == CblasTrans) {

      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          BASE temp_real = 0.0;
          BASE temp_imag = 0.0;
          for (k = 0; k < K; k++) {
            const BASE Fki_real = CONST_REAL(F, ldf * k + i);
            const BASE Fki_imag = conjF * CONST_IMAG(F, ldf * k + i);
            const BASE Gjk_real = CONST_REAL(G, ldg * j + k);
            const BASE Gjk_imag = conjG * CONST_IMAG(G, ldg * j + k);

            temp_real += Fki_real * Gjk_real - Fki_imag * Gjk_imag;
            temp_imag += Fki_real * Gjk_imag + Fki_imag * Gjk_real;
          }
          REAL(C, ldc * i + j) += alpha_real * temp_real - alpha_imag * temp_imag;
          IMAG(C, ldc * i + j) += alpha_real * temp_imag + alpha_imag * temp_real;
        }
      }

    } else {
      BLAS_ERROR("unrecognized operation");
    }
  }
}
