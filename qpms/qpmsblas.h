// If included together with <cblas.h>, this must be include _afterwards_ because of the typedefs!
#ifndef QPMSBLAS_H
#define QPMSBLAS_H
#define QPMS_BLAS_INDEX_T long long int

#ifndef CBLAS_H
typedef enum {CblasRowMajor=101, CblasColMajor=102} CBLAS_LAYOUT;
typedef enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113} CBLAS_TRANSPOSE;
typedef enum {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
typedef enum {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
typedef enum {CblasLeft=141, CblasRight=142} CBLAS_SIDE;
#endif

void qpms_zgemm(CBLAS_LAYOUT Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
		const QPMS_BLAS_INDEX_T M, const QPMS_BLAS_INDEX_T N, const QPMS_BLAS_INDEX_T K,
		const _Complex double *alpha, const _Complex double *A, const QPMS_BLAS_INDEX_T lda,
		const _Complex double *B, const QPMS_BLAS_INDEX_T ldb,
		const _Complex double *beta, _Complex double *C, const QPMS_BLAS_INDEX_T ldc);

#endif //QPMSBLAS_H
