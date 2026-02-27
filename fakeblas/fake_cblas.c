#include <stdio.h>
#include "cblas.h"

// Используем blasint, если он определен в cblas.h, иначе int
#ifndef blasint
#define blasint int
#endif

/* --- GEMV --- */
void cblas_sgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, 
                 const blasint M, const blasint N, const float alpha, const float *A, 
                 const blasint lda, const float *X, const blasint incX, const float beta, 
                 float *Y, const blasint incY) {
    // ОШИБКА: Ничего не делаем, тест завалится
}

void cblas_dgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, 
                 const blasint M, const blasint N, const double alpha, const double *A, 
                 const blasint lda, const double *X, const blasint incX, const double beta, 
                 double *Y, const blasint incY) {
    // ОШИБКА: Заполняем результат мусором
    if(M > 0) Y[0] = -999.0;
}

void cblas_cgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, 
                 const blasint M, const blasint N, const void *alpha, const void *A, 
                 const blasint lda, const void *X, const blasint incX, const void *beta, 
                 void *Y, const blasint incY) {}

void cblas_zgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, 
                 const blasint M, const blasint N, const void *alpha, const void *A, 
                 const blasint lda, const void *X, const blasint incX, const void *beta, 
                 void *Y, const blasint incY) {}

/* --- SYMV / HEMV --- */
void cblas_ssymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const blasint N, const float alpha, const float *A, const blasint lda,
                 const float *X, const blasint incX, const float beta, float *Y, const blasint incY) {}

void cblas_dsymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const blasint N, const double alpha, const double *A, const blasint lda,
                 const double *X, const blasint incX, const double beta, double *Y, const blasint incY) {}

void cblas_chemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const blasint N, const void *alpha, const void *A, const blasint lda,
                 const void *X, const blasint incX, const void *beta, void *Y, const blasint incY) {}

// ИСПРАВЛЕНО: Добавлен пропущенный аргумент lda
void cblas_zhemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const blasint N, const void *alpha, const void *A, const blasint lda,
                 const void *X, const blasint incX, const void *beta, void *Y, const blasint incY) {
    // Пустая реализация приведет к FAIL в тестах
}

/* --- TRMV / TRSV --- */
void cblas_strmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const float *A, const blasint lda, float *X, const blasint incX) {}

void cblas_dtrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const double *A, const blasint lda, double *X, const blasint incX) {}

void cblas_ctrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const void *A, const blasint lda, void *X, const blasint incX) {}

void cblas_ztrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const void *A, const blasint lda, void *X, const blasint incX) {}

void cblas_strsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const float *A, const blasint lda, float *X, const blasint incX) {}

void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const double *A, const blasint lda, double *X, const blasint incX) {}

/* --- GER --- */
void cblas_sger(const enum CBLAS_ORDER order, const blasint M, const blasint N,
                const float alpha, const float *X, const blasint incX,
                const float *Y, const blasint incY, float *A, const blasint lda) {}

void cblas_dger(const enum CBLAS_ORDER order, const blasint M, const blasint N,
                const double alpha, const double *X, const blasint incX,
                const double *Y, const blasint incY, double *A, const blasint lda) {}

void cblas_cgeru(const enum CBLAS_ORDER order, const blasint M, const blasint N,
                 const void *alpha, const void *X, const blasint incX,
                 const void *Y, const blasint incY, void *A, const blasint lda) {}

void cblas_zgeru(const enum CBLAS_ORDER order, const blasint M, const blasint N,
                 const void *alpha, const void *X, const blasint incX,
                 const void *Y, const blasint incY, void *A, const blasint lda) {}

void cblas_cgerc(const enum CBLAS_ORDER order, const blasint M, const blasint N,
                 const void *alpha, const void *X, const blasint incX,
                 const void *Y, const blasint incY, void *A, const blasint lda) {}

void cblas_zgerc(const enum CBLAS_ORDER order, const blasint M, const blasint N,
                 const void *alpha, const void *X, const blasint incX,
                 const void *Y, const blasint incY, void *A, const blasint lda) {}

/* --- SYR / HER --- */
void cblas_ssyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const blasint N, const float alpha, const float *X,
                const blasint incX, float *A, const blasint lda) {}

void cblas_dsyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const blasint N, const double alpha, const double *X,
                const blasint incX, double *A, const blasint lda) {}

void cblas_cher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const blasint N, const float alpha, const void *X, const blasint incX,
                void *A, const blasint lda) {}

void cblas_zher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const blasint N, const double alpha, const void *X, const blasint incX,
                void *A, const blasint lda) {}

/* --- SYR2 / HER2 --- */
void cblas_ssyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const blasint N, const float alpha, const float *X,
                const blasint incX, const float *Y, const blasint incY, float *A,
                const blasint lda) {}

void cblas_dsyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const blasint N, const double alpha, const double *X,
                const blasint incX, const double *Y, const int incY, double *A,
                const blasint lda) {}

void cblas_cher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N,
                const void *alpha, const void *X, const blasint incX,
                const void *Y, const blasint incY, void *A, const blasint lda) {}

void cblas_zher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N,
                const void *alpha, const void *X, const blasint incX,
                const void *Y, const blasint incY, void *A, const blasint lda) {}