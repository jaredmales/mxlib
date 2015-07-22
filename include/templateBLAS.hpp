/** \file templateBLAS.hpp
  * \brief Declares and defines templatized wrappers for the BLAS
  * \ingroup template_blas
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#include <complex>

extern "C"
{
#include <cblas.h>
}


#ifndef __templateBLAS_hpp__
#define __templateBLAS_hpp__


/// Template Wrapper for cblas xGEMM
/** 
  *
  * 
  */
template<typename dataT>
void gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
          const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
          const int K, const dataT alpha, const dataT *A,
          const int lda, const dataT *B, const int ldb,
          const dataT beta, dataT *C, const int ldc)
{
   //Throw exception here!
   return; //No BLAS for this time.
}



template<>
inline
void gemm<float>(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc)
{
   cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template<>
inline
void gemm<double>(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                  const int K, const double alpha, const double *A,
                  const int lda, const double *B, const int ldb,
                  const double beta, double *C, const int ldc)
{
   cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


template<>
inline
void gemm<std::complex<float> >(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                                const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                                const int K, const std::complex<float> alpha, const std::complex<float> *A,
                                const int lda, const std::complex<float> *B, const int ldb,
                                const std::complex<float> beta, std::complex<float> *C, const int ldc)
{
   cblas_cgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

template<>
inline
void gemm<std::complex<double> >(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                                 const int K, const std::complex<double> alpha, const std::complex<double> *A,
                                 const int lda, const std::complex<double> *B, const int ldb,
                                 const std::complex<double> beta, std::complex<double> *C, const int ldc)
{
   cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc);
}



/// Template Wrapper for cblas xSYRK
/** 
  *
  * 
  */
template<typename dataT>
void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
          const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
          const dataT alpha, const dataT *A, const int lda,
          const dataT beta, dataT *C, const int ldc)
{
   //Throw exception here!
   return; //No BLAS for this time.
}

template<>
inline
void syrk<float>(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc)
{
   cblas_ssyrk( Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template<>
inline
void syrk<double>(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double beta, double *C, const int ldc)
{
   cblas_dsyrk( Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template<>
inline
void syrk<std::complex<float> >(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                                const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                                const std::complex<float> alpha, const std::complex<float> *A, const int lda,
                                const std::complex<float> beta, std::complex<float> *C, const int ldc)
{
   cblas_csyrk( Order, Uplo, Trans, N, K, &alpha, A, lda, &beta, C, ldc);
}

template<>
inline
void syrk<std::complex<double> >(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                                 const std::complex<double> alpha, const std::complex<double> *A, const int lda,
                                 const std::complex<double> beta, std::complex<double> *C, const int ldc)
{
   cblas_zsyrk( Order, Uplo, Trans, N, K, &alpha, A, lda, &beta, C, ldc);
}

#endif //__templateBLAS_hpp__


