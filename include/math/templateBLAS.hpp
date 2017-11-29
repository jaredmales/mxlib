/** \file templateBLAS.hpp
  * \brief Declares and defines templatized wrappers for the BLAS
  * \ingroup template_blas
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#include <complex>

extern "C"
{
#ifdef MXLIB_MKL

   #include <mkl.h>
#elif defined(__APPLE__)
#include <vecLib/cblas.h>
#else

   #include <cblas.h>

#endif
}


#ifndef templateBLAS_hpp
#define templateBLAS_hpp


namespace mx
{
namespace math
{
   
/// Template wrapper for cblas xSCAL
/**
  * \tparam dataT the data type of the alpha, X, and Y
  * 
  * \ingroup template_blas 
  */
template<typename dataT>
void scal( const int N,
           const dataT alpha,
           dataT * X,
           const int incX )
{
   //static_assert(0, "templateBLAS: no scal wrapper defined for type dataT");
   return; //No BLAS for this type.
}

template<>
inline
void scal<float>( const int N,
                  const float alpha,
                  float * X,
                  const int incX )
{
   cblas_sscal(N, alpha, X, incX);
}

template<>
inline
void scal<double>( const int N,
                  const double alpha,
                  double * X,
                  const int incX )
{
   cblas_dscal(N, alpha, X, incX);
}

template<>
inline
void scal<std::complex<float> >( const int N,
                  const std::complex<float>  alpha,
                  std::complex<float>  * X,
                  const int incX )
{
   cblas_cscal(N, &alpha, X, incX);
}


template<>
inline
void scal<std::complex<double> >( const int N,
                  const std::complex<double>  alpha,
                  std::complex<double>  * X,
                  const int incX )
{
   cblas_zscal(N, &alpha, X, incX);
}

///Implementation of the Hadamard (element-wise) product of two vectors 
/** Computes the the Hadamard or element-wise product: X <- alpha*X*Y
  * 
  * \param N [in] the length of the two vectors
  * \param alpha [in] scalar to multiply each element by
  * \param Y [in] vector to perform element-wise multiplication with
  * \param incY [in] in-memory increment or stride for Y
  * \param X [in/out] vector which is multiplied by alpha and element-wise multiplied by Y
  * \param incX [in] in-memeory increment or stride for X
  * 
  * \tparam dataT the data type of the alpha, X, and Y
  * 
  * \ingroup template_blas
  */
template<typename dataT>
inline
void hadp_impl( const int N,
                dataT * __restrict__ Y,
                dataT * __restrict__ X)
{
   dataT *x = (dataT *) __builtin_assume_aligned(X, 16);
   dataT *y = (dataT *) __builtin_assume_aligned(Y, 16);
   
   #pragma omp simd
   for(int i=0; i<N; i++)
   {
      x[i] *= y[i];
   }
}

/// Template wrapper for cblas-extension xHADP
/** Computes the the Hadamard or element-wise product: X \<- alpha*X*Y
  * 
  * \param N [in] the length of the two vectors
  * \param alpha [in] scalar to multiply each element by
  * \param Y [in] vector to perform element-wise multiplication with
  * \param incY [in] in-memory increment or stride for Y
  * \param X [in/out] vector which is multiplied by alpha and element-wise multiplied by Y
  * \param incX [in] in-memeory increment or stride for X
  * 
  * \tparam dataT the data type of the alpha, X, and Y
  * 
  * \ingroup template_blas
  */
template<typename dataT>
inline
void hadp( const int N,
           dataT * Y,
           dataT * X )
{
   hadp_impl(N, Y, X);
}

///Implementation of the Hadamard (element-wise) division of two vectors 
/** Computes the the Hadamard or element-wise product: X \<- alpha*X/Y
  * 
  * \param N [in] the length of the two vectors
  * \param alpha [in] scalar to multiply each element by
  * \param Y [in] vector to perform element-wise division with
  * \param incY [in] in-memory increment or stride for Y
  * \param X [in/out] vector which is multiplied by alpha and element-wise divided by Y
  * \param incX [in] in-memeory increment or stride for X
  * 
  * \tparam dataT the data type of the alpha, X, and Y
  * 
  * \ingroup template_blas
  */
template<typename dataT>
void hadd_impl( const int N,
                const dataT alpha,
                const dataT * Y,
                const int incY,
                dataT * X,
                const int incX )
{
   #pragma omp parallel for
   for(int i=0; i<N; ++i)
   {
      X[i*incX] = alpha*X[i*incX]/Y[i*incY];
   }
}

/// Template wrapper for cblas-extension xHADD
/** Computes the the Hadamard or element-wise product: X <- alpha*X/Y
  * 
  * \param N [in] the length of the two vectors
  * \param alpha [in] scalar to multiply each element by
  * \param Y [in] vector to perform element-wise division with
  * \param incY [in] in-memory increment or stride for Y
  * \param X [in/out] vector which is multiplied by alpha and element-wise divided by Y
  * \param incX [in] in-memeory increment or stride for X
  * 
  * \tparam dataT the data type of the alpha, X, and Y
  * 
  * \ingroup template_blas
  */
template<typename dataT>
void hadd( const int N,
           const dataT alpha,
           const dataT * Y,
           const int incY,
           dataT * X,
           const int incX )
{
   hadd_impl(N, alpha, Y, incY, X, incX);
}

/// Template Wrapper for cblas xGEMM
/** 
  *
  * \ingroup template_blas
  */
template<typename dataT>
/*void gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
          const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
          const int K, const dataT alpha, const dataT *A,
          const int lda, const dataT *B, const int ldb,
          const dataT beta, dataT *C, const int ldc)*/
void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
          const CBLAS_TRANSPOSE TransB, const int M, const int N,
          const int K, const dataT alpha, const dataT *A,
          const int lda, const dataT *B, const int ldb,
          const dataT beta, dataT *C, const int ldc)
{
   //static_assert(0, "templateBLAS: no gemm wrapper defined for type dataT");
   return; //No BLAS for this type.
}



template<>
inline
/*void gemm<float>(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc)*/
void gemm<float>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc)
{
   cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template<>
inline
/*void gemm<double>(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                  const int K, const double alpha, const double *A,
                  const int lda, const double *B, const int ldb,
                  const double beta, double *C, const int ldc)*/
void gemm<double>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                  const CBLAS_TRANSPOSE TransB, const int M, const int N,
                  const int K, const double alpha, const double *A,
                  const int lda, const double *B, const int ldb,
                  const double beta, double *C, const int ldc)
{
   cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


template<>
inline
/*void gemm<std::complex<float> >(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                                const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                                const int K, const std::complex<float> alpha, const std::complex<float> *A,
                                const int lda, const std::complex<float> *B, const int ldb,
                                const std::complex<float> beta, std::complex<float> *C, const int ldc)*/
void gemm<std::complex<float> >(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                                const CBLAS_TRANSPOSE TransB, const int M, const int N,
                                const int K, const std::complex<float> alpha, const std::complex<float> *A,
                                const int lda, const std::complex<float> *B, const int ldb,
                                const std::complex<float> beta, std::complex<float> *C, const int ldc)
{
   cblas_cgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

template<>
inline
/*void gemm<std::complex<double> >(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                                 const int K, const std::complex<double> alpha, const std::complex<double> *A,
                                 const int lda, const std::complex<double> *B, const int ldb,
                                 const std::complex<double> beta, std::complex<double> *C, const int ldc)*/
void gemm<std::complex<double> >(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                                 const CBLAS_TRANSPOSE TransB, const int M, const int N,
                                 const int K, const std::complex<double> alpha, const std::complex<double> *A,
                                 const int lda, const std::complex<double> *B, const int ldb,
                                 const std::complex<double> beta, std::complex<double> *C, const int ldc)
{
   cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc);
}



/// Template Wrapper for cblas xSYRK
/** 
  *
  * \ingroup template_blas
  */
template<typename dataT>
/*void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
          const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
          const dataT alpha, const dataT *A, const int lda,
          const dataT beta, dataT *C, const int ldc)*/
void syrk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
          const CBLAS_TRANSPOSE Trans, const int N, const int K,
          const dataT alpha, const dataT *A, const int lda,
          const dataT beta, dataT *C, const int ldc)
{
   //static_assert(0, "templateBLAS: no syrk wrapper defined for type dataT");
   return; //No BLAS for this time.
}

template<>
inline
/*void syrk<float>(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc)*/
void syrk<float>(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc)
{
   cblas_ssyrk( Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template<>
inline
/*void syrk<double>(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double beta, double *C, const int ldc)*/
void syrk<double>(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double beta, double *C, const int ldc)
{
   cblas_dsyrk( Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template<>
inline
/*void syrk<std::complex<float> >(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                                const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                                const std::complex<float> alpha, const std::complex<float> *A, const int lda,
                                const std::complex<float> beta, std::complex<float> *C, const int ldc)*/
void syrk<std::complex<float> >(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                                const CBLAS_TRANSPOSE Trans, const int N, const int K,
                                const std::complex<float> alpha, const std::complex<float> *A, const int lda,
                                const std::complex<float> beta, std::complex<float> *C, const int ldc)
{
   cblas_csyrk( Order, Uplo, Trans, N, K, &alpha, A, lda, &beta, C, ldc);
}

template<>
inline
/*void syrk<std::complex<double> >(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                                 const std::complex<double> alpha, const std::complex<double> *A, const int lda,
                                 const std::complex<double> beta, std::complex<double> *C, const int ldc)*/
void syrk<std::complex<double> >(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                                 const CBLAS_TRANSPOSE Trans, const int N, const int K,
                                 const std::complex<double> alpha, const std::complex<double> *A, const int lda,
                                 const std::complex<double> beta, std::complex<double> *C, const int ldc)
{
   cblas_zsyrk( Order, Uplo, Trans, N, K, &alpha, A, lda, &beta, C, ldc);
}

} //namespace math
} //namespace mx

#endif //templateBLAS_hpp


