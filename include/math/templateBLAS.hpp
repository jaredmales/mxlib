/** \file templateBLAS.hpp
  * \brief Declares and defines templatized wrappers for the BLAS
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#include <complex>

extern "C"
{
#ifdef MXLIB_MKL

   #include <mkl.h>

#else

   #include <gsl/gsl_cblas.h>

#endif
}


#ifndef math_templateBLAS_hpp
#define math_templateBLAS_hpp


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
           const dataT & alpha,
           dataT * X,
           const int incX )
{
   //static_assert(0, "templateBLAS: no scal wrapper defined for type dataT");
   return; //No BLAS for this type.
}

template<>
void scal<float>( const int N,
                  const float & alpha,
                  float * X,
                  const int incX );

template<>
void scal<double>( const int N,
                   const double & alpha,
                   double * X,
                   const int incX 
                 );

template<>
void scal<std::complex<float> >( const int N,
                                 const std::complex<float>  & alpha,
                                 std::complex<float>  * X,
                                 const int incX 
                               );

template<>
void scal<std::complex<double> >( const int N,
                                  const std::complex<double> & alpha,
                                  std::complex<double>  * X,
                                  const int incX 
                                );

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
           const int incX 
         )
{
   hadd_impl(N, alpha, Y, incY, X, incX);
}

/// Template Wrapper for cblas xGEMM
/** 
  *
  * \ingroup template_blas
  */
template<typename dataT>
void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
          const CBLAS_TRANSPOSE TransB, const int M, const int N,
          const int K, const dataT & alpha, const dataT *A,
          const int lda, const dataT *B, const int ldb,
          const dataT & beta, dataT *C, const int ldc)
{
   //static_assert(0, "templateBLAS: no gemm wrapper defined for type dataT");
   return; //No BLAS for this type.
}



template<>
void gemm<float>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float & alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float & beta, float *C, const int ldc
                );

template<>
void gemm<double>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                  const CBLAS_TRANSPOSE TransB, const int M, const int N,
                  const int K, const double & alpha, const double *A,
                  const int lda, const double *B, const int ldb,
                  const double & beta, double *C, const int ldc
                 );

template<>
void gemm<std::complex<float> >(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                                const CBLAS_TRANSPOSE TransB, const int M, const int N,
                                const int K, const std::complex<float> & alpha, const std::complex<float> *A,
                                const int lda, const std::complex<float> *B, const int ldb,
                                const std::complex<float> & beta, std::complex<float> *C, const int ldc
                               );

template<>
void gemm<std::complex<double> >(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                                 const CBLAS_TRANSPOSE TransB, const int M, const int N,
                                 const int K, const std::complex<double> & alpha, const std::complex<double> *A,
                                 const int lda, const std::complex<double> *B, const int ldb,
                                 const std::complex<double> & beta, std::complex<double> *C, const int ldc
                                );

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
          const dataT & alpha, const dataT *A, const int lda,
          const dataT & beta, dataT *C, const int ldc)
{
   //static_assert(0, "templateBLAS: no syrk wrapper defined for type dataT");
   return; //No BLAS for this time.
}

template<>
void syrk<float>(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float & alpha, const float *A, const int lda,
                 const float & beta, float *C, const int ldc);

template<>
void syrk<double>(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double & alpha, const double *A, const int lda,
                  const double & beta, double *C, const int ldc
                 );

template<>
void syrk<std::complex<float> >(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                                const CBLAS_TRANSPOSE Trans, const int N, const int K,
                                const std::complex<float> & alpha, const std::complex<float> *A, const int lda,
                                const std::complex<float> & beta, std::complex<float> *C, const int ldc
                               );

template<>
void syrk<std::complex<double> >(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                                 const CBLAS_TRANSPOSE Trans, const int N, const int K,
                                 const std::complex<double> & alpha, const std::complex<double> *A, const int lda,
                                 const std::complex<double> & beta, std::complex<double> *C, const int ldc
                                );

} //namespace math
} //namespace mx

#endif //math_templateBLAS_hpp


