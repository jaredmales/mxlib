/** \file templateBLAS.cpp
 * \brief Implementation of templatized wrappers for the BLAS
 * \ingroup gen_math_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
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

#include "math/templateBLAS.hpp"

namespace mx
{
namespace math
{

template <>
void scal<float>( const int N, const float &alpha, float *X, const int incX )
{
    cblas_sscal( N, alpha, X, incX );
}

template <>
void scal<double>( const int N, const double &alpha, double *X, const int incX )
{
    cblas_dscal( N, alpha, X, incX );
}

template <>
void scal<std::complex<float>>( const int N, const std::complex<float> &alpha, std::complex<float> *X, const int incX )
{
    cblas_cscal( N, &alpha, X, incX );
}

template <>
void scal<std::complex<double>>( const int N,
                                 const std::complex<double> &alpha,
                                 std::complex<double> *X,
                                 const int incX )
{
    cblas_zscal( N, &alpha, X, incX );
}

template <>
void gemm<float>( const CBLAS_ORDER Order,
                  const CBLAS_TRANSPOSE TransA,
                  const CBLAS_TRANSPOSE TransB,
                  const int M,
                  const int N,
                  const int K,
                  const float &alpha,
                  const float *A,
                  const int lda,
                  const float *B,
                  const int ldb,
                  const float &beta,
                  float *C,
                  const int ldc )
{
    cblas_sgemm( Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc );
}

template <>
void gemm<double>( const CBLAS_ORDER Order,
                   const CBLAS_TRANSPOSE TransA,
                   const CBLAS_TRANSPOSE TransB,
                   const int M,
                   const int N,
                   const int K,
                   const double &alpha,
                   const double *A,
                   const int lda,
                   const double *B,
                   const int ldb,
                   const double &beta,
                   double *C,
                   const int ldc )
{
    cblas_dgemm( Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc );
}

template <>
void gemm<std::complex<float>>( const CBLAS_ORDER Order,
                                const CBLAS_TRANSPOSE TransA,
                                const CBLAS_TRANSPOSE TransB,
                                const int M,
                                const int N,
                                const int K,
                                const std::complex<float> &alpha,
                                const std::complex<float> *A,
                                const int lda,
                                const std::complex<float> *B,
                                const int ldb,
                                const std::complex<float> &beta,
                                std::complex<float> *C,
                                const int ldc )
{
    cblas_cgemm( Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc );
}

template <>
void gemm<std::complex<double>>( const CBLAS_ORDER Order,
                                 const CBLAS_TRANSPOSE TransA,
                                 const CBLAS_TRANSPOSE TransB,
                                 const int M,
                                 const int N,
                                 const int K,
                                 const std::complex<double> &alpha,
                                 const std::complex<double> *A,
                                 const int lda,
                                 const std::complex<double> *B,
                                 const int ldb,
                                 const std::complex<double> &beta,
                                 std::complex<double> *C,
                                 const int ldc )
{
    cblas_zgemm( Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc );
}

template <>
void syrk<float>( const CBLAS_ORDER Order,
                  const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans,
                  const int N,
                  const int K,
                  const float &alpha,
                  const float *A,
                  const int lda,
                  const float &beta,
                  float *C,
                  const int ldc )
{
    cblas_ssyrk( Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc );
}

template <>
void syrk<double>( const CBLAS_ORDER Order,
                   const CBLAS_UPLO Uplo,
                   const CBLAS_TRANSPOSE Trans,
                   const int N,
                   const int K,
                   const double &alpha,
                   const double *A,
                   const int lda,
                   const double &beta,
                   double *C,
                   const int ldc )
{
    cblas_dsyrk( Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc );
}

template <>
void syrk<std::complex<float>>( const CBLAS_ORDER Order,
                                const CBLAS_UPLO Uplo,
                                const CBLAS_TRANSPOSE Trans,
                                const int N,
                                const int K,
                                const std::complex<float> &alpha,
                                const std::complex<float> *A,
                                const int lda,
                                const std::complex<float> &beta,
                                std::complex<float> *C,
                                const int ldc )
{
    cblas_csyrk( Order, Uplo, Trans, N, K, &alpha, A, lda, &beta, C, ldc );
}

template <>
void syrk<std::complex<double>>( const CBLAS_ORDER Order,
                                 const CBLAS_UPLO Uplo,
                                 const CBLAS_TRANSPOSE Trans,
                                 const int N,
                                 const int K,
                                 const std::complex<double> &alpha,
                                 const std::complex<double> *A,
                                 const int lda,
                                 const std::complex<double> &beta,
                                 std::complex<double> *C,
                                 const int ldc )
{
    cblas_zsyrk( Order, Uplo, Trans, N, K, &alpha, A, lda, &beta, C, ldc );
}

} // namespace math
} // namespace mx
