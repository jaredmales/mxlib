/** \file templateCublas.hpp
  * \author Jared R. Males
  * \brief A template interface to cuBlas
  * \ingroup cuda_files
  *
  */

//***********************************************************************//
// Copyright 2019,2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_templateCublas_hpp
#define math_templateCublas_hpp

#include <cuda_runtime.h>
#include <cublas_v2.h>

namespace mx
{
namespace cuda
{

/// Multiplies a vector by a scalar, overwriting the vector with the result.
/** Implements
  * \f[
  * \vec{x} = \alpha \vec{x}
  * \f]
  * 
  * Specializations are provided for float, double, complex-float, and complex-double
  *
  * \tparam floatT a floating-point type, either float, double, complex-float, or complex-double
  * 
  * \test Scenario: scaling a vector with cublas \ref test_math_templateCublas_scal "[test doc]"
  * 
  * \ingroup cublas
  */ 
template<typename floatT>
cublasStatus_t cublasTscal( cublasHandle_t handle, ///< [in] The cublas context handle
                            int n,                 ///< [in] Number of elements in the vector
                            const floatT *alpha,   ///< [in] The scalar
                            floatT *x,             ///< [in/out] The vector of length n
                            int incx               ///< [in] The stride of the vector
                          );

/// Multiplies a vector by a scalar, adding it to a second vector which is overwritten by the result.
/** Implements  
  * \f[
  * \vec{y} = \alpha \vec{x} + \vec{y}
  * \f]
  * 
  * Specializations are provided for float, double, complex-float, and complex-double
  *
  * \tparam floatT a floating-point type, either float, double, complex-float, or complex-double
  * 
  * \test Scenario: scaling and accumulating a vector with cublas \ref test_math_templateCublas_axpy "[test doc]"
  * 
  * \ingroup cublas
  */ 
template<typename floatT>
cublasStatus_t cublasTaxpy( cublasHandle_t handle, ///< [in] handle to the cuBLAS library context.
                            int n,                 ///< [in] scalar used for multiplication.
                            const floatT *alpha,   ///< [in] number of elements in the vector x and y
                            const floatT *x,       ///< [in] vector with n elements. 
                            int incx,              ///< [in] stride between consecutive elements of x
                            floatT *y,             ///< [in/out] vector with n elements. 
                            int incy               ///< [in] stride between consecutive elements of y
                          );

//----------------------------------------------------
// Element-wise (Hadamard) products of vectors


/// Calculates the element-wise product of two vectors, storing the result in the first.
/** Calculates
  * \f$
  * x = x * y
  * \f$
  * element by element, a.k.a. the Hadamard product.
  * 
  * Specializations are provided for:
  * - float,float
  * - complex-float, float
  * - complex-float, complex-float
  * - double, double
  * - complex-double, double
  * - complex-double, complex-double
  * 
  * \test Scenario: multiplying two vectors element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
  * 
  * \ingroup cublas
  */
template<typename dataT1, typename dataT2>
cudaError_t elementwiseXxY( dataT1 * x, ///< [in/out] device pointer for the 1st vector.  Is replaced with the product of the two vectors
                            dataT2 * y, ///< [in] device pointer for the 2nd vector.
                            int size    ///< [in] the number of elements in the vectors.
                          );

//----------------------------------------------------
// Tgemv

/// Perform a matrix-vector multiplication.
/** Implements  
  * \f[
  * \vec{y} = \alpha \mathbf{A} \vec{x} + \beta \vec{y}
  * \f]
  * 
  * Specializations are provided for float, double, complex-float, and complex-double
  *
  * \tparam floatT a floating-point type, either float, double, complex-float, or complex-double
  * 
  * \tests Scenario: multiplying a vectory by a matrix \ref test_math_templateCublas_cublasTgemv_inc "[code doc]"
  * 
  * \ingroup cublas
  */ 
template<typename floatT>
cublasStatus_t cublasTgemv( cublasHandle_t handle,   ///< [in] handle to the cuBLAS library context.
                            cublasOperation_t trans, ///< [in] operation on a, CUBLAS_OP_N for none, and CUBLAS_OP_T for transpose
                            int m,                   ///< [in] rows in matrix A. 
                            int n,                   ///< [in] columns in matrix A. 
                            const floatT *alpha,     ///< [in] scalar used for multiplication of A
                            const floatT *A,         ///< [in] array of dimension lda x n with lda >= max(1,m). The leading m by n part of the array A is multiplied by alpha and x.  Unchanged. 
                            int lda,                 ///< [in] leading dimension of A. lda must be at least max(1,m). 
                            const floatT *x,         ///< [in] vector of at least (1+(n-1)*abs(incx)) elements if transa==CUBLAS_OP_N and at least (1+(m-1)*abs(incx)) elements otherwise. 
                            int incx,                ///< [in] stride of x. 
                            const floatT *beta,      ///< [in] scalar used for multiplication of y, if beta==0 then y does not need to be initialized. 
                            floatT *y,               ///< [in/out] vector of at least (1+(m-1)*abs(incy)) elements if transa==CUBLAS_OP_N and at least (1+(n-1)*abs(incy)) elements otherwise. 
                            int incy                 ///< [in] stride of y
                          );

/// Perform a matrix-vector multiplication for stride-less arrays
/** Implements  
  * \f[
  * \vec{y} = \alpha \mathbf{A} \vec{x} + \beta \vec{y}
  * \f]
  * 
  * Specializations are provided for float, double, complex-float, and complex-double
  *
  * \overload 
  * This version assumes stride is 1 in all arrays.
  * 
  * \tparam floatT a floating-point type, either float, double, complex-float, or complex-double
  * 
  * \ingroup cublas
  */ 
template<typename floatT>
cublasStatus_t cublasTgemv( cublasHandle_t handle,   ///< [in] handle to the cuBLAS library context.
                            cublasOperation_t trans, ///< [in] operation on a, CUBLAS_OP_N for none, and CUBLAS_OP_T for transpose
                            int m,                   ///< [in] rows in matrix A. 
                            int n,                   ///< [in] columns in matrix A. 
                            const floatT *alpha,     ///< [in] scalar used for multiplication of A
                            const floatT *A,         ///< [in] array of dimension m x n.  Unchanged. 
                            const floatT *x,         ///< [in] vector of at least (1+(n-1)*abs(incx)) elements if transa==CUBLAS_OP_N and at least (1+(m-1)*abs(incx)) elements otherwise. 
                            const floatT *beta,      ///< [in] scalar used for multiplication of y, if beta==0 then y does not need to be initialized. 
                            floatT *y                ///< [in/out] vector of at least (1+(m-1)*abs(incy)) elements if transa==CUBLAS_OP_N and at least (1+(n-1)*abs(incy)) elements otherwise. 
                          );

template<>
cublasStatus_t cublasTgemv<float>( cublasHandle_t handle,   
                                   cublasOperation_t trans, 
                                   int m,                   
                                   int n,                   
                                   const float *alpha,     
                                   const float *A,          
                                   int lda,                 
                                   const float *x,         
                                   int incx,                
                                   const float *beta,      
                                   float *y,                
                                   int incy                 
                                 );

template<>
cublasStatus_t cublasTgemv<double>( cublasHandle_t handle,   
                                    cublasOperation_t trans, 
                                    int m,                   
                                    int n,                   
                                    const double *alpha,     
                                    const double *A,          
                                    int lda,                 
                                    const double *x,         
                                    int incx,                
                                    const double *beta,      
                                    double *y,                
                                    int incy                 
                                  );

template<>
cublasStatus_t cublasTgemv<float>( cublasHandle_t handle,   
                                   cublasOperation_t trans, 
                                   int m,                   
                                   int n,                   
                                   const float *alpha,     
                                   const float *A,          
                                   const float *x,         
                                   const float *beta,      
                                   float *y                
                                 );

template<>
cublasStatus_t cublasTgemv<double>( cublasHandle_t handle,   
                                    cublasOperation_t trans, 
                                    int m,                   
                                    int n,                   
                                    const double *alpha,     
                                    const double *A,          
                                    const double *x, 
                                    const double *beta,      
                                    double *y                
                                  );

}//namespace cuda 
}//namespace mx
#endif // math_templateCublas_hpp
