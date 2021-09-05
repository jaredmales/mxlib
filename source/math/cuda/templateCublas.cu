/** \file templateCublas.cu
  * \author Jared R. Males
  * \brief Implementation of the template interface to cuBlas
  * \ingroup cuda_files
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#include "math/cuda/templateCublas.hpp"
#include <iostream>

namespace mx
{
namespace cuda
{

//----------------------------------------------------
// Tscal

template<>
cublasStatus_t cublasTscal<float>( cublasHandle_t handle, 
                                   int n,
                                   const float *alpha,
                                   float *x, 
                                   int incx
                                 )
{
   return ::cublasSscal( handle, n, alpha, x, incx);
}

template<>
cublasStatus_t cublasTscal<double>( cublasHandle_t handle, 
                                    int n,
                                    const double *alpha,
                                    double *x, 
                                    int incx
                                  )
{
   return ::cublasDscal( handle, n, alpha, x, incx);
}

template<>
cublasStatus_t cublasTscal<cuComplex>( cublasHandle_t handle, 
                                       int n,
                                       const cuComplex *alpha,
                                       cuComplex *x, 
                                       int incx
                                     )
{
   return ::cublasCscal( handle, n, alpha, x, incx);
}

template<>
cublasStatus_t cublasTscal<cuDoubleComplex>( cublasHandle_t handle, 
                                             int n,
                                             const cuDoubleComplex *alpha,
                                             cuDoubleComplex *x, 
                                             int incx
                                           )
{
   return ::cublasZscal( handle, n, alpha, x, incx);
}

//----------------------------------------------------
// Taxpy

template<>
cublasStatus_t cublasTaxpy<float>( cublasHandle_t handle, 
                                   int n,
                                   const float *alpha,
                                   const float *x, 
                                   int incx,
                                   float *y, 
                                   int incy
                                 )
{
   return ::cublasSaxpy(handle, n, alpha, x, incx, y, incy);
}

template<>
cublasStatus_t cublasTaxpy<double>( cublasHandle_t handle, 
                                    int n,
                                    const double *alpha,
                                    const double *x, 
                                    int incx,
                                    double *y, 
                                    int incy
                                  )
{
   return ::cublasDaxpy(handle, n, alpha, x, incx, y, incy);
}
 
template<>
cublasStatus_t cublasTaxpy<cuComplex>( cublasHandle_t handle, 
                                       int n,
                                       const cuComplex *alpha,
                                       const cuComplex *x, 
                                       int incx,
                                       cuComplex *y, 
                                       int incy
                                     )
{
   return ::cublasCaxpy(handle, n, alpha, x, incx, y, incy);
}

template<>
cublasStatus_t cublasTaxpy<cuDoubleComplex>( cublasHandle_t handle, 
                                             int n,
                                             const cuDoubleComplex *alpha,
                                             const cuDoubleComplex *x, 
                                             int incx,
                                             cuDoubleComplex *y, 
                                             int incy
                                           )
{
   return ::cublasZaxpy(handle, n, alpha, x, incx, y, incy);
}

//----------------------------------------------------
// Element-wise (Hadamard) products of vectors

   
// \test Scenario: multiplying two vector element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
template<typename dataT1, typename dataT2>
__device__
dataT1 elementMul( dataT1 & a, 
                   dataT2 & b
                 )
{
    return a*b;
}

// complex-float by complex-float multiplication
// \test Scenario: multiplying two vector element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
template<>
__device__ 
cuComplex elementMul<cuComplex, cuComplex>( cuComplex & a, 
                                            cuComplex & b
                                          )
{
    cuComplex c;
    
    ((float*) &c)[0] = ((float*) &a)[0] * ((float*) &b)[0] - ((float*) &a)[1] * ((float*) &b)[1];
    ((float*) &c)[1] = ((float*) &a)[0] * ((float*) &b)[1] + ((float*) &a)[1] * ((float*) &b)[0];
    return c;

    
}

// complex-float by scalar multiplication
// \test Scenario: multiplying two vector element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
template<>
__device__
cuComplex elementMul<cuComplex, float>( cuComplex & a, 
                                        float & b
                                      )
{
    cuComplex c;
    
    ((float*) &c)[0] = ((float*) &a)[0] * b; 
    ((float*) &c)[1] = ((float*) &a)[1] * b; 
    return c;

    
}


// complex-double by complex-double multiplication
// \test Scenario: multiplying two vector element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
template<>
__device__ 
cuDoubleComplex elementMul<cuDoubleComplex, cuDoubleComplex>( cuDoubleComplex & a, 
                                                              cuDoubleComplex & b
                                                            )
{
    cuDoubleComplex c;
    
    ((double*) &c)[0] = ((double*) &a)[0] * ((double*) &b)[0] - ((double*) &a)[1] * ((double*) &b)[1];
    ((double*) &c)[1] = ((double*) &a)[0] * ((double*) &b)[1] + ((double*) &a)[1] * ((double*) &b)[0];
    return c;

    
}

// complex-double by real-double multiplication
// \test Scenario: multiplying two vector element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
template<>
__device__
cuDoubleComplex elementMul<cuDoubleComplex, double>( cuDoubleComplex & a, 
                                                     double & b
                                                   )
{
    cuDoubleComplex c;
    
    ((double*) &c)[0] = ((double*) &a)[0] * b; 
    ((double*) &c)[1] = ((double*) &a)[1] * b;
    
    return c;

    
}

// \test Scenario: multiplying two vector element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
template<typename dataT1, typename dataT2>
__global__ 
void elwiseMul(dataT1 *a, dataT2 *b, int size)
{   
   #ifdef __CUDACC__

   const int numThreads = blockDim.x * gridDim.x;
   const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

   for (int i = threadID; i < size; i += numThreads)
   {
       a[i] = elementMul<dataT1, dataT2>( a[i],  b[i]);
   }
    
   #endif //__CUDACC__
}

// Calculates the element-wise product of two vectors, storing the result in the first.
/* Calculates x = x * y element by element, a.k.a. the Hadamard product.
 * \test Scenario: multiplying two vector element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
 */
template<typename dataT1, typename dataT2>
cudaError_t elementwiseXxY_impl( dataT1 * x,
                                 dataT2 * y,
                                 int size
                               )
{

   cudaError_t rv = cudaSuccess;

   #ifdef __CUDACC__
   rv = cudaGetLastError();
   elwiseMul<dataT1,dataT2><<<(size+255)/256, 256>>>( x, y, size);
   rv = cudaGetLastError();
   #endif

   return rv;
}

template<>
cudaError_t elementwiseXxY<float,float>( float * x,
                                         float * y,
                                         int size
                                       )
{
   return elementwiseXxY_impl<float,float>(x,y,size);
}

template<>
cudaError_t elementwiseXxY<double,double>( double * x,
                                           double * y,
                                           int size
                                         )
{
   return elementwiseXxY_impl<double,double>(x,y,size);
}

template<>
cudaError_t elementwiseXxY<cuComplex,float>( cuComplex * x,
                                             float * y,
                                             int size
                                           )
{
   return elementwiseXxY_impl<cuComplex,float>(x,y,size);
}

template<>
cudaError_t elementwiseXxY<cuComplex,cuComplex>( cuComplex * x,
                                                 cuComplex * y,
                                                 int size
                                               )
{
   return elementwiseXxY_impl<cuComplex,cuComplex>(x,y,size);
}

template<>
cudaError_t elementwiseXxY<cuDoubleComplex,double>( cuDoubleComplex * x,
                                                    double * y,
                                                    int size
                                                  )
{
   return elementwiseXxY_impl<cuDoubleComplex,double>(x,y,size);
}

template<>
cudaError_t elementwiseXxY<cuDoubleComplex,cuDoubleComplex>( cuDoubleComplex * x,
                                                             cuDoubleComplex * y,
                                                             int size
                                                           )
{
   return elementwiseXxY_impl<cuDoubleComplex,cuDoubleComplex>(x,y,size);
}

//----------------------------------------------------
// Tgemv

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
                                 )
{
   return ::cublasSgemv(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

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
                                  )
{
   return ::cublasDgemv(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

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
                                 )
{
   return ::cublasSgemv(handle, trans, m, n, alpha, A, m, x, 1, beta, y, 1);
}

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
                                  )
{
   return ::cublasDgemv(handle, trans, m, n, alpha, A, m, x, 1, beta, y, 1);
}

}//namespace cuda 
}//namespace mx
