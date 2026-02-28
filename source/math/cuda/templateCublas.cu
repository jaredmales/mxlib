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

template <>
cublasStatus_t cublasTscal<float>( cublasHandle_t handle, int n, const float *alpha, float *x, int incx )
{
    return ::cublasSscal( handle, n, alpha, x, incx );
}

template <>
cublasStatus_t cublasTscal<double>( cublasHandle_t handle, int n, const double *alpha, double *x, int incx )
{
    return ::cublasDscal( handle, n, alpha, x, incx );
}

template <>
cublasStatus_t cublasTscal<cuComplex>( cublasHandle_t handle, int n, const cuComplex *alpha, cuComplex *x, int incx )
{
    return ::cublasCscal( handle, n, alpha, x, incx );
}

template <>
cublasStatus_t
cublasTscal<cuDoubleComplex>( cublasHandle_t handle, int n, const cuDoubleComplex *alpha, cuDoubleComplex *x, int incx )
{
    return ::cublasZscal( handle, n, alpha, x, incx );
}

//----------------------------------------------------
// Taxpy

template <>
cublasStatus_t
cublasTaxpy<float>( cublasHandle_t handle, int n, const float *alpha, const float *x, int incx, float *y, int incy )
{
    return ::cublasSaxpy( handle, n, alpha, x, incx, y, incy );
}

template <>
cublasStatus_t
cublasTaxpy<double>( cublasHandle_t handle, int n, const double *alpha, const double *x, int incx, double *y, int incy )
{
    return ::cublasDaxpy( handle, n, alpha, x, incx, y, incy );
}

template <>
cublasStatus_t cublasTaxpy<cuComplex>(
    cublasHandle_t handle, int n, const cuComplex *alpha, const cuComplex *x, int incx, cuComplex *y, int incy )
{
    return ::cublasCaxpy( handle, n, alpha, x, incx, y, incy );
}

template <>
cublasStatus_t cublasTaxpy<cuDoubleComplex>( cublasHandle_t handle,
                                             int n,
                                             const cuDoubleComplex *alpha,
                                             const cuDoubleComplex *x,
                                             int incx,
                                             cuDoubleComplex *y,
                                             int incy )
{
    return ::cublasZaxpy( handle, n, alpha, x, incx, y, incy );
}

//----------------------------------------------------
// Addition, with overloads for complex

template <typename dataT1, typename dataT2>
__device__ dataT1 elementAdd( dataT1 &a, dataT2 &b )
{
    return a + b;
}

// complex-float by complex-float addition
template<>
__device__ cuComplex elementAdd<cuComplex, cuComplex>( cuComplex &a, cuComplex &b )
{
    cuComplex c;

    ( (float *)&c )[0] = ( (float *)&a )[0] + ( (float *)&b )[0];
    ( (float *)&c )[1] = ( (float *)&a )[1] * ( (float *)&b )[1];

    return c;
}

// complex-float by scalar addition
template <>
__device__ cuComplex elementAdd<cuComplex, float>( cuComplex &a, float &b )
{
    cuComplex c;

    ( (float *)&c )[0] = ( (float *)&a )[0] + b;
    ( (float *)&c )[1] = ( (float *)&a )[1];

    return c;
}

// complex-double by complex-double addition
template <>
__device__ cuDoubleComplex elementAdd<cuDoubleComplex, cuDoubleComplex>( cuDoubleComplex &a, cuDoubleComplex &b )
{
    cuDoubleComplex c;

    ( (double *)&c )[0] = ( (double *)&a )[0] + ( (double *)&b )[0];
    ( (double *)&c )[1] = ( (double *)&a )[1] + ( (double *)&b )[1];

    return c;
}

// complex-double by real-double addition
template<>
__device__ cuDoubleComplex elementAdd<cuDoubleComplex, double>( cuDoubleComplex &a, double &b )
{
    cuDoubleComplex c;

    ( (double *)&c )[0] = ( (double *)&a )[0] + b;
    ( (double *)&c )[1] = ( (double *)&a )[1];

    return c;
}

//----------------------------------------------------
// Multiplication with overloads for complex

template <typename dataT1, typename dataT2>
__device__ dataT1 elementMul( dataT1 &a, dataT2 &b )
{
    return a * b;
}

// complex-float by complex-float multiplication
template <>
__device__ cuComplex elementMul<cuComplex, cuComplex>( cuComplex &a, cuComplex &b )
{
    cuComplex c;

    ( (float *)&c )[0] = ( (float *)&a )[0] * ( (float *)&b )[0] - ( (float *)&a )[1] * ( (float *)&b )[1];
    ( (float *)&c )[1] = ( (float *)&a )[0] * ( (float *)&b )[1] + ( (float *)&a )[1] * ( (float *)&b )[0];
    return c;
}

// complex-float by scalar multiplication
template <>
__device__ cuComplex elementMul<cuComplex, float>( cuComplex &a, float &b )
{
    cuComplex c;

    ( (float *)&c )[0] = ( (float *)&a )[0] * b;
    ( (float *)&c )[1] = ( (float *)&a )[1] * b;
    return c;
}

// complex-double by complex-double multiplication
template <>
__device__ cuDoubleComplex elementMul<cuDoubleComplex, cuDoubleComplex>( cuDoubleComplex &a, cuDoubleComplex &b )
{
    cuDoubleComplex c;

    ( (double *)&c )[0] = ( (double *)&a )[0] * ( (double *)&b )[0] - ( (double *)&a )[1] * ( (double *)&b )[1];
    ( (double *)&c )[1] = ( (double *)&a )[0] * ( (double *)&b )[1] + ( (double *)&a )[1] * ( (double *)&b )[0];
    return c;
}

// complex-double by real-double multiplication
template <>
__device__ cuDoubleComplex elementMul<cuDoubleComplex, double>( cuDoubleComplex &a, double &b )
{
    cuDoubleComplex c;

    ( (double *)&c )[0] = ( (double *)&a )[0] * b;
    ( (double *)&c )[1] = ( (double *)&a )[1] * b;

    return c;
}

//----------------------------------------------------
// Element-wise (Hadamard) products of vectors

// Kernel for hadamard product
template <typename dataT0, typename dataT1, typename dataT2>
__global__ void elwiseMul( dataT0 *c, dataT1 *a, dataT2 *b, int size )
{
#ifdef __CUDACC__

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for( int i = threadID; i < size; i += numThreads )
    {
        c[i] = elementMul<dataT1, dataT2>( a[i], b[i] );
    }

#endif //__CUDACC__
}

// Kernel for hadamard product with accumulation
template <typename dataT0, typename dataT1, typename dataT2>
__global__ void elwiseMulAccum( dataT0 *c, dataT1 *a, dataT2 *b, int size )
{
#ifdef __CUDACC__

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for( int i = threadID; i < size; i += numThreads )
    {
        dataT1 ctmp = elementMul<dataT1, dataT2>( a[i], b[i] );

        c[i] = elementAdd<dataT0, dataT2>( c[i] , ctmp );
    }

#endif //__CUDACC__
}

// Calculates the element-wise product of two vectors, storing the result in the first.
/* Calculates z = x * y element by element, a.k.a. the Hadamard product.
 */
template <typename dataT0, typename dataT1, typename dataT2>
cudaError_t elementwiseXxY_impl( dataT0 *z, dataT1 *x, dataT2 *y, int size )
{

    cudaError_t rv = cudaSuccess;

#ifdef __CUDACC__
    rv = cudaGetLastError();
    elwiseMul<dataT0, dataT1, dataT2><<<( size + 255 ) / 256, 256>>>( z, x, y, size );
    rv = cudaGetLastError();
#endif

    return rv;
}

// Calculates the element-wise product of two vectors, accumulating the result in the first.
/* Calculates z += x * y element by element, a.k.a. the Hadamard product.
 *
 */
template <typename dataT0, typename dataT1, typename dataT2>
cudaError_t elementwiseXxYAccum_impl( dataT0 *z, dataT1 *x, dataT2 *y, int size )
{

    cudaError_t rv = cudaSuccess;

#ifdef __CUDACC__
    rv = cudaGetLastError();
    elwiseMulAccum<dataT0, dataT1, dataT2><<<( size + 255 ) / 256, 256>>>( z, x, y, size );
    rv = cudaGetLastError();
#endif

    return rv;
}

template <>
cudaError_t elementwiseXxY<float, float>( float *x, float *y, int size )
{
    return elementwiseXxY_impl<float, float, float>( x, x, y, size );
}

template <>
cudaError_t elementwiseXxY<double, double>( double *x, double *y, int size )
{
    return elementwiseXxY_impl<double, double, double>( x, x, y, size );
}

template <>
cudaError_t elementwiseXxY<cuComplex, float>( cuComplex *x, float *y, int size )
{
    return elementwiseXxY_impl<cuComplex, cuComplex, float>( x, x, y, size );
}

template <>
cudaError_t elementwiseXxY<cuComplex, cuComplex>( cuComplex *x, cuComplex *y, int size )
{
    return elementwiseXxY_impl<cuComplex, cuComplex, cuComplex>( x, x, y, size );
}

template <>
cudaError_t elementwiseXxY<cuDoubleComplex, double>( cuDoubleComplex *x, double *y, int size )
{
    return elementwiseXxY_impl<cuDoubleComplex, cuDoubleComplex, double>( x, x, y, size );
}

template <>
cudaError_t elementwiseXxY<cuDoubleComplex, cuDoubleComplex>( cuDoubleComplex *x, cuDoubleComplex *y, int size )
{
    return elementwiseXxY_impl<cuDoubleComplex, cuDoubleComplex, cuDoubleComplex>( x, x, y, size );
}

template <>
cudaError_t elementwiseXxY<float, float, float>( float *z, float *x, float *y, int size )
{
    return elementwiseXxY_impl<float, float, float>( z, x, y, size );
}

template <>
cudaError_t elementwiseXxY<double, double, double>( double *z, double *x, double *y, int size )
{
    return elementwiseXxY_impl<double, double, double>( z, x, y, size );
}


template <>
cudaError_t elementwiseXxY<cuComplex, cuComplex, float>( cuComplex *z, cuComplex *x, float *y, int size )
{
    return elementwiseXxY_impl<cuComplex, cuComplex, float>( z, x, y, size );
}

template <>
cudaError_t elementwiseXxY<cuComplex, cuComplex, cuComplex>( cuComplex *z, cuComplex *x, cuComplex *y, int size )
{
    return elementwiseXxY_impl<cuComplex, cuComplex, cuComplex>( z, x, y, size );
}

template <>
cudaError_t
elementwiseXxY<cuDoubleComplex, cuDoubleComplex, double>( cuDoubleComplex *z, cuDoubleComplex *x, double *y, int size )
{
    return elementwiseXxY_impl<cuDoubleComplex, cuDoubleComplex, double>( z, x, y, size );
}

template <>
cudaError_t elementwiseXxY<cuDoubleComplex, cuDoubleComplex, cuDoubleComplex>( cuDoubleComplex *z,
                                                                               cuDoubleComplex *x,
                                                                               cuDoubleComplex *y,
                                                                               int size )
{
    return elementwiseXxY_impl<cuDoubleComplex, cuDoubleComplex, cuDoubleComplex>( z, x, y, size );
}

template <>
cudaError_t elementwiseXxYAccum<float, float, float>( float *z, float *x, float *y, int size )
{
    return elementwiseXxYAccum_impl<float, float, float>( z, x, y, size );
}

template <>
cudaError_t elementwiseXxYAccum<double, double, double>( double *z, double *x, double *y, int size )
{
    return elementwiseXxYAccum_impl<double, double, double>( z, x, y, size );
}

/*
template <>
cudaError_t elementwiseXxYAccum<cuComplex, cuComplex, float>( cuComplex *z, cuComplex *x, float *y, int size )
{
    return elementwiseXxYAccum_impl<cuComplex, cuComplex, float>( z, x, y, size );
}

template <>
cudaError_t elementwiseXxYAccum<cuComplex, cuComplex, cuComplex>( cuComplex *z, cuComplex *x, cuComplex *y, int size )
{
    return elementwiseXxYAccum_impl<cuComplex, cuComplex, cuComplex>( z, x, y, size );
}

template <>
cudaError_t elementwiseXxYAccum<cuDoubleComplex, cuDoubleComplex, double>( cuDoubleComplex *z,
                                                                           cuDoubleComplex *x,
                                                                           double *y,
                                                                           int size )
{
    return elementwiseXxYAccum_impl<cuDoubleComplex, cuDoubleComplex, double>( z, x, y, size );
}

template <>
cudaError_t elementwiseXxYAccum<cuDoubleComplex, cuDoubleComplex, cuDoubleComplex>( cuDoubleComplex *z,
                                                                                    cuDoubleComplex *x,
                                                                                    cuDoubleComplex *y,
                                                                                    int size )
{
    return elementwiseXxYAccum_impl<cuDoubleComplex, cuDoubleComplex, cuDoubleComplex>( z, x, y, size );
}

*/

//-----------------------------------

#if 0

template <typename realImageT, typename complexImageT>
__global__ void extractIntensityImageAccum_impl( realImageT &im,    /**<[in] */
                                                 int imX0,          /**<[in] The x-coord of the starting image pixel*/
                                                 int imXsz,         /**<[in] The number of x image pixels*/
                                                 int imY0,          /**<[in] The y-coord of the starting image pixel */
                                                 int imYsz,         /**<[in] The number of y image pixels */
                                                 complexImageT &wf, /**<[in] */
                                                 int wfX0, /**<[in] The x-coord of the starting wavefront pixel*/
                                                 int wfY0  /**<[in] The y-coord of the starting wavefront pixel*/
)
{
    #ifdef __CUDACC__

    int im_rows = im.cols();

    int wf_rows = wf.cols();

    typename realImageT::Scalar *im_data;
    typename complexImageT::Scalar *wf_data;

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for( int j = threadID; j < imXsz; j += numThreads )
    {
        im_data = &im.data()[imX0 + ( imY0 + j ) * im_rows];
        wf_data = &wf.data()[wfX0 + ( wfY0 + j ) * wf_rows];
        for( int i = 0; i < imYsz; ++i )
        {
            im_data[i] += norm( wf_data[i] );
        }
    }

    #endif //__CUDACC__
}

template <typename realImageT, typename complexImageT>
cudaError_t extractIntensityImageAccum( realImageT &im,    /**<[in] */
                                        int imX0,          /**<[in] The x-coord of the starting image pixel*/
                                        int imXsz,         /**<[in] The number of x image pixels*/
                                        int imY0,          /**<[in] The y-coord of the starting image pixel */
                                        int imYsz,         /**<[in] The number of y image pixels */
                                        complexImageT &wf, /**<[in] */
                                        int wfX0,          /**<[in] The x-coord of the starting wavefront pixel*/
                                        int wfY0           /**<[in] The y-coord of the starting wavefront pixel*/
)
{

    cudaError_t rv = cudaSuccess;

    #ifdef __CUDACC__
    rv = cudaGetLastError();
    extractIntensityImageAccum_impl<realImageT, complexImageT>
        <<<( imXsz + 255 ) / 256, 256>>>( im, imX0, imXsz, imY0, imYsz, wf, wfX0, wfY0 );
    rv = cudaGetLastError();
    #endif

    return rv;
}
#endif

//----------------------------------------------------
// Tgemv

template <>
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
                                   int incy )
{
    return ::cublasSgemv( handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy );
}

template <>
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
                                    int incy )
{
    return ::cublasDgemv( handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy );
}

template <>
cublasStatus_t cublasTgemv<float>( cublasHandle_t handle,
                                   cublasOperation_t trans,
                                   int m,
                                   int n,
                                   const float *alpha,
                                   const float *A,
                                   const float *x,
                                   const float *beta,
                                   float *y )
{
    return ::cublasSgemv( handle, trans, m, n, alpha, A, m, x, 1, beta, y, 1 );
}

template <>
cublasStatus_t cublasTgemv<double>( cublasHandle_t handle,
                                    cublasOperation_t trans,
                                    int m,
                                    int n,
                                    const double *alpha,
                                    const double *A,
                                    const double *x,
                                    const double *beta,
                                    double *y )
{
    return ::cublasDgemv( handle, trans, m, n, alpha, A, m, x, 1, beta, y, 1 );
}

} // namespace cuda
} // namespace mx
