

#ifndef templateCublas_hpp
#define templateCublas_hpp

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
  */ 
template<typename floatT>
cublasStatus_t cublasTscal( cublasHandle_t handle, ///< [in] The cublas context handle
                            int n,                 ///< [in] Number of elements in the vector
                            const floatT *alpha,   ///< [in] The scalar
                            floatT *x,             ///< [in/out] The vector of length n
                            int incx               ///< [in] The stride of the vector
                          );

template<>
inline
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
inline
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
inline
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
inline
cublasStatus_t cublasTscal<cuDoubleComplex>( cublasHandle_t handle, 
                                             int n,
                                             const cuDoubleComplex *alpha,
                                             cuDoubleComplex *x, 
                                             int incx
                                           )
{
   return ::cublasZscal( handle, n, alpha, x, incx);
}

/// Multiplies a vector by a scalar, adding it to a second vector which is overwritten by the result.
/** Implements  
  * \f[
  * \vec{y} = \alpha \vec{x} + \vec{y}
  * \f]
  * 
  * Specializations are provided for float, double, complex-float, and complex-double
  *
  * \tparam floatT a floating-point type, either float, double, complex-float, or complex-double
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

template<>
inline
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
inline
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
inline
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
inline
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

/// Perform a matrix-vector multiplication.
/** Implements  
  * \f[
  * \vec{y} = \alpha \bm{A} \vec{x} + \beta \vec{y}
  * \f]
  * 
  * Specializations are provided for float, double, complex-float, and complex-double
  *
  * \tparam floatT a floating-point type, either float, double, complex-float, or complex-double
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
   

template<>
inline
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
inline
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

/// Perform a matrix-vector multiplication for stride-less arrays
/** Implements  
  * \f[
  * \vec{y} = \alpha \bm{A} \vec{x} + \beta \vec{y}
  * \f]
  * 
  * Specializations are provided for float, double, complex-float, and complex-double
  *
  * \tparam floatT a floating-point type, either float, double, complex-float, or complex-double
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
inline
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
inline
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
#endif // templateCublas_hpp
