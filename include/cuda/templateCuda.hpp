

#ifndef templateCuda_hpp
#define templateCuda_hpp

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <curand.h>


namespace mx
{
namespace cuda
{

template<typename realT>
struct complex;

template<>
struct complex<float>
{
   typedef cuComplex cudaType;
};

template<>
struct complex<double>
{
   typedef cuDoubleComplex cudaType;
};

template<typename floatT>
cublasStatus_t cublasTscal( cublasHandle_t handle, 
                            int n,
                            const floatT *alpha,
                            floatT *x, 
                            int incx
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
   return cublasSscal( handle, n, alpha, x, incx);
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
   return cublasDscal( handle, n, alpha, x, incx);
}

template<typename floatT>
cublasStatus_t cublasTaxpy( cublasHandle_t handle, ///< [in] handle to the cuBLAS library context.
                            int n,                 ///< [in] scalar used for multiplication.
                            const floatT *alpha,   ///< [in] number of elements in the vector x and y
                            const floatT *x,       ///< [in] vector with n elements. 
                            int incx,              ///< [in] stride between consecutive elements of x
                            floatT *y,             ///< [in] vector with n elements. 
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
   return cublasSaxpy(handle, n, alpha, x, incx, y, incy);
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
   return cublasDaxpy(handle, n, alpha, x, incx, y, incy);
}
   
template<typename dataT>
static __device__ __host__ inline 
dataT elementMul( dataT & a, 
                  dataT & b
                )
{
    return a*b;
}


// Complex multiplication
template<>
__device__ __host__ inline 
cuComplex elementMul<cuComplex>( cuComplex & a, 
                                 cuComplex & b
                               )
{
    cuComplex c;
    
    ((float*) &c)[0] = ((float*) &a)[0] * ((float*) &b)[0] - ((float*) &a)[1] * ((float*) &b)[1];
    ((float*) &c)[1] = ((float*) &a)[0] * ((float*) &b)[1] + ((float*) &a)[1] * ((float*) &b)[0];
    return c;

    
}

// Complex multiplication
template<>
__device__ __host__ inline 
cuDoubleComplex elementMul<cuDoubleComplex>( cuDoubleComplex & a, 
                                             cuDoubleComplex & b
                                           )
{
    cuDoubleComplex c;
    
    ((double*) &c)[0] = ((double*) &a)[0] * ((double*) &b)[0] - ((double*) &a)[1] * ((double*) &b)[1];
    ((double*) &c)[1] = ((double*) &a)[0] * ((double*) &b)[1] + ((double*) &a)[1] * ((double*) &b)[0];
    return c;

    
}

// Complex pointwise multiplicationtemplate<>template<>

#if 0
template<typename dataT>
static __global__ 
void scalarMul(dataT *a, dataT * b, int size)
{
   #ifdef __CUDACC__

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = threadID; i < size; i += numThreads)
    {
        a[i] = elementMul<cuComplex>( a[i],  b[0]);
    }
    
    #endif //__CUDACC__
}
#endif

template<typename dataT>
static __global__ 
void pointwiseMul(dataT *a, dataT *b, int size)
{   
   #ifdef __CUDACC__

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = threadID; i < size; i += numThreads)
    {
        a[i] = elementMul<dataT>( a[i],  b[i]);
    }
    
    #endif //__CUDACC__
}


template<typename inputT, typename outputT>
cufftResult cufftPlan2d( cufftHandle *plan, 
                         int nx, 
                         int ny
                       );

template<>
inline
cufftResult cufftPlan2d<std::complex<float>, std::complex<float>>( cufftHandle *plan, 
                                                                  int nx, 
                                                                  int ny
                                                                )
{
   return ::cufftPlan2d(plan, nx, ny, CUFFT_C2C);
}

template<>
inline
cufftResult cufftPlan2d<cuComplex, cuComplex>( cufftHandle *plan, 
                                               int nx, 
                                               int ny
                                             )
{
   return ::cufftPlan2d(plan, nx, ny, CUFFT_C2C);
}

template<>
inline
cufftResult cufftPlan2d<std::complex<double>, std::complex<double>>( cufftHandle *plan, 
                                                                    int nx, 
                                                                    int ny
                                                                  )
{
   return ::cufftPlan2d(plan, nx, ny, CUFFT_Z2Z);
}

template<>
inline
cufftResult cufftPlan2d<cuDoubleComplex, cuDoubleComplex>( cufftHandle *plan, 
                                                           int nx, 
                                                           int ny
                                                         )
{
   return ::cufftPlan2d(plan, nx, ny, CUFFT_Z2Z);
}

template<typename inputT, typename outputT>
cufftResult cufftExec( cufftHandle plan, 
                       inputT *idata, 
                       inputT *odata, 
                       int direction
                     );

template<>
inline
cufftResult cufftExec<std::complex<float>, std::complex<float>>( cufftHandle plan, 
                                                                 std::complex<float> *idata, 
                                                                 std::complex<float> *odata, 
                                                                 int direction
                                                               )
{
   return cufftExecC2C(plan, (cuComplex *) idata, (cuComplex *) odata, direction);
}

template<>
inline
cufftResult cufftExec<cuComplex, cuComplex>( cufftHandle plan, 
                                             cuComplex *idata, 
                                             cuComplex *odata, 
                                             int direction
                                           )
{
   return cufftExecC2C(plan, idata, odata, direction);
}

template<>
inline
cufftResult cufftExec<std::complex<double>, std::complex<double>>( cufftHandle plan, 
                                                                   std::complex<double> *idata, 
                                                                   std::complex<double> *odata, 
                                                                   int direction
                                                                 )
{
   return cufftExecZ2Z(plan, (cuDoubleComplex *) idata, (cuDoubleComplex *) odata, direction);
}

template<>
inline
cufftResult cufftExec<cuDoubleComplex, cuDoubleComplex>( cufftHandle plan, 
                                                         cuDoubleComplex *idata, 
                                                         cuDoubleComplex *odata, 
                                                         int direction
                                                       )
{
   return cufftExecZ2Z(plan, idata, odata, direction);
}

template<typename realT>
curandStatus_t curandGenerateNormal( curandGenerator_t generator, 
                                     realT *outputPtr, 
                                     size_t n, 
                                     realT mean, 
                                     realT stddev
                                   );

template<>
inline
curandStatus_t curandGenerateNormal<float>( curandGenerator_t generator, 
                                            float *outputPtr, 
                                            size_t n, 
                                            float mean, 
                                            float stddev
                                          )
{
   return curandGenerateNormal(generator, outputPtr, n, mean, stddev);
}

template<>
inline
curandStatus_t curandGenerateNormal<double>( curandGenerator_t generator, 
                                             double *outputPtr, 
                                             size_t n, 
                                             double mean, 
                                             double stddev
                                           )
{
   return curandGenerateNormalDouble(generator, outputPtr, n, mean, stddev);
}

}//namespace cuda 
}//namespace mx
#endif // templateCudaPtr_hpp
