

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


template<typename cppType>
struct cpp2cudaType
{
   typedef cppType cudaType;
};

template<>
struct cpp2cudaType<std::complex<float>>
{
   typedef complex<float>::cudaType cudaType;
};

template<>
struct cpp2cudaType<std::complex<double>>
{
   typedef complex<double>::cudaType cudaType;
};


template<typename dataT1, typename dataT2>
static __device__ __host__ inline 
dataT1 elementMul( dataT1 & a, 
                   dataT2 & b
                 )
{
    return a*b;
}


// complex-float by complex-float multiplication
template<>
__device__ __host__ inline 
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
template<>
__device__ __host__ inline 
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

// complex-double by real-double multiplication
template<>
__device__ __host__ inline 
cuDoubleComplex elementMul<cuDoubleComplex, double>( cuDoubleComplex & a, 
                                                     double & b
                                                   )
{
    cuDoubleComplex c;
    
    ((double*) &c)[0] = ((double*) &a)[0] * b; 
    ((double*) &c)[1] = ((double*) &a)[1] * b;
    
    return c;

    
}

template<typename dataT1, typename dataT2>
static __global__ 
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

/// Calculates the element-wise product of two vectors, storing the result in the first.
/** Calculates
  \f$
  x = x * y
  \f$
  * element by element, a.k.a. the Hadamard product.
  */
template<typename dataT1, typename dataT2>
void elementwiseXxY( dataT1 * x,
                     dataT2 * y,
                     int size
                   )
{
   #ifdef __CUDACC__
   elwiseMul<dataT1,dataT2><<<(size+255)/256, 256>>>( x, y, size);
   #endif
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
   return ::cufftExecC2C(plan, (cuComplex *) idata, (cuComplex *) odata, direction);
}

template<>
inline
cufftResult cufftExec<cuComplex, cuComplex>( cufftHandle plan, 
                                             cuComplex *idata, 
                                             cuComplex *odata, 
                                             int direction
                                           )
{
   return ::cufftExecC2C(plan, idata, odata, direction);
}

template<>
inline
cufftResult cufftExec<std::complex<double>, std::complex<double>>( cufftHandle plan, 
                                                                   std::complex<double> *idata, 
                                                                   std::complex<double> *odata, 
                                                                   int direction
                                                                 )
{
   return ::cufftExecZ2Z(plan, (cuDoubleComplex *) idata, (cuDoubleComplex *) odata, direction);
}

template<>
inline
cufftResult cufftExec<cuDoubleComplex, cuDoubleComplex>( cufftHandle plan, 
                                                         cuDoubleComplex *idata, 
                                                         cuDoubleComplex *odata, 
                                                         int direction
                                                       )
{
   return ::cufftExecZ2Z(plan, idata, odata, direction);
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
   return ::curandGenerateNormal(generator, outputPtr, n, mean, stddev);
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
   return ::curandGenerateNormalDouble(generator, outputPtr, n, mean, stddev);
}

}//namespace cuda 
}//namespace mx
#endif // templateCudaPtr_hpp
