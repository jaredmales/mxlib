

#ifndef templateCuda_hpp
#define templateCuda_hpp

#include <cuda_runtime.h>
#include <cuComplex.h>

namespace mx
{
namespace cuda
{
 
   
   
   
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

// Complex pointwise multiplicationtemplate<>template<>


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

template<typename dataT>
static __global__ 
void pointwiseMul(dataT *a, dataT *b, int size)
{   
   #ifdef __CUDACC__

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = threadID; i < size; i += numThreads)
    {
        a[i] = elementMul<cuComplex>( a[i],  b[i]);
    }
    
    #endif //__CUDACC__
}

}//namespace cuda 
}//namespace mx
#endif // templateCudaPtr_hpp
