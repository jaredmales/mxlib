

#ifndef templateCuda_hpp
#define templateCuda_hpp

#include <cuda_runtime.h>

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
typename std::complex<float> elementMul<std::complex<float>>( std::complex<float> & a, 
                                                              std::complex<float> & b
                                                            )
{
    std::complex<float> c;
    
    ((float*) &c)[0] = ((float*) &a)[0] * ((float*) &b)[0] - ((float*) &a)[1] * ((float*) &b)[1];
    ((float*) &c)[1] = ((float*) &a)[0] * ((float*) &b)[1] + ((float*) &a)[1] * ((float*) &b)[0];
    return c;
}

// Complex pointwise multiplicationtemplate<>template<>


template<typename dataT>
static __global__ 
void pointwiseMul(dataT *a, dataT *b, int size)
{
    //typedef typename cudaType<dataT>::vectorType vT;
   
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = threadID; i < size; i += numThreads)
    {
        //((vT*)a)[i] = elementMul<std::complex<float>>( ((vT*)a)[i],  ((vT*)b)[i]);
        a[i] = elementMul<std::complex<float>>( a[i],  b[i]);
    }
}

}//namespace cuda 
}//namespace mx
#endif // templateCudaPtr_hpp
