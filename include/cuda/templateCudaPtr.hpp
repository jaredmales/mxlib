

#ifndef templateCudaPtr_hpp
#define templateCudaPtr_hpp

#include <iostream>

#include <cuda_runtime.h>

namespace mx
{
namespace cuda
{
   
template<typename T>
struct cudaPtr
{
   ///The host data type.
   //typedef typename cudaType<T>::hostType hostPtrT;
   typedef T hostPtrT;
   
   ///The device data type
   //typedef typename cudaType<T>::deviceType devicePtrT;
   typedef T devicePtrT;
   
   ///The device pointer
   devicePtrT * m_devicePtr {nullptr};

   ///The allocated size
   size_t m_size {0};
   
   ///Destructor, frees memory if allocated.
   ~cudaPtr();
   
   ///Resize the memory allocation
   /** If no size change, this is a no-op.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     * 
     */
   int resize( size_t sz /**< [in] the new size */);
   
   ///Free the memory allocation
   /** 
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     * 
     */
   int free();
   
   ///Copy from the host to the device, after allocation.
   /**
     * The device pointer must be allocated.
     * 
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     * 
     */ 
   int upload( hostPtrT * src /**< [in] The host location */);
   
   ///Copy from the host to the device with allocation.
   /**
     * The device pointer will be re-allocated as needed.
     * 
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     * 
     */
   int upload( hostPtrT * src, ///< [in] The host location
               size_t sz ///< [in] The size of the array
             );
   
   ///Copy from the device to the host.
   int download( hostPtrT * dest /**< [in] The host location, allocated.*/ );
   
   ///Conversion operator, accesses the device pointer for use in Cuda functions.
   operator devicePtrT *()
   {
      return m_devicePtr;
   }
   
};

template<typename T>
cudaPtr<T>::~cudaPtr()
{
   free();
}

template<typename T>
int cudaPtr<T>::resize( size_t sz )
{
   if( m_size == sz ) return 0;
   
   m_size = sz;
   
   int rv = cudaMalloc((void **)&m_devicePtr, sz*sizeof(devicePtrT));
   
   if(rv != cudaSuccess)
   {
      std::cerr << "Cuda Malloc Error \n";
      return rv;
   }
   
   return 0;
   
}

template<typename T>
int cudaPtr<T>::free()
{
   if(m_devicePtr)
   {
      int rv = cudaFree(m_devicePtr);
      
      if(rv != cudaSuccess)
      {
         std::cerr << "Cuda Free Error \n";
         return rv;
      }  
   }
   
   m_devicePtr = 0;
   m_size = 0;
   
   return 0;
}

template<typename T>
int cudaPtr<T>::upload( hostPtrT * src )
{
    // Copy host memory to device
   int rv = cudaMemcpy( m_devicePtr, src, m_size*sizeof(devicePtrT), cudaMemcpyHostToDevice);
   
   if(rv != cudaSuccess)
   {
      std::cerr << "Cuda Memcpy error \n";
      return rv;
   }
   
   return 0;
}

template<typename T>
int cudaPtr<T>::upload( hostPtrT * src,
                        size_t sz
                      )
{
   int rv;
   
   rv = resize(sz);
   
   if(rv) return rv;
   
   return upload(src);
   
   return 0;
}

template<typename T>
int cudaPtr<T>::download( hostPtrT * dest )
{
    // Copy device memory to host
   int rv = cudaMemcpy( dest,  m_devicePtr, m_size*sizeof(devicePtrT), cudaMemcpyDeviceToHost);
   
   if(rv != cudaSuccess)
   {
      std::cerr << "Cuda Memcpy error \n";
      return rv;
   }
   
   return 0;
}

}//namespace cuda 
}//namespace mx
#endif // templateCudaPtr_hpp
