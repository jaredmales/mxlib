/** \file cudaPtr.hpp
  * \author Jared R. Males
  * \brief A wrapper for cuda device pointers
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

#ifndef math_cudaPtr_hpp
#define math_cudaPtr_hpp

#include <iostream>

#include <cuda_runtime.h>

#include "templateCuda.hpp"

namespace mx
{
namespace cuda
{
   
/// A smart-pointer wrapper for cuda device pointers.
/**
  * \ingroup cuda 
  */ 
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
   
   size_t size()
   {
      return m_size;
   }
   
   ///Resize the memory allocation, in 1D
   /** If no size change, this is a no-op.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     * 
     * \test Scenario: scaling a vector with cublas \ref test_math_templateCublas_scal "[test doc]"
     * \test Scenario: scaling and accumulating a vector with cublas \ref test_math_templateCublas_axpy "[test doc]"
     * \test Scenario: multiplying two vectors element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
     */
   int resize( size_t sz /**< [in] the new size */);
   
   ///Resize the memory allocation, in 2D
   /** If no size change, this is a no-op.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     * 
     */
   int resize( size_t x_sz, ///< [in] the new x size,
               size_t y_sz  ///< [in] the new y size
             );
   
   ///Resize the memory allocation, in 3D
   /** If no size change, this is a no-op.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     * 
     */
   int resize( size_t x_sz, ///< [in] the new x size,
               size_t y_sz, ///< [in] the new y size,
               size_t z_sz  ///< [in] the new z size
             );
   
   /// Initialize the array bytes to 0.
   /** Just a wrapper to cudaMemset.
     *
     */  
   cudaError_t initialize();
   
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
     * \test Scenario: multiplying two vectors element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
     */ 
   int upload( const hostPtrT * src /**< [in] The host location */);
   
   ///Copy from the host to the device with allocation.
   /**
     * The device pointer will be re-allocated as needed.
     * 
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     * 
     * \test Scenario: scaling a vector with cublas \ref test_math_templateCublas_scal "[test doc]"
     * \test Scenario: scaling and accumulating a vector with cublas \ref test_math_templateCublas_axpy "[test doc]"
     * \test Scenario: multiplying two vectors element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
     */
   int upload( const hostPtrT * src, ///< [in] The host location
               size_t sz ///< [in] The size of the array
             );
   
   ///Copy from the device to the host.
   /**
     *
     * \test Scenario: scaling a vector with cublas \ref test_math_templateCublas_scal "[test doc]"
     * \test Scenario: scaling and accumulating a vector with cublas \ref test_math_templateCublas_axpy "[test doc]"
     * \test Scenario: multiplying two vectors element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
     */ 
   int download( hostPtrT * dest /**< [in] The host location, allocated.*/ );
   
   ///Conversion operator, accesses the device pointer for use in Cuda functions.
   /**
     * \test Scenario: scaling and accumulating a vector with cublas \ref test_math_templateCublas_axpy "[test doc]"
     * \test Scenario: multiplying two vectors element by element \ref test_math_templateCublas_elementwiseXxY "[test doc]"
     */ 
   typename cpp2cudaType<devicePtrT>::cudaType* operator()()
   {
      return (typename cpp2cudaType<devicePtrT>::cudaType*) m_devicePtr;
   }
   
   /// Conversion operator, accesses the device pointer for use in Cuda functions.
   /**
     *
     * \test Scenario: scaling a vector with cublas \ref test_math_templateCublas_scal "[test doc]"
     */ 
   const typename cpp2cudaType<devicePtrT>::cudaType* operator()() const
   {
      return (typename cpp2cudaType<devicePtrT>::cudaType*) m_devicePtr;
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
   
   cudaError_t rv = cudaMalloc((void **)&m_devicePtr, sz*sizeof(devicePtrT));
   
   if(rv != cudaSuccess)
   {
      std::cerr << "Error from cudaMalloc: ";
      printf("[%s] %s\n", cudaGetErrorName(rv), cudaGetErrorString(rv));
   }
   
   return 0;
   
}

template<typename T>
int cudaPtr<T>::resize( size_t x_sz,
                        size_t y_sz
                      )
{
   return resize(x_sz*y_sz);
}

template<typename T>
int cudaPtr<T>::resize( size_t x_sz,
                        size_t y_sz,
                        size_t z_sz
                      )
{
   return resize(x_sz*y_sz*z_sz);
}

template<typename T>
cudaError_t cudaPtr<T>::initialize()
{
   return ::cudaMemset(m_devicePtr, 0, m_size*sizeof(devicePtrT));
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
int cudaPtr<T>::upload( const hostPtrT * src )
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
int cudaPtr<T>::upload( const hostPtrT * src,
                        size_t sz
                      )
{
   int rv;
   
   rv = resize(sz);
   
   if(rv) return rv;
   
   return upload(src);
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
#endif // math_cudaPtr_hpp
