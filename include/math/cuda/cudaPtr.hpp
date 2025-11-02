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
#include <cstdint>

#include <cuda_runtime.h>

#include "../../mxlib.hpp"

#include "templateCuda.hpp"

namespace mx
{
namespace cuda
{

/// A smart-pointer wrapper for cuda device pointers.
/**
 * \ingroup cuda
 */
template <typename T, class verboseT = verbose::d>
struct cudaPtr
{

    /// The host data type.
    // typedef typename cudaType<T>::hostType hostPtrT;
    typedef T hostPtrT;

    /// The device data type
    // typedef typename cudaType<T>::deviceType devicePtrT;
    typedef T devicePtrT;

    /// The device pointer
    devicePtrT *m_devicePtr{ nullptr };

    /// The allocated size
    size_t m_size{ 0 };

    /// The number of rows set on allocation
    uint32_t m_rows {0};

    /// The number of columns set on allocation
    uint32_t m_cols {0};

    /// The number of planes set on allocation
    uint32_t m_planes {0};

    /// Destructor, frees memory if allocated.
    ~cudaPtr();

    size_t size()
    {
        return m_size;
    }

    uint32_t rows()
    {
        return m_rows;
    }

    uint32_t cols()
    {
        return m_cols;
    }

    uint32_t planes()
    {
        return m_planes;
    }

    private:
    /// Resize the memory allocation, in 1D
    /** If no size change, this is a no-op.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t resizeImpl( size_t sz /**< [in] the new size */ );

    public:
    /// Resize the memory allocation, in 1D
    /** If no size change, this is a no-op.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t resize( size_t sz /**< [in] the new size */ );

    /// Resize the memory allocation, in 2D
    /** If no size change, this is a no-op.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t resize( uint32_t x_sz, ///< [in] the new x size,
                uint32_t y_sz  ///< [in] the new y size
    );

    /// Resize the memory allocation, in 3D
    /** If no size change, this is a no-op.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t resize( uint32_t x_sz, ///< [in] the new x size,
                uint32_t y_sz, ///< [in] the new y size,
                uint32_t z_sz  ///< [in] the new z size
    );

    /// Initialize the array bytes to 0.
    /** Same as setZero, just a wrapper to cudaMemset.
     *
     */
    error_t initialize();

    /// Initialize the array bytes to 0.
    /** Same as initialize, just a wrapper to cudaMemset.
     *
     */
    error_t setZero();

    /// Free the memory allocation
    /**
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t free();

    /// Copy from the host to the device, after allocation.
    /**
     * The device pointer must be allocated.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t upload( const hostPtrT *src /**< [in] The host location */ );

    /// Copy from the host to the device with 1D allocation.
    /**
     * The device pointer will be re-allocated as needed.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t upload( const hostPtrT *src, ///< [in] The host location
                size_t x_sz            ///< [in] The x size of the array
    );

    /// Copy from the host to the device with 2D allocation.
    /**
     * The device pointer will be re-allocated as needed.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t upload( const hostPtrT *src, ///< [in] The host location
                uint32_t x_sz,            ///< [in] The x size of the array
                uint32_t y_sz            ///< [in] The x size of the array
    );

    /// Copy from the host to the device with #D allocation.
    /**
     * The device pointer will be re-allocated as needed.
     *
     * \returns 0 on success.
     * \returns a cuda error code otherwise.
     *
     */
    error_t upload( const hostPtrT *src, ///< [in] The host location
                uint32_t x_sz,            ///< [in] The x size of the array
                uint32_t y_sz,            ///< [in] The x size of the array
                uint32_t z_sz            ///< [in] The x size of the array
    );

    /// Copy from the device to the host.
    /**
     *
     */
    error_t download( hostPtrT *dest /**< [in] The host location, allocated.*/ );

    /// Accesses the device pointer for use in Cuda functions.
    /**
     */
    typename cpp2cudaType<devicePtrT>::cudaType *data()
    {
        return reinterpret_cast<typename cpp2cudaType<devicePtrT>::cudaType *>(m_devicePtr);
    }

    /// Conversion operator, accesses the device pointer for use in Cuda functions.
    /**
     */
    typename cpp2cudaType<devicePtrT>::cudaType *operator()()
    {
       return reinterpret_cast<typename cpp2cudaType<devicePtrT>::cudaType *>(m_devicePtr);
    }

    /// Conversion operator, accesses the device pointer for use in Cfuda functions.
    /**
     *
     */
    const typename cpp2cudaType<devicePtrT>::cudaType *operator()() const
    {
        return reinterpret_cast<typename cpp2cudaType<devicePtrT>::cudaType *>(m_devicePtr);
    }
};

template <typename T, class verboseT>
cudaPtr<T, verboseT>::~cudaPtr()
{
    free();
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::resizeImpl( size_t sz )
{
    if( m_size == sz )
    {
        return error_t::noerror;
    }

    m_size = sz;

    cudaError_t rv = cudaMalloc( (void **)&m_devicePtr, sz * sizeof( devicePtrT ) );

    if( rv != cudaSuccess )
    {
        return internal::mxlib_error_report<verboseT>(cudaError2error_t(rv), "cudaMalloc");
    }

    return error_t::noerror;
}


template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::resize( size_t sz )
{
    m_rows =sz;
    m_cols = 1;
    m_planes = 1;
    return resizeImpl(sz);
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::resize( uint32_t x_sz, uint32_t y_sz )
{
    m_rows = x_sz;
    m_cols = y_sz;
    m_planes = 1;
    return resizeImpl( x_sz * y_sz );
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::resize( uint32_t x_sz, uint32_t y_sz, uint32_t z_sz )
{
    m_rows = x_sz;
    m_cols = y_sz;
    m_planes = z_sz;
    return resize( x_sz * y_sz * z_sz );
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::initialize()
{
    cudaError_t rv = ::cudaMemset( m_devicePtr, 0, m_size * sizeof( devicePtrT ) );

    if(rv != cudaSuccess)
    {
        return internal::mxlib_error_report<verboseT>(cudaError2error_t(rv), "cudaMemset");
    }

    return error_t::noerror;
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::setZero()
{
    return initialize();
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::free()
{
    if( m_devicePtr )
    {
        cudaError_t rv = cudaFree( m_devicePtr );

        if( rv != cudaSuccess )
        {
            return internal::mxlib_error_report<verboseT>(cudaError2error_t(rv), "cudaFree");
        }
        return error_t::noerror;
    }

    m_devicePtr = 0;
    m_size = 0;

    return error_t::noerror;
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::upload( const hostPtrT *src )
{
    // Copy host memory to device
    cudaError_t rv = cudaMemcpy( m_devicePtr, src, m_size * sizeof( devicePtrT ), cudaMemcpyHostToDevice );

    if( rv != cudaSuccess )
    {
        return internal::mxlib_error_report<verboseT>(cudaError2error_t(rv), "cudaMemcpy");
    }

    return error_t::noerror;
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::upload( const hostPtrT *src, size_t x_sz )
{
    error_t rv;

    rv = resize( x_sz );

    if( !!rv )
    {
        return internal::mxlib_error_report<verboseT>(rv);
    }

    return upload( src );
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::upload( const hostPtrT *src, uint32_t x_sz, uint32_t y_sz )
{
    error_t rv;

    rv = resize( x_sz, y_sz );

    if( !!rv )
    {
        return internal::mxlib_error_report<verboseT>(rv);
    }

    return upload( src );
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::upload( const hostPtrT *src, uint32_t x_sz, uint32_t y_sz, uint32_t z_sz )
{
    error_t rv;

    rv = resize( x_sz, y_sz, z_sz );

    if( !!rv )
    {
        return internal::mxlib_error_report<verboseT>(rv);
    }

    return upload( src );
}

template <typename T, class verboseT>
error_t cudaPtr<T, verboseT>::download( hostPtrT *dest )
{
    // Copy device memory to host
    cudaError_t rv = cudaMemcpy( dest, m_devicePtr, m_size * sizeof( devicePtrT ), cudaMemcpyDeviceToHost );

    if( rv != cudaSuccess )
    {
        return internal::mxlib_error_report<verboseT>(cudaError2error_t(rv), "cudaMemcpy");
    }

    return error_t::noerror;
}

} // namespace cuda
} // namespace mx
#endif // math_cudaPtr_hpp
