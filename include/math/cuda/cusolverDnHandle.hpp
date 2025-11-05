/** \file cusolverDnHandle.hpp
 * \author Jared R. Males
 * \brief Management of a cusolverDn handle
 * \ingroup cuda_files
 *
 */

//***********************************************************************//
// Copyright 2025 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_cusolverDnHandle_hpp
#define math_cusolverDnHandle_hpp

#ifdef MXLIB_CUDA

#include <cuda_runtime.h>
#include <cusolverDn.h>

namespace mx
{
namespace cuda
{

/// Management of a cusolverDn handle
/** RAII management of a cusolverDn handle.
 *
 * The handle is not created automatically on default construction, e.g. in case
 * it is desired to do so in a critical block scope.
 *
 * The handle is destroyed on the call to this class's destructor.
 *
 * \todo throw exceptions in cuda::cusolverDNHandle
 */
struct cusolverDnHandle
{

  private:
    cusolverDnHandle_t m_handle{ NULL };

  public:
    /// Default c'tor
    /** Creates (allocates) the handle and sets the stream to nullptr
     */
    cusolverDnHandle()
    {
        create( nullptr );
    }

    /// Constructor with option to create / not create the handle
    explicit cusolverDnHandle( bool create /**< [in] if true the handle is created. if false it is not created. */ )
    {
        if( create )
        {
            this->create();
        }
    }

    /// Constructor which creates the handle and sets the stream
    explicit cusolverDnHandle( cudaStream_t stream /**< [in] cuda stream to associate with this handle. */ )
    {
        create( stream );
    }

    /// Destructor
    ~cusolverDnHandle()
    {
        if( m_handle )
        {
            cusolverDnDestroy( m_handle );
        }
    }

    /// Create (allocate) the handle.
    void create()
    {
        cusolverStatus_t csec = cusolverDnCreate( &m_handle );
        if( csec != CUSOLVER_STATUS_SUCCESS )
        {
            std::cerr << __FILE__ << " " << __LINE__ << " " << csec << "\n";
            exit( -1 );
        }
    }

    /// Create (allocate) the handle.
    void create( cudaStream_t stream )
    {
        cusolverStatus_t csec = cusolverDnCreate( &m_handle );
        if( csec != CUSOLVER_STATUS_SUCCESS )
        {
            std::cerr << __FILE__ << " " << __LINE__ << " " << csec << "\n";
            exit( -1 );
        }
        setStream( stream );
    }

    /// Create (allocate) the handle.
    void setStream( cudaStream_t stream )
    {
        if( m_handle == NULL )
        {
            std::cerr << __FILE__ << " " << __LINE__ << " cusolverDnHandle::setStream m_handle not set";
            exit( -1 );
        }

        cusolverStatus_t csec = cusolverDnSetStream( m_handle, stream );
        if( csec != CUSOLVER_STATUS_SUCCESS )
        {
            std::cerr << __FILE__ << " " << __LINE__ << "\n";
            exit( -1 );
        }
    }

    /// Get the handle for use in calls to cusolverDN routines
    /**
     * \returns the cusolverDn handle
     */
    cusolverDnHandle_t operator()()
    {
        return m_handle;
    }

    /// Conversion operator, allows objects of this class to be used as if they are the handle
    /**
     * \returns the cusolverDn handle
     */
    operator cusolverDnHandle_t()
    {
        return m_handle;
    }
};

} // namespace cuda
} // namespace mx

#endif // #MXLIB_CUDA

#endif // math_cusolverDnHandle_hpp
