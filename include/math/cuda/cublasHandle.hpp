/** \file cublasHandle.hpp
 * \author Jared R. Males
 * \brief Management of a cublas handle
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

#ifndef math_cublasHandle_hpp
#define math_cublasHandle_hpp

#include <cuda_runtime.h>
#include <cublas_v2.h>

namespace mx
{
namespace cuda
{

/// Management of a cublas handle
/** RAII management of a cublas handle.
 *
 * The handle is not created automatically on default construction, e.g. in case
 * it is desired to do so in a critical block scope.
 *
 * The handle is destroyed on the call to this class's destructor.
 *
 * \todo throw exceptions in cuda::cublasHandle
 */
struct cublasHandle
{

  private:
    cublasHandle_t m_handle{ NULL };

  public:
    /// Default c'tor
    /** Creates the handle.
     */
    cublasHandle()
    {
        create();
    }

    /// Constructor with option to create / not create the handle
    explicit cublasHandle( bool create /**< [in] if true the handle is created. if false it is not created. */ )
    {
        if( create )
        {
            this->create();
        }
    }

    /// Destructor
    ~cublasHandle()
    {
        if( m_handle )
        {
            cublasDestroy( m_handle );
        }
    }

    /// Create (allocate) the handle.
    void create()
    {
        cublasStatus_t cbec = cublasCreate( &m_handle );
        if( cbec != CUBLAS_STATUS_SUCCESS )
        {
            std::cerr << __FILE__ << " " << __LINE__ << " " << cbec << "\n";
            exit( -1 );
        }
    }

    /// Get the handle for use in calls to cublas routines
    /**
     * \returns the cublas handle
     */
    cublasHandle_t operator()()
    {
        return m_handle;
    }

    /// Conversion operator, allows objects of this class to be used as if they are the handle
    /**
     * \returns the cublas handle
     */
    operator cublasHandle_t()
    {
        return m_handle;
    }
};

} // namespace cuda
} // namespace mx
#endif // math_cublasHandle_hpp
