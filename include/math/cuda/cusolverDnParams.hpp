/** \file cusolverDnParams.hpp
 * \author Jared R. Males
 * \brief Management of a cusolverDnParams structure
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

#ifndef math_cusolverDnParams_hpp
#define math_cusolverDnParams_hpp

#ifdef MXLIB_CUDA

#include <cuda_runtime.h>
#include <cusolverDn.h>

namespace mx
{
namespace cuda
{

/// Management of a cusolverDnParams structure
/** RAII management of a cusolverDnParams structure.
 *
 * The handle is not created automatically on default construction, e.g. in case
 * it is desired to do so in a critical block scope.
 *
 * The handle is destroyed on the call to this class's destructor.
 *
 * \todo throw exceptions in cuda::cusolverDnParams
 */
struct cusolverDnParams
{

  private:
    cusolverDnParams_t m_handle{ NULL };

  public:
    /// Default c'tor
    /** Creates the structure
     */
    cusolverDnParams()
    {
        create();
    }

    /// Constructor with option to create /not-create the structure
    explicit cusolverDnParams( bool create /**< [in] if true the structure is created. if false it is not created. */ )
    {
        if( create )
        {
            this->create();
        }
    }

    /// Destructor
    ~cusolverDnParams()
    {
        if( m_handle )
        {
            cusolverDnDestroyParams( m_handle );
        }
    }

    /// Create (allocate) the structure.
    void create()
    {
        cusolverStatus_t csec = cusolverDnCreateParams( &m_handle );
        if( csec != CUSOLVER_STATUS_SUCCESS )
        {
            std::cerr << __FILE__ << " " << __LINE__ << " " << csec << "\n";
            exit( -1 );
        }
    }

    /// Get the structure for use in calls to cusolverDN routines
    /**
     * \returns the cusolverDnParams structure
     */
    cusolverDnParams_t operator()()
    {
        return m_handle;
    }

    /// Conversion operator, allows objects of this class to be used as if they are the structure
    /**
     * \returns the cusolverDnParams structure
     */
    operator cusolverDnParams_t()
    {
        return m_handle;
    }
};

} // namespace cuda
} // namespace mx

#endif //MXLIB_CUDA

#endif // math_cusolverDnParams_hpp
