/** \file cudaTestUtils.hpp
 * \brief CUDA runtime checks shared by device-dependent tests.
 */

#ifndef mxlib_tests_cudaTestUtils_hpp
#define mxlib_tests_cudaTestUtils_hpp

#include <cuda_runtime.h>

/// \cond
namespace mxlibTest
{

/// Return whether the CUDA runtime can see at least one device.
inline bool cudaDeviceAvailable()
{
    int deviceCount = 0;
    const cudaError_t status = cudaGetDeviceCount( &deviceCount );
    return status == cudaSuccess && deviceCount > 0;
}

} // namespace mxlibTest
/// \endcond

#endif // mxlib_tests_cudaTestUtils_hpp
