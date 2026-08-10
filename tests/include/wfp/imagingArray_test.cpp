/** \file imagingArray_test.cpp
 * \brief Placeholder tests for APIs declared in include/wfp/imagingArray.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/wfp/imagingArray.hpp"

/** \cond */
/// Compile the CPU imaging-array implementation with its FFTW allocator.
template struct mx::wfp::fftwAllocator<std::complex<double>>;
template class mx::wfp::imagingArray<std::complex<double>, mx::wfp::fftwAllocator<std::complex<double>>, 0>;
/** \endcond */

/** \defgroup imagingArray_unit_tests imagingArray Unit Tests
 * \ingroup wfp_unit_tests
 */

namespace unitTest::placeholder::wfp_imagingArray_test
{

/** \brief Verifies that APIs declared in include/wfp/imagingArray.hpp remain available to the unit-test target.
 *
 * \ingroup imagingArray_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/wfp/imagingArray.hpp.
 */
TEST_CASE( "imagingArray API has a test placeholder", "[wfp::imagingArray][placeholder]" )
{
    SUCCEED( "imagingArray behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::wfp_imagingArray_test
