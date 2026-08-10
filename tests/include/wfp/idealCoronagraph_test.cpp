/** \file idealCoronagraph_test.cpp
 * \brief Placeholder tests for APIs declared in include/wfp/idealCoronagraph.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/wfp/idealCoronagraph.hpp"

/** \cond */
/// Compile the ideal-coronagraph implementation for a representative scalar type.
template struct mx::wfp::idealCoronagraph<double>;
/** \endcond */

/** \defgroup idealCoronagraph_unit_tests idealCoronagraph Unit Tests
 * \ingroup wfp_unit_tests
 */

namespace unitTest::placeholder::wfp_idealCoronagraph_test
{

/** \brief Verifies that APIs declared in include/wfp/idealCoronagraph.hpp remain available to the unit-test target.
 *
 * \ingroup idealCoronagraph_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/wfp/idealCoronagraph.hpp.
 */
TEST_CASE( "idealCoronagraph API has a test placeholder", "[wfp::idealCoronagraph][placeholder]" )
{
    SUCCEED( "idealCoronagraph behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::wfp_idealCoronagraph_test
