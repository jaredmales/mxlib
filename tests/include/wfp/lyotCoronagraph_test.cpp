/** \file lyotCoronagraph_test.cpp
 * \brief Placeholder tests for APIs declared in include/wfp/lyotCoronagraph.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/wfp/lyotCoronagraph.hpp"

/** \cond */
/// Compile the Lyot-coronagraph implementation for representative wavefront and mask scalar types.
template struct mx::wfp::lyotCoronagraph<double, double>;
/** \endcond */

/** \defgroup lyotCoronagraph_unit_tests lyotCoronagraph Unit Tests
 * \ingroup wfp_unit_tests
 */

namespace unitTest::placeholder::wfp_lyotCoronagraph_test
{

/** \brief Verifies that APIs declared in include/wfp/lyotCoronagraph.hpp remain available to the unit-test target.
 *
 * \ingroup lyotCoronagraph_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/wfp/lyotCoronagraph.hpp.
 */
TEST_CASE( "lyotCoronagraph API has a test placeholder", "[wfp::lyotCoronagraph][placeholder]" )
{
    SUCCEED( "lyotCoronagraph behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::wfp_lyotCoronagraph_test
