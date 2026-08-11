/** \file vectorUtils_test.cpp
 * \brief Tests vector statistics used by HCI preprocessing.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/math/vectorUtils.hpp"

#include <vector>

namespace unitTest::math_vectorUtils_test
{

/** \brief Verifies the sample-variance overload used by HCI time-series normalization.
 *
 * \ingroup vectorUtils_unit_tests
 */
TEST_CASE( "vectorVariance computes sample variance", "[math::vectorUtils::vectorVariance]" )
{
    const std::vector<double> values{ 1.0, 3.0 };

    REQUIRE( mx::math::vectorVariance( values ) == Approx( 2.0 ) );
    REQUIRE( mx::math::vectorVariance( values, 0.0 ) == Approx( 10.0 ) );
}

} // namespace unitTest::math_vectorUtils_test
