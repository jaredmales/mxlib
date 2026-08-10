/** \file ftTypes_test.cpp
 * \brief Tests Fourier-transform type traits.
 */
#include "../../../catch2/catch.hpp"

#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../../include/math/ft/ftTypes.hpp"

/// Test the ftTypes
/** This is just a compilation check.
 *
 * \ingroup ftTypes_unit_tests
 */
TEST_CASE( "Test the ftTypes", "[math::ft]" )
{
    mx::math::ft::dir F = mx::math::ft::dir::forward;

    REQUIRE(F == mx::math::ft::dir::forward);
    REQUIRE(F != mx::math::ft::dir::backward);

    mx::math::ft::dir B = mx::math::ft::dir::backward;

    REQUIRE(B != mx::math::ft::dir::forward);
    REQUIRE(B == mx::math::ft::dir::backward);
}
