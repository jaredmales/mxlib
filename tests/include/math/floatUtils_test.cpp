/** \file floatUtils_test.cpp
 * \brief Tests of floating-point classification utilities.
 */

#include "../../catch2/catch.hpp"

#include "../../../include/math/floatUtils.hpp"

#include <limits>

/// Verify floating-point classification for finite and nonfinite values.
/** Exercises the NaN and finite-value classifiers for float, double, and long double.
 * The same test applies to standard and fast-math builds.
 */
/**
 * \ingroup floatUtils_unit_tests
 */
TEST_CASE( "Floating-point classification handles finite and nonfinite values", "[math::floatUtils]" )
{
    REQUIRE( mx::math::isFinite( std::numeric_limits<float>::max() ) );
    REQUIRE_FALSE( mx::math::isFinite( std::numeric_limits<float>::infinity() ) );
    REQUIRE( mx::math::isNan( std::numeric_limits<float>::quiet_NaN() ) );
    REQUIRE_FALSE( mx::math::isNan( std::numeric_limits<float>::infinity() ) );

    REQUIRE( mx::math::isFinite( std::numeric_limits<double>::max() ) );
    REQUIRE_FALSE( mx::math::isFinite( std::numeric_limits<double>::quiet_NaN() ) );
    REQUIRE( mx::math::isNan( std::numeric_limits<double>::quiet_NaN() ) );
    REQUIRE_FALSE( mx::math::isNan( std::numeric_limits<double>::infinity() ) );

    REQUIRE( mx::math::isFinite( std::numeric_limits<long double>::max() ) );
    REQUIRE_FALSE( mx::math::isFinite( std::numeric_limits<long double>::infinity() ) );
    REQUIRE( mx::math::isNan( std::numeric_limits<long double>::quiet_NaN() ) );
    REQUIRE_FALSE( mx::math::isNan( std::numeric_limits<long double>::infinity() ) );
}
