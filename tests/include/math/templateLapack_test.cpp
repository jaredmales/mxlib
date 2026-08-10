/** \file templateLapack_test.cpp
 *  \brief Tests the typed LAPACK wrappers.
 */

#include "../../catch2/catch.hpp"

#include "../../../include/math/templateLapack.hpp"

/// Verify typed lamch wrappers for supported real precisions.
/**
 * \ingroup templateLapack_unit_tests
 */
TEST_CASE( "getting lamch values", "[templateLapack]" )
{
    GIVEN( "a precision" )
    {

        WHEN( "precision is float" )
        {
            // This is just a compilation test.
            float ch = mx::math::lamch<float>( 'E' );
            REQUIRE( ch > 0 );
        }
        WHEN( "precision is double" )
        {
            // This is just a compilation test.
            double ch = mx::math::lamch<double>( 'E' );
            REQUIRE( ch > 0 );
        }
    }
}
