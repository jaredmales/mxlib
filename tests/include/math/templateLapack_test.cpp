#include "../../catch2/catch.hpp"

#include "../../../include/math/templateLapack.hpp"/// getting lamch values


/// getting lamch values
/**
 * \ingroup templateLapack_unit_tests
 */
TEST_CASE( "getting lamch values", "[templateLapack]" )
{
    SECTION( "a precision" )
    {

        SECTION( "precision is float" )
        {
            // This is just a compilation test.
            float ch = mx::math::lamch<float>( 'E' );
            REQUIRE( ch > 0 );
        }
        SECTION( "precision is double" )
        {
            // This is just a compilation test.
            double ch = mx::math::lamch<double>( 'E' );
            REQUIRE( ch > 0 );
        }
    }
}
