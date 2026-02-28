#include "../../catch2/catch.hpp"

#include "../../../include/math/templateBLAS.hpp"/// testing scal


/// testing scal
/**
 * \ingroup templateBLAS_unit_tests
 */
TEST_CASE( "testing scal", "[templateBLAS]" )
{
    SECTION( "a precision" )
    {

        SECTION( "precision is float" )
        {
            int N = 2;
            float alpha = 2.;
            float x[] = { 1., 2. };
            int incX = 1;
            mx::math::scal( N, alpha, x, incX );
            REQUIRE_THAT( x[0], Catch::Matchers::WithinAbs( 2.0, 1e-6 ) );
            REQUIRE_THAT( x[1], Catch::Matchers::WithinAbs( 4.0, 1e-6 ) );
        }

        SECTION( "precision is double" )
        {
            int N = 2;
            double alpha = 2.;
            double x[] = { 1., 2. };
            int incX = 1;
            mx::math::scal( N, alpha, x, incX );
            REQUIRE_THAT( x[0], Catch::Matchers::WithinAbs( 2.0, 1e-12 ) );
            REQUIRE_THAT( x[1], Catch::Matchers::WithinAbs( 4.0, 1e-12 ) );
        }
    }
}
