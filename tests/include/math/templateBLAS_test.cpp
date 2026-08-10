/** \file templateBLAS_test.cpp
 *  \brief Tests the typed BLAS wrappers.
 */

#include "../../catch2/catch.hpp"

#include "../../../include/math/templateBLAS.hpp"

/// Verify typed scal wrappers for supported real and complex precisions.
/**
 * \ingroup templateBLAS_unit_tests
 */
TEST_CASE( "testing scal", "[templateBLAS]" )
{
    GIVEN( "a precision" )
    {

        WHEN( "precision is float" )
        {
            int N = 2;
            float alpha = 2.;
            float x[] = { 1., 2. };
            int incX = 1;
            mx::math::scal( N, alpha, x, incX );
            REQUIRE_THAT( x[0], Catch::Matchers::WithinAbs( 2.0, 1e-6 ) );
            REQUIRE_THAT( x[1], Catch::Matchers::WithinAbs( 4.0, 1e-6 ) );
        }

        WHEN( "precision is double" )
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
