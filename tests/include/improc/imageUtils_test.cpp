/** \file imageUtils_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"
#include "../../../include/improc/imageUtils.hpp"
#include "../../../include/improc/eigenCube.hpp"/// Verify center of light calculation


/// Verify center of light calculation
/**
 * \ingroup imageUtils_unit_tests
 */
TEST_CASE( "Verify center of light calculation", "[improc::imageCenterOfLight]" )
{
    SECTION( "a Gaussian" )
    {
        SECTION( "geometric center" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 64, 64 );

            mx::math::func::gaussian2D<double>( im.data(), im.rows(), im.cols(), 0., 1.0, 31.5, 31.5, 2 );

            double x, y;
            mx::improc::imageCenterOfLight( x, y, im );

            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 31.5, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 31.5, 1e-8 ) );
        }
        SECTION( "geometric quarter" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 64, 64 );

            mx::math::func::gaussian2D<double>( im.data(), im.rows(), im.cols(), 0., 1.0, 15.5, 15.5, 2 );

            double x, y;
            mx::improc::imageCenterOfLight( x, y, im );

            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 15.5, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 15.5, 1e-8 ) );
        }
    }
}
