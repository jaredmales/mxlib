/** \file imageUtils_test.cpp
 * \brief Tests of image-processing utilities.
 */
#include "../../catch2/catch.hpp"

#include <Eigen/Dense>

#include <limits>
#include <vector>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"
#include "../../../include/improc/imageUtils.hpp"
#include "../../../include/improc/eigenCube.hpp"

/// Verify image invalid-pixel classification for the sentinel and nonfinite values.
/** Preserves the established invalid-pixel contract. */
TEST_CASE( "Image invalid-pixel detection handles sentinel and nonfinite values", "[improc::isInvalidPixel]" )
{
    REQUIRE_FALSE( mx::improc::isInvalidPixel( 0.0F ) );
    REQUIRE( mx::improc::isInvalidPixel( std::numeric_limits<float>::quiet_NaN() ) );
    REQUIRE( mx::improc::isInvalidPixel( std::numeric_limits<float>::infinity() ) );
    REQUIRE( mx::improc::isInvalidPixel( -std::numeric_limits<float>::infinity() ) );
    REQUIRE( mx::improc::isInvalidPixel( mx::improc::invalidNumber<float>() ) );
}

/** Scenario: centroiding Gaussians with center of light
 *
 * Verify center of light calculation
 *
 * \anchor tests_improc_imageUtils_imageCenterOfLight
 */
SCENARIO( "Verify center of light calculation", "[improc::imageCenterOfLight]" )
{
    GIVEN( "a Gaussian" )
    {
        WHEN( "geometric center" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 64, 64 );

            mx::math::func::gaussian2D<double>( im.data(), im.rows(), im.cols(), 0., 1.0, 31.5, 31.5, 2 );

            double x, y;
            mx::improc::imageCenterOfLight( x, y, im );

            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 31.5, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 31.5, 1e-8 ) );
        }
        WHEN( "geometric quarter" )
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
