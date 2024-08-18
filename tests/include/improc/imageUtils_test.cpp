/** \file imageUtils_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"
#include "../../../include/improc/imageUtils.hpp"
#include "../../../include/improc/eigenCube.hpp"

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

            REQUIRE( fabs( x - 31.5 ) < 1e-8 );
            REQUIRE( fabs( y - 31.5 ) < 1e-8 );
        }
        WHEN( "geometric quarter" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 64, 64 );

            mx::math::func::gaussian2D<double>( im.data(), im.rows(), im.cols(), 0., 1.0, 15.5, 15.5, 2 );

            double x, y;
            mx::improc::imageCenterOfLight( x, y, im );

            REQUIRE( fabs( x - 15.5 ) < 1e-8 );
            REQUIRE( fabs( y - 15.5 ) < 1e-8 );
        }
    }
}
