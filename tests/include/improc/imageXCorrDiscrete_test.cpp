/** \file imageXCorrDiscrete_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"
#include "../../../include/improc/imageXCorrDiscrete.hpp"
#include "../../../include/improc/eigenCube.hpp"/// Verify X-Corr with center of light calculation


/// Verify X-Corr with center of light calculation
/**
 * \ingroup imageXCorrDiscrete_unit_tests
 */
TEST_CASE( "Verify X-Corr with center of light calculation", "[improc::imageXCorrDiscrete]" )
{
    SECTION( "two Gaussians" )
    {
        SECTION( "1 at geometric center" )
        {
            mx::improc::eigenImage<double> im0, im2;
            im0.resize( 64, 64 );
            im2.resize( 64, 64 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, 31.5, 31.5, 2 );
            mx::math::func::gaussian2D<double>( im2.data(), im2.rows(), im2.cols(), 0., 1.0, 31.5 + 4, 31.5 + 4, 2 );

            double x, y;
            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf( 5 );

            mx::improc::eigenImage<double> refIm = im0.block( 10, 10, im0.rows() - 11, im0.cols() - 11 );
            xcf.refIm( refIm );
            xcf.m_peakMethod = mx::improc::xcorrPeakMethod::centroid;

            xcf( x, y, im2 );

            std::cerr << x << " " << y << "\n";

            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 2, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 2, 1e-8 ) );
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
