/** \file imageMasks_test.cpp
 * \brief Tests image-mask generation.
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/improc/imageMasks.hpp"
#include "../../../include/improc/eigenImage.hpp"

/** Masking wedges in an image
 *
 * Verify wedge masking, including that all pixels are masked for continuous rotations of the wedge
 *
 */
/**
 * \ingroup imageMasks_unit_tests
 */
TEST_CASE( "Masking wedges in an image", "[improc::imageMasks::maskWedge]" )
{
    GIVEN( "a single wedge" )
    {
        WHEN( "geometric center, 0-90 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 45.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        WHEN( "geometric center, 90-180 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 135.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        WHEN( "geometric center, 180-270 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 225.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        WHEN( "geometric center, 270-360 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 315.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        WHEN( "geometric center, 45-135 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 90.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        WHEN( "geometric center, 135-225 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 180.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        WHEN( "geometric center, 225-315 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 270.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        WHEN( "geometric center, 315-45 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 0.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        WHEN( "geometric center, 3 wedges of 120 degrees" )
        {
            mx::improc::eigenImage<double> im60;
            im60.resize( 1024, 1024 );
            im60.setZero();

            double xcen = 0.5 * ( im60.rows() - 1 );
            double ycen = 0.5 * ( im60.cols() - 1 );

            mx::improc::maskWedge( im60, xcen, ycen, 60.0, 60.0, 1 );

            mx::improc::eigenImage<double> im180;
            im180.resize( 1024, 1024 );
            im180.setZero();

            xcen = 0.5 * ( im180.rows() - 1 );
            ycen = 0.5 * ( im180.cols() - 1 );

            mx::improc::maskWedge( im180, xcen, ycen, 180.0, 60.0, 1 );

            mx::improc::eigenImage<double> im300;
            im300.resize( 1024, 1024 );
            im300.setZero();

            xcen = 0.5 * ( im300.rows() - 1 );
            ycen = 0.5 * ( im300.cols() - 1 );

            mx::improc::maskWedge( im300, xcen, ycen, 300.0, 60.0, 1 );

            REQUIRE( im60.sum() + im180.sum() + im300.sum() == 1024 * 1024 );
        }
    }
}

/** \brief Verifies radius and degree-angle images for a centered rectangular grid.
 *
 * \ingroup imageMasks_unit_tests
 */
TEST_CASE( "radAngImage produces HCI radius and angle coordinates", "[improc::imageMasks::radAngImage]" )
{
    mx::improc::eigenImage<double> radius( 3, 5 );
    mx::improc::eigenImage<double> angle;

    mx::improc::radAngImage<mx::math::degreesT<double>>( radius, angle, 1.0, 2.0 );

    REQUIRE( angle.rows() == 3 );
    REQUIRE( angle.cols() == 5 );
    REQUIRE( radius( 1, 2 ) == Approx( 0.0 ) );
    REQUIRE( angle( 1, 2 ) == Approx( 0.0 ) );
    REQUIRE( radius( 2, 2 ) == Approx( 1.0 ) );
    REQUIRE( angle( 2, 2 ) == Approx( 0.0 ) );
    REQUIRE( radius( 1, 3 ) == Approx( 1.0 ) );
    REQUIRE( angle( 1, 3 ) == Approx( 90.0 ) );
    REQUIRE( radius( 0, 2 ) == Approx( 1.0 ) );
    REQUIRE( angle( 0, 2 ) == Approx( 180.0 ) );
    REQUIRE( radius( 1, 1 ) == Approx( 1.0 ) );
    REQUIRE( angle( 1, 1 ) == Approx( 270.0 ) );

    mx::improc::radAngImage<mx::math::degreesT<double>>( radius, angle, 1.0, 2.0, 2.0 );
    REQUIRE( radius( 2, 2 ) == Approx( 2.0 ) );
}

/** \brief Verifies HCI annulus indices for radial, angular, image-boundary, and mask constraints.
 *
 * \ingroup imageMasks_unit_tests
 */
TEST_CASE( "annulusIndices selects bounded and masked HCI regions", "[improc::imageMasks::annulusIndices]" )
{
    using angleT = mx::math::degreesT<double>;

    mx::improc::eigenImage<double> radius( 5, 5 );
    mx::improc::eigenImage<double> angle;
    mx::improc::radAngImage<angleT>( radius, angle, 2.0, 2.0 );

    const std::vector<size_t> full =
        mx::improc::annulusIndices<angleT>( radius, angle, 2.0, 2.0, 1.0, 2.0, 0.0, 360.0 );
    REQUIRE( full == std::vector<size_t>{ 6, 7, 8, 11, 13, 16, 17, 18 } );

    mx::improc::eigenImage<double> mask( 5, 5 );
    mask.setOnes();
    mask( 1, 1 ) = 0;
    const std::vector<size_t> masked =
        mx::improc::annulusIndices<angleT>( radius, angle, 2.0, 2.0, 1.0, 2.0, 0.0, 360.0, &mask );
    REQUIRE( masked == std::vector<size_t>{ 7, 8, 11, 13, 16, 17, 18 } );

    const std::vector<size_t> clipped =
        mx::improc::annulusIndices<angleT>( radius, angle, 2.0, 2.0, 0.0, 10.0, -360.0, 360.0 );
    REQUIRE( clipped.size() == 25 );

    const std::vector<size_t> half =
        mx::improc::annulusIndices<angleT>( radius, angle, 2.0, 2.0, 0.0, 10.0, 0.0, 180.0 );
    REQUIRE_FALSE( half.empty() );
    REQUIRE( half.size() < clipped.size() );

    const std::vector<size_t> wrapped =
        mx::improc::annulusIndices<angleT>( radius, angle, 2.0, 2.0, 0.0, 10.0, 300.0, 60.0 );
    REQUIRE_FALSE( wrapped.empty() );
    REQUIRE( wrapped.size() < clipped.size() );
}
