/** \file imageMasks_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/improc/imageMasks.hpp"
#include "../../../include/improc/eigenImage.hpp"/// Masking wedges in an image


/// Masking wedges in an image
/**
 * \ingroup imageMasks_unit_tests
 */
TEST_CASE( "Masking wedges in an image", "[improc::imageMasks::maskWedge]" )
{
    SECTION( "a single wedge" )
    {
        SECTION( "geometric center, 0-90 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 45.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        SECTION( "geometric center, 90-180 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 135.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        SECTION( "geometric center, 180-270 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 225.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        SECTION( "geometric center, 270-360 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 315.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        SECTION( "geometric center, 45-135 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 90.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        SECTION( "geometric center, 135-225 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 180.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        SECTION( "geometric center, 225-315 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 270.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        SECTION( "geometric center, 315-45 degrees" )
        {
            mx::improc::eigenImage<double> im;
            im.resize( 1024, 1024 );
            im.setZero();

            double xcen = 0.5 * ( im.rows() - 1 );
            double ycen = 0.5 * ( im.cols() - 1 );

            mx::improc::maskWedge( im, xcen, ycen, 0.0, 45.0, 1 );

            REQUIRE( im.sum() == 512 * 512 );
        }
        SECTION( "geometric center, 3 wedges of 120 degrees" )
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
