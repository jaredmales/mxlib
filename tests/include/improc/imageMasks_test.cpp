/** \file imageMasks_test.cpp
 * \brief Tests image-mask generation.
 */
#include "../../catch2/catch.hpp"

#include <algorithm>
#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/improc/imageMasks.hpp"
#include "../../../include/improc/eigenImage.hpp"

/** \cond */
/// Emit the binary-mask rotation specialization for LCOV accounting.
template void mx::improc::rotateMask<mx::improc::eigenImage<double>, mx::improc::cubicConvolTransform<double>>(
    mx::improc::eigenImage<double> &, mx::improc::eigenImage<double> &, double );
/** \endcond */

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

/** \brief Verifies mx::improc::rotateMask resizes, rotates, and thresholds binary masks.
 *
 * \ingroup imageMasks_unit_tests
 */
TEST_CASE( "rotateMask preserves binary masks through interpolation", "[improc::rotateMask]" )
{
    mx::improc::eigenImage<double> mask( 9, 9 );
    mask.setZero();
    mask( 3, 4 ) = 0.6;
    mask( 4, 3 ) = 0.49;

    mx::improc::eigenImage<double> rotated( 1, 1 );
    mx::improc::rotateMask( rotated, mask, 0.0 );
    REQUIRE( rotated.rows() == mask.rows() );
    REQUIRE( rotated.cols() == mask.cols() );
    REQUIRE( rotated( 3, 4 ) == 1.0 );
    REQUIRE( rotated( 4, 3 ) == 0.0 );

    mx::improc::rotateMask( rotated, mask, mx::math::half_pi<double>() );
    REQUIRE( rotated( 4, 3 ) == 1.0 );
    REQUIRE( rotated( 3, 4 ) == 0.0 );
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

/** \brief Verifies annulusCoords and annulusIndices forward radial pixel buffers to their shared worker.
 *
 * \ingroup imageMasks_unit_tests
 */
TEST_CASE( "Annulus coordinate wrappers apply radial pixel buffers", "[improc::imageMasks][annulusCoords][pixbuf]" )
{
    using angleT = mx::math::degreesT<double>;

    mx::improc::eigenImage<double> radius( 7, 7 );
    mx::improc::eigenImage<double> angle;
    mx::improc::radAngImage<angleT>( radius, angle, 3.0, 3.0 );
    mx::improc::eigenImage<double> *noMask = nullptr;

    const auto coordinates = mx::improc::annulusCoords<angleT>( radius, angle, 3.0, 3.0, 1.5, 2.5, 0.0, 360.0 );
    const auto bufferedCoordinates =
        mx::improc::annulusCoords<angleT>( radius, angle, 3.0, 3.0, 1.5, 2.5, 0.0, 360.0, noMask, 0.5 );

    REQUIRE( coordinates.size() == 12 );
    REQUIRE( bufferedCoordinates.size() == 24 );

    const auto hasCoordinate = []( const std::vector<std::vector<int>> &region, int row, int column )
    { return std::find( region.begin(), region.end(), std::vector<int>{ row, column } ) != region.end(); };

    REQUIRE_FALSE( hasCoordinate( coordinates, 2, 3 ) );
    REQUIRE( hasCoordinate( bufferedCoordinates, 2, 3 ) );
    REQUIRE_FALSE( hasCoordinate( coordinates, 1, 1 ) );
    REQUIRE( hasCoordinate( bufferedCoordinates, 1, 1 ) );
    REQUIRE_FALSE( hasCoordinate( bufferedCoordinates, 0, 3 ) );

    mx::improc::eigenImage<double> mask( 7, 7 );
    mask.setOnes();
    mask( 2, 3 ) = 0;
    const auto maskedBufferedCoordinates =
        mx::improc::annulusCoords<angleT>( radius, angle, 3.0, 3.0, 1.5, 2.5, 0.0, 360.0, &mask, 0.5 );
    REQUIRE_FALSE( hasCoordinate( maskedBufferedCoordinates, 2, 3 ) );

    const auto indices = mx::improc::annulusIndices<angleT>( radius, angle, 3.0, 3.0, 1.5, 2.5, 0.0, 360.0 );
    const auto bufferedIndices =
        mx::improc::annulusIndices<angleT>( radius, angle, 3.0, 3.0, 1.5, 2.5, 0.0, 360.0, noMask, 0.5 );

    REQUIRE( indices.size() == coordinates.size() );
    REQUIRE( bufferedIndices.size() == bufferedCoordinates.size() );
    REQUIRE( std::find( indices.begin(), indices.end(), 23 ) == indices.end() );
    REQUIRE( std::find( bufferedIndices.begin(), bufferedIndices.end(), 23 ) != bufferedIndices.end() );
    REQUIRE( std::find( indices.begin(), indices.end(), 8 ) == indices.end() );
    REQUIRE( std::find( bufferedIndices.begin(), bufferedIndices.end(), 8 ) != bufferedIndices.end() );

    mx::improc::eigenImage<double> fractionalRadius( 10, 10 );
    mx::improc::eigenImage<double> fractionalAngle;
    mx::improc::radAngImage<angleT>( fractionalRadius, fractionalAngle, 4.5, 4.5 );

    const auto fractionalCoordinates = mx::improc::annulusCoords<
        angleT>( fractionalRadius, fractionalAngle, 4.5, 4.5, 3.0, 3.5, 0.0, 360.0, noMask, 0.5 );
    REQUIRE( hasCoordinate( fractionalCoordinates, 1, 4 ) );
    REQUIRE( hasCoordinate( fractionalCoordinates, 8, 4 ) );

    const auto fractionalIndices = mx::improc::annulusIndices<
        angleT>( fractionalRadius, fractionalAngle, 4.5, 4.5, 3.0, 3.5, 0.0, 360.0, noMask, 0.5 );
    REQUIRE( std::find( fractionalIndices.begin(), fractionalIndices.end(), 41 ) != fractionalIndices.end() );
    REQUIRE( std::find( fractionalIndices.begin(), fractionalIndices.end(), 48 ) != fractionalIndices.end() );

    const auto emptyBufferedCoordinates =
        mx::improc::annulusCoords<angleT>( radius, angle, 3.0, 3.0, 0.0, 0.5, 0.0, 360.0, noMask, -0.5 );
    REQUIRE( emptyBufferedCoordinates.empty() );
}

/** \brief Verifies cutImageRegion and insertImageRegion preserve indexed pixel order and resize policy.
 *
 * \ingroup imageMasks_unit_tests
 */
TEST_CASE( "Image regions are cut and inserted by linear index", "[improc::imageMasks][imageRegion]" )
{
    mx::improc::eigenImage<double> image( 3, 4 );
    for( Eigen::Index index = 0; index < image.size(); ++index )
    {
        image( index ) = static_cast<double>( 10 + index );
    }
    const std::vector<size_t> indices{ 0, 5, 11 };

    mx::improc::eigenImage<double> cut;
    mx::improc::cutImageRegion( cut, image, indices );
    REQUIRE( cut.rows() == 3 );
    REQUIRE( cut.cols() == 1 );
    REQUIRE( cut( 0 ) == Approx( 10 ) );
    REQUIRE( cut( 1 ) == Approx( 15 ) );
    REQUIRE( cut( 2 ) == Approx( 21 ) );

    mx::improc::eigenImage<double> preallocated( 3, 1 );
    preallocated.setConstant( -1 );
    mx::improc::cutImageRegion( preallocated, image, indices, false );
    REQUIRE( preallocated.isApprox( cut ) );

    mx::improc::eigenImage<double> inserted( image.size(), 1 );
    inserted.setConstant( -1 );
    auto insertedView = inserted.block( 0, 0, inserted.rows(), 1 );
    mx::improc::insertImageRegion( insertedView, cut, indices );
    REQUIRE( inserted( 0 ) == Approx( 10 ) );
    REQUIRE( inserted( 5 ) == Approx( 15 ) );
    REQUIRE( inserted( 11 ) == Approx( 21 ) );
    REQUIRE( inserted( 1 ) == Approx( -1 ) );

    mx::improc::cutImageRegion( cut, image, {} );
    REQUIRE( cut.rows() == 0 );
    REQUIRE( cut.cols() == 1 );
    mx::improc::insertImageRegion( insertedView, cut, {} );
    REQUIRE( inserted( 0 ) == Approx( 10 ) );
}
