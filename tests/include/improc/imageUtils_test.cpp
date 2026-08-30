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

/// Verify image invalid-pixel generation and classification for floating-point values.
/** Exercises the default NaN contract or a configured finite override, plus nonfinite classification. */
/**
 * \ingroup imageUtils_unit_tests
 */
TEST_CASE( "Image invalid-pixel detection handles sentinel and nonfinite values", "[improc::isInvalidPixel]" )
{
    REQUIRE_FALSE( mx::improc::isInvalidPixel( 0.0F ) );
    REQUIRE( mx::improc::isInvalidPixel( std::numeric_limits<float>::quiet_NaN() ) );
    REQUIRE( mx::improc::isInvalidPixel( std::numeric_limits<double>::quiet_NaN() ) );
    REQUIRE( mx::improc::isInvalidPixel( std::numeric_limits<float>::infinity() ) );
    REQUIRE( mx::improc::isInvalidPixel( -std::numeric_limits<float>::infinity() ) );
    REQUIRE( mx::improc::isInvalidPixel( mx::improc::invalidNumber<float>() ) );
    REQUIRE( mx::improc::isInvalidPixel( mx::improc::invalidNumber<double>() ) );

#ifdef MXLIB_INVALID_NUMBER_VALUE
    REQUIRE( mx::improc::invalidNumber<float>() == static_cast<float>( MXLIB_INVALID_NUMBER_VALUE ) );
    REQUIRE( mx::improc::invalidNumber<double>() == static_cast<double>( MXLIB_INVALID_NUMBER_VALUE ) );
#else
    REQUIRE( mx::math::isNan( mx::improc::invalidNumber<float>() ) );
    REQUIRE( mx::math::isNan( mx::improc::invalidNumber<double>() ) );
#endif
}

/** centroiding Gaussians with center of light
 *
 * Verify center of light calculation
 *
 */
/**
 * \ingroup imageUtils_unit_tests
 */
TEST_CASE( "Verify center of light calculation", "[improc::imageCenterOfLight]" )
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

/** \brief Verifies masked and unmasked image medians with internal and caller-owned workspaces.
 *
 * \ingroup imageUtils_unit_tests
 */
TEST_CASE( "imageMedian selects valid image samples", "[improc::imageMedian]" )
{
    mx::improc::eigenImage<double> image( 2, 3 );
    image << 1.0, 9.0, 3.0, 7.0, 5.0, 11.0;

    REQUIRE( mx::improc::imageMedian( image ) == Approx( 6.0 ) );

    std::vector<double> work;
    REQUIRE( mx::improc::imageMedian( image, &work ) == Approx( 6.0 ) );
    REQUIRE( work.size() == static_cast<size_t>( image.size() ) );

    Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> mask( 2, 3 );
    mask << 1, 0, 0, 0, 1, 1;

    REQUIRE( mx::improc::imageMedian( image, &mask, &work ) == Approx( 5.0 ) );
    REQUIRE( work.size() == 3 );
    REQUIRE( mx::improc::imageMedian( image, &mask ) == Approx( 5.0 ) );
}

/// Verify zeroNaNCube replaces invalid cube pixels and records their locations.
/** Exercises both mx::improc::zeroNaNCube overloads, including mask allocation and finite pixels. */
/**
 * \ingroup imageUtils_unit_tests
 */
TEST_CASE( "zeroNaNCube replaces invalid pixels and reports a mask", "[improc::zeroNaNCube]" )
{
    mx::improc::eigenCube<float> cube( 2, 3, 2 );
    cube.cube().setOnes();
    cube.image( 0 )( 0, 1 ) = std::numeric_limits<float>::quiet_NaN();
    cube.image( 1 )( 1, 2 ) = std::numeric_limits<float>::infinity();

    mx::improc::eigenCube<float> mask;
    mx::improc::zeroNaNCube( cube, &mask );

    REQUIRE( cube.image( 0 )( 0, 1 ) == 0.0F );
    REQUIRE( cube.image( 1 )( 1, 2 ) == 0.0F );
    REQUIRE( cube.image( 0 )( 1, 1 ) == 1.0F );
    REQUIRE( mask.rows() == cube.rows() );
    REQUIRE( mask.cols() == cube.cols() );
    REQUIRE( mask.planes() == cube.planes() );
    REQUIRE( mask.image( 0 )( 0, 1 ) == 1.0F );
    REQUIRE( mask.image( 1 )( 1, 2 ) == 1.0F );
    REQUIRE( mask.image( 0 )( 1, 1 ) == 0.0F );

    mx::improc::zeroNaNCube( cube );
    REQUIRE( cube.image( 0 )( 1, 1 ) == 1.0F );
}
