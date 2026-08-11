/** \file imageTransforms_test.cpp
 * \brief Tests image transformations.
 */
#include "../../catch2/catch.hpp"

#include <Eigen/Dense>
#include <numbers>
#include <vector>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"

#include "../../../include/improc/imageUtils.hpp"
#include "../../../include/improc/eigenCube.hpp"

#include "../../../include/improc/imageTransforms.hpp"

/** \cond
 * Explicit instantiations compile-check the double-precision bilinear and cubic-convolution transform policies and
 * emit their header-defined methods for coverage accounting in this test translation unit.
 */
template struct mx::improc::bilinearTransform<double>;
template struct mx::improc::cubicConvolTransform<double>;
/** \endcond */

namespace
{
double imageMSE( const mx::improc::eigenImage<double> &a, const mx::improc::eigenImage<double> &b )
{
    return ( a - b ).square().mean();
}
} // namespace

/** Verify direction and accuracy of various image shifts
 *
 * Tests image shifts by fractional pixels.
 *
 */
/**
 * \ingroup imageTransforms_unit_tests
 */
TEST_CASE( "Verify direction and accuracy of various image shifts", "[improc::imageTransforms]" )
{
    GIVEN( "a Gaussian image" )
    {
        WHEN( "shifting" )
        {
            mx::improc::eigenImage<double> im, shift, ref;

            im.resize( 256, 256 );
            shift.resize( im.rows(), im.cols() );
            ref.resize( im.rows(), im.cols() );

            // Use sigma = 8 to get a well oversampled image, making shifts more accurate
            mx::math::func::gaussian2D<double>( im.data(), im.rows(), im.cols(), 0., 1.0, 127.5, 127.5, 8 );

            mx::improc::imageShift( shift, im, -0.5, -0.5, mx::improc::cubicConvolTransform<double>() );
            mx::math::func::gaussian2D<double>( ref.data(), ref.rows(), ref.cols(), 0., 1.0, 127.0, 127.0, 8 );
            REQUIRE_THAT( imageMSE( shift, ref ), Catch::Matchers::WithinAbs( 0.0, 1e-5 ) );

            mx::improc::imageShift( shift, im, +0.5, +0.5, mx::improc::cubicConvolTransform<double>() );
            mx::math::func::gaussian2D<double>( ref.data(), ref.rows(), ref.cols(), 0., 1.0, 128.0, 128.0, 8 );
            REQUIRE_THAT( imageMSE( shift, ref ), Catch::Matchers::WithinAbs( 0.0, 1e-5 ) );

            mx::improc::imageShift( shift, im, +1.0, +1.0, mx::improc::cubicConvolTransform<double>() );
            mx::math::func::gaussian2D<double>( ref.data(), ref.rows(), ref.cols(), 0., 1.0, 128.5, 128.5, 8 );
            REQUIRE_THAT( imageMSE( shift, ref ), Catch::Matchers::WithinAbs( 0.0, 1e-5 ) );

            mx::improc::imageShift( shift, im, +0.5, -0.5, mx::improc::cubicConvolTransform<double>() );
            mx::math::func::gaussian2D<double>( ref.data(), ref.rows(), ref.cols(), 0., 1.0, 128.0, 127.0, 8 );
            REQUIRE_THAT( imageMSE( shift, ref ), Catch::Matchers::WithinAbs( 0.0, 1e-5 ) );

            mx::improc::imageShift( shift, im, -0.3, +0.7, mx::improc::cubicConvolTransform<double>() );
            mx::math::func::gaussian2D<double>( ref.data(), ref.rows(), ref.cols(), 0., 1.0, 127.2, 128.2, 8 );
            REQUIRE_THAT( imageMSE( shift, ref ), Catch::Matchers::WithinAbs( 0.0, 1e-5 ) );

            mx::improc::imageShift( shift, im, 1.3, -0.7, mx::improc::cubicConvolTransform<double>() );
            mx::math::func::gaussian2D<double>( ref.data(), ref.rows(), ref.cols(), 0., 1.0, 128.8, 126.8, 8 );
            REQUIRE_THAT( imageMSE( shift, ref ), Catch::Matchers::WithinAbs( 0.0, 1e-5 ) );
        }
    }
}

/** \brief Verifies cubic-convolution image rotation direction, flux, allocation, and edge handling.
 *
 * \ingroup imageTransforms_unit_tests
 */
TEST_CASE( "imageRotate rotates an impulse counterclockwise", "[improc::imageRotate]" )
{
    mx::improc::eigenImage<double> image = mx::improc::eigenImage<double>::Zero( 33, 33 );
    image( 10, 16 ) = 1.0;

    mx::improc::eigenImage<double> rotated;
    mx::improc::imageRotate( rotated,
                             image,
                             0.5 * std::numbers::pi_v<double>,
                             mx::improc::cubicConvolTransform<double>() );

    REQUIRE( rotated.rows() == image.rows() );
    REQUIRE( rotated.cols() == image.cols() );
    REQUIRE( rotated( 16, 10 ) == Approx( 1.0 ) );
    REQUIRE( rotated.sum() == Approx( 1.0 ).margin( 1e-12 ) );
    REQUIRE( rotated( 0, 0 ) == 0.0 );
}
