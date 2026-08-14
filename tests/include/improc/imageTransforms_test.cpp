/** \file imageTransforms_test.cpp
 * \brief Tests image transformations.
 */
#include "../../catch2/catch.hpp"

#include <array>
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

/** \brief Verifies the default single-precision cubic-convolution kernel and its complete image footprint.
 *
 * \ingroup imageTransforms_unit_tests
 */
TEST_CASE( "cubicConvolTransform produces normalized four-pixel kernels",
           "[improc::cubicConvolTransform][improc::imageTransforms]" )
{
    using transformT = mx::improc::cubicConvolTransform<float>;
    constexpr Eigen::Index transformWidth = transformT::width;
    constexpr Eigen::Index leftBuffer = transformT::lbuff;

    transformT transform;

    REQUIRE( transform.cubic == -0.5F );
    REQUIRE( transformWidth == 4 );
    REQUIRE( leftBuffer == 1 );

    REQUIRE( transform.cubicConvolKernel( 0.0F ) == 1.0F );
    REQUIRE( transform.cubicConvolKernel( 0.25F ) == 0.8671875F );
    REQUIRE( transform.cubicConvolKernel( 0.5F ) == 0.5625F );
    REQUIRE( transform.cubicConvolKernel( 0.75F ) == 0.2265625F );
    REQUIRE( transform.cubicConvolKernel( 1.0F ) == 0.0F );
    REQUIRE( transform.cubicConvolKernel( 1.25F ) == -0.0703125F );
    REQUIRE( transform.cubicConvolKernel( 1.5F ) == -0.0625F );
    REQUIRE( transform.cubicConvolKernel( 1.75F ) == -0.0234375F );
    REQUIRE( transform.cubicConvolKernel( 2.0F ) == 0.0F );

    const std::array<std::array<float, 2>, 4> phases{
        { { 0.0F, 0.0F }, { 0.25F, 0.75F }, { 0.5F, 0.5F }, { 0.75F, 0.25F } } };

    mx::improc::eigenImage<float> kernel( transformWidth, transformWidth );
    for( const auto &phase : phases )
    {
        transform( kernel, phase[0], phase[1] );
        REQUIRE( kernel.sum() == 1.0F );
    }

    const std::array<float, 4> quarterWeights{ -0.0703125F, 0.8671875F, 0.2265625F, -0.0234375F };
    const std::array<float, 4> threeQuarterWeights{ -0.0234375F, 0.2265625F, 0.8671875F, -0.0703125F };

    transform( kernel, 0.25F, 0.75F );
    for( Eigen::Index row = 0; row < kernel.rows(); ++row )
    {
        for( Eigen::Index col = 0; col < kernel.cols(); ++col )
        {
            REQUIRE( kernel( row, col ) == quarterWeights[row] * threeQuarterWeights[col] );
        }
    }

    const mx::improc::eigenImage<float> constantImage = mx::improc::eigenImage<float>::Constant( 8, 8, 3.25F );
    const Eigen::Index firstAnchor = leftBuffer;
    const Eigen::Index lastAnchor = constantImage.rows() - transformWidth + leftBuffer;
    const std::array<Eigen::Index, 2> anchors{ firstAnchor, lastAnchor };

    REQUIRE( firstAnchor - leftBuffer == 0 );
    REQUIRE( lastAnchor - leftBuffer + transformWidth == constantImage.rows() );
    REQUIRE( firstAnchor - 1 - leftBuffer < 0 );
    REQUIRE( lastAnchor + 1 - leftBuffer + transformWidth > constantImage.rows() );

    for( Eigen::Index rowAnchor : anchors )
    {
        for( Eigen::Index colAnchor : anchors )
        {
            transform( kernel, 0.25F, 0.75F );
            const float sample = ( constantImage.block( rowAnchor - leftBuffer,
                                                        colAnchor - leftBuffer,
                                                        transformWidth,
                                                        transformWidth ) *
                                   kernel )
                                     .sum();
            REQUIRE( sample == 3.25F );
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
