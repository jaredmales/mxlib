/** \file imageFilters_test.cpp
 * \brief Tests of image filtering and smoothing utilities.
 */
#include "../../catch2/catch.hpp"

#include <cmath>
#include <limits>
#include <vector>
#include <Eigen/Dense>

#include "../../../include/improc/imageFilters.hpp"
#include "../../../include/improc/imageUtils.hpp"

/** \cond
 * Explicit instantiations compile-check the default Gaussian and azimuthal box kernels with the canonical
 * double-precision Eigen image and emit their header-defined methods for coverage accounting in this test
 * translation unit.
 */
template struct mx::improc::gaussKernel<mx::improc::eigenImage<double>>;
template struct mx::improc::azBoxKernel<mx::improc::eigenImage<double>>;
/** \endcond */

namespace unitTest
{
namespace improcTest
{
namespace imageFiltersTest
{

/** \cond */
template <typename arrayT>
bool arrayAllFinite( const arrayT &array )
{
    for( int column = 0; column < array.cols(); ++column )
    {
        for( int row = 0; row < array.rows(); ++row )
        {
            if( !mx::math::isFinite( array( row, column ) ) )
            {
                return false;
            }
        }
    }
    return true;
}
/** \endcond */

/// Verify precalcKernel returns the same Gaussian and azimuthal kernels as their production generators.
/** Exercises mx::improc::precalcKernel::setKernel against mx::improc::gaussKernel::setKernel and
 * mx::improc::azBoxKernel::setKernel.
 */
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "precalcKernel caches production kernel values", "[improc::precalcKernel]" )
{
    SECTION( "with Gaussian kernel" )
    {
        typedef mx::improc::gaussKernel<mx::improc::eigenImage<float>> kernelT;

        kernelT kernel( 2 );
        mx::improc::precalcKernel<kernelT> pck( kernel, 8, 8, 3.5, 3.5 );

        REQUIRE( pck.m_kernels.size() == 8 * 8 );

        mx::improc::eigenImage<float> generatedKernel, cachedKernel;
        for( int column = 0; column < 8; ++column )
        {
            for( int row = 0; row < 8; ++row )
            {
                REQUIRE( kernel.setKernel( row - 3.5, column - 3.5, generatedKernel ) == mx::error_t::noerror );
                REQUIRE( pck.setKernel( row - 3.5, column - 3.5, cachedKernel ) == mx::error_t::noerror );
                REQUIRE( ( generatedKernel == cachedKernel ).all() );
            }
        }
    }

    SECTION( "with azBoxKernel kernel" )
    {
        typedef mx::improc::azBoxKernel<mx::improc::eigenImage<float>, 1, mx::verbose::o> kernelT;

        kernelT kernel( 3, 5 );
        mx::improc::precalcKernel<kernelT> pck( kernel, 8, 8, 3.5, 3.5 );

        REQUIRE( pck.m_kernels.size() == 8 * 8 );

        mx::improc::eigenImage<float> generatedKernel, cachedKernel;
        for( int column = 0; column < 8; ++column )
        {
            for( int row = 0; row < 8; ++row )
            {
                REQUIRE( kernel.setKernel( row - 3.5, column - 3.5, generatedKernel ) == mx::error_t::noerror );
                REQUIRE( pck.setKernel( row - 3.5, column - 3.5, cachedKernel ) == mx::error_t::noerror );
                REQUIRE( ( generatedKernel == cachedKernel ).all() );
            }
        }
    }
}

/// Verify azBoxKernel produces finite normalized kernels at image centers and neighboring pixels.
/** Exercises mx::improc::azBoxKernel::setKernel for exact, half-pixel, radial, and azimuthal center-relative
 * coordinates. The exact image center consistently uses the positive x-axis as its radial orientation.
 */
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "azBoxKernel is finite at centers and neighboring pixels", "[improc::azBoxKernel]" )
{
    auto requireFiniteNormalized = []( const mx::improc::eigenImage<double> &kernel )
    {
        REQUIRE( kernel.size() > 0 );
        REQUIRE( arrayAllFinite( kernel ) );
        REQUIRE( kernel.sum() == Approx( 1.0 ).margin( 1e-12 ) );
    };

    SECTION( "exact center uses the positive x-axis and configured half-widths" )
    {
        typedef mx::improc::azBoxKernel<mx::improc::eigenImage<double>, 2, mx::verbose::o> kernelT;
        kernelT kernelGenerator( 3.0, 1.0 );
        mx::improc::eigenImage<double> kernel;

        REQUIRE( kernelGenerator.setKernel( 0.0, 0.0, kernel ) == mx::error_t::noerror );
        requireFiniteNormalized( kernel );
        REQUIRE( kernelGenerator.maxWidth() == 4 );
        REQUIRE( kernel.rows() == 8 );
        REQUIRE( kernel.cols() == 4 );
    }

    SECTION( "an odd kernel includes its finite center sample" )
    {
        typedef mx::improc::azBoxKernel<mx::improc::eigenImage<double>, 1, mx::verbose::o> kernelT;
        kernelT kernelGenerator( 2.0, 2.0 );
        mx::improc::eigenImage<double> centerKernel, radialNeighbor, azimuthalNeighbor, halfPixelKernel;

        REQUIRE( kernelGenerator.setKernel( 0.0, 0.0, centerKernel ) == mx::error_t::noerror );
        REQUIRE( kernelGenerator.setKernel( 1.0, 0.0, radialNeighbor ) == mx::error_t::noerror );
        REQUIRE( kernelGenerator.setKernel( 0.0, 1.0, azimuthalNeighbor ) == mx::error_t::noerror );
        REQUIRE( kernelGenerator.setKernel( 0.5, 0.5, halfPixelKernel ) == mx::error_t::noerror );

        requireFiniteNormalized( centerKernel );
        requireFiniteNormalized( radialNeighbor );
        requireFiniteNormalized( azimuthalNeighbor );
        requireFiniteNormalized( halfPixelKernel );
        REQUIRE( kernelGenerator.maxWidth() == 1 );
        REQUIRE( centerKernel.rows() % 2 == 1 );
        REQUIRE( centerKernel.cols() % 2 == 1 );
        REQUIRE( centerKernel( centerKernel.rows() / 2, centerKernel.cols() / 2 ) > 0 );
        REQUIRE( radialNeighbor( radialNeighbor.rows() / 2, radialNeighbor.cols() / 2 ) > 0 );
    }

    SECTION( "radial and azimuthal neighbors rotate an anisotropic kernel" )
    {
        typedef mx::improc::azBoxKernel<mx::improc::eigenImage<double>, 2, mx::verbose::o> kernelT;
        kernelT kernelGenerator( 3.0, 1.0 );
        mx::improc::eigenImage<double> radialKernel, azimuthalKernel;

        REQUIRE( kernelGenerator.setKernel( 2.0, 0.0, radialKernel ) == mx::error_t::noerror );
        REQUIRE( kernelGenerator.setKernel( 0.0, 2.0, azimuthalKernel ) == mx::error_t::noerror );

        requireFiniteNormalized( radialKernel );
        requireFiniteNormalized( azimuthalKernel );
        REQUIRE( radialKernel.rows() == azimuthalKernel.cols() );
        REQUIRE( radialKernel.cols() == azimuthalKernel.rows() );
    }
}

/// Verify azBoxKernel honors maximum azimuth and validates degenerate configuration values.
/** Exercises mx::improc::azBoxKernel::setKernel normalization, zero-width identity behavior, absolute half-widths,
 * finite-value validation, and maximum-angle limiting.
 */
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "azBoxKernel validates widths coordinates and maxAz", "[improc::azBoxKernel]" )
{
    typedef mx::improc::azBoxKernel<mx::improc::eigenImage<double>, 2, mx::verbose::o> kernelT;
    mx::improc::eigenImage<double> kernel;

    SECTION( "zero half-widths produce an identity kernel" )
    {
        kernelT identityKernel( 0.0, 0.0 );
        REQUIRE( identityKernel.maxWidth() == 0 );
        REQUIRE( identityKernel.setKernel( 0.0, 0.0, kernel ) == mx::error_t::noerror );
        REQUIRE( kernel.rows() == 1 );
        REQUIRE( kernel.cols() == 1 );
        REQUIRE( kernel( 0, 0 ) == 1.0 );
        REQUIRE( identityKernel.setKernel( 4.0, -3.0, kernel ) == mx::error_t::noerror );
        REQUIRE( kernel( 0, 0 ) == 1.0 );
    }

    SECTION( "subpixel half-widths support a one-pixel odd kernel" )
    {
        typedef mx::improc::azBoxKernel<mx::improc::eigenImage<double>, 1, mx::verbose::o> subpixelKernelT;
        subpixelKernelT subpixelKernel( 0.25, 0.25 );

        REQUIRE( subpixelKernel.maxWidth() == 0 );
        REQUIRE( subpixelKernel.setKernel( 0.0, 0.0, kernel ) == mx::error_t::noerror );
        REQUIRE( kernel.rows() == 1 );
        REQUIRE( kernel.cols() == 1 );
        REQUIRE( kernel( 0, 0 ) == 1.0 );
    }

    SECTION( "negative half-widths retain their absolute-value contract" )
    {
        kernelT positiveWidths( 2.0, 4.0 );
        kernelT negativeWidths( -2.0, -4.0 );
        mx::improc::eigenImage<double> positiveKernel, negativeKernel;

        REQUIRE( positiveWidths.setKernel( 3.0, 1.0, positiveKernel ) == mx::error_t::noerror );
        REQUIRE( negativeWidths.setKernel( 3.0, 1.0, negativeKernel ) == mx::error_t::noerror );
        REQUIRE( ( positiveKernel == negativeKernel ).all() );
    }

    SECTION( "nonfinite widths maxAz and coordinates are rejected" )
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double infinity = std::numeric_limits<double>::infinity();
        kernelT invalidRadialWidth( nan, 1.0 );
        kernelT invalidAzimuthalWidth( 1.0, infinity );
        kernelT invalidMaxAz( 1.0, 1.0, nan );
        kernelT validKernel( 1.0, 1.0 );

        REQUIRE( invalidRadialWidth.setKernel( 0.0, 0.0, kernel ) == mx::error_t::invalidconfig );
        REQUIRE( kernel.size() == 0 );
        REQUIRE( invalidAzimuthalWidth.setKernel( 0.0, 0.0, kernel ) == mx::error_t::invalidconfig );
        REQUIRE( invalidMaxAz.setKernel( 0.0, 0.0, kernel ) == mx::error_t::invalidconfig );
        REQUIRE( validKernel.setKernel( nan, 0.0, kernel ) == mx::error_t::invalidarg );
        REQUIRE( validKernel.setKernel( 0.0, infinity, kernel ) == mx::error_t::invalidarg );
    }

    SECTION( "maxAz limits support around the positive x-axis at image center" )
    {
        kernelT unlimitedKernel( 2.0, 5.0 );
        kernelT equivalentUnlimitedKernel( 2.0, 5.0, 180.0 );
        kernelT limitedKernel( 2.0, 5.0, 30.0 );
        kernelT negativeLimitKernel( 2.0, 5.0, -30.0 );
        mx::improc::eigenImage<double> unlimited, equivalentUnlimited, limited, negativeLimit;

        REQUIRE( unlimitedKernel.setKernel( 0.0, 0.0, unlimited ) == mx::error_t::noerror );
        REQUIRE( equivalentUnlimitedKernel.setKernel( 0.0, 0.0, equivalentUnlimited ) == mx::error_t::noerror );
        REQUIRE( limitedKernel.setKernel( 0.0, 0.0, limited ) == mx::error_t::noerror );
        REQUIRE( negativeLimitKernel.setKernel( 0.0, 0.0, negativeLimit ) == mx::error_t::noerror );

        REQUIRE( arrayAllFinite( unlimited ) );
        REQUIRE( arrayAllFinite( limited ) );
        REQUIRE( unlimited.sum() == Approx( 1.0 ).margin( 1e-12 ) );
        REQUIRE( limited.sum() == Approx( 1.0 ).margin( 1e-12 ) );
        REQUIRE( ( unlimited == equivalentUnlimited ).all() );
        REQUIRE( ( limited == negativeLimit ).all() );
        REQUIRE( ( limited != 0.0 ).count() < ( unlimited != 0.0 ).count() );

        const double xCenter = 0.5 * ( limited.rows() - 1.0 );
        const double yCenter = 0.5 * ( limited.cols() - 1.0 );
        const double maxAngle = mx::math::dtor( 30.0 );
        for( int column = 0; column < limited.cols(); ++column )
        {
            for( int row = 0; row < limited.rows(); ++row )
            {
                if( limited( row, column ) > 0 )
                {
                    REQUIRE( std::fabs( std::atan2( column - yCenter, row - xCenter ) ) <= maxAngle + 1e-12 );
                }
            }
        }
    }

    SECTION( "finite dimension overflow is rejected" )
    {
        kernelT oversizedKernel( std::numeric_limits<double>::max(), 1.0 );

        REQUIRE( oversizedKernel.maxWidth() == 0 );
        REQUIRE( oversizedKernel.setKernel( 1.0, 0.0, kernel ) == mx::error_t::invalidconfig );
        REQUIRE( kernel.size() == 0 );
    }

    SECTION( "inconsistent cached width bounds are rejected" )
    {
        kernelT widthKernel( 3.0, 1.0 );
        widthKernel.m_maxWidth = 0;
        REQUIRE( widthKernel.setKernel( 1.0, 0.0, kernel ) == mx::error_t::sizeerr );

        kernelT heightKernel( 1.0, 3.0 );
        heightKernel.m_maxWidth = 2;
        REQUIRE( heightKernel.setKernel( 1.0, 0.0, kernel ) == mx::error_t::sizeerr );
    }

    SECTION( "a kernel with no accepted samples is rejected" )
    {
        kernelT emptyKernel( 1.0, 1.0 );
        emptyKernel.m_radWidth = -1.0;
        REQUIRE( emptyKernel.setKernel( 1.0, 0.0, kernel ) == mx::error_t::invalidconfig );
        REQUIRE( kernel.size() == 0 );
    }

    SECTION( "maxAz is evaluated away from the image center" )
    {
        kernelT limitedKernel( 2.0, 5.0, 30.0 );
        REQUIRE( limitedKernel.setKernel( 3.0, 2.0, kernel ) == mx::error_t::noerror );
        REQUIRE( arrayAllFinite( kernel ) );
        REQUIRE( kernel.sum() == Approx( 1.0 ).margin( 1e-12 ) );
    }
}

/// Verify precalcKernel handles odd, even, and non-square image centers without invalid cached kernels.
/** Exercises mx::improc::precalcKernel construction and mx::improc::precalcKernel::setKernel at exact and half-pixel
 * image centers.
 */
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "precalcKernel supports odd even and non-square image centers", "[improc::precalcKernel]" )
{
    typedef mx::improc::azBoxKernel<mx::improc::eigenImage<double>, 1, mx::verbose::o> kernelT;
    kernelT kernelGenerator( 2.0, 2.0 );
    mx::improc::precalcKernel<kernelT> oddKernelCache( kernelGenerator, 5, 7, 2.0, 3.0 );
    mx::improc::precalcKernel<kernelT> evenKernelCache( kernelGenerator, 4, 6, 1.5, 2.5 );

    REQUIRE( oddKernelCache.m_kernels.size() == 5 * 7 );
    REQUIRE( evenKernelCache.m_kernels.size() == 4 * 6 );

    for( const auto &cachedKernel : oddKernelCache.m_kernels )
    {
        REQUIRE( cachedKernel.size() > 0 );
        REQUIRE( arrayAllFinite( cachedKernel ) );
        REQUIRE( cachedKernel.sum() == Approx( 1.0 ).margin( 1e-12 ) );
    }
    for( const auto &cachedKernel : evenKernelCache.m_kernels )
    {
        REQUIRE( cachedKernel.size() > 0 );
        REQUIRE( arrayAllFinite( cachedKernel ) );
        REQUIRE( cachedKernel.sum() == Approx( 1.0 ).margin( 1e-12 ) );
    }

    mx::improc::eigenImage<double> centerKernel;
    REQUIRE( oddKernelCache.setKernel( 0.0, 0.0, centerKernel ) == mx::error_t::noerror );
    REQUIRE( arrayAllFinite( centerKernel ) );
    REQUIRE( evenKernelCache.setKernel( 0.5, 0.5, centerKernel ) == mx::error_t::noerror );
    REQUIRE( arrayAllFinite( centerKernel ) );
}

/// Verify precalcKernel propagates generation errors and rejects invalid lookup coordinates.
/** Exercises mx::improc::precalcKernel construction error propagation and mx::improc::precalcKernel::setKernel
 * coordinate bounds.
 */
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "precalcKernel propagates errors and validates lookup bounds", "[improc::precalcKernel]" )
{
    typedef mx::improc::azBoxKernel<mx::improc::eigenImage<double>, 1, mx::verbose::o> kernelT;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    try
    {
        kernelT invalidKernel( nan, 1.0 );
        mx::improc::precalcKernel<kernelT> invalidCache( invalidKernel, 3, 3, 1.0, 1.0 );
        static_cast<void>( invalidCache );
        FAIL( "precalcKernel ignored a setKernel failure" );
    }
    catch( const mx::exception<mx::verbose::o> &error )
    {
        REQUIRE( error.code() == mx::error_t::invalidconfig );
    }

    kernelT validKernel( 2.0, 2.0 );
    mx::improc::precalcKernel<kernelT> cache( validKernel, 3, 5, 1.0, 2.0 );
    mx::improc::eigenImage<double> kernel;

    REQUIRE( cache.setKernel( -1.0, -2.0, kernel ) == mx::error_t::noerror );
    REQUIRE( cache.setKernel( 1.0, 2.0, kernel ) == mx::error_t::noerror );
    REQUIRE( cache.setKernel( -2.0, 0.0, kernel ) == mx::error_t::invalidarg );
    REQUIRE( cache.setKernel( 2.0, 0.0, kernel ) == mx::error_t::invalidarg );
    REQUIRE( cache.setKernel( 0.0, -3.0, kernel ) == mx::error_t::invalidarg );
    REQUIRE( cache.setKernel( 0.0, 3.0, kernel ) == mx::error_t::invalidarg );
    REQUIRE( cache.setKernel( 0.25, 0.0, kernel ) == mx::error_t::invalidarg );
    REQUIRE( cache.setKernel( nan, 0.0, kernel ) == mx::error_t::invalidarg );

    mx::improc::precalcKernel<kernelT> emptyCache( validKernel, 0, 0, 0.0, 0.0 );
    REQUIRE( emptyCache.setKernel( 0.0, 0.0, kernel ) == mx::error_t::invalidconfig );
}

/** \brief Verifies the radial-profile subtraction used by HCI preprocessing.
 *
 * Exercises mx::improc::radprofim with its centered-radius overload and subtraction enabled, followed by the
 * mx::improc::zeroNaNs cleanup used by HCI.
 *
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "radprofim subtracts a centered radial background", "[improc::radprofim][hciReduce]" )
{
    mx::improc::eigenImage<double> image( 9, 9 );
    image.setConstant( 5.0 );
    mx::improc::eigenImage<double> profile;

    mx::improc::radprofim( profile, image, true );
    mx::improc::zeroNaNs( image );

    REQUIRE( profile.rows() == image.rows() );
    REQUIRE( profile.cols() == image.cols() );
    REQUIRE( arrayAllFinite( image ) );
    REQUIRE( image.abs().maxCoeff() < 1e-12 );
}

/** \brief Verifies that masked pixels are excluded from and preserved by radial-profile subtraction.
 *
 * Exercises the masked mx::improc::radprofim overload used when preprocessing must omit invalid pixels.
 *
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "radprofim excludes masked pixels", "[improc::radprofim][hciReduce]" )
{
    mx::improc::eigenImage<double> image( 9, 9 );
    image.setConstant( 5.0 );

    mx::improc::eigenImage<double> radius( 9, 9 );
    mx::improc::radiusImage( radius );

    mx::improc::eigenImage<double> mask( 9, 9 );
    mask.setOnes();
    mask( 4, 4 ) = 0;

    mx::improc::eigenImage<double> profile;
    mx::improc::radprofim( profile, image, radius, &mask, true );

    REQUIRE( profile.rows() == image.rows() );
    REQUIRE( profile.cols() == image.cols() );
    REQUIRE( profile( 4, 4 ) == 0.0 );
    REQUIRE( image( 4, 4 ) == 5.0 );

    mx::improc::zeroNaNs( image );
    REQUIRE( image( 4, 5 ) == Approx( 0.0 ).margin( 1e-12 ) );
}

/** \brief Verifies HCI azimuthal median filtering through a precomputed azimuthal kernel.
 *
 * Exercises mx::improc::medianFilterImage with mx::improc::precalcKernel and mx::improc::azBoxKernel.
 *
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "medianFilterImage handles HCI azimuthal kernels", "[improc::medianFilterImage][hciReduce]" )
{
    typedef mx::improc::azBoxKernel<mx::improc::eigenImage<double>, 1, mx::verbose::o> kernelT;

    kernelT kernelGenerator( 2.0, 2.0 );
    mx::improc::precalcKernel<kernelT> kernelCache( kernelGenerator, 9, 9, 4.0, 4.0 );

    mx::improc::eigenImage<double> image( 9, 9 );
    image.setConstant( 5.0 );
    image( 4, 4 ) = 100.0;

    mx::improc::eigenImage<double> filtered;
    mx::improc::medianFilterImage( filtered, image, kernelCache );

    REQUIRE( filtered.rows() == image.rows() );
    REQUIRE( filtered.cols() == image.cols() );
    REQUIRE( arrayAllFinite( filtered ) );
    REQUIRE( filtered( 4, 4 ) == 5.0 );
    REQUIRE( filtered( 0, 0 ) == 5.0 );
}

/// Verify that filterImage preserves non-square dimensions with a Gaussian kernel.
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "Gaussian filter handles non-square image dimensions", "[improc::imageFilters]" )
{
    mx::improc::eigenImage<float> im;
    mx::improc::eigenImage<float> fim;

    im.resize( 640, 480 );
    for( int r = 0; r < im.rows(); ++r )
    {
        for( int c = 0; c < im.cols(); ++c )
        {
            im( r, c ) = static_cast<float>( ( r + 2 * c ) % 17 );
        }
    }

    mx::improc::gaussKernel<mx::improc::eigenImage<float>, 2> kernel( 10.0f );
    auto err = mx::improc::filterImage( fim, im, kernel );

    REQUIRE( err == mx::error_t::noerror );
    REQUIRE( fim.rows() == im.rows() );
    REQUIRE( fim.cols() == im.cols() );
    REQUIRE( std::isfinite( fim( 0, 0 ) ) );
    REQUIRE( std::isfinite( fim( fim.rows() - 1, fim.cols() - 1 ) ) );
    REQUIRE( std::isfinite( fim( fim.rows() / 2, fim.cols() / 2 ) ) );
}

/// Verify that medianSmooth uses exactly the configured odd or even square window.
/** Exercises both mx::improc::medianSmooth overloads, including the arithmetic mean of the two central samples for an
 * even-sized window and the caller-owned untouched edge pixels.
 */
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "medianSmooth honors exact odd and even full widths", "[improc::medianSmooth]" )
{
    SECTION( "even width uses width squared samples" )
    {
        mx::improc::eigenImage<double> input( 4, 4 );
        mx::improc::eigenImage<double> output( 4, 4 );
        for( int column = 0; column < input.cols(); ++column )
        {
            for( int row = 0; row < input.rows(); ++row )
            {
                input( row, column ) = row + column * input.rows();
            }
        }
        output.setConstant( -99.0 );

        int xMax = 17;
        int yMax = 19;
        double pMax = 23.0;
        const int result = mx::improc::medianSmooth( output, xMax, yMax, pMax, input, 4 );

        REQUIRE( result == 0 );
        REQUIRE( output( 2, 2 ) == Approx( 7.5 ) );
        REQUIRE( xMax == 2 );
        REQUIRE( yMax == 2 );
        REQUIRE( pMax == Approx( 7.5 ) );
        for( int column = 0; column < output.cols(); ++column )
        {
            for( int row = 0; row < output.rows(); ++row )
            {
                if( row != 2 || column != 2 )
                {
                    REQUIRE( output( row, column ) == -99.0 );
                }
            }
        }
    }

    SECTION( "odd width retains its centered footprint" )
    {
        mx::improc::eigenImage<double> input( 3, 3 );
        mx::improc::eigenImage<double> output( 3, 3 );
        for( int column = 0; column < input.cols(); ++column )
        {
            for( int row = 0; row < input.rows(); ++row )
            {
                input( row, column ) = row + column * input.rows();
            }
        }
        output.setConstant( -99.0 );

        REQUIRE( mx::improc::medianSmooth( output, input, 3 ) == 0 );
        REQUIRE( output( 1, 1 ) == Approx( 4.0 ) );
        REQUIRE( output( 0, 0 ) == -99.0 );
        REQUIRE( output( 2, 2 ) == -99.0 );
    }

    SECTION( "unit width is an identity operation" )
    {
        mx::improc::eigenImage<double> input( 2, 3 );
        input << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
        mx::improc::eigenImage<double> output( 2, 3 );
        output.setZero();

        REQUIRE( mx::improc::medianSmooth( output, input, 1 ) == 0 );
        REQUIRE( ( output == input ).all() );
    }
}

/// Verify that medianSmooth rejects invalid dimensions without modifying caller-owned pixels.
/** Exercises mx::improc::medianSmooth validation and deterministic maximum outputs. */
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "medianSmooth validates widths and output dimensions", "[improc::medianSmooth]" )
{
    mx::improc::eigenImage<double> input( 3, 4 );
    mx::improc::eigenImage<double> output( 3, 4 );
    input.setZero();
    output.setConstant( -42.0 );

    int xMax = 17;
    int yMax = 19;
    double pMax = 23.0;

    REQUIRE( mx::improc::medianSmooth( output, xMax, yMax, pMax, input, 0 ) == -1 );
    REQUIRE( xMax == -1 );
    REQUIRE( yMax == -1 );
    REQUIRE( pMax == std::numeric_limits<double>::lowest() );
    REQUIRE( ( output == -42.0 ).all() );

    REQUIRE( mx::improc::medianSmooth( output, xMax, yMax, pMax, input, -1 ) == -1 );
    REQUIRE( mx::improc::medianSmooth( output, xMax, yMax, pMax, input, 4 ) == -1 );
    REQUIRE( ( output == -42.0 ).all() );

    mx::improc::eigenImage<double> wrongSize( 2, 4 );
    wrongSize.setConstant( -42.0 );
    REQUIRE( mx::improc::medianSmooth( wrongSize, xMax, yMax, pMax, input, 3 ) == -1 );
    REQUIRE( ( wrongSize == -42.0 ).all() );
}

/// Verify that meanSmooth uses an exact even-width footprint.
/** Exercises both mx::improc::meanSmooth overloads with the same even-window placement used by medianSmooth. */
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "meanSmooth honors exact even full widths", "[improc::meanSmooth]" )
{
    mx::improc::eigenImage<double> input( 4, 4 );
    mx::improc::eigenImage<double> output( 4, 4 );
    for( int column = 0; column < input.cols(); ++column )
    {
        for( int row = 0; row < input.rows(); ++row )
        {
            input( row, column ) = row + column * input.rows();
        }
    }
    output.setConstant( -99.0 );

    REQUIRE( mx::improc::meanSmooth( output, input, 4 ) == 0 );
    REQUIRE( output( 2, 2 ) == Approx( 7.5 ) );
    REQUIRE( output( 1, 1 ) == -99.0 );

    output.setConstant( -99.0 );
    int xMax = 17;
    int yMax = 19;
    double pMax = 23.0;
    REQUIRE( mx::improc::meanSmooth( output, xMax, yMax, pMax, input, 4, true ) == 0 );
    REQUIRE( output( 2, 2 ) == Approx( 7.5 ) );
    REQUIRE( xMax == 2 );
    REQUIRE( yMax == 2 );
    REQUIRE( pMax == Approx( 7.5 ) );

    SECTION( "reject extrema in an all-positive window" )
    {
        mx::improc::eigenImage<double> positiveInput( 2, 2 );
        positiveInput << 1.0, 2.0, 3.0, 100.0;
        mx::improc::eigenImage<double> positiveOutput( 2, 2 );
        positiveOutput.setConstant( -99.0 );

        REQUIRE( mx::improc::meanSmooth( positiveOutput, positiveInput, 2, true ) == 0 );
        REQUIRE( positiveOutput( 1, 1 ) == Approx( 2.5 ) );
        REQUIRE( positiveOutput( 0, 0 ) == -99.0 );
    }

    SECTION( "reject extrema in an all-negative window" )
    {
        mx::improc::eigenImage<double> negativeInput( 2, 2 );
        negativeInput << -100.0, -3.0, -2.0, -1.0;
        mx::improc::eigenImage<double> negativeOutput( 2, 2 );
        negativeOutput.setConstant( -99.0 );

        xMax = 17;
        yMax = 19;
        pMax = 23.0;
        REQUIRE( mx::improc::meanSmooth( negativeOutput, xMax, yMax, pMax, negativeInput, 2, true ) == 0 );
        REQUIRE( negativeOutput( 1, 1 ) == Approx( -2.5 ) );
        REQUIRE( negativeOutput( 0, 0 ) == -99.0 );
        REQUIRE( xMax == 1 );
        REQUIRE( yMax == 1 );
        REQUIRE( pMax == Approx( -2.5 ) );
    }
}

/// Verify that vectorSmoothMedian uses no more than the configured even width.
/** Exercises mx::math::vectorSmoothMedian at interior and truncated edge locations. */
/**
 * \ingroup vectorUtils_unit_tests
 */
TEST_CASE( "vectorSmoothMedian honors exact even full widths", "[math::vectorSmoothMedian]" )
{
    std::vector<double> input{ 0.0, 1.0, 2.0, 3.0, 100.0 };
    const std::vector<double> original = input;
    std::vector<double> output;

    REQUIRE( mx::math::vectorSmoothMedian( output, input, 4 ) == 0 );
    REQUIRE( input == original );
    REQUIRE( output.size() == input.size() );
    REQUIRE( output[0] == Approx( 0.5 ) );
    REQUIRE( output[1] == Approx( 1.0 ) );
    REQUIRE( output[2] == Approx( 1.5 ) );
    REQUIRE( output[3] == Approx( 2.5 ) );
    REQUIRE( output[4] == Approx( 3.0 ) );

    output.assign( 1, -99.0 );
    REQUIRE( mx::math::vectorSmoothMedian( output, input, 0 ) == -1 );
    REQUIRE( output == std::vector<double>{ -99.0 } );
}

} // namespace imageFiltersTest
} // namespace improcTest
} // namespace unitTest
