/** \file imageXCorrDiscrete_test.cpp
 * \brief Tests discrete image cross-correlation.
 */
#include "../../catch2/catch.hpp"

#include <cmath>
#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"
#include "../../../include/improc/imageXCorrDiscrete.hpp"
#include "../../../include/improc/eigenCube.hpp"

/** \cond
 * Explicit instantiation compile-checks discrete correlation with the canonical double-precision Eigen image and
 * emits its header-defined methods for coverage accounting in this test translation unit.
 */
template class mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>>;
/** \endcond */

/** Recovering Gaussian displacement and centroid
 *
 * Verify imageXCorrDiscrete interpolation and imageCenterOfLight calculation.
 *
 */
/**
 * \ingroup imageXCorrDiscrete_unit_tests
 */
TEST_CASE( "Recover Gaussian displacement and centroid", "[improc::imageXCorrDiscrete]" )
{
    GIVEN( "two Gaussians" )
    {
        WHEN( "the target is shifted from a centered reference" )
        {
            mx::improc::eigenImage<double> im0, im2;
            im0.resize( 64, 64 );
            im2.resize( 64, 64 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, 31.5, 31.5, 2 );
            mx::math::func::gaussian2D<double>( im2.data(), im2.rows(), im2.cols(), 0., 1.0, 31.5 + 4, 31.5 + 4, 2 );

            double x, y;
            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf( 5 );

            mx::improc::eigenImage<double> refIm = im0.block( 10, 10, im0.rows() - 20, im0.cols() - 20 );
            REQUIRE( xcf.refIm( refIm ) == 0 );
            xcf.m_peakMethod = mx::improc::xcorrPeakMethod::centroid;

            REQUIRE( xcf( x, y, im2 ) == 0 );

            int xPeak = -1;
            int yPeak = -1;
            xcf.ccIm().maxCoeff( &xPeak, &yPeak );
            REQUIRE( xPeak == 9 );
            REQUIRE( yPeak == 9 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 4, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 4, 1e-8 ) );
            REQUIRE( std::isfinite( x ) );
            REQUIRE( std::isfinite( y ) );
        }

        WHEN( "the target has an asymmetric signed shift" )
        {
            mx::improc::eigenImage<double> centered( 64, 64 );
            mx::improc::eigenImage<double> shifted( 64, 64 );
            mx::math::func::gaussian2D<double>(
                centered.data(), centered.rows(), centered.cols(), 0., 1.0, 31.5, 31.5, 2 );
            mx::math::func::gaussian2D<double>(
                shifted.data(), shifted.rows(), shifted.cols(), 0., 1.0, 35.5, 28.5, 2 );

            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf( 5 );
            mx::improc::eigenImage<double> reference = centered.block( 10, 10, 44, 44 );
            REQUIRE( xcf.refIm( reference ) == 0 );
            xcf.m_peakMethod = mx::improc::xcorrPeakMethod::centroid;

            double x = 0;
            double y = 0;
            REQUIRE( xcf( x, y, shifted ) == 0 );

            int xPeak = -1;
            int yPeak = -1;
            xcf.ccIm().maxCoeff( &xPeak, &yPeak );
            REQUIRE( xPeak == 9 );
            REQUIRE( yPeak == 2 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 4, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( -3, 1e-8 ) );
            REQUIRE( std::isfinite( x ) );
            REQUIRE( std::isfinite( y ) );
        }

        WHEN( "the target has zero shift" )
        {
            mx::improc::eigenImage<double> target( 64, 64 );
            mx::math::func::gaussian2D<double>(
                target.data(), target.rows(), target.cols(), 0., 1.0, 31.5, 31.5, 2 );

            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf( 5 );
            mx::improc::eigenImage<double> reference = target.block( 10, 10, 44, 44 );
            REQUIRE( xcf.refIm( reference ) == 0 );
            xcf.m_peakMethod = mx::improc::xcorrPeakMethod::centroid;

            double x = 0;
            double y = 0;
            REQUIRE( xcf( x, y, target ) == 0 );

            int xPeak = -1;
            int yPeak = -1;
            xcf.ccIm().maxCoeff( &xPeak, &yPeak );
            REQUIRE( xPeak == 5 );
            REQUIRE( yPeak == 5 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 0, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 0, 1e-8 ) );
            REQUIRE( std::isfinite( x ) );
            REQUIRE( std::isfinite( y ) );
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

        WHEN( "automatic lags span only valid target blocks" )
        {
            mx::improc::eigenImage<double> target( 9, 10 );
            mx::math::func::gaussian2D<double>(
                target.data(), target.rows(), target.cols(), 0., 1.0, 4.0, 4.5, 1.5 );

            mx::improc::eigenImage<double> reference = target.block( 2, 2, 4, 5 );
            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf;
            REQUIRE( xcf.refIm( reference ) == 0 );

            double x, y;
            REQUIRE( xcf( x, y, target ) == 0 );
            REQUIRE( xcf.ccIm().rows() == 5 );
            REQUIRE( xcf.ccIm().cols() == 5 );
            int xPeak = -1;
            int yPeak = -1;
            xcf.ccIm().maxCoeff( &xPeak, &yPeak );
            REQUIRE( xPeak == 2 );
            REQUIRE( yPeak == 2 );
            REQUIRE( std::isfinite( x ) );
            REQUIRE( std::isfinite( y ) );
        }

        WHEN( "interpolation uses the lag range allowed by the image dimensions" )
        {
            mx::improc::eigenImage<double> centered( 64, 64 );
            mx::improc::eigenImage<double> shifted( 64, 64 );
            mx::math::func::gaussian2D<double>(
                centered.data(), centered.rows(), centered.cols(), 0., 1.0, 31.5, 31.5, 2 );
            mx::math::func::gaussian2D<double>(
                shifted.data(), shifted.rows(), shifted.cols(), 0., 1.0, 32.5, 31.5, 2 );

            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf( 5 );
            mx::improc::eigenImage<double> reference = centered.block( 2, 2, 60, 60 );
            REQUIRE( xcf.refIm( reference ) == 0 );
            xcf.m_peakMethod = mx::improc::xcorrPeakMethod::interp;

            double x = 0;
            double y = 0;
            REQUIRE( xcf( x, y, shifted ) == 0 );
            REQUIRE( xcf.ccIm().rows() == 5 );
            REQUIRE( xcf.ccIm().cols() == 5 );
            REQUIRE( xcf.magIm().rows() == 41 );
            REQUIRE( xcf.magIm().cols() == 41 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 1, 0.11 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 0, 0.11 ) );
            REQUIRE( std::isfinite( x ) );
            REQUIRE( std::isfinite( y ) );
        }

        WHEN( "interpolation falls back safely for a three-sample lag range" )
        {
            mx::improc::eigenImage<double> centered( 64, 64 );
            mx::improc::eigenImage<double> shifted( 64, 64 );
            mx::math::func::gaussian2D<double>(
                centered.data(), centered.rows(), centered.cols(), 0., 1.0, 31.5, 31.5, 2 );
            mx::math::func::gaussian2D<double>(
                shifted.data(), shifted.rows(), shifted.cols(), 0., 1.0, 32.5, 31.5, 2 );

            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf( 5 );
            mx::improc::eigenImage<double> reference = centered.block( 1, 1, 62, 62 );
            REQUIRE( xcf.refIm( reference ) == 0 );
            xcf.m_peakMethod = mx::improc::xcorrPeakMethod::interp;

            double x = 0;
            double y = 0;
            REQUIRE( xcf( x, y, centered ) == 0 );
            REQUIRE( xcf.ccIm().rows() == 3 );
            REQUIRE( xcf.ccIm().cols() == 3 );
            REQUIRE( xcf.magIm().size() == 0 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 0, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 0, 1e-8 ) );

            REQUIRE( xcf( x, y, shifted ) == 0 );
            REQUIRE( xcf.magIm().size() == 0 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 1, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 0, 1e-8 ) );
            REQUIRE( std::isfinite( x ) );
            REQUIRE( std::isfinite( y ) );
        }

        WHEN( "interpolation falls back safely for a boundary correlation peak" )
        {
            mx::improc::eigenImage<double> centered( 64, 64 );
            mx::improc::eigenImage<double> shiftedNegative( 64, 64 );
            mx::improc::eigenImage<double> shiftedPositive( 64, 64 );
            mx::math::func::gaussian2D<double>(
                centered.data(), centered.rows(), centered.cols(), 0., 1.0, 31.5, 31.5, 2 );
            mx::math::func::gaussian2D<double>( shiftedNegative.data(),
                                                shiftedNegative.rows(),
                                                shiftedNegative.cols(),
                                                0.,
                                                1.0,
                                                29.5,
                                                31.5,
                                                2 );
            mx::math::func::gaussian2D<double>( shiftedPositive.data(),
                                                shiftedPositive.rows(),
                                                shiftedPositive.cols(),
                                                0.,
                                                1.0,
                                                33.5,
                                                31.5,
                                                2 );

            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf( 5 );
            mx::improc::eigenImage<double> reference = centered.block( 2, 2, 60, 60 );
            REQUIRE( xcf.refIm( reference ) == 0 );
            xcf.m_peakMethod = mx::improc::xcorrPeakMethod::interp;

            double x = 0;
            double y = 0;
            REQUIRE( xcf( x, y, shiftedNegative ) == 0 );
            REQUIRE( xcf.ccIm().rows() == 5 );
            REQUIRE( xcf.ccIm().cols() == 5 );
            REQUIRE( xcf.magIm().size() == 0 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( -2, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 0, 1e-8 ) );

            REQUIRE( xcf( x, y, shiftedPositive ) == 0 );
            REQUIRE( xcf.magIm().size() == 0 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( 2, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 0, 1e-8 ) );
            REQUIRE( std::isfinite( x ) );
            REQUIRE( std::isfinite( y ) );
        }

        WHEN( "interpolation falls back safely for a negative correlation surface" )
        {
            mx::improc::eigenImage<double> target( 7, 7 );
            mx::improc::eigenImage<double> reference( 3, 3 );
            for( int row = 0; row < target.rows(); ++row )
            {
                for( int col = 0; col < target.cols(); ++col )
                {
                    target( row, col ) = row + 0.2 * col + 1.2 * std::sin( row + 0.3 * col );
                }
            }
            for( int row = 0; row < reference.rows(); ++row )
            {
                for( int col = 0; col < reference.cols(); ++col )
                {
                    reference( row, col ) = -( row + 0.2 * col );
                }
            }

            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf( 5 );
            REQUIRE( xcf.refIm( reference ) == 0 );
            xcf.m_peakMethod = mx::improc::xcorrPeakMethod::interp;

            double x = 0;
            double y = 0;
            REQUIRE( xcf( x, y, target ) == 0 );
            REQUIRE( xcf.ccIm().rows() == 5 );
            REQUIRE( xcf.ccIm().cols() == 5 );

            int xPeak = -1;
            int yPeak = -1;
            REQUIRE( xcf.ccIm().maxCoeff( &xPeak, &yPeak ) < 0 );
            REQUIRE( xPeak == 1 );
            REQUIRE( yPeak == 3 );
            REQUIRE( xcf.magIm().size() == 0 );
            REQUIRE_THAT( x, Catch::Matchers::WithinAbs( -1, 1e-8 ) );
            REQUIRE_THAT( y, Catch::Matchers::WithinAbs( 1, 1e-8 ) );
            REQUIRE( std::isfinite( x ) );
            REQUIRE( std::isfinite( y ) );
        }

        WHEN( "a mask differs from the reference in exactly one dimension" )
        {
            mx::improc::eigenImage<double> reference( 5, 6 );
            reference.setOnes();

            mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf;

            mx::improc::eigenImage<double> wrongRows( 4, 6 );
            wrongRows.setOnes();
            REQUIRE( xcf.maskIm( wrongRows ) == 0 );
            REQUIRE( xcf.refIm( reference ) == -1 );

            mx::improc::eigenImage<double> wrongCols( 5, 7 );
            wrongCols.setOnes();
            REQUIRE( xcf.maskIm( wrongCols ) == 0 );
            REQUIRE( xcf.refIm( reference ) == -1 );
        }
    }
}
