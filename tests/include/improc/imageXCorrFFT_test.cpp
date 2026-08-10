/** \file imageXCorrFFT_test.cpp
 * \brief Tests FFT-based image cross-correlation.
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"
#include "../../../include/improc/imageXCorrFFT.hpp"
#include "../../../include/improc/eigenCube.hpp"

/** \cond
 * Explicit instantiation compile-checks FFT correlation with the canonical double-precision Eigen image and emits
 * its header-defined methods for coverage accounting in this test translation unit.
 */
template class mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>>;
/** \endcond */

/** centroiding Gaussians with center of light
 *
 * Verify center of light calculation
 *
 */
/**
 * \ingroup imageXCorrFFT_unit_tests
 */
TEST_CASE( "Image cross-correlation with FFT using center of light", "[improc::imageXCorrFFT]" )
{
    GIVEN( "two Gaussians" )
    {
        WHEN( "ref at geometric center, equal even sizes, shift=(+4,+4)" )
        {
            double xshift = 4;
            double yshift = 4;
            int xsz = 64;
            int ysz = 64;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);
            xcf.peakMethod( mx::improc::xcorrPeakMethod::centerOfLight );
            xcf.maxLag( 20 );
            xcf.refIm( refIm );
            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal odd sizes, shift=(+4,+4)" )
        {
            double xshift = 4;
            double yshift = 4;
            int xsz = 65;
            int ysz = 65;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );

            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::centerOfLight );
            xcf.maxLag( 20 );

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);
            xcf.refIm( refIm );

            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal even sizes, shift=(-3.6,+2.25)" )
        {
            double xshift = -3.6;
            double yshift = 2.25;
            int xsz = 64;
            int ysz = 64;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::centerOfLight );
            xcf.maxLag( 20 );

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);
            xcf.refIm( refIm );

            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal odd sizes, shift=(+1.3,-0.6)" )
        {
            double xshift = 1.3;
            double yshift = -0.6;
            int xsz = 65;
            int ysz = 65;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::centerOfLight );
            xcf.maxLag( 20 );

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);
            xcf.refIm( refIm );

            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
    }
}

/** centroiding by magnification
 *
 * Verify magnification peak finding
 *
 */
/**
 * \ingroup imageXCorrFFT_unit_tests
 */
TEST_CASE( "Image cross-correlation with FFT using magnification peak finding", "[improc::imageXCorrFFT]" )
{
    GIVEN( "two Gaussians" )
    {
        WHEN( "ref at geometric center, equal even sizes, shift=(+4,+4)" )
        {
            double xshift = 4;
            double yshift = 4;
            int xsz = 64;
            int ysz = 64;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::interpPeak );
            xcf.maxLag( 5 );

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);

            xcf.refIm( refIm );
            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal odd sizes, shift=(+4,+4)" )
        {
            double xshift = 4;
            double yshift = 4;
            int xsz = 65;
            int ysz = 65;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);

            xcf.peakMethod( mx::improc::xcorrPeakMethod::interpPeak );
            xcf.maxLag( 5 );

            xcf.refIm( refIm );
            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal even sizes, shift=(-3.6,+2.3)" )
        {
            double xshift = -3.6;
            double yshift = 2.3;
            int xsz = 64;
            int ysz = 64;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::interpPeak );
            xcf.tol( 0.1 );

            mx::improc::eigenImage<double> refIm = im0;
            xcf.refIm( refIm );

            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal odd sizes, shift=(+1.3,-0.6)" )
        {
            double xshift = 1.3;
            double yshift = -0.6;
            int xsz = 65;
            int ysz = 65;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::interpPeak );
            xcf.maxLag( 5 );

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);
            xcf.refIm( refIm );

            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
    }
}

/// Verify imageXCorrFFT displacement recovery with Gaussian peak fitting.
/**
 * \ingroup imageXCorrFFT_unit_tests
 */
TEST_CASE( "Image cross-correlation with FFT using Gaussian peak fit", "[improc::imageXCorrFFT]" )
{
    GIVEN( "two Gaussians" )
    {
        WHEN( "ref at geometric center, equal even sizes, shift=(+4,+4)" )
        {
            double xshift = 4;
            double yshift = 4;
            int xsz = 64;
            int ysz = 64;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;

            mx::improc::eigenImage<double> refIm = im0;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::gaussFit );
            xcf.maxLag( 32 );
            xcf.refIm( refIm );
            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal odd sizes, shift=(+4,+4)" )
        {
            double xshift = 4;
            double yshift = 4;
            int xsz = 65;
            int ysz = 65;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::gaussFit );
            xcf.maxLag( 32 );

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);
            xcf.refIm( refIm );

            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal even sizes, shift=(-3.6,+2.25)" )
        {
            double xshift = -3.6;
            double yshift = 2.25;
            int xsz = 64;
            int ysz = 64;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);
            xcf.peakMethod( mx::improc::xcorrPeakMethod::gaussFit );
            xcf.maxLag( 32 );
            xcf.refIm( refIm );

            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
        WHEN( "ref at geometric center, equal odd sizes, shift=(+1.3,-0.6)" )
        {
            double xshift = 1.3;
            double yshift = -0.6;
            int xsz = 65;
            int ysz = 65;

            mx::improc::eigenImage<double> im0, im2;
            im0.resize( xsz, ysz );
            im2.resize( xsz, ysz );

            double xcen = 0.5 * ( im0.rows() - 1 );
            double ycen = 0.5 * ( im0.cols() - 1 );

            mx::math::func::gaussian2D<double>( im0.data(), im0.rows(), im0.cols(), 0., 1.0, xcen, ycen, 2 );
            mx::math::func::gaussian2D<double>( im2.data(),
                                                im2.rows(),
                                                im2.cols(),
                                                0.,
                                                1.0,
                                                xcen + xshift,
                                                ycen + yshift,
                                                2 );

            double x, y, peak;
            mx::improc::imageXCorrFFT<mx::improc::eigenImage<double>> xcf;
            xcf.peakMethod( mx::improc::xcorrPeakMethod::gaussFit );
            xcf.maxLag( 32 );

            mx::improc::eigenImage<double> refIm = im0; //.block(10,10, im0.rows()-11, im0.cols()-11);
            xcf.refIm( refIm );

            xcf( x, y, peak, im2 );

            REQUIRE( x == Approx( xshift ) );
            REQUIRE( y == Approx( yshift ) );
            REQUIRE( peak > 0 );
        }
    }
}
