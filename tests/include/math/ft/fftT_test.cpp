/** \file fftT_test.cpp
 * \brief Tests CPU FFT transforms.
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include <iostream>

#include "../../../../include/math/ft/fftT.hpp"
#include "../../../../include/improc/eigenImage.hpp"

/// 1D FFT with FFTW
/**
 * \ingroup fftT_unit_tests
 */
TEST_CASE( "1D c2c FFT with FFTW, float", "[math::ft]" )
{
    srand( 1 );

    SECTION( "out-of-place, forward, default constructed, raw interface" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 1> fft;

        fft.plan( 1024 );

        mx::improc::eigenImage<std::complex<float>> in( 1024, 1 ), out( 1024, 1 );
        in.setRandom();

        fft( out.data(), in.data() );

        float sin = in.abs2().sum();
        float sout = out.abs2().sum() / out.rows();

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinRel( sout, 5e-6f ) );
    }

    SECTION( "out-of-place, forward, default constructed, eigen interface" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 1> fft;

        fft.plan( 1024 );

        mx::improc::eigenImage<std::complex<float>> in( 1024, 1 ), out( 1024, 1 );
        in.setRandom();

        fft( out, in );

        float sin = in.abs2().sum();
        float sout = out.abs2().sum() / out.rows();

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinRel( sout, 5e-6f ) );
    }

    SECTION( "out-of-place, forward, plan constructor" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 1> fft( 1024 );

        mx::improc::eigenImage<std::complex<float>> in( 1024, 1 ), out( 1024, 1 );
        in.setRandom();

        fft( out, in );

        float sin = in.abs2().sum();
        float sout = out.abs2().sum() / out.rows();

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinRel( sout, 5e-6f ) );
    }

    SECTION( "out-of-place, backward" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 1> fft( 1024, mx::math::ft::dir::backward );

        mx::improc::eigenImage<std::complex<float>> in( 1024, 1 ), out( 1024, 1 );
        out.setRandom();

        fft( in, out );

        float sin = in.abs2().sum() / in.rows();
        float sout = out.abs2().sum();

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinRel( sout, 5e-6f ) );
    }

    SECTION( "in-place, forward" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 1> fft( 1024, mx::math::ft::dir::forward, true );

        mx::improc::eigenImage<std::complex<float>> in( 1024, 1 );
        in.setRandom();
        mx::improc::eigenImage<std::complex<float>> incheck = in;
        float sin = in.abs2().sum();

        fft( in, in );

        float sout = in.abs2().sum() / in.rows();

        float rmsdiff = ( in - incheck ).abs2().sum();

        // Make sure it isn't ident
        REQUIRE( rmsdiff > 0 );

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinRel( sout, 5e-6f ) );
    }

    SECTION( "in-place, backward" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 1> fft( 1024, mx::math::ft::dir::backward, true );

        mx::improc::eigenImage<std::complex<float>> in( 1024, 1 );
        in.setRandom();

        mx::improc::eigenImage<std::complex<float>> incheck = in;

        float sin = in.abs2().sum();

        fft( in, in );

        float sout = in.abs2().sum() / in.rows();

        float rmsdiff = ( in - incheck ).abs2().sum();

        // Make sure it isn't ident
        REQUIRE( rmsdiff > 0 );

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinRel( sout, 5e-6f ) );
    }
}

/// 2D FFT with FFTW
/**
 * \ingroup fftT_unit_tests
 */
TEST_CASE( "2D c2c FFT with FFTW, float", "[math::ft]" )
{
    srand( 1 );

    SECTION( "out-of-place, forward, default constructed, raw interface" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 2> fft;

        fft.plan( 128, 128 );

        mx::improc::eigenImage<std::complex<float>> in( 128, 128 ), out( 128, 128 );
        in.setRandom();

        fft( out.data(), in.data() );

        float sin = in.abs2().sum() * ( out.rows() * out.cols() );
        float sout = out.abs2().sum();

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinAbs( sout, ( sin ) * ( 1e-3 ) ) );
    }

    SECTION( "out-of-place, forward, default constructed, eigen interface" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 2> fft;

        fft.plan( 128, 128 );

        mx::improc::eigenImage<std::complex<float>> in( 128, 128 ), out( 128, 128 );
        in.setRandom();

        fft( out, in );

        float sin = in.abs2().sum();
        float sout = out.abs2().sum() / ( out.rows() * out.cols() );

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinAbs( sout, ( sin ) * ( 1e-3 ) ) );
    }

    SECTION( "out-of-place, forward, plan constructor" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 2> fft( 128, 128 );

        mx::improc::eigenImage<std::complex<float>> in( 128, 128 ), out( 128, 128 );
        ;
        in.setRandom();

        fft( out, in );

        float sin = in.abs2().sum();
        float sout = out.abs2().sum() / ( in.rows() * out.cols() );

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinAbs( sout, ( sin ) * ( 1e-3 ) ) );
    }

    SECTION( "out-of-place, backward" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 2> fft( 128, 128, mx::math::ft::dir::backward );

        mx::improc::eigenImage<std::complex<float>> in( 128, 128 ), out( 128, 128 );
        out.setRandom();

        fft( in, out );

        float sin = in.abs2().sum() / ( in.rows() * in.cols() );
        float sout = out.abs2().sum();

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinAbs( sout, ( sin ) * ( 1e-3 ) ) );
    }

    SECTION( "in-place, forward" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 2> fft( 128,
                                                                             128,
                                                                             mx::math::ft::dir::forward,
                                                                             true );

        mx::improc::eigenImage<std::complex<float>> in( 128, 128 );
        in.setRandom();
        mx::improc::eigenImage<std::complex<float>> incheck = in;
        float sin = in.abs2().sum();

        fft( in, in );

        float sout = in.abs2().sum() / ( in.rows() * in.cols() );

        float rmsdiff = ( in - incheck ).abs2().sum();

        // Make sure it isn't ident
        REQUIRE( rmsdiff > 0 );

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinAbs( sout, ( sin ) * ( 1e-3 ) ) );
    }

    SECTION( "in-place, backward" )
    {
        mx::math::ft::fftT<std::complex<float>, std::complex<float>, 2> fft( 128,
                                                                             128,
                                                                             mx::math::ft::dir::backward,
                                                                             true );

        mx::improc::eigenImage<std::complex<float>> in( 128, 128 );
        in.setRandom();

        mx::improc::eigenImage<std::complex<float>> incheck = in;

        float sin = in.abs2().sum();

        fft( in, in );

        float sout = in.abs2().sum() / ( in.rows() * in.cols() );

        float rmsdiff = ( in - incheck ).abs2().sum();

        // Make sure it isn't ident
        REQUIRE( rmsdiff > 0 );

        // Test by Parsevals
        REQUIRE_THAT( sin, Catch::Matchers::WithinAbs( sout, ( sin ) * ( 1e-3 ) ) );
    }
}
#ifdef HASQUAD

SECTION( "quad" )
{
}

#endif // HASQUD
