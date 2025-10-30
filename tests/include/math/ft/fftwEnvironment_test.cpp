/** \file fftwEnvironment_test.cpp
 */
#include "../../../catch2/catch.hpp"

#include <stdlib.h>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../../include/math/ft/fftwEnvironment.hpp"
#include "../../../../include/math/ft/fftwTemplates.hpp"
#include "../../../../include/ioutils/fileUtils.hpp"

/// Create and destroy fftwEnvironment, no threads, environment set
/**
 * \ingroup math_ft_fftwEnvironment_test
 */
TEST_CASE( "Create and destroy fftwEnvironment, no threads, environment set", "[math::ft]" )
{
    int rv = setenv( "MXFFTW_WISDOM", "/tmp/.fftw", 1 );

    REQUIRE( rv == 0 );

    SECTION( "float" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<float>();
        REQUIRE( fname == "/tmp/.fftw/fftw_wisdom.float" );
        remove( fname.c_str() );
        mx::ioutils::createDirectories( "/tmp/.fftw" );

        {
            mx::math::ft::fftwEnvironment<float> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    SECTION( "double" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<double>();
        REQUIRE( fname == "/tmp/.fftw/fftw_wisdom.double" );
        remove( fname.c_str() );
        mx::ioutils::createDirectories( "/tmp/.fftw" );

        {
            mx::math::ft::fftwEnvironment<double> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    SECTION( "long double" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<long double>();
        REQUIRE( fname == "/tmp/.fftw/fftw_wisdom.long_double" );
        remove( fname.c_str() );
        mx::ioutils::createDirectories( "/tmp/.fftw" );

        {
            mx::math::ft::fftwEnvironment<long double> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    // clang-format off
    #ifdef HASQUAD // clang-format on

    SECTION( "quad" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<__float128>();
        REQUIRE( fname == "/tmp/.fftw/fftw_wisdom.quad" );
        remove( fname.c_str() );
        mx::ioutils::createDirectories( "/tmp/.fftw" );

        {
            mx::math::ft::fftwEnvironment<__float128> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    // clang-format off
    #endif // clang-format on
}

/// Create and destroy fftwEnvironment, no threads, environment not set
/**
 * \ingroup math_ft_fftwEnvironment_test
 */
TEST_CASE( "Create and destroy fftwEnvironment, no threads, environment not set", "[math::ft]" )
{
    int rv = setenv( "MXFFTW_WISDOM", "", 1 );

    REQUIRE( rv == 0 );

    SECTION( "float" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<float>();
        REQUIRE( fname == "./.fftw_wisdom.float" );
        remove( fname.c_str() );

        {
            mx::math::ft::fftwEnvironment<float> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    SECTION( "double" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<double>();
        REQUIRE( fname == "./.fftw_wisdom.double" );
        remove( fname.c_str() );

        {
            mx::math::ft::fftwEnvironment<double> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    SECTION( "long double" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<long double>();
        REQUIRE( fname == "./.fftw_wisdom.long_double" );
        remove( fname.c_str() );

        {
            mx::math::ft::fftwEnvironment<long double> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    // clang-format off
    #ifdef HASQUAD // clang-format on

    SECTION( "quad" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<__float128>();
        REQUIRE( fname == "./.fftw_wisdom.quad" );
        remove( fname.c_str() );

        {
            mx::math::ft::fftwEnvironment<__float128> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    // clang-format off
    #endif // clang-format on
}

/// Create and destroy fftwEnvironment, with threads, environment set
/**
 * \ingroup math_ft_fftwEnvironment_test
 */
TEST_CASE( "Create and destroy fftwEnvironment, with threads, environment set", "[math::ft]" )
{
    int rv = setenv( "MXFFTW_WISDOM", "/tmp/.fftw", 1 );

    REQUIRE( rv == 0 );

    SECTION( "float" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<float>();
        REQUIRE( fname == "/tmp/.fftw/fftw_wisdom.float" );
        remove( fname.c_str() );
        mx::ioutils::createDirectories( "/tmp/.fftw" );

        {
            mx::math::ft::fftwEnvironment<float, true> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    SECTION( "double" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<double>();
        REQUIRE( fname == "/tmp/.fftw/fftw_wisdom.double" );
        remove( fname.c_str() );
        mx::ioutils::createDirectories( "/tmp/.fftw" );

        {
            mx::math::ft::fftwEnvironment<double, true> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    SECTION( "long double" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<long double>();
        REQUIRE( fname == "/tmp/.fftw/fftw_wisdom.long_double" );
        remove( fname.c_str() );
        mx::ioutils::createDirectories( "/tmp/.fftw" );

        {
            mx::math::ft::fftwEnvironment<long double, true> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    // clang-format off
    #ifdef HASQUAD // clang-format on

    SECTION( "quad" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<__float128>();
        REQUIRE( fname == "/tmp/.fftw/fftw_wisdom.quad" );
        remove( fname.c_str() );
        mx::ioutils::createDirectories( "/tmp/.fftw" );

        {
            mx::math::ft::fftwEnvironment<__float128, true> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    // clang-format off
    #endif // clang-format on
}

/// Create and destroy fftwEnvironment, with threads, environment not set
/**
 * \ingroup math_ft_fftwEnvironment_test
 */
TEST_CASE( "Create and destroy fftwEnvironment, with threads, environment not set", "[math::ft]" )
{
    int rv = setenv( "MXFFTW_WISDOM", "", 1 );

    REQUIRE( rv == 0 );

    SECTION( "float" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<float>();
        REQUIRE( fname == "./.fftw_wisdom.float" );
        remove( fname.c_str() );

        {
            mx::math::ft::fftwEnvironment<float, true> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    SECTION( "double" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<double>();
        REQUIRE( fname == "./.fftw_wisdom.double" );
        remove( fname.c_str() );

        {
            mx::math::ft::fftwEnvironment<double, true> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    SECTION( "long double" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<long double>();
        REQUIRE( fname == "./.fftw_wisdom.long_double" );
        remove( fname.c_str() );

        {
            mx::math::ft::fftwEnvironment<long double, true> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    // clang-format off
    #ifdef HASQUAD // clang-format on

    SECTION( "quad" )
    {
        // first remove any existing file
        std::string fname = mx::math::ft::fftw_wisdom_filename<__float128>();
        REQUIRE( fname == "./.fftw_wisdom.quad" );
        remove( fname.c_str() );

        {
            mx::math::ft::fftwEnvironment<__float128, true> fftwEnv;
        }

        mx::error_t errc;
        bool ex = mx::ioutils::exists<mx::verbose::d>( fname, errc );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( ex == true );
    }

    // clang-format off
    #endif // clang-format on
}
