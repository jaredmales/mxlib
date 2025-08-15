/** \file stringUtils_test.cpp
 */
#include "../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../include/ioutils/stringUtils.hpp"

#include <cmath>

using namespace Catch::Matchers;

using namespace mx::ioutils;
using namespace mx;

/** \test Scenario: Converting strings to numbers
 *
 * \anchor tests_ioutils_stringUtils_stoT
 */
SCENARIO( "Converting strings to numbers", "[ioutils::stringUtils]" )
{
    GIVEN( "a string char" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            char val = stoT<char>( "5" );
            REQUIRE( static_cast<int>( val ) == 5 );
        }

        WHEN( "string valid, negative, w/ error check" )
        {
            mx::error_t errc;
            char val = stoT<char>( "-5", &errc );
            REQUIRE( static_cast<int>( val ) == -5 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            char val = stoT<char>( "128", &errc );
            REQUIRE( static_cast<int>( val ) == std::numeric_limits<char>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "negative overflow, w/ error check" )
        {
            mx::error_t errc;
            char val = stoT<char>( "-129", &errc );
            REQUIRE( static_cast<int>( val ) == std::numeric_limits<char>::lowest() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            char val = stoT<char>( "!", &errc );
            REQUIRE( static_cast<int>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string unsigned char" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            unsigned char val = stoT<unsigned char>( "200" );
            REQUIRE( static_cast<int>( val ) == 200 );
        }

        WHEN( "string valid, w/ error check" )
        {
            mx::error_t errc;
            unsigned char val = stoT<unsigned char>( "100", &errc );
            REQUIRE( static_cast<int>( val ) == 100 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            unsigned char val = stoT<unsigned char>( "256", &errc );
            REQUIRE( static_cast<int>( val ) == std::numeric_limits<unsigned char>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            unsigned char val = stoT<unsigned char>( "*", &errc );
            REQUIRE( static_cast<int>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string short" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            short val = stoT<short>( "500" );
            REQUIRE( static_cast<int>( val ) == 500 );
        }

        WHEN( "string valid, negative, w/ error check" )
        {
            mx::error_t errc;
            short val = stoT<short>( "-1024", &errc );
            REQUIRE( static_cast<int>( val ) == -1024 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            short val = stoT<short>( "47000", &errc );
            REQUIRE( static_cast<int>( val ) == std::numeric_limits<short>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "negative overflow, w/ error check" )
        {
            mx::error_t errc;
            short val = stoT<short>( "-37000", &errc );
            REQUIRE( static_cast<int>( val ) == std::numeric_limits<short>::lowest() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            short val = stoT<short>( "-", &errc );
            REQUIRE( static_cast<int>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string unsigned short" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            unsigned short val = stoT<unsigned short>( "20000" );
            REQUIRE( static_cast<int>( val ) == 20000 );
        }

        WHEN( "string valid, w/ error check" )
        {
            mx::error_t errc;
            unsigned short val = stoT<unsigned short>( "65000", &errc );
            REQUIRE( static_cast<int>( val ) == 65000 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            unsigned short val = stoT<unsigned short>( "70000", &errc );
            REQUIRE( static_cast<int>( val ) == std::numeric_limits<unsigned short>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            unsigned short val = stoT<unsigned short>( "#", &errc );
            REQUIRE( static_cast<int>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string int" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            int val = stoT<int>( "1000000" );
            REQUIRE( static_cast<int>( val ) == 1000000 );
        }

        WHEN( "string valid, negative, w/ error check" )
        {
            mx::error_t errc;
            int val = stoT<int>( "-2000000", &errc );
            REQUIRE( static_cast<int>( val ) == -2000000 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            int val = stoT<int>( "3000000000", &errc );
            REQUIRE( static_cast<int>( val ) == std::numeric_limits<int>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "negative overflow, w/ error check" )
        {
            mx::error_t errc;
            int val = stoT<int>( "-2147483650", &errc );
            REQUIRE( static_cast<int>( val ) == std::numeric_limits<int>::lowest() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            int val = stoT<int>( " w", &errc );
            REQUIRE( static_cast<int>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string unsigned int" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            unsigned int val = stoT<unsigned int>( "4000000000" );
            REQUIRE( static_cast<unsigned int>( val ) == 4000000000 );
        }

        WHEN( "string valid, w/ error check" )
        {
            mx::error_t errc;
            unsigned int val = stoT<unsigned int>( "2", &errc );
            REQUIRE( static_cast<unsigned int>( val ) == 2 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            unsigned int val = stoT<unsigned int>( "6000000000", &errc );
            REQUIRE( static_cast<unsigned int>( val ) == std::numeric_limits<unsigned int>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            unsigned int val = stoT<unsigned int>( "?8", &errc );
            REQUIRE( static_cast<unsigned int>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string long" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            long val = stoT<long>( "1000000" );
            REQUIRE( static_cast<long>( val ) == 1000000 );
        }

        WHEN( "string valid, negative, w/ error check" )
        {
            mx::error_t errc;
            long val = stoT<long>( "-2000000", &errc );
            REQUIRE( static_cast<long>( val ) == -2000000 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            long val = stoT<long>( "9223372036854775808", &errc );
            REQUIRE( static_cast<long>( val ) == std::numeric_limits<long>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "negative overflow, w/ error check" )
        {
            mx::error_t errc;
            long val = stoT<long>( "-9223372036854775809", &errc );
            REQUIRE( static_cast<long>( val ) == std::numeric_limits<long>::lowest() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            long val = stoT<long>( " w8", &errc );
            REQUIRE( static_cast<long>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string unsigned long" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            unsigned long val = stoT<unsigned long>( "400000000000" );
            REQUIRE( static_cast<unsigned long>( val ) == 400000000000 );
        }

        WHEN( "string valid, w/ error check" )
        {
            mx::error_t errc;
            unsigned long val = stoT<unsigned long>( "16", &errc );
            REQUIRE( static_cast<unsigned long>( val ) == 16 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            unsigned long val = stoT<unsigned long>( "18523372036854775808", &errc );
            REQUIRE( static_cast<unsigned long>( val ) == std::numeric_limits<unsigned long>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            unsigned long val = stoT<unsigned long>( "?8", &errc );
            REQUIRE( static_cast<unsigned long>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string long long" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            long long val = stoT<long long>( "1000052" );
            REQUIRE( static_cast<long long>( val ) == 1000052 );
        }

        WHEN( "string valid, negative, w/ error check" )
        {
            mx::error_t errc;
            long long val = stoT<long long>( "-2300000", &errc );
            REQUIRE( static_cast<long long>( val ) == -2300000 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            long long val = stoT<long long>( "9223372036854775808", &errc );
            REQUIRE( static_cast<long long>( val ) == std::numeric_limits<long long>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "negative overflow, w/ error check" )
        {
            mx::error_t errc;
            long long val = stoT<long long>( "-9223372036854775809", &errc );
            REQUIRE( static_cast<long long>( val ) == std::numeric_limits<long long>::lowest() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            long long val = stoT<long long>( "-..8", &errc );
            REQUIRE( static_cast<long long>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string unsigned long long" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            unsigned long long val = stoT<unsigned long long>( "400000000000" );
            REQUIRE( static_cast<unsigned long long>( val ) == 400000000000 );
        }

        WHEN( "string valid, w/ error check" )
        {
            mx::error_t errc;
            unsigned long long val = stoT<unsigned long long>( "16", &errc );
            REQUIRE( static_cast<unsigned long long>( val ) == 16 );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            unsigned long long val = stoT<unsigned long long>( "18523372036854775808", &errc );
            REQUIRE( static_cast<unsigned long long>( val ) == std::numeric_limits<unsigned long long>::max() );
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            unsigned long long val = stoT<unsigned long long>( "++9.2", &errc );
            REQUIRE( static_cast<unsigned long long>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string bool" )
    {
        WHEN( "string valid, true, no error check" )
        {
            bool val = stoT<bool>( "true" );
            REQUIRE( static_cast<bool>( val ) == true );
        }

        WHEN( "string valid, false, w/ error check" )
        {
            mx::error_t errc;
            bool val = stoT<bool>( "false", &errc );
            REQUIRE( static_cast<bool>( val ) == false );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "string valid, t, w/ error check" )
        {
            mx::error_t errc;
            bool val = stoT<bool>( "t", &errc );
            REQUIRE( static_cast<bool>( val ) == true );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "string valid, f, w/ error check" )
        {
            mx::error_t errc;
            bool val = stoT<bool>( "false", &errc );
            REQUIRE( static_cast<bool>( val ) == false );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "string valid, 1, w/ error check" )
        {
            mx::error_t errc;
            bool val = stoT<bool>( "1", &errc );
            REQUIRE( static_cast<bool>( val ) == true );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "string valid, 0, w/ error check" )
        {
            mx::error_t errc;
            bool val = stoT<bool>( "0", &errc );
            REQUIRE( static_cast<bool>( val ) == false );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "string valid, number > 1, w/ error check" )
        {
            mx::error_t errc;
            bool val = stoT<bool>( "237", &errc );
            REQUIRE( static_cast<bool>( val ) == true );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "string valid, -0, w/ error check" )
        {
            mx::error_t errc;
            bool val = stoT<bool>( "-0.2", &errc );
            REQUIRE( static_cast<bool>( val ) == false );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            bool val = stoT<bool>( "Xtrue", &errc );
            REQUIRE( static_cast<bool>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }
    }

    GIVEN( "a string float" )
    {
        WHEN( "string valid, positive, no error check" )
        {
            float val = stoT<float>( "2.567" );
            REQUIRE_THAT( val, WithinRel( static_cast<float>( 2.567 ), std::numeric_limits<float>::epsilon() ) );
        }

        WHEN( "string valid, negative, w/ error check" )
        {
            mx::error_t errc;
            float val = stoT<float>( "-2300000.897", &errc );
            REQUIRE_THAT( val, WithinRel( static_cast<float>( -2300000.897 ), std::numeric_limits<float>::epsilon() ) );
            REQUIRE( errc == mx::error_t::noerror );
        }

        WHEN( "positive overflow, w/ error check" )
        {
            mx::error_t errc;
            float val = stoT<float>( "1e55", &errc );
            // inf comparison won't work
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "negative overflow, w/ error check" )
        {
            mx::error_t errc;
            float val = stoT<float>( "-1e55", &errc );
            // inf comparison won't work
            REQUIRE( errc == mx::error_t::erange );
        }

        WHEN( "invalid string, w/ error check" )
        {
            mx::error_t errc;
            float val = stoT<float>( "r5", &errc );
            REQUIRE( static_cast<float>( val ) == 0 );
            REQUIRE( errc == mx::error_t::invalidarg );
        }

        GIVEN( "a string double" )
        {
            WHEN( "string valid, positive, no error check" )
            {
                double val = stoT<double>( "22.2567" );
                REQUIRE_THAT( val, WithinRel( static_cast<double>( 22.2567 ), std::numeric_limits<double>::epsilon() ) );
            }

            WHEN( "string valid, negative, w/ error check" )
            {
                mx::error_t errc;
                double val = stoT<double>( "-2300000.897987", &errc );
                REQUIRE_THAT(
                    val,
                    WithinRel( static_cast<double>( -2300000.897987 ), std::numeric_limits<double>::epsilon() ) );
                REQUIRE( errc == mx::error_t::noerror );
            }

            WHEN( "positive overflow, w/ error check" )
            {
                mx::error_t errc;
                double val = stoT<double>( "1e400", &errc );
                // inf comparison won't work
                REQUIRE( errc == mx::error_t::erange );
            }

            WHEN( "negative overflow, w/ error check" )
            {
                mx::error_t errc;
                double val = stoT<double>( "-1e400", &errc );
                // inf comparison won't work
                REQUIRE( errc == mx::error_t::erange );
            }

            WHEN( "invalid string, w/ error check" )
            {
                mx::error_t errc;
                double val = stoT<double>( "+d", &errc );
                REQUIRE( static_cast<double>( val ) == 0 );
                REQUIRE( errc == mx::error_t::invalidarg );
            }
        }

        GIVEN( "a string long double" )
        {
            //These long double values don't seem quite right, should be able to use epsilon.
            //new catch2 might work better
            WHEN( "string valid, positive, no error check" )
            {
                long double val = stoT<long double>( "22.2567" );
                REQUIRE(
                    fabs(val - static_cast<long double>( 22.2567 )) < fabs(1e-9* static_cast<long double>(22.2567 ))) ;
            }

            WHEN( "string valid, negative, w/ error check" )
            {
                mx::error_t errc;
                long double val = stoT<long double>( "-2300000.897987", &errc );
                REQUIRE(
                    fabs(val - static_cast<long double>( -2300000.897987 )) < fabs(1e-9* static_cast<long double>( -2300000.897987 ))) ;
                REQUIRE( errc == mx::error_t::noerror );
            }

            WHEN( "positive overflow, w/ error check" )
            {
                mx::error_t errc;
                long double val = stoT<long double>( "1e10000", &errc );
                // inf comparison won't work
                REQUIRE( errc == mx::error_t::erange );
            }

            WHEN( "negative overflow, w/ error check" )
            {
                mx::error_t errc;
                long double val = stoT<long double>( "-1e10000", &errc );
                // inf comparison won't work
                REQUIRE( errc == mx::error_t::erange );
            }

            WHEN( "invalid string, w/ error check" )
            {
                mx::error_t errc;
                long double val = stoT<long double>( "+d", &errc );
                REQUIRE( static_cast<long double>( val ) == 0 );
                REQUIRE( errc == mx::error_t::invalidarg );
            }
        }
    }
}
