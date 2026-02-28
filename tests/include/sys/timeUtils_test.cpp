/** \file timeUtils_test.cpp
 */
#include "../../catch2/catch.hpp"

#include "../../../include/sys/timeUtils.hpp"/// Verify operation of get_curr_time


/// Verify operation of get_curr_time
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify operation of get_curr_time", "[timeutils]" )
{
    SECTION( "getting the time" )
    {
        SECTION( "the same timespec is used" )
        {
            timespec ts;
            double t0 = mx::sys::get_curr_time( ts );
            REQUIRE( t0 > 0 );
            std::cout << std::endl; // Do something to make sure some time passes.
            double t1 = mx::sys::get_curr_time( ts );

            REQUIRE( t1 > t0 );
        }

        SECTION( "no timespec is provided" )
        {
            double t0 = mx::sys::get_curr_time();
            REQUIRE( t0 > 0 );
            std::cout << std::endl; // Do something to make sure some time passes.
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 > t0 );
        }
    }
}/// Verify operation of thread sleep functions


/// Verify operation of thread sleep functions
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify operation of thread sleep functions", "[timeutils]" )
{
    SECTION( "sleeping for 1 second" )
    {
        SECTION( "sleeping in seconds" )
        {
            double t0 = mx::sys::get_curr_time();
            mx::sys::sleep( 1 );
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 >= t0 + 1.0 );
        }
        SECTION( "sleeping in milliseconds" )
        {
            double t0 = mx::sys::get_curr_time();
            mx::sys::milliSleep( 1000 );
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 >= t0 + 1.0 );
        }
        SECTION( "sleeping in microseconds" )
        {
            double t0 = mx::sys::get_curr_time();
            mx::sys::microSleep( 1000000 );
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 >= t0 + 1.0 );
        }
        SECTION( "sleeping in nanoseconds" )
        {
            double t0 = mx::sys::get_curr_time();
            mx::sys::nanoSleep( 1000000000 );
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 >= t0 + 1.0 );
        }
    }
}/// Verify operation of timespecAddNsec


/// Verify operation of timespecAddNsec
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify operation of timespecAddNsec", "[timeutils]" )
{
    SECTION( "a timespec" )
    {
        SECTION( "adding less than 1e9 nanoseconds" )
        {
            timespec ts;
            ts.tv_sec = 1;
            ts.tv_nsec = 0;

            mx::sys::timespecAddNsec( ts, 10 );
            REQUIRE( ts.tv_sec == 1 );
            REQUIRE( ts.tv_nsec == 10 );

            mx::sys::timespecAddNsec( ts, 100 );
            REQUIRE( ts.tv_sec == 1 );
            REQUIRE( ts.tv_nsec == 110 );

            mx::sys::timespecAddNsec( ts, 1000 );
            REQUIRE( ts.tv_sec == 1 );
            REQUIRE( ts.tv_nsec == 1110 );

            mx::sys::timespecAddNsec( ts, 10000 );
            REQUIRE( ts.tv_sec == 1 );
            REQUIRE( ts.tv_nsec == 11110 );

            mx::sys::timespecAddNsec( ts, 100000 );
            REQUIRE( ts.tv_sec == 1 );
            REQUIRE( ts.tv_nsec == 111110 );

            mx::sys::timespecAddNsec( ts, 1000000 );
            REQUIRE( ts.tv_sec == 1 );
            REQUIRE( ts.tv_nsec == 1111110 );

            mx::sys::timespecAddNsec( ts, 10000000 );
            REQUIRE( ts.tv_sec == 1 );
            REQUIRE( ts.tv_nsec == 11111110 );

            mx::sys::timespecAddNsec( ts, 100000000 );
            REQUIRE( ts.tv_sec == 1 );
            REQUIRE( ts.tv_nsec == 111111110 );
        }

        SECTION( "adding more than 1e9 nanoseconds but less than 2e9" )
        {
            timespec ts;
            ts.tv_sec = 1;
            ts.tv_nsec = 0;

            mx::sys::timespecAddNsec( ts, 1000000000 );
            REQUIRE( ts.tv_sec == 2 );
            REQUIRE( ts.tv_nsec == 0 );

            mx::sys::timespecAddNsec( ts, 1000000000 + 10 );
            REQUIRE( ts.tv_sec == 3 );
            REQUIRE( ts.tv_nsec == 10 );

            mx::sys::timespecAddNsec( ts, 1000000000 + 100 );
            REQUIRE( ts.tv_sec == 4 );
            REQUIRE( ts.tv_nsec == 110 );

            mx::sys::timespecAddNsec( ts, 1000000000 + 1000 );
            REQUIRE( ts.tv_sec == 5 );
            REQUIRE( ts.tv_nsec == 1110 );

            mx::sys::timespecAddNsec( ts, 1000000000 + 10000 );
            REQUIRE( ts.tv_sec == 6 );
            REQUIRE( ts.tv_nsec == 11110 );

            mx::sys::timespecAddNsec( ts, 1000000000 + 100000 );
            REQUIRE( ts.tv_sec == 7 );
            REQUIRE( ts.tv_nsec == 111110 );

            mx::sys::timespecAddNsec( ts, 1000000000 + 1000000 );
            REQUIRE( ts.tv_sec == 8 );
            REQUIRE( ts.tv_nsec == 1111110 );

            mx::sys::timespecAddNsec( ts, 1000000000 + 10000000 );
            REQUIRE( ts.tv_sec == 9 );
            REQUIRE( ts.tv_nsec == 11111110 );

            mx::sys::timespecAddNsec( ts, 1000000000 + 100000000 );
            REQUIRE( ts.tv_sec == 10 );
            REQUIRE( ts.tv_nsec == 111111110 );
        }

        SECTION( "adding more than 2e9" )
        {
            timespec ts;
            ts.tv_sec = 1;
            ts.tv_nsec = 0;

            mx::sys::timespecAddNsec( ts, 2000000010 );
            REQUIRE( ts.tv_sec == 3 );
            REQUIRE( ts.tv_nsec == 10 );
        }
    }
}/// Verify parsing of a formatted time string


/// Verify parsing of a formatted time string
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify parsing of a formatted time string", "[timeutils]" )
{
    SECTION( "a valid time string" )
    {
        SECTION( "integer seconds" )
        {
            float hr;
            float mn;
            float sec;

            mx::sys::parse_hms( hr, mn, sec, "1:2:3" );

            REQUIRE( hr == 1 );
            REQUIRE( mn == 2 );
            REQUIRE( sec == 3 );
        }

        SECTION( "floating seconds" )
        {
            float hr;
            float mn;
            float sec;

            mx::sys::parse_hms( hr, mn, sec, "1:2:3.23" );

            REQUIRE( hr == 1 );
            REQUIRE( mn == 2 );
            REQUIRE_THAT( sec, Catch::Matchers::WithinAbs( 3.23, 1e-7 ) );
        }

        SECTION( "negative hour" )
        {
            float hr;
            float mn;
            float sec;

            mx::sys::parse_hms( hr, mn, sec, "-1:2:3.23" );

            REQUIRE( hr == -1 );
            REQUIRE( mn == -2 );
            REQUIRE_THAT( sec, Catch::Matchers::WithinAbs( -3.23, 1e-7 ) );
        }

        SECTION( "0 pads" )
        {
            float hr;
            float mn;
            float sec;

            mx::sys::parse_hms( hr, mn, sec, "01:02:03.23" );

            REQUIRE( hr == 1 );
            REQUIRE( mn == 2 );
            REQUIRE_THAT( sec, Catch::Matchers::WithinAbs( 3.23, 1e-7 ) );
        }
    }
}/// Verify calculation of MJD


/// Verify calculation of MJD
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify calculation of MJD", "[timeutils]" )
{
    SECTION( "a valid Gregorian date" )
    {
        SECTION( "integer seconds" )
        {
            double mjd = mx::sys::Cal2mjd( 2020, 12, 31, 0, 0, 0 );
            REQUIRE( mjd == 59214.0 );
        }
    }

    SECTION( "a valid Gregorian date" )
    {
        SECTION( "floating seconds" )
        {
            double mjd = mx::sys::Cal2mjd( 2020, 12, 31, 0, 0, 10.2357 );
            REQUIRE_THAT( mjd, Catch::Matchers::WithinAbs( 59214.00011846875, 1e-14 ) );
        }
    }
}/// Verify parsing of an ISO 8601 time string


/// Verify parsing of an ISO 8601 time string
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify parsing of an ISO 8601 time string", "[timeutils]" )
{
    SECTION( "a valid ISO 8601 time string" )
    {
        SECTION( "integer seconds" )
        {
            int yr, mon, day, hr, min;
            double sec;

            int rv = mx::sys::ISO8601dateBreakdown( yr, mon, day, hr, min, sec, "2020-12-31T00:00:00" );

            REQUIRE( rv == 0 );
            REQUIRE( yr == 2020 );
            REQUIRE( mon == 12 );
            REQUIRE( day == 31 );
            REQUIRE( hr == 0 );
            REQUIRE( min == 0 );
            REQUIRE( sec == 0 );
        }

        SECTION( "fractional seconds" )
        {
            int yr, mon, day, hr, min;
            double sec;

            int rv = mx::sys::ISO8601dateBreakdown( yr, mon, day, hr, min, sec, "2020-12-31T00:00:10.2357" );

            REQUIRE( rv == 0 );
            REQUIRE( yr == 2020 );
            REQUIRE( mon == 12 );
            REQUIRE( day == 31 );
            REQUIRE( hr == 0 );
            REQUIRE( min == 0 );
            REQUIRE_THAT( sec, Catch::Matchers::WithinAbs( 10.2357, 1e-14 ) );
        }
    }

    SECTION( "an invalid ISO 8601 time string" )
    {
        SECTION( "string too short" )
        {
            int yr, mon, day, hr, min;
            double sec;

            int rv = mx::sys::ISO8601dateBreakdown( yr, mon, day, hr, min, sec, "2020-12-31" );

            REQUIRE( rv == -4 );
        }
    }
}/// Verify conversion of an ISO 8601 time string to MJD


/// Verify conversion of an ISO 8601 time string to MJD
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify conversion of an ISO 8601 time string to MJD", "[timeutils]" )
{
    SECTION( "a valid ISO 8601 time string" )
    {
        SECTION( "integer seconds" )
        {
            double mjd = mx::sys::ISO8601date2mjd( "2020-12-31T00:00:00" );

            REQUIRE( mjd == 59214.0 );
        }

        SECTION( "fractional seconds" )
        {
            double mjd = mx::sys::ISO8601date2mjd( "2020-12-31T00:00:10.2357" );

            REQUIRE_THAT( mjd, Catch::Matchers::WithinAbs( 59214.00011846875, 1e-14 ) );
        }
    }
}
