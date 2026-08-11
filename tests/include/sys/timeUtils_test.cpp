/** \file timeUtils_test.cpp
 * \brief Tests time and date utilities.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/sys/timeUtils.hpp"

/** \cond */
namespace
{

using mx::sys::timeUtilsDetail::operations;

int controlledDatResult = 0;
int controlledCal2jdResult = 0;

tm *failGmtimeR( const time_t *, tm * )
{
    return nullptr;
}

int controlledDat( int, int, int, double, double *dat )
{
    *dat = 37.0;
    return controlledDatResult;
}

int controlledCal2jd( int, int, int, double *djm0, double *djm )
{
    *djm0 = DJM0;
    *djm = 57754.0;
    return controlledCal2jdResult;
}

class operationsGuard
{
  public:
    operationsGuard() : m_saved( mx::sys::timeUtilsDetail::operationsInstance() )
    {
        controlledDatResult = 0;
        controlledCal2jdResult = 0;
    }

    ~operationsGuard()
    {
        mx::sys::timeUtilsDetail::operationsInstance() = m_saved;
    }

  private:
    operations m_saved;
};

} // namespace
/** \endcond */

/** Verify operation of get_curr_time.
 *
 * This only checks for increasing time on subsequent calls.
 *
 */
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify operation of get_curr_time", "[timeutils]" )
{
    GIVEN( "getting the time" )
    {
        WHEN( "the same timespec is used" )
        {
            timespec ts;
            double t0 = mx::sys::get_curr_time( ts );
            REQUIRE( t0 > 0 );
            std::cout << std::endl; // Do something to make sure some time passes.
            double t1 = mx::sys::get_curr_time( ts );

            REQUIRE( t1 > t0 );
        }

        WHEN( "no timespec is provided" )
        {
            double t0 = mx::sys::get_curr_time();
            REQUIRE( t0 > 0 );
            std::cout << std::endl; // Do something to make sure some time passes.
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 > t0 );
        }
    }
}

/** \brief Verifies MJD-to-ISO8601 formatting used for HCI coadd header timestamps.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Formatting MJD timestamps as ISO8601", "[timeutils][hciReduce]" )
{
    REQUIRE( mx::sys::ISO8601DateTimeStrMJD( 51544.0, 0 ) == "2000-01-01T00:00:00.000000000" );
    REQUIRE( mx::sys::ISO8601DateTimeStrMJD( 51544.5, 1 ) == "2000-01-01T12:00:00.000000000Z" );
    REQUIRE( mx::sys::ISO8601DateTimeStrMJD( 51545.0, 2 ) == "2000-01-02T00:00:00.000000000+00:00" );
}

/** Verify operation of thread sleep functions
 *
 * Uses mx::sys::get_curr_time to verify duration of sleep.
 *
 */
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify operation of thread sleep functions", "[timeutils]" )
{
    GIVEN( "sleeping for 1 second" )
    {
        WHEN( "sleeping in seconds" )
        {
            double t0 = mx::sys::get_curr_time();
            mx::sys::sleep( 1 );
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 >= t0 + 1.0 );
        }
        WHEN( "sleeping in milliseconds" )
        {
            double t0 = mx::sys::get_curr_time();
            mx::sys::milliSleep( 1000 );
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 >= t0 + 1.0 );
        }
        WHEN( "sleeping in microseconds" )
        {
            double t0 = mx::sys::get_curr_time();
            mx::sys::microSleep( 1000000 );
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 >= t0 + 1.0 );
        }
        WHEN( "sleeping in nanoseconds" )
        {
            double t0 = mx::sys::get_curr_time();
            mx::sys::nanoSleep( 1000000000 );
            double t1 = mx::sys::get_curr_time();

            REQUIRE( t1 >= t0 + 1.0 );
        }
    }
}

/** Verify operation of timespecAddNsec
 *
 */
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify operation of timespecAddNsec", "[timeutils]" )
{
    GIVEN( "a timespec" )
    {
        WHEN( "adding less than 1e9 nanoseconds" )
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

        WHEN( "adding more than 1e9 nanoseconds but less than 2e9" )
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

        WHEN( "adding more than 2e9" )
        {
            timespec ts;
            ts.tv_sec = 1;
            ts.tv_nsec = 0;

            mx::sys::timespecAddNsec( ts, 2000000010 );
            REQUIRE( ts.tv_sec == 3 );
            REQUIRE( ts.tv_nsec == 10 );
        }

        WHEN( "the existing nanoseconds cause a carry" )
        {
            timespec ts{ 1, 900000000 };

            mx::sys::timespecAddNsec( ts, 200000000 );

            REQUIRE( ts.tv_sec == 2 );
            REQUIRE( ts.tv_nsec == 100000000 );
        }
    }
}

/** Verify parsing of a formatted time string
 *
 *  Tests parsing of a string of format hh:mm:ss.s
 *
 */
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify parsing of a formatted time string", "[timeutils]" )
{
    GIVEN( "a valid time string" )
    {
        WHEN( "integer seconds" )
        {
            float hr;
            float mn;
            float sec;

            mx::sys::parse_hms( hr, mn, sec, "1:2:3" );

            REQUIRE( hr == 1 );
            REQUIRE( mn == 2 );
            REQUIRE( sec == 3 );
        }

        WHEN( "floating seconds" )
        {
            float hr;
            float mn;
            float sec;

            mx::sys::parse_hms( hr, mn, sec, "1:2:3.23" );

            REQUIRE( hr == 1 );
            REQUIRE( mn == 2 );
            REQUIRE_THAT( sec, Catch::Matchers::WithinAbs( 3.23, 1e-7 ) );
        }

        WHEN( "negative hour" )
        {
            float hr;
            float mn;
            float sec;

            mx::sys::parse_hms( hr, mn, sec, "-1:2:3.23" );

            REQUIRE( hr == -1 );
            REQUIRE( mn == -2 );
            REQUIRE_THAT( sec, Catch::Matchers::WithinAbs( -3.23, 1e-7 ) );
        }

        WHEN( "0 pads" )
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
}

/** Verify calculation of MJD
 *
 */
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify calculation of MJD", "[timeutils]" )
{
    GIVEN( "a valid Gregorian date" )
    {
        WHEN( "integer seconds" )
        {
            double mjd = mx::sys::Cal2mjd( 2020, 12, 31, 0, 0, 0 );
            REQUIRE( mjd == 59214.0 );
        }
    }

    GIVEN( "a valid Gregorian date" )
    {
        WHEN( "floating seconds" )
        {
            double mjd = mx::sys::Cal2mjd( 2020, 12, 31, 0, 0, 10.2357 );
            REQUIRE_THAT( mjd, Catch::Matchers::WithinAbs( 59214.00011846875, 1e-14 ) );
        }
    }
}

/** Verify parsing of an ISO 8601 time string
 *
 */
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify parsing of an ISO 8601 time string", "[timeutils]" )
{
    GIVEN( "a valid ISO 8601 time string" )
    {
        WHEN( "integer seconds" )
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

        WHEN( "fractional seconds" )
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

    GIVEN( "an invalid ISO 8601 time string" )
    {
        WHEN( "string too short" )
        {
            int yr, mon, day, hr, min;
            double sec;

            int rv = mx::sys::ISO8601dateBreakdown( yr, mon, day, hr, min, sec, "2020-12-31" );

            REQUIRE( rv == -4 );
        }
    }
}

/** Verify conversion of an ISO 8601 time string to MJD
 *
 */
/**
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Verify conversion of an ISO 8601 time string to MJD", "[timeutils]" )
{
    GIVEN( "a valid ISO 8601 time string" )
    {
        WHEN( "integer seconds" )
        {
            double mjd = mx::sys::ISO8601date2mjd( "2020-12-31T00:00:00" );

            REQUIRE( mjd == 59214.0 );
        }

        WHEN( "fractional seconds" )
        {
            double mjd = mx::sys::ISO8601date2mjd( "2020-12-31T00:00:10.2357" );

            REQUIRE_THAT( mjd, Catch::Matchers::WithinAbs( 59214.00011846875, 1e-14 ) );
        }
    }

    GIVEN( "an invalid ISO 8601 time string" )
    {
        REQUIRE( mx::sys::ISO8601date2mjd( "2020-12-31" ) == -4 );
    }
}

/** \brief Verifies calendar conversion errors are returned to callers.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Reject invalid Gregorian calendar dates", "[timeutils]" )
{
    REQUIRE( mx::sys::Cal2mjd( -5000, 1, 1, 0, 0, 0 ) == -1 );
    REQUIRE( mx::sys::Cal2mjd( 2020, 13, 1, 0, 0, 0 ) == -2 );
    REQUIRE( mx::sys::Cal2mjd( 2020, 2, 31, 0, 0, 0 ) == -3 );
}

/** \brief Verifies ISO8601 formatting for epoch, fractional, and current timestamps.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Format UTC timestamps as ISO8601", "[timeutils]" )
{
    const time_t epoch2000 = 946684800;
    const timespec fractional{ epoch2000, 123456789 };

    REQUIRE( mx::sys::ISO8601DateTimeStr( epoch2000, 0 ) == "2000-01-01T00:00:00" );
    REQUIRE( mx::sys::ISO8601DateTimeStr( epoch2000, 1 ) == "2000-01-01T00:00:00Z" );
    REQUIRE( mx::sys::ISO8601DateTimeStr( epoch2000, 2 ) == "2000-01-01T00:00:00+00:00" );
    REQUIRE( mx::sys::ISO8601DateTimeStr( epoch2000, 3 ) == "2000-01-01T00:00:00" );

    REQUIRE( mx::sys::ISO8601DateTimeStr( fractional, 0 ) == "2000-01-01T00:00:00.123456789" );
    REQUIRE( mx::sys::ISO8601DateTimeStr( fractional, 1 ) == "2000-01-01T00:00:00.123456789Z" );
    REQUIRE( mx::sys::ISO8601DateTimeStr( fractional, 2 ) == "2000-01-01T00:00:00.123456789+00:00" );

    REQUIRE( mx::sys::ISO8601DateTimeStr( 0 ).size() == 19 );
    REQUIRE( mx::sys::ISO8601DateTimeStr( 1 ).size() == 20 );
    REQUIRE( mx::sys::ISO8601DateTimeStr( 2 ).size() == 25 );
}

/** \brief Verifies ISO8601 formatting reports broken-down-time failures safely.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Handle ISO8601 broken-down-time failures", "[timeutils]" )
{
    operationsGuard guard;
    mx::sys::timeUtilsDetail::operationsInstance().gmtimeR = failGmtimeR;

    const time_t epoch = 0;
    const timespec fractional{ epoch, 1 };

    REQUIRE( mx::sys::ISO8601DateTimeStr( epoch ).empty() );
    REQUIRE( mx::sys::ISO8601DateTimeStr( fractional ).empty() );
}

/** \brief Verifies compact UTC timestamp formatting and its failure result.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Format compact UTC timestamps", "[timeutils]" )
{
    timespec ts{ 946684800, 123456789 };
    std::string stamp;

    REQUIRE( mx::sys::timeStamp( stamp, ts ) == 0 );
    REQUIRE( stamp == "20000101000000123456789" );

    operationsGuard guard;
    mx::sys::timeUtilsDetail::operationsInstance().gmtimeR = failGmtimeR;

    REQUIRE( mx::sys::timeStamp( stamp, ts ) == -1 );
}

/** \brief Verifies UTC timespec conversion to split TAI MJD values.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Convert UTC timespecs to TAI MJD", "[timeutils]" )
{
    double djm = 0;
    double djmf = 0;
    tm brokenDown{};

    const timespec epoch2000{ 946684800, 250000000 };
    REQUIRE( mx::sys::timespecUTC2TAIMJD( djm, djmf, epoch2000, &brokenDown ) == 0 );
    REQUIRE( djm == 51544.0 );
    REQUIRE_THAT( djmf, Catch::Matchers::WithinAbs( 32.25 / 86400.0, 1e-15 ) );
    REQUIRE( brokenDown.tm_year == 100 );
    REQUIRE( brokenDown.tm_mon == 0 );
    REQUIRE( brokenDown.tm_mday == 1 );

    const timespec rollover{ 1483315190, 500000000 };
    REQUIRE( mx::sys::timespecUTC2TAIMJD( djm, djmf, rollover, nullptr ) == 0 );
    REQUIRE( djm == 57755.0 );
    REQUIRE_THAT( djmf, Catch::Matchers::WithinAbs( 27.5 / 86400.0, 1e-15 ) );
}

/** \brief Verifies UTC-to-TAI conversion propagates all external operation failures.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Report UTC-to-TAI conversion failures", "[timeutils]" )
{
    operationsGuard guard;
    const timespec ts{ 946684800, 0 };
    tm brokenDown{};
    double djm = 0;
    double djmf = 0;

    mx::sys::timeUtilsDetail::operationsInstance().gmtimeR = failGmtimeR;
    REQUIRE( mx::sys::timespecUTC2TAIMJD( djm, djmf, ts, &brokenDown ) == -10 );

    mx::sys::timeUtilsDetail::resetOperations();
    mx::sys::timeUtilsDetail::operationsInstance().iauDat = controlledDat;
    controlledDatResult = -3;
    REQUIRE( mx::sys::timespecUTC2TAIMJD( djm, djmf, ts, &brokenDown ) == -3 );

    controlledDatResult = 0;
    mx::sys::timeUtilsDetail::operationsInstance().iauCal2jd = controlledCal2jd;
    controlledCal2jdResult = -2;
    REQUIRE( mx::sys::timespecUTC2TAIMJD( djm, djmf, ts, &brokenDown ) == -2 );

    controlledCal2jdResult = 0;
    controlledDatResult = 1;
    REQUIRE( mx::sys::timespecUTC2TAIMJD( djm, djmf, ts, &brokenDown ) == 1 );
}

/** \brief Verifies averaging normalized timespec values.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Average timespec values", "[timeutils]" )
{
    timespec mean = mx::sys::meanTimespec( { 0, 0 }, { 2, 0 } );
    REQUIRE( mean.tv_sec == 1 );
    REQUIRE( mean.tv_nsec == 0 );

    mean = mx::sys::meanTimespec( { 0, 0 }, { 1, 0 } );
    REQUIRE( mean.tv_sec == 0 );
    REQUIRE( mean.tv_nsec == 500000000 );

    mean = mx::sys::meanTimespec( { 0, 900000000 }, { 1, 900000000 } );
    REQUIRE( mean.tv_sec == 1 );
    REQUIRE( mean.tv_nsec == 400000000 );
}

/** \brief Verifies all timespec comparison relations by seconds and nanoseconds.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Compare timespec values", "[timeutils]" )
{
    const timespec early{ 1, 100 };
    const timespec laterNsec{ 1, 200 };
    const timespec laterSec{ 2, 0 };
    const timespec same{ 1, 100 };

    REQUIRE( mx::sys::tscomp::operator<( early, laterNsec ) );
    REQUIRE( mx::sys::tscomp::operator<( early, laterSec ) );
    REQUIRE_FALSE( mx::sys::tscomp::operator<( laterSec, early ) );
    REQUIRE( mx::sys::tscomp::operator>( laterNsec, early ) );
    REQUIRE( mx::sys::tscomp::operator>( laterSec, early ) );
    REQUIRE_FALSE( mx::sys::tscomp::operator>( early, laterSec ) );
    REQUIRE( mx::sys::tscomp::operator==( early, same ) );
    REQUIRE_FALSE( mx::sys::tscomp::operator==( early, laterNsec ) );
    REQUIRE( mx::sys::tscomp::operator<=( early, same ) );
    REQUIRE( mx::sys::tscomp::operator<=( early, laterNsec ) );
    REQUIRE_FALSE( mx::sys::tscomp::operator<=( laterNsec, early ) );
    REQUIRE( mx::sys::tscomp::operator>=( early, same ) );
    REQUIRE( mx::sys::tscomp::operator>=( laterNsec, early ) );
    REQUIRE_FALSE( mx::sys::tscomp::operator>=( early, laterNsec ) );
}

/** \brief Verifies arithmetic addition and subtraction normalize timespec values.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Apply arithmetic offsets to timespec values", "[timeutils]" )
{
    timespec ts{ 1, 100000000 };
    timespec result = mx::sys::tsop::operator+( ts, 0.25 );
    REQUIRE( result.tv_sec == 1 );
    REQUIRE( result.tv_nsec == 350000000 );

    ts = { 1, 900000000 };
    result = mx::sys::tsop::operator+( ts, 0.25 );
    REQUIRE( result.tv_sec == 2 );
    REQUIRE( result.tv_nsec == 150000000 );

    ts = { 1, 900000000 };
    result = mx::sys::tsop::operator-( ts, 0.25 );
    REQUIRE( result.tv_sec == 1 );
    REQUIRE( result.tv_nsec == 650000000 );

    ts = { 1, 100000000 };
    result = mx::sys::tsop::operator-( ts, 0.25 );
    REQUIRE( result.tv_sec == 0 );
    REQUIRE( result.tv_nsec == 850000000 );
}

/** \brief Verifies the generic current-time template with a monotonic clock.
 *
 * \ingroup timeUtils_unit_tests
 */
TEST_CASE( "Read a generic monotonic clock value", "[timeutils]" )
{
    timespec ts{};
    const long double value = mx::sys::get_curr_time<long double, CLOCK_MONOTONIC>( ts );

    REQUIRE( value > 0 );
    REQUIRE( ts.tv_sec > 0 );
}
