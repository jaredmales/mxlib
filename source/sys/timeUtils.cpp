/** \file timeUtils.cpp
 * \brief Definitions for utilities for working with time
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 * \ingroup utils_files
 *
 */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#include "sys/timeUtils.hpp"

namespace mx
{
namespace sys
{

template <>
double get_curr_time<double, CLOCK_REALTIME>( timespec &tsp )
{
    clock_gettime( CLOCK_REALTIME, &tsp );

    return ( (double)tsp.tv_sec ) + ( (double)tsp.tv_nsec ) / 1e9;
}

template <>
double get_curr_time<double, CLOCK_REALTIME>()
{
    struct timespec tsp;
    return get_curr_time<double, CLOCK_REALTIME>( tsp );
}

void sleep( unsigned sec )
{
    std::this_thread::sleep_for( std::chrono::seconds( sec ) );
}

void milliSleep( unsigned msec )
{
    std::this_thread::sleep_for( std::chrono::milliseconds( msec ) );
}

void microSleep( unsigned usec )
{
    std::this_thread::sleep_for( std::chrono::microseconds( usec ) );
}

void nanoSleep( unsigned nsec )
{
    std::this_thread::sleep_for( std::chrono::nanoseconds( nsec ) );
}

void timespecAddNsec( timespec &ts, unsigned nsec )
{
    ts.tv_nsec += nsec % 1000000000;
    ts.tv_sec += nsec / 1000000000;

    if( ts.tv_nsec > 999999999 )
    {
        ts.tv_nsec -= 1000000000;
        ts.tv_sec += 1;
    }
}

double Cal2mjd( int yr, int mon, int day, int hr, int min, double sec )
{
    double djm0;
    double djm;

    int rv = sofa::iauCal2jd( yr, mon, day, &djm0, &djm );

    if( rv < 0 )
        return (double)rv;

    djm0 = djm + ( ( (double)hr ) / ( 24.0 ) + ( (double)min ) / ( 60.0 * 24.0 ) + sec / ( 3600.0 * 24.0 ) );

    return djm0;
}

int ISO8601dateBreakdown( int &yr, int &mon, int &day, int &hr, int &min, double &sec, const std::string &fdate )
{
    if( fdate.length() < 19 )
        return -4;

    yr = atoi( fdate.substr( 0, 4 ).c_str() );
    mon = atoi( fdate.substr( 5, 2 ).c_str() );
    day = atoi( fdate.substr( 8, 2 ).c_str() );

    double _hr, _min;
    parse_hms( _hr, _min, sec, fdate.substr( 11 ) );

    hr = floor( _hr );
    min = floor( _min );

    return 0;
}

double ISO8601date2mjd( const std::string &fdate )
{
    if( fdate.length() < 19 )
        return -4;

    int yr, mon, day, hr, min;
    // double hr, min, sec;
    double sec;

    ISO8601dateBreakdown( yr, mon, day, hr, min, sec, fdate );

    return Cal2mjd( yr, mon, day, hr, min, sec );
}

template <>
std::string ISO8601DateTimeStr<time_t>( const time_t &timeIn, int timeZone )
{
    tm bdt;
    gmtime_r( &timeIn, &bdt );

    char tstr[25];

    strftime( tstr, 25, "%FT%H:%M:%S", &bdt );

    std::string result = tstr;

    if( timeZone == 1 )
        result += "Z";
    if( timeZone == 2 )
        result += "+00:00";

    return result;
}

template <>
std::string ISO8601DateTimeStr<timespec>( const timespec &timeIn, int timeZone )
{
    std::string result = ISO8601DateTimeStr<time_t>( timeIn.tv_sec, 0 );

    char tstr[20];

    snprintf( tstr, 20, ".%09ld", timeIn.tv_nsec );

    result += tstr;

    if( timeZone == 1 )
        result += "Z";
    if( timeZone == 2 )
        result += "+00:00";

    return result;
}

std::string ISO8601DateTimeStr( int timeZone )
{
    return ISO8601DateTimeStr<time_t>( ::time( 0 ), timeZone );
}

std::string ISO8601DateTimeStrMJD( const double &timeIn, int timeZone )
{
    int iy, im, id;
    double fd;

    sofa::iauJd2cal( DJM0, timeIn, &iy, &im, &id, &fd );

    int hr, mn;

    hr = floor( fd * 24.0 );
    fd = ( fd - hr / 24.0 ) * 24.0;

    mn = floor( fd * 60. );

    fd = ( fd - mn / 60.0 ) * 3600.0;

    char tstr[32];

    snprintf( tstr, 32, "%04d-%02d-%02dT%02d:%02d:%012.9f", iy, im, id, hr, mn, fd );

    std::string result = tstr;

    if( timeZone == 1 )
        result += "Z";
    if( timeZone == 2 )
        result += "+00:00";

    return result;
}

int timeStamp( std::string &tstamp, timespec &ts )
{
    tm uttime; // The broken down time.

    time_t t0 = ts.tv_sec;

    if( gmtime_r( &t0, &uttime ) == 0 )
    {
        std::cerr << "Error getting UT time (gmtime_r returned 0). At: " << __FILE__ << " " << __LINE__ << "\n";
        return -1;
    }

    char buffer[48];

    snprintf( buffer,
              sizeof( buffer ),
              "%04i%02i%02i%02i%02i%02i%09i",
              uttime.tm_year + 1900,
              uttime.tm_mon + 1,
              uttime.tm_mday,
              uttime.tm_hour,
              uttime.tm_min,
              uttime.tm_sec,
              static_cast<int>( ts.tv_nsec ) ); // casting in case we switch type of time_ns.

    tstamp = buffer;

    return 0;
}

int timespecUTC2TAIMJD( double &djm, double &djmf, const timespec &tsp, tm *tm0 )
{
    double dat, djm0;
    tm *tmrv;
    int rv1, rv2;

    // Get the broken down time corresponding to tsp0
    tmrv = gmtime_r( &tsp.tv_sec, tm0 );
    if( tmrv == 0 )
        return -10;

    // Then determine deltaAT = TAI-UTC
    rv1 = sofa::iauDat( 1900 + tm0->tm_year, 1 + tm0->tm_mon, tm0->tm_mday, 0.0, &dat );
    if( rv1 < 0 )
        return rv1;

    // And get the MJD
    rv2 = sofa::iauCal2jd( 1900 + tm0->tm_year, 1 + tm0->tm_mon, tm0->tm_mday, &djm0, &djm );
    if( rv2 < 0 )
        return rv2;

    // Finally calculate the day fraction
    djmf = ( (double)tm0->tm_hour ) / 24.0 + ( (double)tm0->tm_min ) / ( 24.0 * 60. ) +
           ( ( (double)tm0->tm_sec ) + ( (double)tsp.tv_nsec / 1e9 ) + dat ) / ( 24.0 * 3600.0 );

    if( djmf >= 1.0 )
    {
        djmf -= 1.0;
        djm += 1.0;
    }

    if( rv1 )
        return rv1;
    return 0;
}

timespec meanTimespec( timespec ts1, timespec ts2 )
{
    double means = ( ts1.tv_sec + ts2.tv_sec ) / 2.0;
    double meanns = ( ts1.tv_nsec + ts2.tv_nsec ) / 2.0;

    ts1.tv_sec = floor( means );
    ts1.tv_nsec = round( meanns );

    if( means != floor( means ) )
    {
        ts1.tv_nsec += 5e8;

        if( ts1.tv_nsec >= 1e9 )
        {
            ts1.tv_sec += 1;
            ts1.tv_nsec -= 1e9;
        }
    }

    return ts1;
}

namespace tscomp
{

bool operator<( timespec const &tsL, timespec const &tsR )
{
    return ( ( ( tsL.tv_sec == tsR.tv_sec ) && ( tsL.tv_nsec < tsR.tv_nsec ) ) || ( tsL.tv_sec < tsR.tv_sec ) );
}

bool operator>( timespec const &tsL, timespec const &tsR )
{
    return ( ( ( tsL.tv_sec == tsR.tv_sec ) && ( tsL.tv_nsec > tsR.tv_nsec ) ) || ( tsL.tv_sec > tsR.tv_sec ) );
}

bool operator==( timespec const &tsL, timespec const &tsR )
{
    return ( ( tsL.tv_sec == tsR.tv_sec ) && ( tsL.tv_nsec == tsR.tv_nsec ) );
}

bool operator<=( timespec const &tsL, timespec const &tsR )
{
    return ( tsL < tsR || tsL == tsR );
}

bool operator>=( timespec const &tsL, timespec const &tsR )
{
    return ( tsL > tsR || tsL == tsR );
}

} // namespace tscomp

} // namespace sys
} // namespace mx
