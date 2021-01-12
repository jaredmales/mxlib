/** \file timeUtils.hpp
  * \brief Utilities for working with time
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  * \ingroup utils_files
  *
  */

//***********************************************************************//
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef timeUtils_hpp
#define timeUtils_hpp

#include <time.h>
#include <sys/time.h>
#include <cmath>

#include <thread>
#include <chrono>

#include <iostream>

#include "../ioutils/stringUtils.hpp"
#include "../astro/sofa.hpp"

namespace mx
{
namespace sys 
{
 
/// Get the current system time in seconds.
/** Uses timespec, so nanosecond resolution is possible.
  *
  * \tparam typeT is the type to return the time in [default=double].
  *               must have cast from integer types, and + and / operators defined.
  * \tparam clk_id is the sys/time.h clock identifier [default=CLOCK_REALTIME]
  *
  * \retval typeT containing the current time in seconds
  *
  * \test Verify operation of get_curr_time. \ref tests_sys_timeUtils_get_curr_time "[test doc]" 
  * 
  * \ingroup timeutils
  */
template<typename typeT=double, clockid_t clk_id=CLOCK_REALTIME>
typeT get_curr_time( timespec & tsp /**<[out] a timespec to populate with the current time */)
{
   clock_gettime(clk_id, &tsp);

   return ((typeT)tsp.tv_sec) + ((typeT)tsp.tv_nsec)/1e9;
}

//Specialization for most common use case.
template<>
double get_curr_time<double, CLOCK_REALTIME>(timespec & tsp);


/// Get the current system time in seconds.
/** Uses timespec, so nanosecond resolution is possible.  This version creates a timespec internally.
  * 
  * \overload
  *
  * \tparam typeT is the type to return the time in [default=double].
  *               must have cast from integer types, and + and / operators defined.
  * \tparam clk_id is the sys/time.h clock identifier [default=CLOCK_REALTIME]
  *
  * \retval typeT containing the current time in seconds
  *
  * \test Verify operation of get_curr_time. \ref tests_sys_timeUtils_get_curr_time "[test doc]" 
  * 
  * \ingroup timeutils
  */
template<typename typeT=double, clockid_t clk_id=CLOCK_REALTIME>
typeT get_curr_time()
{
   struct timespec tsp;
   return get_curr_time<typeT, clk_id>(tsp);
}

//Specialization for most common use case.
template<>
double get_curr_time<double, CLOCK_REALTIME>();


/// Sleep for a specified period in seconds.
/** 
  * \test Verify operation of thread sleep functions. \ref tests_sys_timeUtils_sleep "[test doc]" 
  * 
  * \ingroup timeutils_sleep
  */
void sleep( unsigned sec /**< [in] the number of seconds to sleep. */);

/// Sleep for a specified period in milliseconds.
/** 
  * \test Verify operation of thread sleep functions. \ref tests_sys_timeUtils_sleep "[test doc]" 
  * 
  * \ingroup timeutils_sleep
  */
void milliSleep( unsigned msec /**< [in] the number of milliseconds to sleep. */);

/// Sleep for a specified period in microseconds.
/**
  * \test Verify operation of thread sleep functions. \ref tests_sys_timeUtils_sleep "[test doc]"
  *
  * \ingroup timeutils_sleep 
  */
void microSleep( unsigned usec /**< [in] the number of microseconds to sleep. */);

/// Sleep for a specified period in nanoseconds.
/** 
  * \test Verify operation of thread sleep functions. \ref tests_sys_timeUtils_sleep "[test doc]" 
  * 
  * \ingroup timeutils_sleep
  */
void nanoSleep( unsigned nsec /**< [in] the number of microseconds to sleep. */);

/// Adds a time offset to an existing timespec
/** The offset is specified in nanoseconds, which can be greater than 1e9.
  *
  * \test Verify operation of timespecAddNsec. \ref tests_sys_timeUtils_timespecAddNsec "[test doc]"  
  * 
  * \ingroup timeutils
  */
void timespecAddNsec( timespec & ts, ///< [in/out] the time to add to
                      unsigned nsec  ///< [in] the number of nanoseconds to add to ts.
                    );

/** Parse a string of format hh:mm:ss.s
  * Breaks a time string into constituent parts.  Handles -h by distributing the sign to m and s.
  *
  * \tparam floatT is a floating point type
  *
  * \test Verify parsing of a formatted time string. \ref tests_sys_timeUtils_parse_hms "[test doc]"
  * 
  * \ingroup timeutils
  */
template<typename floatT>
void parse_hms( floatT & h,                ///< [out] the hour component coverted to floatT
                floatT & m,                ///< [out] the minute component converted to floatT
                floatT & s,                ///< [out] the second component converted to floatT
                const std::string & hmsstr ///< [in] a string of format hh:mm:ss.s where ss.s can be of any precision
              )
{
   int st, en;

   int sgn = 1;

   st = 0;
   en = hmsstr.find(':', st);

   h = ioutils::convertFromString<floatT>(hmsstr.substr(st, en-st).c_str());

   //Check for negative
   if(std::signbit(h)) sgn = -1;

   st = en + 1;

   en = hmsstr.find(':', st);

   m = sgn*ioutils::convertFromString<floatT>(hmsstr.substr(st, en-st));

   st = en+1;

   s = sgn*ioutils::convertFromString<floatT>(hmsstr.substr(st, hmsstr.length()-st).c_str());
}

///Converts a Gregorian calendar date into modified Julian date (MJD).
/** Uses the SOFA function iauCal2jd.  This is not a template in floating point
  * because SOFA is always double precision.
  *
  * \retval double containing the MJD
  * \retval <0 on error (-1 = bad year, -2 = bad month, -3 = bad day)
  *
  * \test Verify calculation of MJD. \ref tests_sys_timeUtils_Cal2mjd "[test doc]"
  * 
  * \ingroup timeutils
  */
double Cal2mjd( int yr,    ///< [in] Gregorian calendar year
                int mon,   ///< [in] Gregorian calendar month
                int day,   ///< [in] Gregorian calendar day
                int hr,    ///< [in] Gregorian calendar hour
                int min,   ///< [in] Gregorian calendar minute
                double sec ///< [in] Gregorian calendar second
              );

/// Parse an ISO8601 date of the form "YYYY-MM-DDTHH:MM:SS.S" into the individual components.
/** Parsing is currently only for the exact form above, which is the form in a FITS file.
  * See https://en.wikipedia.org/?title=ISO_8601.
  *
  * \returns 0 on success
  * \returns -4 if fdate is not long enough
  * 
  * \test Verify parsing of an ISO 8601 time string \ref tests_sys_timeUtils_ISO8601dateBreakdown "[test doc]"
  * 
  * \ingroup timeutils
  */
int ISO8601dateBreakdown( int & yr,                 ///< [out] Gregorian calendar year
                          int & mon,                ///< [out] Gregorian calendar month
                          int & day,                ///< [out] Gregorian calendar day
                          int & hr,                 ///< [out] Gregorian calendar hour
                          int & min,                ///< [out] Gregorian calendar minute
                          double & sec,             ///< [out] Gregorian calendar second
                          const std::string & fdate ///< [in] is a standard ISO8601 date string
                        );

///Parse an ISO8601 date of the form "YYYY-MM-DDTHH:MM:SS.S" and return the modified Julian date (MJD)
/** Parsing is currently only for the exact form above, which is the form in a FITS file.
  * See https://en.wikipedia.org/?title=ISO_8601.
  * After parsing calls Cal2mjd.
  *
  * \test Verify conversion of an ISO 8601 time string to MJD. \ref tests_sys_timeUtils_ISO8601date2mjd "[test doc]"
  * 
  * \ingroup timeutils
  */
double ISO8601date2mjd(const std::string & fdate /**<[in] a standard ISO8601 date string*/);

/// Get a date-time string in ISO 8601 format
/** For recognized time types, returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS.SS, with optional timezone designation such as Z or +00:00.
  *
  * \tparam timeT is the time type
  *
  * \retval std::string containing the format date/time
  *
  * \ingroup timeutils
  */
template<typename timeT>
std::string ISO8601DateTimeStr( const timeT &timeIn, ///< [in] the input time
                                int timeZone = 0     ///< [in] [optional] specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.
                              );

/// Get a date-time string in ISO 8601 format for time_t
/** Returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS, with optional timezone designation such as Z or +00:00.
  *
  * \retval std::string containing the format date/time
  *
  * \ingroup timeutils
  */
template<> 
std::string ISO8601DateTimeStr<time_t>( const time_t &timeIn, ///< [in] the input time
                                        int timeZone          ///< [in] [optional] specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.
                                      );

/// Get a date-time string in ISO 8601 format for timespec
/** Returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS.SSSSSSSSS, with optional timezone designation such as Z or +00:00.
  *
  * \retval std::string containing the format date/time
  *
  * \ingroup timeutils
  */
template<>
std::string ISO8601DateTimeStr<timespec>( const timespec &timeIn, ///< [in] the input time
                                          int timeZone            ///< [in] [optional] specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.
                                        );


/// Get a date-time string in ISO 8601 format for the current UTC time
/** Returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS.SS, with optional timezone designation such as Z or +00:00.
  *
  * \overload
  *
  * \retval std::string containing the format date/time
  *
  * \ingroup timeutils
  */

std::string ISO8601DateTimeStr(int timeZone = 0 /**< [in] [optional] specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.*/);

/// Get a date-time string in ISO 8601 format for an MJD
/** Returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS.SSSSSSSSS, with optional timezone designation such as Z or +00:00.
  *
  * \retval std::string containing the format date/time
  *
  * \ingroup timeutils
  */
std::string ISO8601DateTimeStrMJD( const double &timeIn, ///< [in] the input time
                                   int timeZone = 0      ///< [in] specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.
                                 );
/// Get a timestamp string in the form YYYYMMDDHHMMSS.SSSSSSSSS
/** Assumes the input timespec is in UTC.
  *
  * \returns 0 on success
  * \returns -1 on error.
  * 
  * \ingroup timeutils
  */ 
int timeStamp( std::string & tstamp, ///< [out] the string to hold the formatted time 
               timespec & ts         ///< [in] the timespec from which to produce the timestamp string
             );

/// Convert a UTC timespec to TAI modified Julian date
/** Converts a timespec assumed to be in <a href="https://en.wikipedia.org/wiki/Coordinated_Universal_Time">Coordinated
  *  Universal Time (UTC)</a> to a <a href="https://en.wikipedia.org/wiki/Julian_day">
  *  Modified Julian Date (MJD)</a> in
  * <a href="https://en.wikipedia.org/wiki/International_Atomic_Time">International Atomic Time (TAI)</a>.
  *
  * \retval 1 SOFA dubious year [see SOFA documentation for iauDat]
  * \retval 0 success
  * \retval -1 SOFA bad year [see SOFA documentation for iauDat and iauCal2jd]
  * \retval -2 SOFA bad month [see SOFA documentation for iauDat and iauCal2jd]
  * \retval -3 SOFA bad day [see SOFA documentation for iauDat and iauCal2jd]
  * \retval -4 SOFA bad fractional day [see SOFA documentation for iauDat and iauCal2jd]
  * \retval -5 SOFA internal error [see SOFA documentation for iauDat and iauCal2jd]
  * \retval -10 gmtime_r returned error, check errno
  *
  * \ingroup timeutils
  */
int timespecUTC2TAIMJD( double & djm,         ///< [out] the modified Julian day number
                        double & djmf,        ///< [out] the fraction of the day
                        const timespec & tsp, ///< [in] contains the UTC time
                        tm * tm0              ///< [out] [optional] will be filled with the broken down UTC time
                      );

   
/// Calculate the mean time of two times given by timespecs
/**
  * \returns the mean value of the inputs.
  */ 
timespec meanTimespec( timespec ts1, ///< [in] one of the times to average
                       timespec ts2  ///< [in] the other time to average
                     );

namespace tscomp
{
   
/// Timespec comparison operator \< (see caveats)
/** Caveats:
  * - If the inputs are in UTC (or similar scale) this does not account for leap seconds
  * - Assumes that the `tv_nsec` field does not exceed 999999999 nanoseconds
  * 
  * \returns true if tsL is earlier than tsR
  * \returns false otherwise
  * 
  * \ingroup timeutils_tscomp
  */  
bool operator<( timespec const& tsL, ///< [in] the left hand side of the comparison
                timespec const& tsR  ///< [in] the right hand side of the comparison 
              );

/// Timespec comparison operator \> (see caveats)
/** Caveats:
  * - If the inputs are in UTC (or similar scale) this does not account for leap seconds
  * - Assumes that the `tv_nsec` field does not exceed 999999999 nanoseconds
  * 
  * \returns true if tsL is later than tsR
  * \returns false otherwise
  * 
  * \ingroup timeutils_tscomp
  */
bool operator>( timespec const& tsL, ///< [in] the left hand side of the comparison
                timespec const& tsR  ///< [in] the right hand side of the comparison 
              );

/// Timespec comparison operator == (see caveats)
/** Caveats:
  * - If the inputs are in UTC (or similar scale) this does not account for leap seconds
  * - Assumes that the `tv_nsec` field does not exceed 999999999 nanoseconds
  * 
  * \returns true if tsL is exactly the same as tsR
  * \returns false otherwise
  * 
  * \ingroup timeutils_tscomp
  */
bool operator==( timespec const& tsL, ///< [in] the left hand side of the comparison
                 timespec const& tsR  ///< [in] the right hand side of the comparison 
               );

/// Timespec comparison operator \<= (see caveats)
/** Caveats:
  * - If the inputs are in UTC (or similar scale) this does not account for leap seconds
  * - Assumes that the `tv_nsec` field does not exceed 999999999 nanoseconds.  
  * 
  * \returns true if tsL is earlier than or exactly equal to tsR
  * \returns false otherwise
  * 
  * \ingroup timeutils_tscomp
  */
bool operator<=( timespec const& tsL, ///< [in] the left hand side of the comparison
                 timespec const& tsR  ///< [in] the right hand side of the comparison 
               );

/// Timespec comparison operator \>= (see caveats)
/** Caveats:
  * - If the inputs are in UTC (or similar scale) this does not account for leap seconds
  * - Assumes that the `tv_nsec` field does not exceed 999999999 nanoseconds
  * 
  * \returns true if tsL is exactly equal to or is later than tsR
  * \returns false otherwise
  * 
  * \ingroup timeutils_tscomp
  */
bool operator>=( timespec const& tsL, ///< [in] the left hand side of the comparison
                 timespec const& tsR  ///< [in] the right hand side of the comparison 
               );

} //namespace tscomp

namespace tsop
{
   
/// Add an amount of time specified in seconds to a timespec
/**
  * \returns the resulting timespec after addition
  */  
template<typename arithT>
timespec operator+( timespec ts, ///< [in] the timespec to add to
                    arithT add   ///< [in] the seconds to add
                  )
{
   ts.tv_sec += floor(add);
   
   ts.tv_nsec += (add - floor(add))*1e9;
   
   if(ts.tv_nsec >= 1e9)
   {
      ts.tv_sec += 1;
      ts.tv_nsec -= 1e9;
   }
   
   return ts;
}

/// Subtract an amount of time specified in seconds from a timespec
/**
  * \returns the resulting timespec after subtraction
  */
template<typename arithT>
timespec operator-( timespec ts, ///< [in] the timespec to subtract from
                    arithT sub   ///< [in] the seconds to subtract
                  )
{
   ts.tv_sec -= floor(sub);
   
   ts.tv_nsec -= (sub - floor(sub))*1e9;
   
   if(ts.tv_nsec < 0)
   {
      ts.tv_sec -= 1;
      ts.tv_nsec += 1e9;
   }
   
   return ts;
}

} //namespace tsop

} //namespace sys
} //namespace mx

#endif //timeUtils_hpp
