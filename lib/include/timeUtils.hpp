/** \file timeUtils.hpp
  * \brief Utilities for working with time
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup timeutils
  *
  */
  
#ifndef __timeUtils_hpp__
#define __timeUtils_hpp__


#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <sofa.h>

#include "stringUtils.hpp"


namespace mx
{
   
  
/** Get the current system time.
  * Uses timespec, so nanosecond resolution is possible.
  * 
  * \tparam typeT is the type to return the time in [default=double].  
  *               must have cast from integer types, and + and / operators defined.
  * \tparam clk_id is the sys/time.h clock identifier [default=CLOCK_REALTIME]
  * 
  * \retval typeT containing the current time in seconds
  *
  * \ingroup timeutils
  */
template<typename typeT=double, clockid_t clk_id=CLOCK_REALTIME>
typeT get_curr_time()
{
   struct timespec tsp;
   clock_gettime(clk_id, &tsp);
   
   return ((typeT)tsp.tv_sec) + ((typeT)tsp.tv_nsec)/1e9;
}



/** Parse a string of format hh:mm:ss.s
  * Breaks a time string into constituent parts.  Handles -h by distributing the sign to m and s.
  *
  * \param[in] hmsstr  is a string of format hh:mm:ss.s where ss.s can be of any precision
  * \param[out] h  is the hour component coverted to floatT
  * \param[out] m  is the minute component converted to floatT
  * \param[out] s  is the second component converted to floatT
  * 
  * \tparam floatT is a floating point type
  * 
  * \ingroup timeutils
  */ 
template<typename floatT>
void parse_hms(const std::string & hmsstr, floatT & h, floatT & m, floatT &s)
{
   int st, en;

   int sgn = 1;

   st = 0;
   en = hmsstr.find(':', st);
   
   //h = strtold(hmsstr.substr(st, en-st).c_str(), 0);
   h = convertFromString<floatT>(hmsstr.substr(st, en-st).c_str());
   
   //Check for negative
   if(std::signbit(h)) sgn = -1;

   st = en + 1;
   
   en = hmsstr.find(':', st);
   
   //m = sgn*strtold(hmsstr.substr(st, en-st).c_str(), 0);
   m = sgn*convertFromString<floatT>(hmsstr.substr(st, en-st));
   
   st = en+1;
   
   //s = sgn*strtold(hmsstr.substr(st, xmsstr.length()-st).c_str(), 0);
   s = sgn*convertFromString<floatT>(hmsstr.substr(st, hmsstr.length()-st).c_str());
}

///Converts a Gregorian calendar date into modified Julian date (MJD).
/** Uses the SOFA function iauCal2jd.  Note this is not a template in floating point
  * because SOFA is always double precision. 
  *
  * \param [in] yr Gregorian calendar year
  * \param [in] mon Gregorian calendar month
  * \param [in] day Gregorian calendar day
  * \param [in] hr Gregorian calendar hour
  * \param [in] min  Gregorian calendar minute
  * \param [in] sec  Gregorian calendar second
  *
  * \retval double containing the MJD
  * \retval <0 on error (-1 = bad year, -2 = bad month, -3 = bad day)
  * 
  * \ingroup timeutils
  * 
  */ 
inline 
double Cal2mjd(int yr, int mon, int day, int hr, int min, double sec)
{
   double djm0, djm;
   
   int rv = iauCal2jd(yr, mon, day, &djm0, &djm);
   
   if (rv < 0) return (double) rv;
   
   djm += ( (double) hr)/ (24.0);
   djm += ( (double) min)/ (60.0*24.0);
   djm +=  sec / (3600.0*24.0);
   
   return djm;
}
   
///Parse an ISO8601 date of the form "YYYY-MM-DDTHH:MM:SS.S" and return the modified Julian date (MJD)
/** Parsing is currently only for the exact form above, which is the form in a FITS file.
  * See https://en.wikipedia.org/?title=ISO_8601. 
  * After parsing calls Cal2mjd.
  * 
  * \param [in] fdate is a standard ISO8601 date string
  * 
  * \ingroup timeutils
  */
inline 
double ISO8601date2mjd(const std::string & fdate)
{   
   if(fdate.length() < 19) return -4;
   
   int yr, mon, day;
   double hr, min, sec;
   
   yr = atoi(fdate.substr(0,4).c_str());
   mon = atoi(fdate.substr(5,2).c_str());
   day = atoi(fdate.substr(8,2).c_str());
   
   parse_hms(fdate.substr(11), hr, min, sec);
   
   return Cal2mjd(yr,mon,day, floor(hr), floor(min),sec);
   
}

/// Get a date-time string in ISO 8601 format
/** For recognized time types, returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS.SS, with optional timezone designation such as Z or +00:00.
  *
  * \tparam timeT is the time type
  *
  * \param [in] timeIn is the input time
  * \param [in] timeZone specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.
  *
  * \retval std::string containing the format date/time 
  * 
  * \ingroup timeutils
  */ 
template<typename timeT>
std::string ISO8601DateTimeStr(const timeT &timeIn, int timeZone = 0)
{
   return "";
}




/// Get a date-time string in ISO 8601 format for time_t
/** Returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS, with optional timezone designation such as Z or +00:00.
  *
  * \overload
  *
  * \param [in] timeIn is the input time
  * \param [in] timeZone specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.
  *
  * \retval std::string containing the format date/time 
  * 
  * \ingroup timeutils
  */ 
template<> inline
std::string ISO8601DateTimeStr<time_t>(const time_t &timeIn, int timeZone)
{
   tm bdt;
   gmtime_r(&timeIn, &bdt);
   
   char tstr[25];
   
   strftime(tstr, 25, "%FT%H:%M:%S", &bdt);
   
   std::string result = tstr;
   
   if(timeZone == 1) result += "Z";
   if(timeZone == 2) result += "+00:00";
   
   return result;
}

/// Get a date-time string in ISO 8601 format for timespec
/** Returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS.SSSSSSSSS, with optional timezone designation such as Z or +00:00.
  *
  * \overload
  *
  * \param [in] timeIn is the input time
  * \param [in] timeZone specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.
  *
  * \retval std::string containing the format date/time 
  * 
  * \ingroup timeutils
  */ 
template<> inline
std::string ISO8601DateTimeStr<timespec>(const timespec &timeIn, int timeZone)
{
   std::string result = ISO8601DateTimeStr<time_t>(timeIn.tv_sec, 0);
   
   char tstr[20];
   
   snprintf(tstr, 20, ".%09ld", timeIn.tv_nsec);
   
   result += tstr;
   
   if(timeZone == 1) result += "Z";
   if(timeZone == 2) result += "+00:00";
   
   return result;
}

/// Get a date-time string in ISO 8601 format for the current UTC time
/** Returns a string in the ISO 8601 format:
  * YYYY-MM-DDYHH:MM:SS.SS, with optional timezone designation such as Z or +00:00.
  *
  * \overload
  *
  * \param [in] timeZone [optional] specifies whether to include a timezone designation.  0=> none, 1=> letter, 2=>offset.
  *
  * \retval std::string containing the format date/time 
  * 
  * \ingroup timeutils
  */ 
inline 
std::string ISO8601DateTimeStr(int timeZone = 0)
{
   return ISO8601DateTimeStr<time_t>(time(0), timeZone);
}


/// Convert a UTC timespec to TAI modified Julian date
/** Converts a timespec assumed to be in <a href="https://en.wikipedia.org/wiki/Coordinated_Universal_Time">Coordinated
  *  Universal Time (UTC)</a> to a <a href="https://en.wikipedia.org/wiki/Julian_day">
  *  Modified Julian Date (MJD)</a> in 
  * <a href="https://en.wikipedia.org/wiki/International_Atomic_Time">International Atomic Time (TAI)</a>.
  * 
  * \param [out] djm the modified Julian day number
  * \param [out] djmf the fraction of the day
  * \param [in] tsp contains the UTC time
  * \param [out] tm0 [optional] will be filled with the broken down UTC time
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
inline 
int timespecUTC2TAIMJD( double & djm, double & djmf, const timespec & tsp, tm * tm0)
{
   double dat, djm0;
   tm *tmrv;
   int rv1, rv2;
   
   // Get the broken down time corresponding to tsp0
   tmrv = gmtime_r(&tsp.tv_sec, tm0);
   if(tmrv == 0) return -10;
   
   // Then determine deltaAT = TAI-UTC
   rv1 = iauDat(1900+tm0->tm_year, 1+tm0->tm_mon, tm0->tm_mday, 0.0, &dat);
   if(rv1 < 0) return rv1;
   
   // And get the MJD
   rv2 = iauCal2jd(1900+tm0->tm_year, 1+tm0->tm_mon, tm0->tm_mday, &djm0, &djm);
   if(rv2 < 0) return rv2;
   
   // Finally calculate the day fraction
   djmf = ((double) tm0->tm_hour)/24.0 + ((double) tm0->tm_min)/(24.0*60.) + (((double) tm0->tm_sec)+((double) tsp.tv_nsec/1e9) + dat  )/(24.0*3600.0);
    
   if(djmf >= 1.0)
   {
      djmf -= 1.0;
      djm += 1.0;
   }
   
   if(rv1) return rv1;
   return 0;
}

} //namespace mx

#endif //__timeUtils_hpp__
