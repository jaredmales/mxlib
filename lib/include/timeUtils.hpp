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



/** Parse a string of format hh:mm:ss.x
  * Breaks a time string into constituent parts.  Handles -h by distributing the sign to m and s.
  *
  * \param[in] hmsstr  is a string of format hh:mm:ss.x where ss.x can be of any precision
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
double Cal2mjd(int yr, int mon, int day, int hr, int min, double sec);
   
///Parse a FITS date of the form "YYYY-MM-DDTHH:MM:SS.S" and return the modified Julian date (MJD)
/** After parsing calls Cal2mjd.
  * 
  * \param [in] fdate is a standard FITS date string
  * 
  * \ingroup timeutils
  */
double FITSdate2mjd(const std::string & fdate);


} //namespace mx

#endif //__timeUtils_hpp__
