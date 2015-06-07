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


#include "time.h"
#include "sys/time.h"
#include "stringUtils.hpp"

namespace mx
{
   
/** \addtogroup timeutils
  * @{
  */
  
/** Get the current system time.
  * Uses timespec, so nanosecond resolution is possible.
  * 
  * \tparam typeT is the type to return the time in [default=double].  
  *               must have cast from integer types, and + and / operators defined.
  * \tparam clk_id is the sys/time.h clock identifier [default=CLOCK_REALTIME]
  * 
  * \returns the current time in seconds
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
  * \param hmsstr [input] is a string of format hh:mm:ss.x where ss.x can be of any precision
  * \param h [output] is the hour component coverted to floatT
  * \param m [output] is the minute component converted to floatT
  * \param s [output] is the second component converted to floatT
  * 
  * \tparam floatT is a floating point type
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


/// @}

} //namespace mx

#endif //__timeUtils_hpp__
