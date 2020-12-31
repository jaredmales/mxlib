/** \file timeUtils_test.cpp
 */
#include "../../catch2/catch.hpp"


#include "../../../include/sys/timeUtils.hpp"

/** Verify operation of get_curr_time.  This only 
  * checks for increasing time on subsequent calls.
  * 
  *  \anchor tests_sys_timeUtils_get_curr_time
  */
SCENARIO( "getting the current time in seconds", "[timeutils]")
{
   GIVEN("getting the time")
   {
      WHEN("the same timespec is used")
      {
         timespec ts;
         double t0 = mx::sys::get_curr_time(ts);
         REQUIRE(t0 > 0);
         double t1 = mx::sys::get_curr_time(ts);
         
         REQUIRE(t1 > t0);
         
      }
   
      WHEN("no timespec is provided")
      {
         double t0 = mx::sys::get_curr_time();
         REQUIRE(t0 > 0);
         double t1 = mx::sys::get_curr_time();
         
         REQUIRE(t1 > t0);
         
      }
   }
}

/** Verify operation of thread sleep functions.  
  * 
  * Uses mx::sys::get_curr_time to verify duration of sleep.
  * 
  *  \anchor tests_sys_timeUtils_sleep
  */
SCENARIO( "putting a thread to sleep", "[timeutils]" ) 
{
   GIVEN("sleeping for 1 second")
   {
      WHEN("sleeping in seconds")
      {
         double t0 = mx::sys::get_curr_time();
         mx::sys::sleep(1);
         double t1 = mx::sys::get_curr_time();
         
         REQUIRE(t1 >= t0 + 1.0);
      }
      WHEN("sleeping in milliseconds")
      {
         double t0 = mx::sys::get_curr_time();
         mx::sys::milliSleep(1000);
         double t1 = mx::sys::get_curr_time();
         
         REQUIRE(t1 >= t0 + 1.0);
         
      }
      WHEN("sleeping in microseconds")
      {
         double t0 = mx::sys::get_curr_time();
         mx::sys::microSleep(1000000);
         double t1 = mx::sys::get_curr_time();
         
         REQUIRE(t1 >= t0 + 1.0);
      }
      WHEN("sleeping in nanoseconds")
      {
         double t0 = mx::sys::get_curr_time();
         mx::sys::nanoSleep(1000000000);
         double t1 = mx::sys::get_curr_time();
         
         REQUIRE(t1 >= t0 + 1.0);
      }
      
   }   
}

/** Verify operation of timespecAddNsec.  
  * 
  *  \anchor tests_sys_timeUtils_timespecAddNsec
  */
SCENARIO( "adding time to a timespec", "[timeutils]" ) 
{
   GIVEN("a timespec")
   {
      WHEN("adding less than 1e9 nanoseconds")
      {
         timespec ts;
         ts.tv_sec = 1;
         ts.tv_nsec = 0;
         
         mx::sys::timespecAddNsec(ts, 10);
         REQUIRE(ts.tv_sec == 1);
         REQUIRE(ts.tv_nsec == 10);
         
         mx::sys::timespecAddNsec(ts, 100);
         REQUIRE(ts.tv_sec == 1);
         REQUIRE(ts.tv_nsec == 110);
         
         mx::sys::timespecAddNsec(ts, 1000);
         REQUIRE(ts.tv_sec == 1);
         REQUIRE(ts.tv_nsec == 1110);
         
         mx::sys::timespecAddNsec(ts, 10000);
         REQUIRE(ts.tv_sec == 1);
         REQUIRE(ts.tv_nsec == 11110);
         
         mx::sys::timespecAddNsec(ts, 100000);
         REQUIRE(ts.tv_sec == 1);
         REQUIRE(ts.tv_nsec == 111110);
         
         mx::sys::timespecAddNsec(ts, 1000000);
         REQUIRE(ts.tv_sec == 1);
         REQUIRE(ts.tv_nsec == 1111110);
         
         mx::sys::timespecAddNsec(ts, 10000000);
         REQUIRE(ts.tv_sec == 1);
         REQUIRE(ts.tv_nsec == 11111110);
         
         mx::sys::timespecAddNsec(ts, 100000000);
         REQUIRE(ts.tv_sec == 1);
         REQUIRE(ts.tv_nsec == 111111110);
      }
      
      WHEN("adding more than 1e9 nanoseconds but less than 2e9")
      {
         timespec ts;
         ts.tv_sec = 1;
         ts.tv_nsec = 0;
         
         mx::sys::timespecAddNsec(ts, 1000000000);
         REQUIRE(ts.tv_sec == 2);
         REQUIRE(ts.tv_nsec == 0);
         
         mx::sys::timespecAddNsec(ts, 1000000000+10);
         REQUIRE(ts.tv_sec == 3);
         REQUIRE(ts.tv_nsec == 10);
         
         mx::sys::timespecAddNsec(ts, 1000000000+100);
         REQUIRE(ts.tv_sec == 4);
         REQUIRE(ts.tv_nsec == 110);
         
         mx::sys::timespecAddNsec(ts, 1000000000+1000);
         REQUIRE(ts.tv_sec == 5);
         REQUIRE(ts.tv_nsec == 1110);
         
         mx::sys::timespecAddNsec(ts, 1000000000+10000);
         REQUIRE(ts.tv_sec == 6);
         REQUIRE(ts.tv_nsec == 11110);
         
         mx::sys::timespecAddNsec(ts, 1000000000+100000);
         REQUIRE(ts.tv_sec == 7);
         REQUIRE(ts.tv_nsec == 111110);
         
         mx::sys::timespecAddNsec(ts, 1000000000+1000000);
         REQUIRE(ts.tv_sec == 8);
         REQUIRE(ts.tv_nsec == 1111110);
         
         mx::sys::timespecAddNsec(ts, 1000000000+10000000);
         REQUIRE(ts.tv_sec == 9);
         REQUIRE(ts.tv_nsec == 11111110);
         
         mx::sys::timespecAddNsec(ts, 1000000000+100000000);
         REQUIRE(ts.tv_sec == 10);
         REQUIRE(ts.tv_nsec == 111111110);
      }
      
      WHEN("adding more than 2e9")
      {
         timespec ts;
         ts.tv_sec = 1;
         ts.tv_nsec = 0;
         
         mx::sys::timespecAddNsec(ts, 2000000010);
         REQUIRE(ts.tv_sec == 3);
         REQUIRE(ts.tv_nsec == 10);
      }
   }
}

/** Verify parsing of a formatted time string.
  *  Tests parsing of a string of format hh:mm:ss.s  
  * 
  *  \anchor tests_sys_timeUtils_parse_hms
  */
SCENARIO( "parsing a formatted time string", "[timeutils]" ) 
{
   GIVEN("a valid time string")
   {
      WHEN("integer seconds")
      {
         float hr;
         float mn;
         float sec;
         
         mx::sys::parse_hms(hr, mn, sec, "1:2:3");
         
         REQUIRE(hr == 1);
         REQUIRE(mn == 2);
         REQUIRE(sec == 3);
         
      }
      
      WHEN("floating seconds")
      {
         float hr;
         float mn;
         float sec;
         
         mx::sys::parse_hms(hr, mn, sec, "1:2:3.23");
         
         REQUIRE(hr == 1);
         REQUIRE(mn == 2);
         REQUIRE(fabs(sec-3.23) < 1e-7);
         
      }
      
      WHEN("negative hour")
      {
         float hr;
         float mn;
         float sec;
         
         mx::sys::parse_hms(hr, mn, sec, "-1:2:3.23");
         
         REQUIRE(hr == -1);
         REQUIRE(mn == -2);
         REQUIRE(fabs(sec - -3.23) < 1e-7);
         
      }
      
      WHEN("0 pads")
      {
         float hr;
         float mn;
         float sec;
         
         mx::sys::parse_hms(hr, mn, sec, "01:02:03.23");
         
         REQUIRE(hr == 1);
         REQUIRE(mn == 2);
         REQUIRE(fabs(sec-3.23) < 1e-7);
         
      }
   }
}

/** Verify calculation of MJD.
  * 
  *  \anchor tests_sys_timeUtils_Cal2mjd
  */
SCENARIO( "calculating MJD from a Gregorian date", "[timeutils]" ) 
{
   GIVEN("a valid Gregorian date")
   {
      WHEN("integer seconds")
      {
         double mjd = mx::sys::Cal2mjd(2020,12,31,0,0,0);
         REQUIRE(mjd == 59214.0);
      }
   }
   
   GIVEN("a valid Gregorian date")
   {
      WHEN("floating seconds")
      {
         double mjd = mx::sys::Cal2mjd(2020,12,31,0,0,10.2357);
         REQUIRE(fabs(mjd-59214.00011846875) < 1e-14);
      }
   }
}

/** Verify parsing of an ISO 8601 time string.
  * 
  *  \anchor tests_sys_timeUtils_ISO8601dateBreakdown
  */
SCENARIO( "parsing an ISO 8601 time string", "[timeutils]" ) 
{
   GIVEN("a valid ISO 8601 time string")
   {
      WHEN("integer seconds")
      {
         int yr, mon, day, hr, min;
         double sec;
         
         int rv =  mx::sys::ISO8601dateBreakdown(yr, mon, day, hr, min, sec, "2020-12-31T00:00:00");
         
         REQUIRE(rv == 0);
         REQUIRE(yr == 2020);
         REQUIRE(mon == 12);
         REQUIRE(day == 31);
         REQUIRE(hr == 0);
         REQUIRE(min == 0);
         REQUIRE(sec == 0);
         
         
      }
      
      WHEN("fractional seconds")
      {
         int yr, mon, day, hr, min;
         double sec;
         
         int rv =  mx::sys::ISO8601dateBreakdown(yr, mon, day, hr, min, sec, "2020-12-31T00:00:10.2357");
         
         REQUIRE(rv == 0);
         REQUIRE(yr == 2020);
         REQUIRE(mon == 12);
         REQUIRE(day == 31);
         REQUIRE(hr == 0);
         REQUIRE(min == 0);
         REQUIRE(fabs(sec-10.2357)<1e-14);
         
         
      }
   }
   
   GIVEN("an invalid ISO 8601 time string")
   {
      WHEN("string too short")
      {
         int yr, mon, day, hr, min;
         double sec;
         
         int rv =  mx::sys::ISO8601dateBreakdown(yr, mon, day, hr, min, sec, "2020-12-31");
         
         REQUIRE(rv == -4);
      }
   }
}

/** Verify conversion of an ISO 8601 time string to MJD.
  * 
  *  \anchor tests_sys_timeUtils_ISO8601date2mjd
  */
SCENARIO( "coverting an ISO 8601 time string to MJD", "[timeutils]" ) 
{
   GIVEN("a valid ISO 8601 time string")
   {
      WHEN("integer seconds")
      {
         double mjd =  mx::sys::ISO8601date2mjd("2020-12-31T00:00:00");
         
         REQUIRE(mjd == 59214.0);
      }
      
      WHEN("fractional seconds")
      {
         double mjd =  mx::sys::ISO8601date2mjd("2020-12-31T00:00:10.2357");
         
         REQUIRE(fabs(mjd-59214.00011846875) < 1e-14);
      }
   }
}
