#include "../catch2/catch.hpp"


#include "../../include/timeUtils.hpp"

SCENARIO( "putting a thread to sleep", "[timeutils]" ) 
{
   GIVEN("sleeping for 1 second")
   {
      WHEN("sleeping in seconds")
      {
         double t0 = mx::get_curr_time();
         mx::sleep(1);
         double t1 = mx::get_curr_time();
         
         REQUIRE(t1 >= t0 + 1.0);
      }
      WHEN("sleeping in milliseconds")
      {
         double t0 = mx::get_curr_time();
         mx::milliSleep(1000);
         double t1 = mx::get_curr_time();
         
         REQUIRE(t1 >= t0 + 1.0);
         
      }
      WHEN("sleeping in microseconds")
      {
         double t0 = mx::get_curr_time();
         mx::microSleep(1000000);
         double t1 = mx::get_curr_time();
         
         REQUIRE(t1 >= t0 + 1.0);
      }
      WHEN("sleeping in nanoseconds")
      {
         double t0 = mx::get_curr_time();
         mx::nanoSleep(1000000000);
         double t1 = mx::get_curr_time();
         
         REQUIRE(t1 >= t0 + 1.0);
      }
      
   }   
}

