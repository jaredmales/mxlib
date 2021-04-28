/** \file astroDynamics_test.cpp
 */

#include "../../catch2/catch.hpp"


#include "../../../include/astro/astroDynamics.hpp"

/** Scenario: calculating parallactic angles
  *
  * Compares par. ang. calc vs. an actual observation.
  * 
  * \anchor tests_astrodynamics_parang
  */
SCENARIO( "calculating parallactic angles", "[astroDynamics::parang]" ) 
{
   GIVEN("a typical target")
   {
      //Using an actual obervation of beta Pic and FITS headers.
      
      double dec = -51.06575; //beta Pic.
      double lat = -29.015; //LCO
      WHEN("before transit")
      {
         double ha = -0.739167/24.*360.;

         double parang = mx::astro::parAngDeg( ha, dec, lat, true);
         double diff = fabs( parang - (-24.903)); //The result from the TCS.

         REQUIRE(parang < 0);         
         REQUIRE(diff <= 1e-1); //This is just a roughly close enough test -- don't have exact precision of TCS, etc.
      }

      WHEN("after transit")
      {
         double ha = 2.008889/24.*360.;

         double parang = mx::astro::parAngDeg( ha, dec, lat, true);
         double diff = fabs( parang - (57.1241)); //The result from the TCS.

         REQUIRE(parang > 0);         
         REQUIRE(diff <= 1e-1); //This is just a roughly close enough test -- don't have exact precision of TCS, etc.
      }
      
   }
}

