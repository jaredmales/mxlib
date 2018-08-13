#define CATCH_CONFIG_MAIN
#include "../../catch2/catch.hpp"


#include "../../../include/math/templateLapack.hpp"


SCENARIO( "getting lamch values", "[templateLapack]" ) 
{
   GIVEN("a precision")
   {
      
      WHEN("precision is float")
      {
         //This is just a compilation test.
         float ch = mx::math::lamch<float>('E');
         REQUIRE(ch > 0);
      }
      WHEN("precision is double")
      {
         //This is just a compilation test.
         double ch = mx::math::lamch<double>('E');
         REQUIRE(ch > 0);
      }
   }
}

