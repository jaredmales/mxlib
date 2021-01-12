#include "../../catch2/catch.hpp"


#include "../../../include/math/templateBLAS.hpp"


SCENARIO( "testing scal", "[templateBLAS]" ) 
{
   GIVEN("a precision")
   {
      
      WHEN("precision is float")
      {
         int N = 2;
         float alpha = 2.;
         float x[] = {1.,2.};
         int incX = 1;
         mx::math::scal(N, alpha, x, incX);
         REQUIRE(x[0] == 2.0);
         REQUIRE(x[1] == 4.0);
      }

      WHEN("precision is double")
      {
         int N = 2;
         double alpha = 2.;
         double x[] = {1.,2.};
         int incX = 1;
         mx::math::scal(N, alpha, x, incX);
         REQUIRE(x[0] == 2.0);
         REQUIRE(x[1] == 4.0);
      }
   }
}

