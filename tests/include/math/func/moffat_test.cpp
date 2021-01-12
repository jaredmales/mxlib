/** \file moffat_test.cpp
 */
#include "../../../catch2/catch.hpp"

#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../../include/math/func/moffat.hpp"

/** Scenario: compiling 1D Moffat function.
  *
  * 
  * \anchor tests_math_func_moffat1D
  */
SCENARIO( "compiling 1D Moffat function", "[math::func::moffat]" ) 
{
   GIVEN("the 1D Moffat function")
   {
      WHEN("compiling")
      {
         double mv = mx::math::func::moffat<double>(0., 0., 1., 0., 1., 1.);
         REQUIRE(mv == 1.0);
      }
   }
}

/** Scenario: compiling 2D Moffat function.
  *
  * 
  * \anchor tests_math_func_moffat2D
  */
SCENARIO( "compiling 2D Moffat function", "[math::func::moffat]" ) 
{
   GIVEN("the 2D Moffat function")
   {
      WHEN("compiling")
      {
         double mv = mx::math::func::moffat2D<double>(0., 0., 0., 1., 0., 0., 1., 1.);
         REQUIRE(mv == 1.0);
      }
   }
}

/** Scenario: compiling Moffat FWHM
  *
  * 
  * \anchor tests_math_func_moffatFWHM
  */
SCENARIO( "compiling Moffat FWHM", "[math::func::moffat]" ) 
{
   GIVEN("the Moffat FWHM")
   {
      WHEN("compiling")
      {
         double mv = mx::math::func::moffatFWHM<double>(1.,1.);
         REQUIRE(mv == 2.0);
      }
   }
}
