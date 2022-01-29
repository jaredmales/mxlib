/** \file imageTransforms_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"

#include "../../../include/improc/imageUtils.hpp"
#include "../../../include/improc/eigenCube.hpp"

#include "../../../include/improc/imageTransforms.hpp"

#include "../../../include/math/fit/fitGaussian.hpp"

/** Scenario: Verify direction and accuracy of various image shifts
  * 
  * Tests image shifts by fractional pixels.
  * 
  * \anchor tests_improc_imageTransforms_imageShift
  */
SCENARIO( "Verify direction and accuracy of various image shifts", "[improc::imageTransforms]" ) 
{
   GIVEN("a Gaussian image")
   {
      WHEN("shifting")
      {
         mx::improc::eigenImage<double> im, shift;
         mx::math::fit::fitGaussian2Dsym<double> fit;
         double x, y;

         im.resize(256,256);
         shift.resize(im.rows(), im.cols());

         fit.setArray(shift.data(), shift.rows(), shift.cols());
         
         //Use sigma = 8 to get a well oversampled image, making shifts more accurate
         mx::math::func::gaussian2D<double>(im.data(), im.rows(), im.cols(), 0., 1.0, 127.5, 127.5, 8);
         fit.setGuess(0.0, 1.0, 127.5, 127.5, 8.0);

         mx::improc::imageShift(shift, im, -0.5, -0.5, mx::improc::cubicConvolTransform<double>());
         fit.fit();
         REQUIRE(fabs(fit.x0() - 127.0) < 1e-4);//should be much better than this, but this is a test
         REQUIRE(fabs(fit.y0() - 127.0) < 1e-4);

         mx::improc::imageShift(shift, im, +0.5, +0.5, mx::improc::cubicConvolTransform<double>());
         fit.fit();
         REQUIRE(fabs(fit.x0() - 128.0) < 1e-4);
         REQUIRE(fabs(fit.y0() - 128.0) < 1e-4);

         mx::improc::imageShift(shift, im, +1.0, +1.0, mx::improc::cubicConvolTransform<double>());
         fit.fit();
         REQUIRE(fabs(fit.x0() - 128.5) < 1e-4);
         REQUIRE(fabs(fit.y0() - 128.5) < 1e-4);

   
         mx::improc::imageShift(shift, im, +0.5, -0.5, mx::improc::cubicConvolTransform<double>());
         fit.fit();
         REQUIRE(fabs(fit.x0() - 128.0) < 1e-4);
         REQUIRE(fabs(fit.y0() - 127.0) < 1e-4);

         mx::improc::imageShift(shift, im, -0.3, +0.7, mx::improc::cubicConvolTransform<double>());
         fit.fit();
         REQUIRE(fabs(fit.x0() - 127.2) < 1e-3); //non-0.5 pixel shifts are harder
         REQUIRE(fabs(fit.y0() - 128.2) < 1e-3);

         mx::improc::imageShift(shift, im, 1.3, -0.7, mx::improc::cubicConvolTransform<double>());
         fit.fit();
         REQUIRE(fabs(fit.x0() - 128.8) < 1e-3);
         REQUIRE(fabs(fit.y0() - 126.8) < 1e-3);
      }
   }
}

         

