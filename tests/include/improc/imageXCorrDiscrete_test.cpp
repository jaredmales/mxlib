/** \file imageXCorrDiscrete_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/func/gaussian.hpp"
#include "../../../include/improc/imageXCorrDiscrete.hpp"
#include "../../../include/improc/eigenCube.hpp"
#include "../../../include/ioutils/fits/fitsFile.hpp"

/** Scenario: centroiding Gaussians with center of light
  * 
  * Verify center of light calculation
  * 
  * \anchor tests_improc_imageUtils_imageXCorrDiscrete
  */
SCENARIO( "Verify X-Corr with center of light calculation", "[improc::imageXCorrDiscrete]" ) 
{
   GIVEN("two Gaussians")
   {
      WHEN("1 at geometric center")
      {
         mx::improc::eigenImage<double> im0, im2;
         im0.resize(64,64);
         im2.resize(64,64);
         
         mx::math::func::gaussian2D<double>(im0.data(), im0.rows(), im0.cols(), 0., 1.0, 31.5, 31.5, 2);
         mx::math::func::gaussian2D<double>(im2.data(), im2.rows(), im2.cols(), 0., 1.0, 31.5+4, 31.5+4, 2);
         
         double x, y;
         mx::improc::imageXCorrDiscrete<mx::improc::eigenImage<double>> xcf(5);
   
         mx::improc::eigenImage<double> refIm = im0.block(10,10, im0.rows()-11, im0.cols()-11);
         xcf.refIm(refIm);
         xcf.m_peakMethod = mx::improc::xcorrPeakMethod::centroid;
         
   
         xcf(x, y, im2);
   
         mx::fits::fitsFile<double> ff;
         ff.write("refIm.fits", xcf.m_refIm);
         ff.write("ccIm.fits", xcf.m_ccIm);
         
         std::cerr << x << " " << y << "\n";
         
         REQUIRE(fabs(x-2) < 1e-8 );
         REQUIRE(fabs(y-2) < 1e-8 );
      }
      WHEN("geometric quarter")
      {
         mx::improc::eigenImage<double> im;
         im.resize(64,64);
         
         mx::math::func::gaussian2D<double>(im.data(), im.rows(), im.cols(), 0., 1.0, 15.5, 15.5, 2);
         
         double x, y;
         mx::improc::imageCenterOfLight( x, y, im);
         
         REQUIRE(fabs(x-15.5) < 1e-8 );
         REQUIRE(fabs(y-15.5) < 1e-8 );
      }
   }
}

         

