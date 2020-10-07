/** \file psdFilter_test.cpp
 */
#define CATCH_CONFIG_MAIN
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/sigproc/psdFilter.hpp"
#include "../../../include/sigproc/psdUtils.hpp"
#include "../../../include/improc/eigenCube.hpp"
#include "../../../include/math/randomT.hpp"
#include "../../../include/math/vectorUtils.hpp"

/** Verify compilation and initilization of the 3 ranks for psdFilter.
  * 
  * \anchor tests_sigproc_psdFilter_compile
  */
SCENARIO( "compiling psdFilter", "[sigproc::psdFilter]" ) 
{
   GIVEN("a psdFilter, sqrt pointer")
   {
      WHEN("rank==1")
      {
         mx::sigproc::psdFilter<double, 1> psdF;
         
         std::vector<double> psdSqrt(1024,1);
         
         int rv = psdF.psdSqrt(&psdSqrt, 1);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 1024);
         REQUIRE(psdF.cols() == 1);
         REQUIRE(psdF.planes() == 1);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
         
      }
      WHEN("rank==2")
      {
         mx::sigproc::psdFilter<double, 2> psdF;
         
         Eigen::Array<double,-1,-1> psdSqrt;
         psdSqrt.resize(256,256);
         psdSqrt.setConstant(1);
         
         int rv = psdF.psdSqrt(&psdSqrt, 1, 1);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 256);
         REQUIRE(psdF.cols() == 256);
         REQUIRE(psdF.planes() == 1);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
      }
      WHEN("rank==3")
      {
         mx::sigproc::psdFilter<double, 3> psdF;
         
         mx::improc::eigenCube<double> psdSqrt;
         psdSqrt.resize(128,128,256);
         //psdSqrt.setConstant(1);
         
         int rv = psdF.psdSqrt(&psdSqrt, 1, 1, 1);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 128);
         REQUIRE(psdF.cols() == 128);
         REQUIRE(psdF.planes() == 256);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
      }
   }
   
   GIVEN("a psdFilter, sqrt reference")
   {
      WHEN("rank==1")
      {
         mx::sigproc::psdFilter<double, 1> psdF;
         
         std::vector<double> psdSqrt(1024,1);
         
         int rv = psdF.psdSqrt(psdSqrt, 1);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 1024);
         REQUIRE(psdF.cols() == 1);
         REQUIRE(psdF.planes() == 1);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
      }
      WHEN("rank==2")
      {
         mx::sigproc::psdFilter<double, 2> psdF;
         
         Eigen::Array<double,-1,-1> psdSqrt;
         psdSqrt.resize(256,256);
         psdSqrt.setConstant(1);
         
         int rv = psdF.psdSqrt(psdSqrt, 1, 1);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 256);
         REQUIRE(psdF.cols() == 256);
         REQUIRE(psdF.planes() == 1);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
      }
      WHEN("rank==3")
      {
         mx::sigproc::psdFilter<double, 3> psdF;
         
         mx::improc::eigenCube<double> psdSqrt;
         psdSqrt.resize(128,128,256);
         //psdSqrt.setConstant(1);
         
         int rv = psdF.psdSqrt(psdSqrt, 1, 1, 1);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 128);
         REQUIRE(psdF.cols() == 128);
         REQUIRE(psdF.planes() == 256);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
      }
   }
   
   GIVEN("a psdFilter, psd reference")
   {
      WHEN("rank==1")
      {
         mx::sigproc::psdFilter<double, 1> psdF;
         
         std::vector<double> psd(1024,1);
         
         int rv = psdF.psd(psd, 1.0);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 1024);
         REQUIRE(psdF.cols() == 1);
         REQUIRE(psdF.planes() == 1);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
      }
      WHEN("rank==2")
      {
         mx::sigproc::psdFilter<double, 2> psdF;
         
         Eigen::Array<double,-1,-1> psd;
         psd.resize(256,256);
         psd.setConstant(1);
         
         int rv = psdF.psd(psd, 1.0, 1.0);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 256);
         REQUIRE(psdF.cols() == 256);
         REQUIRE(psdF.planes() == 1);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
      }
      WHEN("rank==3")
      {
         mx::sigproc::psdFilter<double, 3> psdF;
         
         mx::improc::eigenCube<double> psd;
         psd.resize(128,128,256);
         //psdSqrt.setConstant(1);
         
         int rv = psdF.psd(psd, 1, 1, 1);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 128);
         REQUIRE(psdF.cols() == 128);
         REQUIRE(psdF.planes() == 256);
         
         psdF.clear();
         REQUIRE(psdF.rows() == 0);
         REQUIRE(psdF.cols() == 0);
         REQUIRE(psdF.planes() == 0);
      }
   }
}

/** Verify filtering and noise normalization
  * Conducts random noise tests, verifying that the resultant rms is within 2% of expected value on average over many trials.
  * Results are usually better than 1%, but 2% makes sure we don't get false failures.
  * 
  * \anchor tests_sigproc_psdFilter_filter
  */
SCENARIO( "filtering with psdFilter", "[sigproc::psdFilter]" ) 
{
   GIVEN("a rank 1 psd")
   {
      WHEN("alpha=-2.5, df Nyquist matched to array size, var=1")
      {
         mx::sigproc::psdFilter<double, 1> psdF;
         
         std::vector<double> f(2049), psd(2049);
         
         for(size_t n=0;n<psd.size();++n) f[n] = n*1.0/4096.;
         
         for(size_t n=1;n<psd.size();++n) psd[n] = pow(f[n], -2.5);
         psd[0] = psd[1];
         
         std::vector<double> f2s, psd2s;
         mx::sigproc::augment1SidedPSDFreq(f2s, f);
         mx::sigproc::augment1SidedPSD(psd2s, psd);
         
         double df = f2s[1]-f2s[0];
         mx::sigproc::normPSD(psd2s, f2s, 1.0, -1e5, 1e5);
         int rv = psdF.psd(psd2s, df);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 2.*psd.size()-2);
         REQUIRE(psdF.cols() == 1);
         REQUIRE(psdF.planes() == 1);
         
         std::vector<double> noise(psdF.rows());
         
         mx::math::normDistT<double> normVar;
         
         double avgRms = 0;
         
         for(int k=0;k<10000;++k)
         {
            for(size_t n=0;n<noise.size();++n) noise[n] = normVar;
            psdF(noise);
            avgRms += (mx::math::vectorVariance(noise,0.0));
         }
         
         avgRms = sqrt(avgRms/10000);
         
         REQUIRE(fabs(avgRms - 1.0) < 0.02);
         
      }
      WHEN("alpha=-1.5, df arbitrary, var = 2.2")
      {
         mx::sigproc::psdFilter<double, 1> psdF;
         
         std::vector<double> f(1025), psd(1025);
         
         for(size_t n=0;n<psd.size();++n) f[n] = n*1.0/7000.;
         
         for(size_t n=1;n<psd.size();++n) psd[n] = pow(f[n], -1.5);
         psd[0] = psd[1];
         
         std::vector<double> f2s, psd2s;
         mx::sigproc::augment1SidedPSDFreq(f2s, f);
         mx::sigproc::augment1SidedPSD(psd2s, psd);
         
         double df = f2s[1]-f2s[0];
         mx::sigproc::normPSD(psd2s, f2s, 2.2,-1e5,1e5);
         
         
         int rv = psdF.psd(psd2s, df);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == 2.*psd.size()-2);
         REQUIRE(psdF.cols() == 1);
         REQUIRE(psdF.planes() == 1);
         
         std::vector<double> noise(psdF.rows());
         
         mx::math::normDistT<double> normVar;
         
         double avgRms = 0;
         
         for(int k=0;k<10000;++k)
         {
            for(size_t n=0;n<noise.size();++n) noise[n] = normVar;
            psdF(noise);
            avgRms += (mx::math::vectorVariance(noise,0.0));
         }
         
         avgRms = sqrt(avgRms/10000);
         
         REQUIRE(fabs(avgRms - sqrt(2.2)) < 0.02*sqrt(2.2));
         
      }
   }
   GIVEN("a rank 2 psd")
   {
      WHEN("alpha=-2.5, dk Nyquist matched to array size, var=1")
      {
         mx::sigproc::psdFilter<double, 2> psdF;
         
         Eigen::Array<double, -1, -1> k, psd;
         
         k.resize(64, 64);
         psd.resize(64, 64);
         
         mx::sigproc::frequency_grid(k, 1./128., 1.0);
         for(int cc=0; cc< psd.cols(); ++cc)
         {
            for(int rr=0; rr<psd.rows(); ++rr)
            {
               if(k(rr,cc) == 0) psd(rr,cc) = 0;
               else psd(rr,cc) = pow(k(rr,cc), -2.5);
            }
         }
         
         double dk = k(0,1) - k(0,0);
         
         mx::sigproc::normPSD(psd, k, 1.0);
         
         int rv = psdF.psd(psd, dk, dk);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == psd.rows());
         REQUIRE(psdF.cols() == psd.cols());
         REQUIRE(psdF.planes() == 1);
         
         Eigen::Array<double, -1, -1> noise(psdF.rows(), psdF.cols());
         
         mx::math::normDistT<double> normVar;
         
         double avgRms = 0;
         
         for(int k=0;k<10000;++k)
         {
            for(int cc=0; cc< psd.cols(); ++cc)
            {
               for(int rr=0; rr<psd.rows(); ++rr)
               {
                  noise(rr,cc) = normVar;
               }
            }
            
            psdF(noise);
            avgRms += noise.square().sum(); //(mx::math::vectorVariance(noise,0.0));
         }
         
         avgRms = sqrt(avgRms/(psd.rows()*psd.cols())/10000);
         
         
         REQUIRE(fabs(avgRms - 1.0) < 0.02);
         
      }
      WHEN("alpha=-1.5, dk arb, var=2.2")
      {
         mx::sigproc::psdFilter<double, 2> psdF;
         
         Eigen::Array<double, -1, -1> k, psd;
         
         k.resize(64, 64);
         psd.resize(64, 64);
         
         mx::sigproc::frequency_grid(k, 1./302., 1.0);
         for(int cc=0; cc< psd.cols(); ++cc)
         {
            for(int rr=0; rr<psd.rows(); ++rr)
            {
               if(k(rr,cc) == 0) psd(rr,cc) = 0;
               else psd(rr,cc) = pow(k(rr,cc), -1.5);
            }
         }
         
         double dk = k(0,1) - k(0,0);
         
         mx::sigproc::normPSD(psd, k, 2.2);
         
         int rv = psdF.psd(psd, dk, dk);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == psd.rows());
         REQUIRE(psdF.cols() == psd.cols());
         REQUIRE(psdF.planes() == 1);
         
         Eigen::Array<double, -1, -1> noise(psdF.rows(), psdF.cols());
         
         mx::math::normDistT<double> normVar;
         
         double avgRms = 0;
         
         for(int k=0;k<10000;++k)
         {
            for(int cc=0; cc< psd.cols(); ++cc)
            {
               for(int rr=0; rr<psd.rows(); ++rr)
               {
                  noise(rr,cc) = normVar;
               }
            }
            
            psdF(noise);
            avgRms += noise.square().sum(); //(mx::math::vectorVariance(noise,0.0));
         }
         
         avgRms = sqrt(avgRms/(psd.rows()*psd.cols())/10000);
         
         
         REQUIRE(fabs(avgRms - sqrt(2.2)) < 0.02*sqrt(2.2));
         
      }
   }
   GIVEN("a rank 3 psd")
   {
      WHEN("k-alpha=-2.5, f-alph=-2.5, dk Nyquist matched to array size, df Nyquist matched to array size, var=1")
      {
         mx::sigproc::psdFilter<double, 3> psdF;
         
         Eigen::Array<double, -1, -1> k, psdk;
         std::vector<double> f, f2s, psd2s;
         
         mx::improc::eigenCube<double> psd;
         
         
         k.resize(32, 32);
         f.resize(33);
         
         
         mx::sigproc::frequency_grid(k, 1./64., 1.0);
         psdk.resize(k.rows(), k.cols());
         for(int cc=0; cc< psdk.cols(); ++cc)
         {
            for(int rr=0; rr<psdk.rows(); ++rr)
            {
               if(k(rr,cc) == 0) psdk(rr,cc) = 0;
               else psdk(rr,cc) = pow(k(rr,cc), -2.5);
            }
         }
         mx::sigproc::normPSD(psdk, k, 1.0);
         
         for(size_t n=0;n<f.size();++n) f[n] = n*1.0/64.;
         mx::sigproc::augment1SidedPSDFreq(f2s, f);
         psd2s.resize(f2s.size());
         for(size_t n=0;n<psd2s.size();++n) psd2s[n] = pow(fabs(f2s[n]), -2.5);
         psd2s[0] = psd2s[1];
         
         
         psd.resize(k.rows(), k.cols(), f2s.size());
         
         
         for(int cc=0; cc< psd.cols(); ++cc)
         {
            for(int rr=0; rr<psd.rows(); ++rr)
            {
               if(k(rr,cc) == 0) psd.pixel(rr,cc).setZero();
               else
               {
                  double p = psdk(rr,cc);
                  mx::sigproc::normPSD(psd2s, f2s, p, -1e5, 1e5);
                  
                  for(int pp=0;pp<psd.planes();++pp) psd.image(pp)(rr,cc)=psd2s[pp];
               }
            }
         }
         
         double dk = k(0,1) - k(0,0);
         double df = f[1]- f[0];
         
         int rv = psdF.psd(psd, dk, dk, df);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == psd.rows());
         REQUIRE(psdF.cols() == psd.cols());
         REQUIRE(psdF.planes() == psd.planes());
         
         mx::improc::eigenCube<double> noise(psdF.rows(), psdF.cols(), psdF.planes());
         
         mx::math::normDistT<double> normVar;
         
         double avgRms = 0;
         
         for(int k=0;k<10000;++k)
         {
            for(int pp=0;pp<psd.planes(); ++pp)
            {
               for(int cc=0; cc< psd.cols(); ++cc)
               {
                  for(int rr=0; rr<psd.rows(); ++rr)
                  {
                     noise.image(pp)(rr,cc) = normVar;
                  }
               }
            }
            
            psdF(noise);
            for(int pp=0;pp<noise.planes();++pp) avgRms += noise.image(pp).square().sum(); 
         }
         
         avgRms = sqrt(avgRms/(psd.rows()*psd.cols()*psd.planes())/10000);
         
         
         REQUIRE(fabs(avgRms - 1.0) < 0.02);
         
      }
      WHEN("k-alpha=-3.5, f-alph=-1.5, dk arb, df arb, var=2")
      {
         mx::sigproc::psdFilter<double, 3> psdF;
         
         Eigen::Array<double, -1, -1> k, psdk;
         std::vector<double> f, f2s, psd2s;
         
         mx::improc::eigenCube<double> psd;
         
         
         k.resize(32, 32);
         f.resize(33);
         
         
         mx::sigproc::frequency_grid(k, 1./640., 1.0);
         psdk.resize(k.rows(), k.cols());
         for(int cc=0; cc< psdk.cols(); ++cc)
         {
            for(int rr=0; rr<psdk.rows(); ++rr)
            {
               if(k(rr,cc) == 0) psdk(rr,cc) = 0;
               else psdk(rr,cc) = pow(k(rr,cc), -3.5);
            }
         }
         mx::sigproc::normPSD(psdk, k, 2.0);
         
         for(size_t n=0;n<f.size();++n) f[n] = n*1.0/78.;
         mx::sigproc::augment1SidedPSDFreq(f2s, f);
         psd2s.resize(f2s.size());
         for(size_t n=0;n<psd2s.size();++n) psd2s[n] = pow(fabs(f2s[n]), -1.5);
         psd2s[0] = psd2s[1];
         
         
         psd.resize(k.rows(), k.cols(), f2s.size());
         
         
         for(int cc=0; cc< psd.cols(); ++cc)
         {
            for(int rr=0; rr<psd.rows(); ++rr)
            {
               if(k(rr,cc) == 0) psd.pixel(rr,cc).setZero();
               else
               {
                  double p = psdk(rr,cc);
                  mx::sigproc::normPSD(psd2s, f2s, p, -1e5, 1e5);
                  
                  for(int pp=0;pp<psd.planes();++pp) psd.image(pp)(rr,cc)=psd2s[pp];
               }
            }
         }
         
         double dk = k(0,1) - k(0,0);
         double df = f[1]- f[0];
         
         int rv = psdF.psd(psd, dk, dk, df);
         
         REQUIRE(rv == 0);
         REQUIRE(psdF.rows() == psd.rows());
         REQUIRE(psdF.cols() == psd.cols());
         REQUIRE(psdF.planes() == psd.planes());
         
         mx::improc::eigenCube<double> noise(psdF.rows(), psdF.cols(), psdF.planes());
         
         mx::math::normDistT<double> normVar;
         
         double avgRms = 0;
         
         for(int k=0;k<10000;++k)
         {
            for(int pp=0;pp<psd.planes(); ++pp)
            {
               for(int cc=0; cc< psd.cols(); ++cc)
               {
                  for(int rr=0; rr<psd.rows(); ++rr)
                  {
                     noise.image(pp)(rr,cc) = normVar;
                  }
               }
            }
            
            psdF(noise);
            for(int pp=0;pp<noise.planes();++pp) avgRms += noise.image(pp).square().sum(); 
         }
         
         avgRms = sqrt(avgRms/(psd.rows()*psd.cols()*psd.planes())/10000);
         
         
         REQUIRE(fabs(avgRms - sqrt(2.0)) < 0.02*sqrt(2));
         
      }
   }
}
