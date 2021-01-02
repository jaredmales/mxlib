/** \file zernike_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/sigproc/zernike.hpp"

/** Scenario: testing noll_nm
  * 
  * Verify calculation of Noll nm values from j.
  * Goes through each of the cases in Table 1 of \cite noll_1976
  * 
  * \anchor tests_sigproc_zernike_noll_nm
  */
SCENARIO( "testing noll_nm", "[sigproc::zernike]" ) 
{
   
   GIVEN("a j value")
   {
      WHEN("j==0")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,0);
         REQUIRE(rv == -1);
      }

      WHEN("j==1")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,1);
         REQUIRE(rv == 0);
         REQUIRE(n == 0);
         REQUIRE(m == 0);
      }
      
      WHEN("j==2")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,2);
         REQUIRE(rv == 0);
         REQUIRE(n == 1);
         REQUIRE(m == 1);
      }
      
      WHEN("j==3")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,3);
         REQUIRE(rv == 0);
         REQUIRE(n == 1);
         REQUIRE(m == -1);
      }
      
      WHEN("j==4")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,4);
         REQUIRE(rv == 0);
         REQUIRE(n == 2);
         REQUIRE(m == 0);
      }
      
      WHEN("j==5")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,5);
         REQUIRE(rv == 0);
         REQUIRE(n == 2);
         REQUIRE(m == -2);
      }
      
      WHEN("j==6")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,6);
         REQUIRE(rv == 0);
         REQUIRE(n == 2);
         REQUIRE(m == 2);
      }
      
      WHEN("j==7")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,7);
         REQUIRE(rv == 0);
         REQUIRE(n == 3);
         REQUIRE(m == -1);
      }
      
      WHEN("j==8")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,8);
         REQUIRE(rv == 0);
         REQUIRE(n == 3);
         REQUIRE(m == 1);
      }
      
      WHEN("j==9")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,9);
         REQUIRE(rv == 0);
         REQUIRE(n == 3);
         REQUIRE(m == -3);
      }
      
      WHEN("j==10")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,10);
         REQUIRE(rv == 0);
         REQUIRE(n == 3);
         REQUIRE(m == 3);
      }
      
      WHEN("j==11")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,11);
         REQUIRE(rv == 0);
         REQUIRE(n == 4);
         REQUIRE(m == 0);
      }
      
      WHEN("j==12")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,12);
         REQUIRE(rv == 0);
         REQUIRE(n == 4);
         REQUIRE(m == 2);
      }
      
      WHEN("j==13")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,13);
         REQUIRE(rv == 0);
         REQUIRE(n == 4);
         REQUIRE(m == -2);
      }
      
      WHEN("j==14")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,14);
         REQUIRE(rv == 0);
         REQUIRE(n == 4);
         REQUIRE(m == 4);
      }
      
      WHEN("j==15")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,15);
         REQUIRE(rv == 0);
         REQUIRE(n == 4);
         REQUIRE(m == -4);
      }
      
      WHEN("j==16")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,16);
         REQUIRE(rv == 0);
         REQUIRE(n == 5);
         REQUIRE(m == 1);
      }
      
      WHEN("j==17")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,17);
         REQUIRE(rv == 0);
         REQUIRE(n == 5);
         REQUIRE(m == -1);
      }
      
      WHEN("j==18")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,18);
         REQUIRE(rv == 0);
         REQUIRE(n == 5);
         REQUIRE(m == 3);
      }
      
      WHEN("j==19")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,19);
         REQUIRE(rv == 0);
         REQUIRE(n == 5);
         REQUIRE(m == -3);
      }
      
      WHEN("j==20")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,20);
         REQUIRE(rv == 0);
         REQUIRE(n == 5);
         REQUIRE(m == 5);
      }
      
      WHEN("j==21")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,21);
         REQUIRE(rv == 0);
         REQUIRE(n == 5);
         REQUIRE(m == -5);
      }
      
      WHEN("j==22")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,22);
         REQUIRE(rv == 0);
         REQUIRE(n == 6);
         REQUIRE(m == 0);
      }
      
      WHEN("j==23")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,23);
         REQUIRE(rv == 0);
         REQUIRE(n == 6);
         REQUIRE(m == -2);
      }
      
      WHEN("j==24")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,24);
         REQUIRE(rv == 0);
         REQUIRE(n == 6);
         REQUIRE(m == 2);
      }
      
      WHEN("j==25")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,25);
         REQUIRE(rv == 0);
         REQUIRE(n == 6);
         REQUIRE(m == -4);
      }
      
      WHEN("j==26")
      {
         int m, n;
         int rv = mx::sigproc::noll_nm(n,m,26);
         REQUIRE(rv == 0);
         REQUIRE(n == 6);
         REQUIRE(m == 4);
      }
   }
}

/** Scenario: testing zernikeQNorm
  * Verify compilation and execution of zernikeQNorm.
  * This does not validate the output. 
  * \anchor tests_sigproc_zernike_zernikeQNorm
  */
SCENARIO( "testing zernikeQNorm", "[sigproc::zernike]" ) 
{
   GIVEN("an array")
   {
      WHEN("j==1")
      {
         Eigen::Array<double, -1, -1> arr, k, phi;
         arr.resize(32,32);
         k.resize(32,32);
         phi.resize(32,32);
         
         for(int i=0;i<32;++i)
         {
            for(int j=0;j<32;++j)
            {
               double kx = i-5;
               double ky = j-15;
               k(i,j) = sqrt(kx*kx + ky*ky);
               phi(i,j) = atan(ky/kx);
            }
         }
         int rv = mx::sigproc::zernikeQNorm( arr, k, phi, 1);
         REQUIRE(rv == 0);

      }
   }
}
