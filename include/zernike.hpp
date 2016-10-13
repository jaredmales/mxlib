
#ifndef __zernike_hpp__
#define __zernike_hpp__

#include <cmath>

#include <boost/math/special_functions/bessel.hpp>

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <boost/math/special_functions/sign.hpp>

// this function returns the Zernike coefficients corresponding to Noll
// index j (http://en.wikipedia.org/wiki/Zernike_polynomials)
void noll_nm( int & n, int & m, int j )
{

   n = ceil(-1.5 + sqrt(0.25 + 2*j) - 1e-10); // 1e-10 is to avoid wrong rounding due to numerical precision

   int  jrem = j - (n * (n+1)/2+1);
   m = (jrem + (jrem % 2) * abs( (n % 2)-1) + fabs( (jrem % 2)-1) * (n %  2)) * (-boost::math::sign( (j % 2)-0.5));

}




template<typename floatT>
floatT zernikeQNorm(int n, int m, floatT k, floatT phi)
{
   
   floatT B;

   if(k < 0.00001)
   {
      B = 1.0;
   }
   else
   {
      B = boost::math::cyl_bessel_j(n+1, 2*pi<floatT>()*k) / (pi<floatT>()*k);
   }

   floatT Q2 = (n+1) * (B * B);

   if (m != 0)
   {
      Q2 = Q2*2; 
    
      if (m > 0) // Even j (see Noll 1976)
      {
        Q2 = Q2 * pow(cos(m*phi), 2);
      }
      else //%Odd j (see Noll 1976)
      {
        Q2 = Q2 * pow(sin(-m*phi), 2);
      }
   }

   return Q2;
}


#endif //__zernike_hpp__

