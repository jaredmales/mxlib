/** \file zernike.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Working with the Zernike polynomials.
  * 
  * \todo the basic zernike polys should be in math::func.
  * 
  * \ingroup signal_processing_files
  * 
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef zernike_hpp
#define zernike_hpp

#include <cmath>

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <boost/math/special_functions/bessel.hpp>

#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/factorials.hpp>


#include "../mxError.hpp"

namespace mx
{
namespace sigproc 
{
   
/** 
  * \ingroup zernike_basis
  * @{
  */


/// Get the Zernike coefficients n,m corrresponding the Noll index j.
/** Calculates the single index j for (n,m) following Noll (1976) \cite noll_1976
  * See also: http://en.wikipedia.org/wiki/Zernike_polynomials
  * 
  * \param[out] n the radial index of the Zernike polynomial
  * \param[out] m the azimuthal index of the Zernnike polynomial  
  * \param[in]  j the Noll index, \f$j > 0\f$
  * 
  * \retval 0 on success
  * \retval -1 on error (j < 1)
  */
int noll_nm( int & n, int & m, int j )
{
   if( j < 1)
   {
      mxError("noll_nm", MXE_INVALIDARG, "The Noll index j cannot be less than 1 in the Zernike polynomials");
      return -1;
   }
   
   
   n = ceil(-1.5 + sqrt(0.25 + 2*j) - 1e-10); // 1e-10 is to avoid wrong rounding due to numerical precision

   int  jrem = j - (n * (n+1)/2+1);
   m = (jrem + (jrem % 2) * abs( (n % 2)-1) + fabs( (jrem % 2)-1) * (n %  2)) * (-boost::math::sign( (j % 2)-0.5));
   
   return 0;
}

/// Calculate the coefficients of a Zernike radial polynomial
/** 
  * \param[out] c is allocated to length \f$ 0.5(n-m)+1\f$ and filled with the coefficients.
  * \param[in] n is the radial index of the Zernike polynomial.
  * \param[in] m is the azimuthal index of the Zernike polynomial.
  * 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
int zernikeRCoeffs( std::vector<realT> & c, int n, int m)
{
   m = abs(m);
   
   if( n < m)
   {
      mxError("zernikeRCoeffs", MXE_INVALIDARG, "n cannot be less than m in the Zernike polynomials");
      return -1;
   }
   
   //If odd, it's 0.
   if( (n-m) % 2 > 0)
   {
      c.resize(1,0);
      return 0;
   }
   
   int ul = 0.5*(n-m) + 1;
   
   c.resize(ul);
   
   for(int k=0; k < ul; ++k)
   {
      c[k] = pow(-1.0, k) * boost::math::factorial<realT>(n - k) / ( boost::math::factorial<realT>(k) * boost::math::factorial<realT>(0.5*(n+m)  - k)* boost::math::factorial<realT>(0.5*(n-m)  - k));
   }
   
   return 0;
}

/// Calculate the value of a Zernike radial polynomial at a given separation.
/** 
  * \param[in] rho is the radial coordinate, \f$ 0 \le \rho \le 1 \f$. 
  * \param[in] n is the radial index of the Zernike polynomial.
  * \param[in] m is the azimuthal index of the Zernike polynomial.
  * \param[in] c is contains the radial polynomial coeeficients, and must be of length \f$ 0.5(n-m)+1\f$.
  * 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernikeR( realT rho, int n, int m, std::vector<realT> & c )
{
   m = abs(m);
   
   //If odd, it's 0.
   if( (n-m) % 2 > 0)
   {
      c.resize(1,0);
      return 0.0;
   }
   
   if( c.size() != 0.5*(n-m) + 1)
   {
      mxError("zernikeR", MXE_INVALIDARG, "c vector has incorrect length for n and m.");
      return -9999;
   }
   
   realT R = 0.0;
   for(size_t k=0; k< c.size(); ++k)
   {
      R += c[k] * pow(rho, n-2*k);
   }
   
   return R;
   
}

/// Calculate the value of a Zernike radial polynomial at a given separation.
/** 
  * \param[in] rho is the radial coordinate, \f$ 0 \le \rho \le 1 \f$. 
  * \param[in] n is the radial index of the Zernike polynomial.
  * \param[in] m is the azimuthal index of the Zernike polynomial.
  * 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernikeR( realT rho, int n, int m)
{
   m = abs(m);
   
   //If odd, it's 0.
   if( (n-m) % 2 > 0)
   {
      return 0.0;
   }
   
   std::vector<realT> c;
   
   if(zernikeRCoeffs(c, n, m) < 0) return -9999;
   
   return zernikeR(rho, n, m, c);
   
}

/// Calculate the value of a Zernike radial polynomial at a given radius and angle.
/** 
  * \param[in] rho is the radial coordinate, \f$ 0 \le \rho \le 1 \f$. 
  * \param[in] phi is the azimuthal angle (in radians)
  * \param[in] n is the radial index of the Zernike polynomial.
  * \param[in] m is the azimuthal index of the Zernike polynomial.
  * \param[in] c is contains the radial polynomial coeeficients, and must be of length \f$ 0.5(n-m)+1\f$.
  * 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernike(realT rho, realT phi, int n, int m, std::vector<realT> & c)
{
   realT azt;
   
   if( n == 0 && m == 0)
   {
      return 1.0;
   }
   
   if( m < 0 )
   {
      azt = root_two<realT>()*sin(-m*phi);
   }
   else if (m > 0)
   {
      azt = root_two<realT>()*cos(m*phi);
   }
   else
   {
      azt = 1.0;
   }
   
   return sqrt((realT) n+1) * zernikeR(rho, n, m, c) * azt;
   
}           

/// Calculate the value of a Zernike radial polynomial at a given radius and angle.
/** 
  * \param[in] rho is the radial coordinate, \f$ 0 \le \rho \le 1 \f$. 
  * \param[in] phi is the azimuthal angle (in radians)
  * \param[in] n is the radial index of the Zernike polynomial.
  * \param[in] m is the azimuthal index of the Zernike polynomial.
  * 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernike(realT rho, realT phi, int n, int m)
{
   
   std::vector<realT> c;

   if( zernikeRCoeffs<realT>(c, n, m) < 0) return -9999;
   
   return  zernike(rho, phi, n, m, c);
}

/// Calculate the value of a Zernike radial polynomial at a given radius and angle.
/** 
  * \param[in] rho is the radial coordinate, \f$ 0 \le \rho \le 1 \f$. 
  * \param[in] phi is the azimuthal angle (in radians)
  * \param[in] j is the Noll index of the Zernike polynomial.
  * 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernike(realT rho, realT phi, int j)
{
   int n, m;
   
   //Get n and m from j
   if(noll_nm(n, m, j) < 0) return -9999;
   
   return zernike(rho, phi, n, m);
}

///Fill in an Eigen-like array with a Zernike polynomial
/** Sets any pixel which is at rad \<= r \< rad+0.5 pixels to rho = 1, to be consistent with mx::circularPupil
  *
  * \param[out] arr is the allocated array with an Eigen-like interface. The rows() and cols() members are used to size the polynomial.
  * \param[in] n is the radial index of the polynomial
  * \param[in] m is the azimuthal index of the polynomial
  * \param[in] xcen is the x coordinate of the desired center of the polynomial, in pixels
  * \param[in] ycen is the y coordinate of the desired center of the polynomial, in pixels
  * \param[in] rad [optional] is the desired radius. If rad \<= 0, then the maximum radius based on dimensions of m is used.
  */
template<typename arrayT, int overscan=2>
int zernike( arrayT & arr, 
             int n, 
             int m, 
             typename arrayT::Scalar xcen, 
             typename arrayT::Scalar ycen,
             typename arrayT::Scalar rad = -1 )
{
   typename arrayT::Scalar x;
   typename arrayT::Scalar y;
   typename arrayT::Scalar r, rho;
   typename arrayT::Scalar phi;
   
   std::vector<typename arrayT::Scalar> c;
   
   if(zernikeRCoeffs(c, n, m) < 0) return -1;
   
   size_t l0 = arr.rows();
   size_t l1 = arr.cols();
   
   if(rad <= 0) rad = 0.5*std::min(l0-1, l1-1);
         
   for(size_t i=0; i < l0; ++i)
   {
      for(size_t j=0; j < l1; ++j)
      {
         x = i - xcen;
         y = j - ycen;
         
         
         r = std::sqrt( x*x + y*y );
         
         //This is to be consistent with mx::circularPupil while still respecting the Zernike rules
         if(r  > rad && r <= rad+(1.0/overscan)) r = rad;
         
         rho =  r / rad;
         
         if(rho <= 1.0)
         {
            phi = std::atan2( y, x);         
            arr(i,j) = zernike(rho, phi, n, m, c);
         }
         else
         {
            arr(i,j) = 0.0;
         }
      }
   }
   return 0;
}

///Fill in an Eigen-like array with a Zernike polynomial
/** Sets any pixel which is at rad \<= r \<= rad+0.5 pixels to rho = 1, to be consistent with mx::circularPupil
  *
  * \param[out] arr is the allocated array with an Eigen-like interface. The rows() and cols() members are used to size the polynomial.
  * \param[in] j is the Noll index of the polynomial
  * \param[in] xcen is the x coordinate of the desired center of the polynomial, in pixels
  * \param[in] ycen is the y coordinate of the desired center of the polynomial, in pixels
  * \param[in] rad [optional] is the desired radius. If rad \<= 0, then the maximum radius based on dimensions of m is used.
  */
template<typename arrayT>
int zernike( arrayT & arr, 
             int j, 
             typename arrayT::Scalar xcen, 
             typename arrayT::Scalar ycen,
             typename arrayT::Scalar rad = -1 )
{
   int n, m;
   
   if(noll_nm(n, m, j) < 0) return -1;
   
   return zernike(arr, n, m, xcen, ycen, rad);
}

///Fill in an Eigen-like array with a Zernike polynomial
/** The geometric center of the array, 0.5*(arr.rows()-1), 0.5*(arr.cols()-1), is used as the center.
  * Sets any pixel which is at rad \<= r \< rad+0.5 pixels to rho = 1, to be consistent with mx::circularPupil
  *
  * \param[out] arr is the allocated array with an Eigen-like interface. The rows() and cols() members are used to size the polynomial.
  * \param[in] j is the Noll index of the polynomial
  * \param[in] rad [optional] is the desired radius. If rad \<= 0, then the maximum radius based on dimensions of m is used.
  */
template<typename arrayT>
int zernike( arrayT & arr, 
             int n, 
             int m, 
             typename arrayT::Scalar rad = -1 )
{
   typename arrayT::Scalar xcen = 0.5*(arr.rows()-1.0);
   typename arrayT::Scalar ycen = 0.5*(arr.cols()-1.0);
   
   return zernike(m, n, m, xcen, ycen, rad);
}

///Fill in an Eigen-like array with a Zernike polynomial
/** The geometric center of the array, 0.5*(arr.rows()-1), 0.5*(arr.cols()-1), is used as the center.
  * Sets any pixel which is at rad \<= r \< rad+0.5 pixels to rho = 1, to be consistent with mx::circularPupil
  *
  * \param[out] arr is the allocated array with an Eigen-like interface. The rows() and cols() members are used to size the polynomial.
  * \param[in] j is the Noll index of the polynomial
  * \param[in] rad [optional] is the desired radius. If rad \<= 0, then the maximum radius based on dimensions of m is used.
  */
template<typename arrayT>
int zernike( arrayT & arr, 
             int j,
             typename arrayT::Scalar rad = -1 )
{
   typename arrayT::Scalar xcen = 0.5*( arr.rows()-1.0);
   typename arrayT::Scalar ycen = 0.5*( arr.cols()-1.0);
   
   return zernike(arr, j, xcen, ycen, rad);
}

template<typename cubeT>
int zernikeBasis( cubeT & cube,
                  typename cubeT::Scalar rad = -1,
                  int minj = 2 )
{
   typename cubeT::imageT im;
   
   im.resize(cube.rows(), cube.cols());
   
   int rv;
   for(int i=0; i < cube.planes(); ++i)
   {
      rv = zernike( im, minj+i, rad);
   
      if(rv < 0) 
      {
         return rv;
      }
      cube.image(i) = im;
   }
   
   return 0;
}

template<typename realT>
realT zernikeQNorm(int n, int m, realT k, realT phi)
{
   
   realT B;

   if(k < 0.00001)
   {
      B = 1.0;
   }
   else
   {
      B = boost::math::cyl_bessel_j(n+1, 2*pi<realT>()*k) / (pi<realT>()*k);
   }

   realT Q2 = (n+1) * (B * B);

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

template<typename realT>
realT zernikeQNorm(int j, realT k, realT phi)
{
   int n, m;
   
   noll_nm(n,m,j);
   
   return zernikeQNorm(n, m, k, phi);
}

///@} signal_processing

} //namespace sigproc 
} //namespace mx

#endif //zernike_hpp

