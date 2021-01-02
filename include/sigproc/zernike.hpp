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
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_zernike_hpp
#define math_zernike_hpp

#include <cmath>

#include <vector>

#include "../math/func/bessel.hpp"
#include "../math/func/factorial.hpp"
#include "../math/func/sign.hpp"
#include "../math/constants.hpp"
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
/** Calculates the values of (n,m) for an index j following Noll (1976) \cite noll_1976
  * See also: http://en.wikipedia.org/wiki/Zernike_polynomials
  * 
  * If j is odd, this returns m <= 0.
  * 
  * 
  * \retval 0 on success
  * \retval -1 on error (j < 1)
  * 
  * \test Scenario: testing noll_nm \ref tests_sigproc_zernike_noll_nm "[test doc]"
  */
int noll_nm( int & n, ///< [out] n the radial index of the Zernike polynomial
             int & m, ///< [out] m the azimuthal index of the Zernnike polynomial.  m < 0 if j odd.
             int j    ///< [in]  j the Noll index, j > 0.
           );

/// Calculate the coefficients of a Zernike radial polynomial
/** 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
int zernikeRCoeffs( std::vector<realT> & c, ///< [out] allocated to length \f$ 0.5(n-m)+1\f$ and filled with the coefficients.
                    int n,                  ///< [in] the radial index of the Zernike polynomial.
                    int m                   ///< [in] the azimuthal index of the Zernike polynomial.
                  )
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
      c[k] = pow(-1.0, k) * math::func::factorial<realT>(n - k) / ( math::func::factorial<realT>(k) * math::func::factorial<realT>(0.5*(n+m)  - k)* math::func::factorial<realT>(0.5*(n-m)  - k));
   }
   
   return 0;
}

//Explicit instantiations:
extern template
int zernikeRCoeffs<float>( std::vector<float> & c, int n, int m);

extern template
int zernikeRCoeffs<double>( std::vector<double> & c, int n, int m);

extern template
int zernikeRCoeffs<long double>( std::vector<long double> & c, int n, int m);

#ifdef HASQUAD
extern template
int zernikeRCoeffs<__float128>( std::vector<__float128> & c, int n, int m);
#endif

/// Calculate the value of a Zernike radial polynomial at a given separation.
/** 
  * 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernikeR( realT rho,             ///< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$.
                int n,                 ///< [in] the radial index of the Zernike polynomial.
                int m,                 ///< [in] the azimuthal index of the Zernike polynomial.
                std::vector<realT> & c ///< [in] contains the radial polynomial coeeficients, and must be of length \f$ 0.5(n-m)+1\f$.
              )
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

extern template
float zernikeR(float rho, int n, int m, std::vector<float> & c);

extern template
double zernikeR(double rho, int n, int m, std::vector<double> & c);

extern template
long double zernikeR(long double rho, int n, int m, std::vector<long double> & c);

#ifdef HASQUAD
extern template
__float128 zernikeR(__float128 rho, int n, int m, std::vector<__float128> & c);
#endif

/// Calculate the value of a Zernike radial polynomial at a given separation.
/** 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernikeR( realT rho, ///< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$. 
                int n,     ///< [in] the radial index of the Zernike polynomial.
                int m      ///< [in] the azimuthal index of the Zernike polynomial.
              )
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

extern template
float zernikeR( float rho, int n, int m);

extern template
double zernikeR( double rho, int n, int m);

extern template
long double zernikeR( long double rho, int n, int m);

#ifdef HASQUAD
extern template
__float128 zernikeR( __float128 rho, int n, int m);
#endif

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
realT zernike( realT rho, 
               realT phi, 
               int n, 
               int m, 
               std::vector<realT> & c
             )
{
   realT azt;
   
   if( n == 0 && m == 0)
   {
      return 1.0;
   }
   
   if( m < 0 )
   {
      azt = math::root_two<realT>()*sin(-m*phi);
   }
   else if (m > 0)
   {
      azt = math::root_two<realT>()*cos(m*phi);
   }
   else
   {
      azt = 1.0;
   }
   
   return sqrt((realT) n+1) * zernikeR(rho, n, m, c) * azt;
   
}           

extern template
float zernike(float rho, float phi, int n, int m, std::vector<float> & c);

extern template
double zernike(double rho, double phi, int n, int m, std::vector<double> & c);

extern template
long double zernike(long double rho, long double phi, int n, int m, std::vector<long double> & c);

#ifdef HASQUAD
extern template
__float128 zernike(__float128 rho, __float128 phi, int n, int m, std::vector<__float128> & c);
#endif

/// Calculate the value of a Zernike radial polynomial at a given radius and angle.
/** 
  * 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernike( realT rho, ///< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$.
               realT phi, ///< [in] the azimuthal angle (in radians)
               int n,     ///< [in] the radial index of the Zernike polynomial.
               int m      ///< [in] the azimuthal index of the Zernike polynomial.
             )
{
   
   std::vector<realT> c;

   if( zernikeRCoeffs<realT>(c, n, m) < 0) return -9999;
   
   return  zernike(rho, phi, n, m, c);
}

extern template
float zernike( float rho, float phi, int n, int m);

extern template
double zernike( double rho, double phi, int n, int m);

extern template
long double zernike( long double rho, long double phi, int n, int m);

#ifdef HASQUAD
extern template
__float128 zernike( __float128 rho, __float128 phi, int n, int m);
#endif

/// Calculate the value of a Zernike radial polynomial at a given radius and angle.
/** 
  * \retval -9999 indicates a possible error
  * \retval R the value of the Zernike radial polynomial otherwise
  * 
  * \tparam realT is a real floating type
  */ 
template<typename realT>
realT zernike( realT rho,  ///< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$. 
               realT phi,  ///< [in] the azimuthal angle (in radians)
               int j       ///< [in] the Noll index of the Zernike polynomial.
             )
{
   int n, m;
   
   //Get n and m from j
   if(noll_nm(n, m, j) < 0) return -9999;
   
   return zernike(rho, phi, n, m);
}

extern template
float zernike(float rho, float phi, int j);

extern template
double zernike(double rho, double phi, int j);

extern template
long double zernike(long double rho, long double phi, int j);

#ifdef HASQUAD
extern template
__float128 zernike(__float128 rho, __float128 phi, int j);
#endif

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

///Fill in an Eigencube-like array with Zernike polynomials in Noll order
/** The cube is pre-allocated to set the image size and the number of modes.
  *
  * \returns 0 on success
  * \returns -1 on error
  */  
template<typename cubeT>
int zernikeBasis( cubeT & cube,                    ///< [in/out] the pre-allocated cube which will be filled with the Zernike basis
                  typename cubeT::Scalar rad = -1, ///< [in] [optional] the radius of the aperture.  If -1 then the full image size is used.
                  int minj = 2                     ///< [in] [optional] the minimum j value to include.  The default is j=2, which skips piston (j=1). 
                )          
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

///Calculate the square-normed Fourier transform of a Zernike polynomial at position (k,phi)
/** Implements Equation (8) of Noll (1976) \cite noll_1976.
  * 
  * \todo need a more robust jinc_n function for n > 1
  * 
  * \test Scenario: testing zernikeQNorm \ref tests_sigproc_zernike_zernikeQNorm "[test doc]" 
  * 
  * \returns the value of |Q(k,phi)|^2
  * 
  * \tparam realT is the floating point type used for arithmetic
  */
template<typename realT>
realT zernikeQNorm( realT k,   ///< [in] the radial coordinate of normalized spatial frequency
                    realT phi, ///< [in] the azimuthal coordinate of normalized spatial frequency
                    int n,     ///< [in] the Zernike polynomial n
                    int m      ///< [in] the Zernike polynomial m
                  )
{
   
   realT B;

   //sloppy implementation of jinc_n for k ~ 0
   if(k < 0.00001)
   {
      if( n == 0 ) B = 1.0;
      else B = 0.0;
   }
   else
   {
      B = math::func::bessel_j(n+1, math::two_pi<realT>()*k) / (math::pi<realT>()*k);
   }

   realT Q2 = (n+1) * (B * B);

   if (m > 0 ) // Even j (see Noll 1976)
   {
      Q2 = 2*Q2 * pow(cos(m*phi), 2);
   }
   else if( m < 0 ) //%Odd j (see Noll 1976)
   {
      Q2 = 2*Q2 * pow(sin(-m*phi), 2);
   }
  
   return Q2;
}

extern template
float zernikeQNorm<float>(float k, float phi, int n, int m);

extern template
double zernikeQNorm<double>(double k, double phi, int n, int m);

extern template
long double zernikeQNorm<long double>(long double k, long double phi, int n, int m);

#ifdef HASQUAD
extern template
__float128 zernikeQNorm<__float128>(__float128 k, __float128 phi, int n, int m);
#endif


///Calculate the square-normed Fourier transform of a Zernike polynomial at position (k,phi)
/** Implements Equation (8) of Noll (1976) \cite noll_1976.
  * 
  * \returns the value of |Q(k,phi)|^2
  * 
  * \tparam realT is the floating point type used for arithmetic
  * 
  */
template<typename realT>
realT zernikeQNorm( realT k,   ///< [in] the radial coordinate of normalized spatial frequency
                    realT phi, ///< [in] the azimuthal coordinate of normalized spatial frequency
                    int j      ///< [in] the Zernike polynomial index j (Noll convention)
                  )
{
   int n, m;
   
   noll_nm(n,m,j);
   
   return zernikeQNorm(k, phi, n, m);
}

/// Fill in an Eigen-like array with the square-normed Fourier transform of a Zernike polynomial
/** The array is filled in with the values of |Q(k,phi)|^2 according to Equation (8) of Noll (1976) \cite noll_1976.
  * 
  * \test Scenario: testing zernikeQNorm \ref tests_sigproc_zernike_zernikeQNorm "[test doc]" 
  *  
  * \returns 0 on success
  * \returns -1 on error
  * 
  * \tparam arrayT is the Eigen-like array type.  Arithmetic will be done in arrayT::Scalar.
  */
template<typename arrayT>
int zernikeQNorm( arrayT & arr, ///< [out] the allocated array. The rows() and cols() members are used to size the transform.
                  arrayT & k,   ///< [in] the normalized spatial frequency magnitude at each pixel
                  arrayT & phi, ///< [in] the spatial frequency angle at each pixel
                  int j         ///< [in] the polynomial index in the Noll convention \cite noll_1976
                )
{
   if(arr.rows() != k.rows() || arr.cols() != k.cols())
   {
      mxError("zernikeQNorm", MXE_INVALIDARG, "output array and input k are not the same size");
      return -1;
   }
   
   if(arr.rows() != phi.rows() || arr.cols() != phi.cols())
   {
      mxError("zernikeQNorm", MXE_INVALIDARG, "output array and input phi are not the same size");
      return -1;
   }
   
   int n,m;
   if(noll_nm(n,m,j) < 0) return -1; //noll_nm will explain error
   
   for(size_t i=0; i < arr.rows(); ++i)
   {
      for(size_t j=0; j < arr.cols(); ++j)
      {
         arr(i,j) = zernikeQNorm(k(i,j), phi(i,j), n, m);
      }
   }
   return 0;
}

///@} signal_processing

} //namespace sigproc 
} //namespace mx

#endif //math_zernike_hpp

