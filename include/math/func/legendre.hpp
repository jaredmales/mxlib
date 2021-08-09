/** \file legendre.hpp
  * \brief Declares and defines Legendre polynomials.
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#ifndef mx_legendre_hpp
#define mx_legendre_hpp

#include <type_traits>

#ifdef MX_INCLUDE_BOOST
#include <boost/math/special_functions/legendre.hpp>
#endif

namespace mx
{
namespace math
{
namespace func
{
   
/// Legendre Polynomials
/** See https://en.wikipedia.org/wiki/Legendre_polynomials
  * and https://www.boost.org/doc/libs/1_75_0/libs/math/doc/html/math_toolkit/sf_poly/legendre.html
  * 
  * \returns the value of the n-th Legendre polynomial at x.
  * 
  * \ingroup functions
  */ 
template<typename T>
T legendre_p( int n, ///< [in] the order of the Legendre polynomial, n>=0.
              T x    ///< [in] the argument, -1 <= x <= 1.
            )
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::legendre_p<T>(n, x);
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "legendre_p<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
float legendre_p<float>( int n,
                         float x 
                       );

template<>
double legendre_p<double>( int n, 
                           double x
                         );

template<>
long double legendre_p<long double>( int n, 
                                     long double x
                                   );

#ifdef HASQUAD
template<>
__float128 legendre_p<__float128>( int n, 
                                   __float128 x
                                 );
#endif


/// The orthonormal Legendre polynomials
/** A version of the Legendre polynomials which are orthonormal on the interval
  * -1 <= x <= 1.
  * 
  * \returns the value of the n-th orthonormal Legendre polynomial at x.
  */ 
template<typename T>
T orthoNormalLegendre( int n,
                       T x
                     )
{
   return sqrt((static_cast<T>(2*n+1))/static_cast<T>(2)) * legendre_p<T>(n,x);
}

} //namespace func
} //namespace math
} //namespace mx

#endif //mx_legendre_hpp
