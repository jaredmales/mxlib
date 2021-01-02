/** \file bessel.hpp
  * \brief Declares and defines Bessel functions of the first kind.
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef mx_bessel_hpp
#define mx_bessel_hpp

#include <type_traits>

#ifdef MX_INCLUDE_BOOST
#include <boost/math/special_functions/bessel.hpp>
#endif

namespace mx
{
namespace math
{
namespace func
{
   
/// Bessel Functions of the First Kind.
/**
  * \ingroup functions
  */ 
template<typename T1, typename T2>
T2 bessel_j( T1 v, ///< [in] 
             T2 x  ///< [in]
           )
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::cyl_bessel_j<T1, T2>(v,x);
#else
   static_assert(std::is_fundamental<T1>::value || !std::is_fundamental<T1>::value, "bessel_j<T1,T2> not specialized for type T1 and/or T2, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
float bessel_j<float, float>( float v, 
                              float x
                            );

template<>
float bessel_j<int, float>( int v,
                            float x
                          );

template<>
double bessel_j<double, double>( double v, 
                                 double x
                               );

template<>
double bessel_j<int, double>( int v,
                              double x
                            );

template<>
long double bessel_j<long double, long double>( long double v, 
                                 long double x
                               );

template<>
long double bessel_j<int, long double>( int v,
                              long double x
                            );

#ifdef HASQUAD
template<>
__float128 bessel_j<__float128, __float128>( __float128 v, 
                                             __float128 x
                                           );

template<>
__float128 bessel_j<int, __float128>( int v,
                                      __float128 x
                                    );
#endif

}
}
}

#endif //mx_bessel_hpp
