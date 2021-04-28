/** \file factorial.hpp
  * \brief Declares and defines the factorial function.
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

#ifndef mx_factorial_hpp
#define mx_factorial_hpp

#include <type_traits>

#ifdef MX_INCLUDE_BOOST
#include <boost/math/special_functions/factorials.hpp>
#endif

namespace mx
{
namespace math
{
namespace func
{
   
/// The factorial function.
/**
  * \ingroup functions
  */ 
template<typename T>
T factorial( T x  /**< [in] the argument */)
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::factorial<T>(x);
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "factorial<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
float factorial<float>( float x);

template<>
double factorial<double>( double x);

template<>
long double factorial<long double>( long double x);

#ifdef HASQUAD
template<>
__float128 factorial<__float128>( __float128 x);
#endif

}
}
}

#endif //mx_factorial_hpp
