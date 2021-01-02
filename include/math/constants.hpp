/** \file math.hpp
  * \author Jared R. Males
  * \brief Definitions of constants
  * \ingroup gen_math_files
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


#ifndef math_constants_hpp
#define math_constants_hpp

#include <type_traits>

#ifdef MX_INCLUDE_BOOST
#include <boost/math/constants/constants.hpp>
#endif

namespace mx
{
namespace math
{
   
/// Get the value of pi
/** Wrapper for boost constant.  Specializations provided for float, double, and long double.
  *
  * \ingroup genconstants 
  */ 
template<typename T>
constexpr T pi() 
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::constants::pi<T>();
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "pi<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

#define MX_INTERNAL_PI_100 (3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679)
#define MX_INTERNAL_ROOT2_100 (1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727)

template<>
constexpr float pi<float>()
{
   return static_cast<float>(MX_INTERNAL_PI_100);
}

template<>
constexpr double pi<double>()
{
   return static_cast<double>(MX_INTERNAL_PI_100);
}


template<>
constexpr long double pi<long double>()
{
   return static_cast<long double>(MX_INTERNAL_PI_100);
}

#ifdef HASQUAD
template<>
constexpr __float128 pi<__float128>()
{
   return static_cast<__float128>(MX_INTERNAL_PI_100);
}
#endif

/// Get the value of 2pi
/** Specializations provided for float, double, and long double.
  *
  * \ingroup genconstants 
  */ 
template<typename T>
constexpr T two_pi() 
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::constants::two_pi<T>();
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "two_pi<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
constexpr float two_pi<float>()
{
   return static_cast<float>(2*MX_INTERNAL_PI_100);
}

template<>
constexpr double two_pi<double>()
{
   return static_cast<double>(2*MX_INTERNAL_PI_100);
}

template<>
constexpr long double two_pi<long double>()
{
   return static_cast<long double>(2*MX_INTERNAL_PI_100);
}

#ifdef HASQUAD
template<>
constexpr __float128 two_pi<__float128>()
{
   return static_cast<__float128>(2*MX_INTERNAL_PI_100);
}
#endif

/// Get the value of pi/2
/** Wrapper for boost constant. Specializations provided for float, double, and long double.
  *
  * \ingroup genconstants 
  */ 
template<typename T>
constexpr T half_pi() 
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::constants::half_pi<T>();
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "half_pi<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
constexpr float half_pi<float>()
{
   return static_cast<float>(MX_INTERNAL_PI_100)/static_cast<float>(2);
}

template<>
constexpr double half_pi<double>()
{
   return static_cast<double>(MX_INTERNAL_PI_100)/static_cast<double>(2);
}

template<>
constexpr long double half_pi<long double>()
{
   return static_cast<long double>(MX_INTERNAL_PI_100)/static_cast<long double>(2);
}

/// Get the value of 180/pi
/** Wrapper for boost constant.  Specializations provided for float, double, and long double.
  *
  * \ingroup genconstants 
  */ 
template<typename T>
constexpr T rad2deg() 
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::constants::radian<T>();
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "rad2deg<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
constexpr float rad2deg<float>()
{
   return static_cast<float>(180)/static_cast<float>(MX_INTERNAL_PI_100);
}

template<>
constexpr double rad2deg<double>()
{
   return static_cast<double>(180)/static_cast<double>(MX_INTERNAL_PI_100);
}

template<>
constexpr long double rad2deg<long double>()
{
   return static_cast<long double>(180)/static_cast<long double>(MX_INTERNAL_PI_100);
}

/// Get the value of sqrt(2)
/** Wrapper for boost constant. Specializations provided for float, double, and long double.
  *
  * \ingroup genconstants 
  */ 
template<typename T>
constexpr T root_two() 
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::constants::root_two<T>();
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "root_two<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

#define MX_INTERNAL_ROOT2_100 (1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727)

template<>
constexpr float root_two<float>()
{
   return static_cast<float>(MX_INTERNAL_ROOT2_100);
}

template<>
constexpr double root_two<double>()
{
   return static_cast<double>(MX_INTERNAL_ROOT2_100);
}

template<>
constexpr long double root_two<long double>()
{
   return static_cast<long double>(MX_INTERNAL_ROOT2_100);
}

/// Get the value of 1/3
/** Wrapper for boost constant.  Specializations provided for float, double, and long double.
  *
  * \ingroup genconstants 
  */ 
template<typename T>
constexpr T third() 
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::constants::third<T>();
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "third<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
constexpr float third<float>()
{
   return static_cast<float>(1)/static_cast<float>(3);
}

template<>
constexpr double third<double>()
{
   return static_cast<double>(1)/static_cast<double>(3);
}

template<>
constexpr long double third<long double>()
{
   return static_cast<long double>(1)/static_cast<long double>(3);
}


///Return 5/3 in the specified precision
/** This constant is used frequently in adaptive optics analysis. 
  * \ingroup genconstants
  */
template<typename floatT>
constexpr floatT five_thirds()
{
   return static_cast<floatT>(5)/static_cast<floatT>(3);
}


///Return 5/6 in the specified precision
/** This constant is used frequently in adaptive optics analysis.
  * \ingroup genconstants
  *
  */
template<typename floatT>
constexpr floatT five_sixths()
{
   return static_cast<floatT>(5)/static_cast<floatT>(6);
}

///Return 11/3 in the specified precision
/** This constant is used frequently in adaptive optics analysis.
  * \ingroup genconstants
  * 
  */
template<typename floatT>
constexpr floatT eleven_thirds()
{
   return static_cast<floatT>(11)/static_cast<floatT>(3);
}

///Return 11/6 in the specified precision
/** This constant is used frequently in adaptive optics analysis.
  * \ingroup genconstants
  * 
  */
template<typename floatT>
constexpr floatT eleven_sixths()
{
   return static_cast<floatT>(11)/static_cast<floatT>(6);
}

///Return 6/5 in the specified precision
/** This constant is used frequently in adaptive optics analysis.
  * \ingroup genconstants
  *
  */
template<typename floatT>
constexpr floatT six_fifths()
{
   return static_cast<floatT>(6)/static_cast<floatT>(5);
}

///Return 3/5 in the specified precision
/** This constant is used frequently in adaptive optics analysis.  
  * \ingroup genconstants
  *  
  */
template<typename floatT>
constexpr floatT three_fifths()
{
   return static_cast<floatT>(3)/static_cast<floatT>(5);
}

///Return 17/3 in the specified precision
/** This constant is used frequently in adaptive optics analysis.
  * \ingroup genconstants
  *  
  */
template<typename floatT>
constexpr floatT seventeen_thirds()
{
   return static_cast<floatT>(17)/static_cast<floatT>(3);
}

} //namespace math
} //namespace mx

#endif //math_constants_hpp
