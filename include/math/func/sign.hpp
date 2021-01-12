/** \file sign.hpp
  * \brief Declares and defines the sign function.
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */


#ifndef mx_sign_hpp
#define mx_sign_hpp

#include <type_traits>

#ifdef MX_INCLUDE_BOOST
#include <boost/math/special_functions/sign.hpp>
#endif

namespace mx
{
namespace math
{
namespace func
{
   
/// The sign function.
/**
  * \ingroup functions
  */ 
template<typename T>
T sign( T x  /**< [in] the argument */)
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::sign<T>(x);
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "sign<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
float sign<float>( float x);

template<>
double sign<double>( double x);

template<>
long double sign<long double>( long double x);

}
}
}

#endif //mx_sign_hpp
