


#ifndef mx_precision_hpp
#define mx_precision_hpp

#include <type_traits>

#ifdef MX_INCLUDE_BOOST
#include <boost/math/tools/precision.hpp>
#endif

namespace mx
{
namespace math
{
namespace func
{
   
/// Get the sqrt(epsilon) where epsilon is machine precision
/** Wrapper for boost function
  *
  * \todo this isn't, but should be constexpr -- need constr sqrt.
  */ 
template<typename T>
T root_epsilon() 
{
#ifdef MX_INCLUDE_BOOST
   return boost::math::tools::root_epsilon<T>();
#else
   static_assert(std::is_fundamental<T>::value || !std::is_fundamental<T>::value, "root_epsilon<T> not specialized for type T, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
#endif
}

template<>
float root_epsilon<float>();

template<>
double root_epsilon<double>();


}
}
}

#endif //mx_precision_hpp
