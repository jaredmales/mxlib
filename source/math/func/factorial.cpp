
#include "math/func/factorial.hpp"

#include <boost/math/special_functions/factorials.hpp>

namespace mx
{
namespace math
{
namespace func
{

template<>
float factorial<float>( float x )
{
   return boost::math::factorial<float>(x);
}

template<>
double factorial<double>( double x )
{
   return boost::math::factorial<double>(x);
}

template<>
long double factorial<long double>( long double x )
{
   return boost::math::factorial<long double>(x);
}


} //namespace mx
} //namespace math
} //namespace func
