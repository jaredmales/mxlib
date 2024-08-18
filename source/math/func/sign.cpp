
#include "math/func/sign.hpp"

#include <boost/math/special_functions/sign.hpp>

namespace mx
{
namespace math
{
namespace func
{

template <>
float sign<float>( float x )
{
    return boost::math::sign<float>( x );
}

template <>
double sign<double>( double x )
{
    return boost::math::sign<double>( x );
}

template <>
long double sign<long double>( long double x )
{
    return boost::math::sign<long double>( x );
}

} // namespace func
} // namespace math
} // namespace mx
