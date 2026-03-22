
#include "math/func/factorial.hpp"

#include <boost/math/special_functions/factorials.hpp>

#ifdef MXLIB_HAS_QUAD
#include <quadmath.h>
#endif

namespace mx
{
namespace math
{
namespace func
{

#ifdef MXLIB_USE_LONGDOUBLE
template <>
float factorial<float>( float x )
{
    return boost::math::factorial<float>( x );
}

template <>
double factorial<double>( double x )
{
    return boost::math::factorial<double>( x );
}

template <>
long double factorial<long double>( long double x )
{
    return boost::math::factorial<long double>( x );
#endif
}

#ifdef MXLIB_HAS_QUAD
template <>
__float128 factorial<__float128>( __float128 x )
{
    return floorq( tgammaq( x + static_cast<__float128>( 1 ) ) + static_cast<__float128>( 0.5 ) );
}
#endif

} // namespace func
} // namespace math
} // namespace mx
