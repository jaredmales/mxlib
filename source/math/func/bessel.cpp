
#include "math/func/bessel.hpp"

#include <boost/math/special_functions/bessel.hpp>

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
float bessel_j<float, float>( float v, float x )
{
    return boost::math::cyl_bessel_j<float, float>( v, x );
}

template <>
float bessel_j<int, float>( int v, float x )
{
    return boost::math::cyl_bessel_j<int, float>( v, x );
}

template <>
double bessel_j<double, double>( double v, double x )
{
    return boost::math::cyl_bessel_j<double, double>( v, x );
}

template <>
double bessel_j<int, double>( int v, double x )
{
    return boost::math::cyl_bessel_j<int, double>( v, x );
}

template <>
long double bessel_j<long double, long double>( long double v, long double x )
{
    return boost::math::cyl_bessel_j<long double, long double>( v, x );
#endif
}

#ifdef MXLIB_USE_LONGDOUBLE
template <>
long double bessel_j<int, long double>( int v, long double x )
{
    return boost::math::cyl_bessel_j<int, long double>( v, x );
#endif
}

#ifdef MXLIB_HAS_QUAD
template <>
__float128 bessel_j<__float128, __float128>( __float128 v, __float128 x )
{
    if( floorq( v ) == v )
    {
        return bessel_j<int, __float128>( static_cast<int>( v ), x );
    }

    return static_cast<__float128>(
        boost::math::cyl_bessel_j<long double, long double>( static_cast<long double>( v ), static_cast<long double>( x ) ) );
}

template <>
__float128 bessel_j<int, __float128>( int v, __float128 x )
{
    if( v == 0 )
    {
        return j0q( x );
    }

    if( v == 1 )
    {
        return j1q( x );
    }

    return jnq( v, x );
}
#endif

} // namespace func
} // namespace math
} // namespace mx
