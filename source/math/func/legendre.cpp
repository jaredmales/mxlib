
#include "math/func/legendre.hpp"

#include <boost/math/special_functions/legendre.hpp>

namespace mx
{
namespace math
{
namespace func
{

#ifdef MXLIB_USE_LONGDOUBLE
template <>
float legendre_p<float>( int n, float x )
{
    return boost::math::legendre_p<float>( n, x );
}

template <>
double legendre_p<double>( int n, double x )
{
    return boost::math::legendre_p<double>( n, x );
}

template <>
long double legendre_p<long double>( int n, long double x )
{
    return boost::math::legendre_p<long double>( n, x );
#endif
}

#ifdef MXLIB_HAS_QUAD
template <>
__float128 legendre_p<__float128>( int n, __float128 x )
{
    if( n == 0 )
    {
        return static_cast<__float128>( 1 );
    }

    if( n == 1 )
    {
        return x;
    }

    __float128 pnm2 = static_cast<__float128>( 1 );
    __float128 pnm1 = x;

    for( int i = 2; i <= n; ++i )
    {
        __float128 pn = ( static_cast<__float128>( 2 * i - 1 ) * x * pnm1 - static_cast<__float128>( i - 1 ) * pnm2 ) /
                        static_cast<__float128>( i );
        pnm2 = pnm1;
        pnm1 = pn;
    }

    return pnm1;
}
#endif

} // namespace func
} // namespace math
} // namespace mx
