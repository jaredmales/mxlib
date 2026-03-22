
#include "math/func/gamma.hpp"

#include <boost/math/special_functions/gamma.hpp>

namespace mx
{
namespace math
{
namespace func
{

#ifdef MXLIB_USE_LONGDOUBLE
template <>
float tgamma<float>( float x )
{
    return boost::math::tgamma<float>( x );
}

template <>
double tgamma<double>( double x )
{
    return boost::math::tgamma<double>( x );
}

template <>
long double tgamma<long double>( long double x )
{
    return boost::math::tgamma<long double>( x );
#endif
}

#ifdef MXLIB_HAS_QUAD
template <>
__float128 tgamma<__float128>( __float128 x )
{
    return boost::math::tgamma<__float128>( x );
}
#endif

} // namespace func
} // namespace math
} // namespace mx
