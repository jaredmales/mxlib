
#include "math/func/legendre.hpp"

#include <boost/math/special_functions/legendre.hpp>

namespace mx
{
namespace math
{
namespace func
{

template<>
float legendre_p<float>( int n, 
                         float x
                       )
{
   return boost::math::legendre_p<float>(n,x);
}

template<>
double legendre_p<double>( int n, 
                           double x
                         )
{
   return boost::math::legendre_p<double>(n,x);
}

template<>
long double legendre_p<long double>( int n, 
                                     long double x
                                   )
{
   return boost::math::legendre_p<long double>(n,x);
}

#ifdef HASQUAD
template<>
__float128 legendre_p<__float128>( int n, 
                                   __float128 x
                                 )
{
   return boost::math::legendre_p<__float128>(n,x);
}
#endif

} //namespace mx
} //namespace math
} //namespace func
