
#include "math/func/bessel.hpp"

#include <boost/math/special_functions/bessel.hpp>

namespace mx
{
namespace math
{
namespace func
{

template<>
float bessel_j<float, float>( float v, 
                              float x
                            )
{
   return boost::math::cyl_bessel_j<float, float>(v,x);
}

template<>
float bessel_j<int, float>( int v,
                            float x
                          )
{
   return boost::math::cyl_bessel_j<int, float>(v,x);
}

template<>
double bessel_j<double, double>( double v, 
                                 double x
                               )
{
   return boost::math::cyl_bessel_j<double, double>(v,x);
}

template<>
double bessel_j<int, double>( int v,
                              double x
                            )
{
   return boost::math::cyl_bessel_j<int, double>(v,x);
}

template<>
long double bessel_j<long double, long double>( long double v, 
                                                long double x
                                              )
{
   return boost::math::cyl_bessel_j<long double, long double>(v,x);
}

template<>
long double bessel_j<int, long double>( int v,
                                        long double x
                                      )
{
   return boost::math::cyl_bessel_j<int, long double>(v,x);
}

#ifdef HASQUAD
template<>
__float128 bessel_j<__float128, __float128>( __float128 v, 
                                             __float128 x
                                           )
{
   return boost::math::cyl_bessel_j<__float128, __float128>(v,x);
}

template<>
__float128 bessel_j<int, __float128>( int v,
                                      __float128 x
                                    )
{
   return boost::math::cyl_bessel_j<int, __float128>(v,x);
}
#endif

} //namespace mx
} //namespace math
} //namespace func
