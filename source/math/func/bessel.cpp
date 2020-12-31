
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


} //namespace mx
} //namespace math
} //namespace func
