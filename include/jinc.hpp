/** \file jinc.hpp
  * \brief Declares and defines the Jinc and Jinc2 functions
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __jinc_hpp__
#define __jinc_hpp__

#include <boost/math/special_functions/bessel.hpp>

namespace mx
{
   
///Returns the the Jinc function
/** The Jinc function is defined here as
  * \f[ 
  *     Ji(x) = \frac{J_1(x)}{x} 
  * \f]
  * where \f$ J_1 \f$ is the cylindrical bessel function of the first kind of order 1.
  * 
  * \param x is the argument
  * 
  * \returns the value of Ji(x)
  * 
  * \todo handle small arguments of Ji(x) with Taylor terms as is done by Boost for Sinc 
  * 
  * \ingroup gen_math
  */
template<typename floatT> 
floatT jinc(floatT x)
{
   if(fabs(x) < .00001) return 0.5;
   else return boost::math::cyl_bessel_j(1,x)/(x);
}

///Returns the the Jinc2 function
/** The Jinc2 function is defined here as
  * \f[ 
  *     Ji2(x) = \frac{J_2(x)}{x} 
  * \f]
  * where \f$ J_2 \f$ is the cylindrical bessel function of the first kind of order 2.
  * 
  * \param x is the argument
  * 
  * \returns the value of Ji2(x)
  * 
  * \todo handle small arguments of Ji2(x) with Taylor terms as is done by Boost for Sinc 
  * 
  * \ingroup gen_math
  */
template<typename floatT> 
floatT jinc2(floatT x)
{
   if(fabs(x) < 1e-5) return 0.0;
   else return boost::math::cyl_bessel_j(2,x)/(x);
}

} //namespace mx

#endif //__jinc_hpp__
