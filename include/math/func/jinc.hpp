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
   
namespace math
{

namespace func 
{
   
///Returns the the Jinc function
/** The Jinc function is defined here as
  * \f[ 
  *     Ji(x) = \frac{J_1(x)}{x} 
  * \f]
  * where \f$ J_1 \f$ is the cylindrical bessel function of the first kind of order 1.
  * 
  * Follows the technique in boost sinc_pi, using the Taylor series for small arguments. If x
  * is smaller than \f$ \epsilon \f$, then it returns 1/2.  If x is larger than \f$ \epsilon \f$ but smaler than \f$ \sqrt{\epsilon} \f$, then
  * this function returns
  * \f[
  *    Ji(x) \approx \frac{1}{2} - \frac{x^2}{16}.
  * \f]
  * 
  * \param x is the argument
  * 
  * \returns the value of Ji(x)
  * 
  * \tparam T is an floating point type 
  * 
  * \ingroup functions
  */
template<typename T> 
T jinc( const T & x)
{
   T const taylor_0_bound = boost::math::tools::epsilon<T>();
   T const taylor_2_bound = boost::math::tools::root_epsilon<T>();
   
   if( abs(x) > taylor_2_bound)
   {
      return boost::math::cyl_bessel_j<T>(1,x)/(x);
   }
   else
   {
      // approximation by taylor series in x at 0
      T  result = static_cast<T>(0.5);
      if (abs(x) >= taylor_0_bound)
      {
         T x2 = x*x;

         // approximation by taylor series in x at 0 up to order 2
         result -= x2/static_cast<T>(16);
      }
      
      return result;
   }
            
}

///Returns the Jinc2 function
/** The Jinc2 function is defined here as
  * \f[ 
  *     Ji_2(x) = \frac{J_2(x)}{x} 
  * \f]
  * where \f$ J_2 \f$ is the cylindrical bessel function of the first kind of order 2.
  * 
  * If x is smaller than \f$ \sqrt{\epsilon} \f$, returns 0.
  * 
  * \param x is the argument
  * 
  * \returns the value of Ji2(x)
  * 
  * \ingroup functions
  */
template<typename T> 
T jinc2(const T & x)
{
   T const taylor_2_bound = boost::math::tools::root_epsilon<T>();
      
   if(abs(x) > taylor_2_bound)
   {
      return boost::math::cyl_bessel_j(2,x)/(x);
   }
   else
   {
      return static_cast<T>(0);
   }
}

} //namespace func 
} //namespace math
} //namespace mx

#endif //__jinc_hpp__
