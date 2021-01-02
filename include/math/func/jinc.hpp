/** \file jinc.hpp
  * \brief Declares and defines the Jinc and Jinc2 functions
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef math_func_jinc_hpp
#define math_func_jinc_hpp

#include <limits>

#include "bessel.hpp"
#include "precision.hpp"

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
  * is smaller than \f$ \epsilon \f$, then it returns 1/2.  If x is larger than \f$ \epsilon \f$ but smaller than \f$ \sqrt{\epsilon} \f$, then
  * this function returns
  * \f[
  *    Ji(x) \approx \frac{1}{2} - \frac{x^2}{16}.
  * \f]
  * 
  * \returns the value of Ji(x)
  * 
  * \tparam T is an floating point type 
  * 
  * \ingroup functions
  */
template<typename T> 
T jinc( const T & x /**< [in] the argument */)
{
   T const taylor_0_bound = std::numeric_limits<T>::epsilon();
   T const taylor_2_bound = root_epsilon<T>();
   
   if( fabs(x) > taylor_2_bound)
   {
      return bessel_j<int, T>(1,x)/(x);
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
  * \returns the value of Ji2(x)
  * 
  * \ingroup functions
  */
template<typename T> 
T jinc2(const T & x /**< [in] the argument */)
{
   T const taylor_2_bound = root_epsilon<T>();
      
   if( fabs(x) > taylor_2_bound )
   {
      return bessel_j<int, T>(2,x)/(x);
   }
   else
   {
      return static_cast<T>(0);
   }
}

} //namespace func 
} //namespace math
} //namespace mx

#endif //math_func_jinc_hpp
