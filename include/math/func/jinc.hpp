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
#include <cmath>

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
   
   if( std::fabs(x) > taylor_2_bound)
   {
      return bessel_j<int, T>(1,x)/(x);
   }
   else
   {
      // approximation by taylor series in x at 0
      T  result = static_cast<T>(0.5);
      if (std::fabs(x) >= taylor_0_bound)
      {
         T x2 = x*x;

         // approximation by taylor series in x at 0 up to order 2
         result -= x2/static_cast<T>(16);
      }
      
      return result;
   }
            
}

extern template
float jinc<float>(const float & x);

extern template
double jinc<double>(const double & x);

extern template
long double jinc<long double>(const long double & x);

#ifdef HASQUAD
extern template
__float128 jinc<__float128>(const __float128 & x);
#endif

///Returns the JincN function
/** The JincN function is defined here as
  * \f[ 
  *     Ji_N(x) = \frac{J_N(x)}{x} 
  * \f]
  * where \f$ J_N \f$ is the cylindrical bessel function of the first kind of order N, \f$ N \ge 1 \f$.
  * 
  * If \f$ N == 1 \f$ this returns \ref mx::math::func::jinc() "jinc(x)".
  * 
  * Otherwise, if x is smaller than \f$ \sqrt{\epsilon} \f$, returns 0.
  * 
  * \returns the value of JiN(x)
  * 
  * \ingroup functions
  */
template<typename T1, typename T2> 
T2 jincN( const T1 & v, ///< [in] the Bessel function order
          const T2 & x  ///< [in] the argument 
        )
{
   if(v == 1) return jinc(x);

   T2 const taylor_2_bound = root_epsilon<T2>();
      
   if( std::fabs(x) > taylor_2_bound )
   {
      return bessel_j<T1, T2>(v,x)/(x);
   }
   else
   {
      return static_cast<T2>(0);
   }
}

extern template
float jincN<float, float>( const float & v,
                           const float & x
                         );

extern template
float jincN<int, float>( const int & v,
                         const float & x
                       );

extern template
double jincN<double, double>( const double & v,
                              const double & x
                            );

extern template
double jincN<int, double>( const int & v,
                           const double & x
                         );


extern template
long double jincN<long double, long double>( const long double & v,
                                             const long double & x
                                           );

extern template
long double jincN<int, long double>( const int & v,
                                     const long double & x
                                   );

#ifdef HASQUAD
extern template

__float128 jincN<__float128, __float128>( const __float128 & v,
                                          const __float128 & x
                                        );

__float128 jincN<int, __float128>( const int & v,
                                   const __float128 & x
                                 );
#endif

} //namespace func 
} //namespace math
} //namespace mx

#endif //math_func_jinc_hpp
