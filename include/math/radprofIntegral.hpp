/** \file radprofIntegral.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief 2D integration of a radial profile
  * \ingroup gen_math_files
  * 
  */

//***********************************************************************//
// Copyright 2022 Jared R. Males (jaredmales@gmail.com)
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
#ifndef radprofIntegral_hpp
#define radprofIntegral_hpp

#include <cmath>
#include <vector>
#include <type_traits>

#include "../mxException.hpp"

#include "constants.hpp"

namespace mx
{
namespace math 
{

/// Calculate the 2D integral given a radial profile.
/** Calculates the integral
  * \f[
    2\pi \int_{min}^{max} p(r)  r dr
    \f]
  * using the <a href="https://en.wikipedia.org/wiki/Trapezoidal_rule">Trapezoid rule</a>.
  *
  * The limits of integration are normally `r[0]` to `r.back()` inclusive.  If the parameter \p inczero
  * is set to true, then the region from `0` to `r[0]` is included, using the value of `p[0]`, where the 
  * the integral is there estimate as \f$ p[0] \pi r[0]^2 \f$.
  * 
  * \todo the trapezoid rule may not be strictly valid as implemented.  Is there a weighting to apply based on r?
  * 
  * \returns the value of the integral 
  * 
  * \throws mxException if the sizes of the input vectors don't match.
  * 
  * \ingroup integration
  */
template< typename vectorScT, typename vectorT>
typename vectorT::value_type radprofIntegral( vectorScT r, /// [in] the r 
                                              vectorT p,
                                              bool inczero = false
                                            )
{
   typedef typename vectorT::value_type floatT;

   static_assert( std::is_floating_point<typename std::remove_cv<floatT>::type>::value, "radprofIntegral: function must be floating point");
   
   if(r.size() != p.size())
   {
      mxThrowException(err::sizeerr, "radprofIntegral", "vectors must have same size");
   }

   if(r.size() < 2)
   {
      mxThrowException(err::sizeerr, "radprofIntegral", "must be at least 2 elements in radial profile");
   }

   floatT s = 0;
   
   if(inczero) s = p[0]*pow(r[0],2);

   size_t n = 0;
   s += 0.5*p[n]*(pow(r[n+1],2) - pow(r[n],2));//r[0] * ( r[1] - r[0] );

   for(n=1; n < r.size()-1; ++n)
   {
      s += p[n]*(pow(r[n],2) - pow(r[n-1],2));//r[n] * ( r[n] - r[n-1]);
   }

   s += 0.5*p[n]*(pow(r[n],2) - pow(r[n-1],2));

   //size_t n = r.size()-1;
   //s += 2*p[n]*r[n] * ( r[n] - r[n-1] );

   return s*pi<floatT>();
}

} //namespace math 
} //namespace mx

#endif //radprofIntegral_hpp
