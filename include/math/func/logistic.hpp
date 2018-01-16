/** \file logistic.hpp
  * \brief The logistic function
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup gen_math_files
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

#ifndef __logistic_hpp__
#define __logistic_hpp__

namespace mx
{

namespace math 
{
   
namespace func 
{
   


///Return the logistic function parameter for a specified rise time
/** The logistic function is
  * \f[
  *  f(t) = \frac{1}{1 + e^{-a(t-t_0)}}
  * \f]
  * The parameter \f$ a \f$ controls how fast the function rises.  Here it is specified by the
  * value \f$ f(t_{1/2}) = x\f$, where \f$ 0 < x < 1 \f$.
  * 
  * \param x [input] the value at which the rist time is specified.
  * \param thalf [input] half the rise time, or the time after 0 when f(t) = x.
  * 
  * \returns the value of a
  * 
  * \tparam floatT is the floating point type of the arguments and the returned value.
  * 
  * \ingroup functions
  */
template<typename floatT>
floatT logistic_param(floatT x, floatT thalf)
{
   return -log( 1.0/x - 1)/thalf;
}

///Return the value of the logistic function
/** The logistic function is
  * \f[
  *  f(t) = \frac{1}{1 + e^{-a(t-t_0)}}
  * \f]
  * 
  * \param t [input] the argument.
  * \param t0 [input] [optional] the center of the curve, defaults to 0.
  * \param a [input] [optional] the exponent parameter, defaults to 1.
  * 
  * \returns the value of the logistic function at t
  * 
  * \tparam floatT is the floating point type of the arguments and the returned value.
  * 
  * \ingroup functions
  */
template<typename floatT>
floatT logistic(floatT t, floatT t0 = 0, floatT a = 1)
{
   return 1.0/(1.0 + exp( -a*(t-t0)));
}



} //namespace func 
} //namespace math 
} //namespace mx

#endif //__logistic_hpp__

