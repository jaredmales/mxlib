/** \file sigmoid.hpp
  * \brief Tools for working with sigmoid curves
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup gen_math_files
  *
  */

#ifndef __sigmoid_hpp__
#define __sigmoid_hpp__

namespace mx
{

/** \ingroup gen_math
  * @{
  */

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
  */
template<typename floatT>
floatT logistic(floatT t, floatT t0 = 0, floatT a = 1)
{
   return 1.0/(1.0 + exp( -a*(t-t0)));
}


/// @}

} //namespace mx

#endif //__sigmoid_hpp__

