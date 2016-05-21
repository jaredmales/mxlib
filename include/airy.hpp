/** \file airy.hpp
  * \author Jared R. Males
  * \brief Utilities for working with the Airy pattern
  * \ingroup imaging
  *
  */

#ifndef __airy_hpp__
#define __airy_hpp__

#include <boost/math/special_functions/bessel.hpp> 

using namespace boost::math;
using namespace boost::math::constants;

namespace mx
{
   
///The classical Airy pattern
/** Returns the intensity distribution of the Airy pattern at a given \f$ \lambda/D \f$
  *
  * \param x [in] is the separation in units of  \f$ \lambda/D \f$.
  * 
  * \tparam floatT is the floating point type used for arithmetic.
  */
template<typename floatT>
floatT airy(floatT x)
{
   if(x == 0) return 1.0;
   
   return pow(2.*cyl_bessel_j(1, pi<floatT>()*x)/(pi<floatT>()*x),2);
}

///The centrally obscured Airy pattern
/** Returns the intensity distribution of the centrally obscured Airy pattern at a given \f$ \lambda/D \f$
  *
  * \param x [in] is the separation in units of  \f$ \lambda/D \f$.
  * \param eps [in] is the ratio of the circular central obscuration diameter to the diameter.
  * 
  * \tparam floatT is the floating point type used for arithmetic.
  */
template<typename floatT>
floatT airy(floatT x, floatT eps)
{
   return (1./pow(1.-eps*eps,2))*pow(2.*cyl_bessel_j(1, pi<floatT>()*x)/(pi<floatT>()*x)-eps*2.*cyl_besel_j(1, eps*pi<floatT>()*x)/(pi<floatT>()*x),2);

}

template<typename floatT>
floatT seeingHalo(floatT x, floatT fwhm)
{
   return (0.488/(fwhm*fwhm))*pow(1. + (11./6.)*pow(x/fwhm,2), -11./6.);
}

}//namespace mx

#endif //__airy_hpp__

