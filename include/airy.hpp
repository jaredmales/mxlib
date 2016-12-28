/** \file airy.hpp
  * \author Jared R. Males
  * \brief Utilities for working with the Airy pattern
  * \ingroup imaging_files
  *
  */

#ifndef __airy_hpp__
#define __airy_hpp__

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;


#include "jinc.hpp"

namespace mx
{
   
///The classical Airy pattern
/** Returns the intensity distribution of the Airy pattern at a given \f$ \lambda/D \f$
  *
  * \param [in] x is the separation in units of \f$ \lambda/D \f$.
  * 
  * \tparam floatT is the floating point type used for arithmetic.
  */
template<typename floatT>
floatT airy(floatT x)
{
   return pow(2*jinc(pi<floatT>()*x),2); 
}

///The centrally obscured Airy pattern
/** Returns the intensity distribution of the centrally obscured Airy pattern at a given \f$ \lambda/D \f$
  *
  * \param [in] x is the separation in units of  \f$ \lambda/D \f$.
  * \param [in] eps is the ratio of the circular central obscuration diameter to the diameter.
  * 
  * \tparam floatT is the floating point type used for arithmetic.
  */
template<typename floatT>
floatT airy(floatT x, floatT eps)
{
   return (1./pow(1.-eps*eps,2))*pow(2.*jinc(2*pi<floatT>()*x)-eps*2.*jinc(eps*pi<floatT>()*x),2);
}

///Seeing Halo
/** \todo document seeingHalo
  */
template<typename floatT>
floatT seeingHalo(floatT x, floatT fwhm)
{
   return (0.488/(fwhm*fwhm))*pow(1. + (11./6.)*pow(x/fwhm,2), -11./6.);
}

}//namespace mx

#endif //__airy_hpp__

