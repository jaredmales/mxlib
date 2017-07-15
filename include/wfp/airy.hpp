/** \file airy.hpp
  * \author Jared R. Males
  * \brief Utilities related to the Airy pattern point spread function.
  * \ingroup imaging_files
  *
  */

#ifndef __airy_hpp__
#define __airy_hpp__

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;


#include "../math/func/jinc.hpp"

namespace mx
{
   
namespace wfp
{
   
///The classical Airy pattern
/** Returns the intensity distribution of the Airy pattern at a given \f$ \lambda/D \f$
  *
  * \tparam floatT is the floating point type used for arithmetic.
  * 
  * \ingroup psfs
  */
template<typename floatT>
floatT airy( floatT x /**< [in] is the separation in units of \f$ \lambda/D \f$. */)
{
   return pow(2*math::func::jinc(pi<floatT>()*x),2); 
}

///The centrally obscured Airy pattern
/** Returns the intensity distribution of the centrally obscured Airy pattern at a given \f$ \lambda/D \f$
  * 
  * \tparam floatT is the floating point type used for arithmetic.
  * 
  * \ingroup psfs
  */
template<typename floatT>
floatT airy( floatT x,  ///< [in] is the separation in units of  \f$ \lambda/D \f$.
             floatT eps ///< [in] is the ratio of the circular central obscuration diameter to the diameter.
           )
{
   return (1./pow(1.-eps*eps,2))*pow(2.*math::func::jinc(2*pi<floatT>()*x)-eps*2.*jinc(eps*pi<floatT>()*x),2);
}

///Seeing Halo Profile
/** A Moffat profile due to Roddier (1981)\cite roddier_1981, see also Racine et al (1999)\cite racine_1999, which
  * can be used to describe the seeing limited point-spread function of a telescope imaging through turbulence, or the halo in partially corrected imaging. 
  *  
  * \returns the value of the profile at x, in units of fractional flux per unit area.
  * 
  * \ingroup psfs
  */
template<typename floatT>
floatT seeingHalo( floatT x, ///< [in] the separation in the same units as fwhm.
                   floatT fwhm ///< [in] the fwhm in arbitrary units.  Note that this defines the area units of the density.
                 )
{
   return (0.488/(fwhm*fwhm))*pow(1. + (11./6.)*pow(x/fwhm,2), -11./6.);
}

}//namespace wfp
}//namespace mx

#endif //__airy_hpp__

