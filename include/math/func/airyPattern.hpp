/** \file airyPattern.hpp
  * \author Jared R. Males
  * \brief Utilities related to the Airy pattern point spread function.
  * \ingroup gen_math_files
  *
  */

#ifndef airyPattern_hpp
#define airyPattern_hpp

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;


#include "jinc.hpp"

namespace mx
{
   
namespace math
{

namespace func
{
   
///The classical Airy pattern
/** Returns the intensity distribution of the Airy pattern at a given \f$ \lambda/D \f$
  * 
  * References: \cite born_and_wolf, \cite mahajan_1986, https://en.wikipedia.org/wiki/Airy_disk.
  * 
  * \tparam realT is the floating point type used for arithmetic.
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */
template<typename realT>
realT airyPattern( realT x /**< [in] is the separation in units of \f$ \lambda/D \f$. */)
{
   return pow(2*jinc(pi<realT>()*x),2); 
}

///The centrally obscured Airy pattern
/** Returns the intensity distribution of the centrally obscured Airy pattern at a given \f$ \lambda/D \f$
  * 
  * References: \cite born_and_wolf, \cite mahajan_1986, https://en.wikipedia.org/wiki/Airy_disk.
  * 
  * \tparam realT is the floating point type used for arithmetic.
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */
template<typename realT>
realT airyPattern( realT x,  ///< [in] is the separation in units of  \f$ \lambda/D \f$.
                   realT eps ///< [in] is the ratio of the circular central obscuration diameter to the diameter.
                 )
{
   return (1./pow(1.-eps*eps,2))*pow(2.*jinc( pi<realT>()*x)-eps*eps*2.*jinc(eps*pi<realT>()*x),2);
}


///The general centrally obscured Airy pattern, with arbitrary center and platescale.
/** Returns the intensity distribution of the centrally obscured Airy pattern at a given \f$ \lambda/D \f$.  This version allows
  * for an arbitrary center, scaling, and platescale.
  * 
  * References: \cite born_and_wolf, \cite mahajan_1986, https://en.wikipedia.org/wiki/Airy_disk.
  * 
  * \tparam realT is the floating point type used for arithmetic.
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */
template<typename realT>
realT airyPattern( realT x,  ///< [in] is the x-coordinate in units of pixels
                   realT y,  ///< [in] is the y-coordinate in units of pixels
                   realT A0, ///< [in] constant value added to the Airy pattern
                   realT A, ///< [in] peak scale of the Airy pattern.
                   realT x0, ///< [in] is the x-center in units of pixels
                   realT y0, ///< [in] is the y-center in units of pixels
                   realT ps, ///< [in] the platescale in \f$ (\lambda/D)/pixel  \f$
                   realT eps ///< [in] is the ratio of the circular central obscuration diameter to the diameter.
                 )
{
   realT r = sqrt( pow(x-x0,2) + pow(y-y0,2)) * ps;
   
   return A0 + A*airyPattern(r, eps);
}


///Seeing Halo Profile
/** A Moffat profile due to Roddier (1981)\cite roddier_1981, see also Racine et al (1999)\cite racine_1999, which
  * can be used to describe the seeing limited point-spread function of a telescope imaging through turbulence, or the halo in partially corrected imaging. 
  *  
  * \returns the value of the profile at x, in units of fractional flux per unit area.
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */
template<typename realT>
realT seeingHalo( realT x, ///< [in] the separation in the same units as fwhm.
                   realT fwhm ///< [in] the fwhm in arbitrary units.  Note that this defines the area units of the density.
                 )
{
   return (0.488/(fwhm*fwhm))*pow(1. + (11./6.)*pow(x/fwhm,2), -11./6.);
}

} //namespace func
}//namespace math
}//namespace mx

#endif //airyPattern_hpp

