/** \file aoConstants.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Calculate and provide constants related to adaptive optics.
  * \ingroup mxAO_files
  * 
  */

#ifndef mx_AO_analysis_aoConstants_hpp
#define mx_AO_analysis_aoConstants_hpp

#include <cmath>


#include "../../math/constants.hpp"
#include "../../math/func/gamma.hpp"



namespace mx
{
namespace AO
{
namespace constants
{
   

///Calculate the AO constants
template<typename floatT>
void calcConstants(floatT & a_sf, floatT & a_psd)
{
   floatT pi = math::pi<floatT>();

   floatT gam1 = math::func::tgamma<floatT>( static_cast<floatT>(6)/ static_cast<floatT>(5) );

   floatT gam2 = math::func::tgamma<floatT>( static_cast<floatT>(11)/ static_cast<floatT>(6) );

   a_sf = pow( static_cast<floatT>(24)/static_cast<floatT>(5) * gam1,  static_cast<floatT>(5) / static_cast<floatT>(6) );
   
   a_psd = a_sf*gam1*gam1*sin(  static_cast<floatT>(5) * pi /  static_cast<floatT>(6) )/ pow(pi, static_cast<floatT>(11)/static_cast<floatT>(3) );
   
   a_sf = a_sf*static_cast<floatT>(2);
   
}

///The scaling constant for the Kolmorogov optical phase structure function 
/** The optical phase structure function for Kolmogorov turbulence is
  * \f[
  * \mathcal{D}_\phi(r) = a_{SF} \left(\frac{r}{r_0}\right)^{5/3}
  * \f]
  * 
  * This constant, given to 50 digits, was caculated using calcConstants with boost::multiprecision::cpp_dec_float_100 (100 digits of precision)
  * 
  * \returns 6.883877.... cast to type floatT
  * 
  * \tparam floatT is the type to cast the value to.
  * \ingroup aoConstants
  * 
  */
template<typename floatT>
constexpr floatT a_SF()
{
   return static_cast<floatT>(6.8838771822938116152935575630969803178936813057678);
}

///The scaling constant for the Kolmorogov optical phase power spectral density 
/** The optical phase PSD for Kolmogorov turbulence is
  * \f[
  * \mathcal{P}_\phi(k) = a_{PSD} \left(\frac{1}{r_0^{5/3} k^{11/3}}\right)
  * \f]
  * 
  * This constant, given to 50 digits, was caculated using calcConstants with boost::multiprecision::cpp_dec_float_100 (100 digits of precision)
  * 
  * \returns 0.02181... cast to type floatT
  * 
  * \tparam floatT is the type to cast the value to
  * 
  * \ingroup aoConstants
  */
template<typename floatT>
constexpr floatT a_PSD()
{
   return static_cast<floatT>(0.0218139977034218241674821945866523430205216037234);
}



} //namespace constants
} //namespace AO
} //namespace mx

#endif //mx_AO_analysis_aoConstants_hpp
