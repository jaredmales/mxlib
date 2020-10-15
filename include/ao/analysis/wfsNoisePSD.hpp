/** \file wfsNoisePSD.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines a function to calculate the measurement noise PSD
  * \ingroup mxAO_files
  * 
  */

#ifndef __wfsNoisePSD_hpp__
#define __wfsNoisePSD_hpp__

namespace mx
{
namespace AO
{
namespace analysis
{
   
///Populate a vector with the PSD of measurement noise given WFS parameters.
/** 
  * \tparam realT is the real floating point type for calculations
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename realT>
void wfsNoisePSD( std::vector<realT> &PSD,  ///< [out] A pre-allocated vector which will be filled with the PSD value.
                  realT beta_p_k,           ///< [in] The WFS \f$ \beta_p \f$ parameter (see Guyon, 2005 \cite{guyon_2005}).
                  realT Fg,                 ///< [in] Total photons/sec used by the WFS
                  realT tau,                ///< [in] WFS Integration time [sec]
                  realT npx,                ///< [in] Number of pixesl used by the WFS
                  realT Fb,                 ///< [in] Background flux, in photons/sec/pixel.
                  realT ron                 ///< [in] Readout noise, in photons/pixel/read.
                )
{

   realT snr2 = pow(Fg*tau,2) / (Fg*tau + npx*Fb*tau + npx*ron*ron);
   
   realT psd = 2*pow(beta_p_k,2)/snr2 * tau;
      
   for(size_t i=0; i< PSD.size(); ++i)
   {
      PSD[i] = psd;
   }
}
 
} //namespace analysis
} //namespace AO
} //namespace mx

#endif //__wfsNoisePSD_hpp__
   
