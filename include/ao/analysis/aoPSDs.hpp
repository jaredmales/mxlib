/** \file aoPSDs.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Spatial power spectra used in adaptive optics.
  * \ingroup mxAO_files
  * 
  */

//***********************************************************************//
// Copyright 2016-2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef aoPSDs_hpp
#define aoPSDs_hpp

#include <string>

#include "../../mxError.hpp"
#include "../../math/constants.hpp"
#include "../../math/func/jinc.hpp"

#include "aoConstants.hpp"

#include "aoAtmosphere.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{
   
namespace PSDComponent
{
   ///Enum to specify which component of the PSD to calcualte
   enum { phase,    ///< The phase or OPD
          amplitude,  ///< The amplitude
          dispPhase,  ///< The phase component of dispersive anisoplanatism
          dispAmplitude  ///< The amplitude component of dispersive anisoplanatism
        };
        
   std::string compName( int cc );
}
   
///Manage calculations using the von Karman spatial power spectrum.
/** This very general PSD has the form
  * 
  * \f[ 
    \mathcal{P}(k) = \frac{\beta}{ ( k^2 + k_0^2) ^{\alpha/2}} 
    \f]
  *  
  * Where \f$ k \f$ is spatial frequency, \f$ \beta \f$ is a normalization constant, and outer scale \f$ L_0 \f$ is included via \f$ k_0 = 1/L_0 \f$.
  * For \f$ L_0 \rightarrow \infty \f$ this becomes a simple power law with index \f$ \alpha \f$.
  * 
  * For atmospheric turbulence \f$ \alpha = 11/3 \f$ and the normalization is
  \f[
  \beta = \frac{\mathcal{A}_p}{r_0^{5/3}}
  \f]
  * where \f$ A_p  = 0.0218...\f$  (see \ref mx::AO::constants::a_PSD ).
  *
  * \ingroup mxAOAnalytic
  */ 
template<typename realT>
struct vonKarmanSpectrum
{
   
protected:
   bool m_subPiston {true};///< flag controlling whether piston is subtracted from the PSD.  Default is true.
   bool m_subTipTilt {false}; ///< flag controlling whether tip and tilt are subtracted from the PSD.  Default is false.
   
   bool m_scintillation {false}; ///< flag controlling whether or not scintillation is included
   int m_component {PSDComponent::phase}; ///< If m_scintillation is true, this controls whether phase (0), amplitude (1), or dispersive contrast (2) is returned.
   
   realT m_D {1.0}; ///< Diameter used for piston and tip/tilt subtraction, in m. Default is 1 m.

   const char * m_id = "von Karman";
   
   //realT m_alpha {eleven_thirds<realT>()}; ///< The power-law index, 11/3 for Kolmogorov
   
public:   
   
   ///Default Constructor
   vonKarmanSpectrum();
   
   ///Constructor specifying the parameters
   /**
     */ 
   vonKarmanSpectrum( bool subP, ///< [in] is the value of m_subPiston.
                      bool subT, ///< [in] is the value of m_subTipTilt.
                      realT D ///< [in] is the value of m_D.
                    );
   
   ///Get the value of m_subPiston
   /**
     * \returns m_subPiston
     */ 
   bool subPiston();
   
   ///Set the value of m_subPiston
   /**
     */ 
   void subPiston( bool sp /**< [in] is the new value of m_subPiston */);
   
   ///Get the value of m_subTipTilt
   /**
     * \returns the value of m_subTipTilt.
     */ 
   bool subTipTilt();
   
   ///Set the value of the m_subTipTilt flag.
   /**
     */ 
   void subTipTilt(bool st /**< [in] the new value of m_subTipTilt.*/);
   
   ///Get the value of m_scintillation
   /**
     * \returns m_scintillation
     */ 
   bool scintillation();
   
   ///Set the value of m_scintillation
   /**
     */ 
   void scintillation(bool sc /**< [in] the new value of m_scintillation*/);
   
   ///Get the value of m_component
   /**
     * \returns m_component
     */ 
   int component();
   
   ///Set the value of m_component
   /**
     */ 
   int component( int cc /**< [in] the new value of m_component */);
   
   ///Get the value of the diameter m_D.
   /**
     * \returns the current value of m_D, the diameter in m.
     */ 
   realT D();
   
   ///Set the aperture diameter
   /**
     */
   void D(realT nd /**< [in] the new diameter in m */);
   
   ///Get the value of the PSD index alpha.
   /**
     * \returns the current value of m_alpha, the PSD index.
     */
  // realT alpha();
   
   ///Set the PSD index alpha
   /**
     */
   void alpha(realT na /**< [in] the new PSD index*/ );
   
   ///Get the value of the PSD at spatial frequency k and a zenith distance.
   /**
     * \returns the von Karman PSD at the specified spatial frequency.
     * \returns -1 if an error occurs.
     * 
     */ 
   template< class psdParamsT >
   realT operator()( psdParamsT & par, ///< [in] gives the PSD parameters.
                     realT k,          ///< [in] is the spatial frequency in m^-1.
                     realT sec_zeta    ///< [in] is the secant of the zenith distance.
                   );
   
   /// Get the value of the PSD at spatial frequency k and wavelength lambda, and a zenith distance, with a WFS at a different wavelength
   /**
     * 
     * \returns the von Karman PSD at the specified spatial frequency for the specified wavelength.
     * \returns -1 if an error occurs.
     */    
   template< class psdParamsT >
   realT operator()( psdParamsT & atm, ///< [in] gives the PSD parameters.
                     realT k,          ///< [in] is the spatial frequency in m^-1.
                     realT lambda,     ///< [in] is the observation wavelength in m
                     realT lambda_wfs, ///< [in] is the wavefront measurement wavelength in m
                     realT secZeta     ///< [in] is the secant of the zenith distance
                   );
   
   /// Get the fitting error for an actuator spacing d.
   /**
     * \todo generatlize for different alpha and beta
     * 
     * \param atm gives the atmosphere parameters r_0 and L_0
     * \param d is the actuator spacing in m
     */    
   realT fittingError(aoAtmosphere<realT> &atm, realT d);

   template<typename iosT>
   iosT & dumpPSD(iosT & ios);

};

template< typename realT>
vonKarmanSpectrum<realT>::vonKarmanSpectrum()
{
}

template< typename realT>
vonKarmanSpectrum<realT>::vonKarmanSpectrum( bool subP, // [in] is the value of m_subPiston.
                   bool subT, // [in] is the value of m_subTipTilt.
                   realT D // [in] is the value of m_D.
                 )
{
   m_subPiston = subP;
   m_subTipTilt = subT;
   m_D = D;
}

template< typename realT>
bool vonKarmanSpectrum<realT>::subPiston()
{
   return m_subPiston;
}


template< typename realT>
void vonKarmanSpectrum<realT>::subPiston( bool sp /* [in] is the new value of m_subPiston */)
{
   m_subPiston = sp;
}

template< typename realT>
bool vonKarmanSpectrum<realT>::subTipTilt()
{
   return m_subTipTilt;
}

template< typename realT>
void vonKarmanSpectrum<realT>::subTipTilt(bool st /* [in] the new value of m_subTipTilt */)
{
   m_subTipTilt = st;
}

template< typename realT>
bool vonKarmanSpectrum<realT>::scintillation()
{
   return m_scintillation;
}

template< typename realT>
void vonKarmanSpectrum<realT>::scintillation(bool sc /* [in] the new value of m_scintillation */)
{
   m_scintillation = sc;
}

template< typename realT>
int vonKarmanSpectrum<realT>::component()
{
   return m_component;
}

template< typename realT>
int vonKarmanSpectrum<realT>::component( int cc /* [in] the new value of m_component */)
{
   if( cc != PSDComponent::phase && cc != PSDComponent::amplitude && cc != PSDComponent::dispPhase && cc != PSDComponent::dispAmplitude )
   {
      mxError("vonKarmanSpectrum::component", MXE_INVALIDARG, "Unknown component");
      return -1;
   }
   
   m_component = cc;
   return 0;         
}

template< typename realT>
realT vonKarmanSpectrum<realT>::D()
{
   return m_D;
}

template< typename realT>
void vonKarmanSpectrum<realT>::D(realT nd /**< [in] the new diameter in m */)
{
   m_D = nd;
}

template< typename realT>
template< class psdParamsT >
realT vonKarmanSpectrum<realT>::operator()( psdParamsT & par, // [in] gives the PSD parameters.
                                            realT k, // [in] is the spatial frequency in m^-1.
                                            realT sec_zeta // [in] is the secant of the zenith distance.
                                          )
{
   realT k02;
   
   ///\todo this needs to handle layers with different L_0
   if(par.L_0(0) > 0)
   {
      k02 = (1)/(par.L_0(0)*par.L_0(0));
   }
   else k02 = 0;

   if(k02 == 0 && k == 0)
   {
      return 0;
   }
   
   realT Ppiston, Ptiptilt;

   if( (m_subPiston || m_subTipTilt) )
   {
      if (m_D == 0)
      {
         mxError("aoAtmosphere", MXE_PARAMNOTSET, "Diameter D not set for Piston and/or TT subtraction.");
         return -1;
      }
      if(m_subPiston)
      {
         Ppiston = pow(2*math::func::jinc(math::pi<realT>()*k*m_D), 2);
      }
      else Ppiston = 0;

      if(m_subTipTilt)
      {
         Ptiptilt = pow(4*math::func::jincN(2, math::pi<realT>()*k*m_D), 2);
      }
      else Ptiptilt = 0;
   }
   else
   {
      Ppiston = 0;
      Ptiptilt = 0;
   }
   
   return par.beta()*pow(k*k+k02, -1*par.alpha()/2) * (1.0-Ppiston - Ptiptilt)*sec_zeta;
}

template< typename realT>
template< class psdParamsT >
realT vonKarmanSpectrum<realT>::operator()( psdParamsT & par, // [in] gives the PSD parameters.
                                            realT k, // [in] is the spatial frequency in m^-1.
                                            realT lambda, // [in] is the observation wavelength in m
                                            realT lambda_wfs, // [in] is the wavefront measurement wavelength in m
                                            realT secZeta // [in] is the secant of the zenith distance
                                          )
{
   realT psd = operator()(par, k, secZeta)* pow( par.lam_0()/lambda, 2);
   
   if(psd < 0) return -1;
   
   if(m_scintillation)
   {
      if(m_component == PSDComponent::phase)
      {
         psd *= (par.X(k, lambda, secZeta));
      }
      else if (m_component == PSDComponent::amplitude)
      {
         psd *= (par.Y(k, lambda, secZeta));
      }
      else if (m_component == PSDComponent::dispPhase)
      {
         psd *= (par.X_Z(k, lambda, lambda_wfs, secZeta));
      }
      else if (m_component == PSDComponent::dispAmplitude)
      {
         mxError("vonKarmanSpectrum::operator()", MXE_NOTIMPL, "Dispersive-aniso amplitude not implemented");
         return 0;
      }
      else
      {
         mxError("vonKarmanSpectrum::operator()", MXE_INVALIDARG, "Invalid component specified");
         return 0;
      }
   }
   
   return psd;
   
}

template< typename realT>
realT vonKarmanSpectrum<realT>::fittingError(aoAtmosphere<realT> &atm, realT d)
{
   realT k0;
   if(atm.L_0(0) > 0)
   {
      k0 = 1/ (atm.L_0(0)*atm.L_0(0));
   }
   else k0 = 0;

   return (math::pi<realT>() * math::six_fifths<realT>())* constants::a_PSD<realT>()/ pow(atm.r_0(), math::five_thirds<realT>()) * (1./pow( pow(0.5/d,2) + k0, math::five_sixths<realT>()));
}

template< typename realT>
template<typename iosT>
iosT & vonKarmanSpectrum<realT>::dumpPSD(iosT & ios)
{
   ios << "# PSD Parameters:" << '\n';
   ios << "#    ID = " << m_id << '\n';
   ios << "#    D = " << m_D  << '\n';
   ios << "#    subPiston = " << std::boolalpha << m_subPiston << '\n';
   ios << "#    subTipTilt = " << std::boolalpha << m_subTipTilt  << '\n';
   ios << "#    Scintillation = " << std::boolalpha << m_scintillation << '\n';
   ios << "#    Component = " << PSDComponent::compName(m_component) << '\n';
   return ios;
}
   
extern template
struct vonKarmanSpectrum<float>;

extern template
struct vonKarmanSpectrum<double>;

extern template
struct vonKarmanSpectrum<long double>;

#ifdef HASQUAD
extern template
struct vonKarmanSpectrum<__float128>;
#endif

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //aoPSDs_hpp
