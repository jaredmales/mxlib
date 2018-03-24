/** \file aoPSDs.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Spatial power spectra used in adaptive optics.
  * \ingroup mxAO_files
  * 
  */

#ifndef aoPSDs_hpp
#define aoPSDs_hpp


#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <mx/mxError.hpp>
#include <mx/math/func/jinc.hpp>

#include "aoConstants.hpp"
using namespace mx::AO::constants;

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
        
   std::string compName( int cc )
   {
      if( cc == phase ) return "phase";
      if( cc == amplitude ) return "amplitude";
      if( cc == dispPhase ) return "dispPhase";
      if( cc == dispAmplitude ) return "dispAmplitude";
      
      return "unknown";
   }
}

///Manage calculations using the von Karman spatial power spectrum.
/** \ingroup mxAOAnalytic
  */ 
template<typename realT>
struct vonKarmanSpectrum
{
   
protected:
   bool m_subPiston;///< flag controlling whether piston is subtracted from the PSD.  Default is true.
   bool m_subTipTilt; ///< flag controlling whether tip and tilt are subtracted from the PSD.  Default is false.
   
   bool m_scintillation; ///< flag controlling whether or not scintillation is included
   int m_component; ///< If m_scintillation is true, this controls whether phase (0), amplitude (1), or dispersive contrast (2) is returned.
   
   realT m_D; ///< Diameter used for piston and tip/tilt subtraction, in m. Default is 1 m.

   const char * m_id = "von Karman";
   
   realT m_alpha {eleven_thirds<realT>()}; ///< The power-law index, 11/3 for Kolmogorov
   
public:   
   
   ///Default Constructor
   vonKarmanSpectrum()
   {
      m_subPiston = true;
      m_subTipTilt = false;

      m_scintillation = false;
      m_component = PSDComponent::phase;
            
      m_D = 0;
   }
   
   ///Constructor specifying the parameters
   /**
     */ 
   vonKarmanSpectrum( bool subP, ///< [in] is the value of m_subPiston.
                      bool subT, ///< [in] is the value of m_subTipTilt.
                      realT D ///< [in] is the value of m_D.
                    )
   {
      m_subPiston = subP;
      m_subTipTilt = subT;
      m_D = D;
   }
   
   ///Get the value of m_subPiston
   /**
     * \returns m_subPiston
     */ 
   bool subPiston()
   {
      return m_subPiston;
   }
   
   ///Set the value of m_subPiston
   /**
     */ 
   void subPiston( bool sp /**< [in] is the new value of m_subPiston */)
   {
      m_subPiston = sp;
   }
   
   ///Get the value of m_subTipTilt
   /**
     * \returns the value of m_subTipTilt.
     */ 
   bool subTipTilt()
   {
      return m_subTipTilt;
   }
   
   ///Set the value of the m_subTipTilt flag.
   /**
     * \param st is the new value of m_subTipTilt.
     */ 
   void subTipTilt(bool st)
   {
      m_subTipTilt = st;
   }
   
   ///Get the value of m_scintillation
   /**
     * \returns m_scintillation
     */ 
   bool scintillation()
   {
      return m_scintillation;
   }
   
   ///Set the value of m_scintillation
   /**
     * \param sc is the new value of m_scintillation
     */ 
   void scintillation(bool sc)
   {
      m_scintillation = sc;
   }
   
   ///Get the value of m_component
   /**
     * \returns m_component
     */ 
   int component()
   {
      return m_component;
   }
   
   ///Set the value of m_component
   /**
     */ 
   int component( int cc /**< [in] the new value of m_component */)
   {
      if( cc != PSDComponent::phase && cc != PSDComponent::amplitude && cc != PSDComponent::dispPhase && cc != PSDComponent::dispAmplitude )
      {
         mxError("vonKarmanSpectrum::component", MXE_INVALIDARG, "Unknown component");
         return -1;
      }
      
      m_component = cc;
      return 0;         
   }

   ///Get the value of the diameter m_D.
   /**
     * \returns the current value of m_D, the diameter in m.
     */ 
   realT D()
   {
      return m_D;
   }
   
   ///Set the aperture diameter
   /**
     */
   void D(realT nd /**< [in] the new diameter in m */)
   {
      m_D = nd;
   }
   
   ///Get the value of the PSD index alpha.
   /**
     * \returns the current value of m_alpha, the PSD index.
     */
   realT alpha()
   {
      return m_alpha;
   }
   
   ///Set the PSD index alpha
   /**
     */
   void alpha(realT na /**< [in] the new PSD index*/ )
   {
      m_alpha = na;
   }
   
   ///Get the value of the PSD at spatial frequency k and a zenith distance.
   /**
     * \returns the von Karman PSD at the specified spatial frequency.
     * \returns -1 if an error occurs.
     * 
     */ 
   realT operator()( aoAtmosphere<realT> & atm, ///< [in] gives the atmosphere parameters r_0 and L_0.
                     realT k,                   ///< [in] is the spatial frequency in m^-1.
                     realT sec_zeta             ///< [in] is the secant of the zenith distance.
                   )
   {
      realT k02;
      
      if(atm.L_0() > 0)
      {
         k02 = (1)/(atm.L_0()*atm.L_0());
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
            Ppiston = pow(2*math::func::jinc(pi<realT>()*k*m_D), 2);
         }
         else Ppiston = 0;
 
         if(m_subTipTilt)
         {
            Ptiptilt = pow(4*math::func::jinc2(pi<realT>()*k*m_D), 2);
         }
         else Ptiptilt = 0;
      }
      else
      {
         Ppiston = 0;
         Ptiptilt = 0;
      }
      
      return constants::a_PSD<realT>()*pow(atm.r_0(), -five_thirds<realT>())*pow(k*k+k02, -1*m_alpha/2) * (1.0-Ppiston - Ptiptilt)*sec_zeta;
   }
   
   /// Get the value of the PSD at spatial frequency k and wavelength lambda, and a zenith distance, with a WFS at a different wavelength
   /**
     * 
     * \returns the von Karman PSD at the specified spatial frequency for the specified wavelength.
     * \returns -1 if an error occurs.
     */    
   realT operator()( aoAtmosphere<realT> & atm, ///< [in] gives the atmosphere parameters r_0 and L_0.
                     realT k,                   ///< [in] is the spatial frequency in m^-1.
                     realT lambda,              ///< [in] is the observation wavelength in m
                     realT lambda_wfs,         ///< [in] is the wavefront measurement wavelength in m
                     realT secZeta             ///< [in] is the secant of the zenith distance
                   )
   {
      realT psd = operator()(atm, k, secZeta)* pow( atm.lam_0()/lambda, 2);
      
      if(psd < 0) return -1;
      
      if(m_scintillation)
      {
         if(m_component == PSDComponent::phase)
         {
            psd *= atm.X(k, lambda, secZeta);
         }
         else if (m_component == PSDComponent::amplitude)
         {
            psd *= atm.Y(k, lambda, secZeta);
         }
         else if (m_component == PSDComponent::dispPhase)
         {
            psd *= atm.X_Z(k, lambda, lambda_wfs, secZeta);
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
   
   /// Get the fitting error for an actuator spacing d.
   /**
     * \param atm gives the atmosphere parameters r_0 and L_0
     * \param d is the actuator spacing in m
     */    
   realT fittingError(aoAtmosphere<realT> &atm, realT d)
   {
      realT k0;
      if(atm.L_0() > 0)
      {
         k0 = 1/ (atm.L_0()*atm.L_0());
      }
      else k0 = 0;

      return (pi<realT>() * six_fifths<realT>())* a_PSD<realT>()/ pow(atm.r_0(), five_thirds<realT>()) * (1./pow( pow(0.5/d,2) + k0, five_sixths<realT>()));
   }

   template<typename iosT>
   iosT & dumpPSD(iosT & ios)
   {
      ios << "# PSD Parameters:" << '\n';
      ios << "#    ID = " << m_id << '\n';
      ios << "#    alpha = " << m_alpha << '\n';
      ios << "#    D = " << m_D  << '\n';
      ios << "#    subPiston = " << std::boolalpha << m_subPiston << '\n';
      ios << "#    subTipTilt = " << std::boolalpha << m_subTipTilt  << '\n';
      ios << "#    Scintillation = " << std::boolalpha << m_scintillation << '\n';
      ios << "#    Component = " << PSDComponent::compName(m_component) << '\n';
      return ios;
   }

};

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //aoPSDs_hpp
