/** \file aoPSDs.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Spatial power spectra used in adaptive optics.
  * \ingroup mxAO_files
  * 
  */

#ifndef __aoPSDs_hpp__
#define __aoPSDs_hpp__


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
   
///Manage calculations using the von Karman spatial power spectrum.
/** \ingroup mxAOAnalytic
  */ 
template<typename realT>
struct vonKarmanSpectrum
{
   
protected:
   bool _subPiston;///< flag controlling whether piston is subtracted from the PSD.  Default is true.
   bool _subTipTilt; ///< flag controlling whether tip and tilt are subtracted from the PSD.  Default is false.
   
   bool _scintillation; ///< flag controlling whether or not scintillation is included
   int _component; ///< If _scintillation is true, this controls whether phase (0), amplitude (1), or dispersive contrast (2) is returned.
   
   realT _D; ///< Diameter used for piston and tip/tilt subtraction, in m. Default is 1 m.

   const char * _id = "von Karman";
   
public:   
   
   ///Default Constructor
   vonKarmanSpectrum()
   {
      _subPiston = true;
      _subTipTilt = false;

      _scintillation = false;
      _component = 0;
            
      _D = 0;
   }
   
   ///Constructor specifying the parameters
   /**
     */ 
   vonKarmanSpectrum( bool subP, ///< [in] is the value of _subPiston.
                      bool subT, ///< [in] is the value of _subTipTilt.
                      realT D ///< [in] is the value of _D.
                    )
   {
      _subPiston = subP;
      _subTipTilt = subT;
      _D = D;
   }
   
   ///Get the value of _subPiston
   /**
     * \returns _subPiston
     */ 
   bool subPiston()
   {
      return _subPiston;
   }
   
   ///Set the value of _subPiston
   /**
     */ 
   void subPiston( bool sp /**< [in] is the new value of _subPiston */)
   {
      _subPiston = sp;
   }
   
   ///Get the value of _subTipTilt
   /**
     * \returns the value of _subTipTilt.
     */ 
   bool subTipTilt()
   {
      return _subTipTilt;
   }
   
   ///Set the value of the _subTipTilt flag.
   /**
     * \param st is the new value of _subTipTilt.
     */ 
   void subTipTilt(bool st)
   {
      _subTipTilt = st;
   }
   
   ///Get the value of _scintillation
   /**
     * \returns _scintillation
     */ 
   bool scintillation()
   {
      return _scintillation;
   }
   
   ///Set the value of _scintillation
   /**
     * \param sc is the new value of _scintillation
     */ 
   void scintillation(bool sc)
   {
      _scintillation = sc;
   }
   
   ///Get the value of _component
   /**
     * \returns _component
     */ 
   bool component()
   {
      return _component;
   }
   
   ///Set the value of _component
   /**
     * \param cc is the new value of _component
     */ 
   void component(int cc)
   {
      _component = cc;
   }

   ///Get the value of the diameter _D.
   /**
     * \returns _D, the diameter in m.
     */ 
   realT D()
   {
      return _D;
   }
   
   ///Set the aperture diameter
   /**
     * \pararm nd is the new diameter in m. 
     */
   void D(realT nd)
   {
      _D = nd;
   }
   
   ///Get the value of the PSD at spatial frequency k and a zenith distance.
   /**
     *   
     * \returns the von Karman PSD at the specified spatial frequency.
     * \returns -1 if an error occurs.
     * 
     */ 
   realT operator()( aoAtmosphere<realT> & atm, ///< [in] gives the atmosphere parameters r_0 and L_0.
                     realT k, ///< [in] is the spatial frequency in m^-1.
                     int n,
                     realT sec_zeta ///< [in] is the secant of the zenith distance.
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
   
      if( (_subPiston || _subTipTilt) )
      {
         if (_D == 0)
         {
            mxError("aoAtmosphere", MXE_PARAMNOTSET, "Diameter D not set for Piston and/or TT subtraction.");
            return -1;
         }
         
         if(_subPiston)
         {
            Ppiston = pow(2*math::func::jinc(pi<realT>()*k*_D), 2);
         }
         else Ppiston = 0;
 
         if(_subTipTilt)
         {
            Ptiptilt = pow(4*math::func::jinc2(pi<realT>()*k*_D), 2);
         }
         else Ptiptilt = 0;
      }
      else
      {
         Ppiston = 0;
         Ptiptilt = 0;
      }
      
      return constants::a_PSD<realT>()*pow(atm.r_0(), -five_thirds<realT>())*pow(k*k+k02, -eleven_sixths<realT>()) * (1.0-Ppiston - Ptiptilt)*sec_zeta;
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
                     int n,
                      realT lambda_wfs,      ///< [in] is the wavefront measurement wavelength in m
                      realT secZeta             ///< [in] is the secant of the zenith distance
                    )
   {
      realT psd = operator()(atm, k, 0, secZeta)* pow( atm.lam_0()/lambda, 2);
      
      if(psd < 0) return -1;
      
      if(_scintillation)
      {
         if(_component == 0)
         {
            psd *= atm.X(k, lambda, secZeta);
         }
         else if (_component == 1)
         {
            psd *= atm.Y(k, lambda, secZeta);
         }
         else
         {
            psd *= atm.X_Z(k, lambda, lambda_wfs, secZeta);
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

      realT lamc = 1.0;//pow( atm.lam_0()/(2*pi<realT>()),2);
      
      return (pi<realT>() * six_fifths<realT>())* a_PSD<realT>()/ pow(atm.r_0(), five_thirds<realT>()) * (1./pow( pow(0.5/d,2) + k0, five_sixths<realT>()))*lamc;
   }

   template<typename iosT>
   iosT & dumpPSD(iosT & ios)
   {
      ios << "# PSD Parameters:" << '\n';
      ios << "#    ID = " << _id << '\n';
      ios << "#    D = " << _D  << '\n';
      ios << "#    subPiston = " << std::boolalpha << _subPiston << '\n';
      ios << "#    subTipTilt = " << std::boolalpha << _subTipTilt  << '\n';
      ios << "#    Scintillation = " << std::boolalpha << _scintillation << '\n';
      ios << "#    Component = " << _component << '\n';
      return ios;
   }

};

} //namespace AO
} //namespace mx

#endif //__aoPSDs_hpp__
