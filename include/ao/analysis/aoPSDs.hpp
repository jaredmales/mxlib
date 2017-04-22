/** \file aoPSDs.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Spatial power spectra used in adaptive optics.
  * \ingroup mxAOAnalytic_files
  * 
  */

#ifndef __aoPSDs_hpp__
#define __aoPSDs_hpp__


#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <mx/mxError.hpp>
#include <mx/jinc.hpp>

#include "aoConstants.hpp"
using namespace mx::AO::constants;

#include "aoAtmosphere.hpp"

namespace mx
{
   
namespace AO
{
   
///Manage calculations using the von Karman spatial power spectrum.
template<typename floatT>
struct vonKarmanSpectrum
{
   
protected:
   bool _subPiston;///< flag controlling whether piston is subtracted from the PSD.  Default is false.
   bool _subTipTilt; ///< flag controlling whether tip and tilt are subtracted from the PSD.  Default is false.
   
   bool _scintillation; ///< flag controlling whether or not scintillation is included
   bool _component; ///< If _scintillation is true, this controls whether phase (false) or amplitude (true) is returned.
   
   floatT _D; ///< Diameter used for piston and tip/tilt subtraction, in m. Default is 1 m.

   const char * _id = "von Karman";
   
public:   
   
   ///Default Constructor
   vonKarmanSpectrum()
   {
      _subPiston = false;
      _subTipTilt = false;

      _scintillation = false;
      _component = false;
            
      _D = 0;
   }
   
   ///Constructor specifying the parameters
   /**
     * \param subP is the value of _subPiston.
     * \param subT is the value of _subTipTilt.
     * \param D is the value of _D.
     */ 
   vonKarmanSpectrum(bool subP, bool subT, floatT D)
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
     * \param sp is the new value of _subPiston
     */ 
   void subPiston(bool sp)
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
   void component(bool cc)
   {
      _component = cc;
   }

   ///Get the value of the diameter _D.
   /**
     * \returns _D, the diameter in m.
     */ 
   floatT D()
   {
      return _D;
   }
   
   ///Set the aperture diameter
   /**
     * \pararm nd is the new diameter in m. 
     */
   void D(floatT nd)
   {
      _D = nd;
   }
   
   ///Get the value of the PSD at spatial frequency k.
   /**
     * \param atm gives the atmosphere parameters r_0 and L_0.
     * \param k is the spatial frequency is m^-1.
     *   
     * \returns the von Karman PSD at the specified spatial frequency.
     * \returns -1 if an error occurs.
     * 
     */ 
   floatT operator()(aoAtmosphere<floatT> & atm, floatT k)
   {
      floatT k02;
      
      if(atm.L_0() > 0)
      {
         k02 = (1)/(atm.L_0()*atm.L_0());
      }
      else k02 = 0;
 
      if(k02 == 0 && k == 0)
      {
         return 0;
      }
      
      floatT Ppiston, Ptiptilt;
   
      if( (_subPiston || _subTipTilt) )
      {
         if (_D == 0)
         {
            mxError("aoAtmosphere", MXE_PARAMNOTSET, "Diameter D not set for Piston and/or TT subtraction.");
            return -1;
         }
         
         if(_subPiston)
         {
            Ppiston = pow(2*mx::jinc(pi<floatT>()*k*_D), 2);
         }
         else Ppiston = 0;
 
         if(_subTipTilt)
         {
            Ptiptilt = pow(4*mx::jinc2(pi<floatT>()*k*_D), 2);
         }
         else Ptiptilt = 0;
      }
      else
      {
         Ppiston = 0;
         Ptiptilt = 0;
      }
      
      return constants::a_PSD<floatT>()*pow(atm.r_0(), -five_thirds<floatT>())*pow(k*k+k02, -eleven_sixths<floatT>()) * (1.0-Ppiston - Ptiptilt);//*exp(-k*k*l0*l0);
   }
   
   /// Get the value of the PSD at spatial frequency k and wavelength lambda
   /**
     * \param atm gives the atmosphere parameters r_0 and L_0.
     * \param k is the spatial frequency in m^-1.
     * \param lambda is the wavelength in m
     * 
     * \returns the von Karman PSD at the specified spatial frequency for the specified wavelength.
     * \returns -1 if an error occurs.
     */    
   floatT operator()(aoAtmosphere<floatT> & atm, floatT k, floatT lambda)
   {
      floatT psd = operator()(atm, k)* pow( atm.lam_0()/lambda, 2);
      
      if(psd < 0) return -1;
      
      if(_scintillation)
      {
         if(_component == false)
         {
            psd *= atm.X(k, lambda);
         }
         else
         {
            psd *= atm.Y(k, lambda);
         }
      }
      
      return psd;
      
   }
   
   /// Get the fitting error for an actuator spacing d.
   /**
     * \param atm gives the atmosphere parameters r_0 and L_0
     * \param d is the actuator spacing in m
     */    
   floatT fittingError(aoAtmosphere<floatT> &atm, floatT d)
   {
      floatT k0;
      if(atm.L_0() > 0)
      {
         k0 = 1/ (atm.L_0()*atm.L_0());
      }
      else k0 = 0;

      floatT lamc = 1.0;//pow( atm.lam_0()/(2*pi<floatT>()),2);
      
      return (pi<floatT>() * six_fifths<floatT>())* a_PSD<floatT>()/ pow(atm.r_0(), five_thirds<floatT>()) * (1./pow( pow(0.5/d,2) + k0, five_sixths<floatT>()))*lamc;
   }

   template<typename iosT>
   iosT & dumpPSD(iosT & ios)
   {
      ios << "# PSD Parameters:" << '\n';
      ios << "#    ID = " << _id << '\n';
      ios << "#    D = " << _D  << '\n';
      ios << "#    subPiston = " << _subPiston << '\n';
      ios << "#    subTipTilt = " << _subTipTilt  << '\n';
      ios << "#    Scintillation = " << _scintillation << '\n';
      ios << "#    Component = " << _component << '\n';
      return ios;
   }

};

} //namespace AO
} //namespace mx

#endif //__aoPSDs_hpp__