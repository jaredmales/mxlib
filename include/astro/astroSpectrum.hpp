/** \file astroSpectrum.hpp
  * \author Jared R. Males
  * \brief A class for working with astronomical spectra.
  * \ingroup astrophot
  *
  */

#ifndef __mx_astro_astroSpectrum_hpp__
#define __mx_astro_astroSpectrum_hpp__

#include <cmath>

#include "../environment.hpp"
#include "../readColumns.hpp"
#include "../gslInterpolation.hpp"

namespace mx
{
namespace astro 
{

///Class to manage an astronomical spectrum.
/** Details of the spectra, including units, file-location, and how they are read from disk, are specified by template parameter _spectrumT.
  * Spectra are loaded and immediately interpolated onto a wavelength grid.  This is to facilitate manipulation of spectra, i.e. for multiplying by filters, etc.
  */
template<typename _spectrumT, bool freq=false>
struct astroSpectrum
{
   typedef _spectrumT spectrumT;
   typedef typename spectrumT::units units;
   typedef typename units::realT realT;
   
   std::string _dataDir; ///< The directory containing the spectrum
      
   std::string _name; ///< The name of the spectrum

   std::vector<realT> _spectrum; ///< Contains the spectrum after it is set.
   
   ///Default c'tor
   astroSpectrum(){}
   
   ///Constructor specifying name, the enivronment will be queried for data directory.
   astroSpectrum( const std::string & name /**< [in] The name of the spectrum */)
   {
      setParameters(name);
   }
   
   ///Constructor specifying name and data directory.
   astroSpectrum( const std::string & name, ///< [in] The name of the spectrum 
                  const std::string & dataDir ///< [in] The directory containing the spectrum
                )
   {
      setParameters(name, dataDir);
   }
   
   void setParameters( const std::string & name /**< [in] The name of the spectrum */)
   {
      _name = name;
      
      if(spectrumT::dataDirEnvVar)
      {
         _dataDir = getEnv(spectrumT::dataDirEnvVar);
      }
   }
   
   void setParameters( const std::string & name,  ///< [in] The name of the spectrum
                       const std::string & dataDir ///< [in] The directory containing the spectrum
                     )
   {
      _name = name;
      _dataDir = dataDir;
   }
   
   ///Load the spectrum and interpolate it on a wavelength scale.  
   template<typename gridT>
   int setSpectrum( gridT & lambda )
   {
      std::vector<realT> rawLambda;
      std::vector<realT> rawSpectrum;
      
      if(_dataDir == "") _dataDir = ".";
      
      std::string path = _dataDir + "/" + spectrumT::fileName(_name);
   
      spectrumT::readSpectrum(rawLambda, rawSpectrum, path);
   
      for(int i=0; i < rawLambda.size(); ++i)
      {
         rawLambda[i] /= spectrumT::wavelengthUnits/(units::length);
         rawSpectrum[i] /= spectrumT::fluxUnits/(units::energy/(units::time * units::length * units::length * units::length));
      }
      
      mx::gsl_interpolate(gsl_interp_linear, rawLambda, rawSpectrum, lambda, _spectrum);
      
      for(int i=0; i < lambda.size(); ++i)
      {
         if( !std::isnormal(_spectrum[i])) _spectrum[i] = 0;
      }
      
      return 0;
   }
    
   realT & operator[](int i)
   {
      return _spectrum[i];
   }
   
   const realT operator[](int i) const
   {
      return _spectrum[i];
   } 
    
   template<typename compSpectrumT>
   std::vector<realT> operator*( const compSpectrumT & spec )
   {
      std::vector<realT> outVec( _spectrum.size() );
      
      for(int i=0; i< _spectrum.size(); ++i)
      {
         outVec[i] = _spectrum[i] * spec[i];
      }
      
      return outVec;
   }
};




} //namespace astro

} //namespace mx





#endif //__mx_astro_astroSpectrum_hpp__
