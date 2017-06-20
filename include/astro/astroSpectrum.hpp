/** \file astroSpectrum.hpp
  * \author Jared R. Males
  * \brief A class for working with astronomical spectra.
  * \ingroup astrophot
  *
  */

#ifndef __mx_astro_astroSpectrum_hpp__
#define __mx_astro_astroSpectrum_hpp__

#include "../environment.hpp"
#include "../readColumns.hpp"
#include "../gslInterpolation.hpp"

namespace mx
{
namespace astro 
{


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
      _dataDir = getEnv(spectrumT::dataDirEnvVar);
   }
   
   void setParameters( const std::string & name,  ///< [in] The name of the spectrum
                       const std::string & dataDir ///< [in] The directory containing the spectrum
                     )
   {
      _name = name;
      _dataDir = dataDir;
   }
   
   template<typename gridT>
   int setSpectrum( gridT & grid )
   {
      std::vector<realT> rawGrid;
      std::vector<realT> rawSpectrum;
      
      std::string path = _dataDir + "/" + spectrumT::fileName(_name);
   
      spectrumT::readSpectrum(rawGrid, rawSpectrum, path);
   
      for(int i=0; i < rawGrid.size(); ++i)
      {
         rawGrid[i] /= spectrumT::wavelengthUnits/(units::length);
         rawSpectrum[i] /= spectrumT::fluxUnits/(units::energy/(units::time * units::length * units::length * units::length));
      }
      
      _spectrum.clear();
      
      mx::gsl_interpolate(gsl_interp_linear, rawGrid, rawSpectrum, grid, _spectrum);
      
      return 0;
   }
    
   realT operator[](int i)
   {
      return _spectrum[i];
   } 
    
};




} //namespace astro

} //namespace mx





#endif //__mx_astro_astroSpectrum_hpp__
