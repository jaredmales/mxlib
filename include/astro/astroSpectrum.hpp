/** \file astroSpectrum.hpp
  * \author Jared R. Males
  * \brief A class for working with astronomical spectra.
  * \ingroup astrophot
  *
  */

#ifndef mx_astro_astroSpectrum_hpp
#define mx_astro_astroSpectrum_hpp

#include <cmath>

#include "../environment.hpp"
#include "../readColumns.hpp"
#include "../gslInterpolation.hpp"

namespace mx
{
namespace astro 
{

///Base spectrum class which provides manipulation and characterization functionality.
/**
  * \ingroup astrophot
  */ 
template<typename realT>
struct baseSpectrum
{
   std::vector<realT> _spectrum; ///< Contains the spectrum after it is set.
   
   ///Access a single point in the spectrum, specified by its vector index.
   /**
     * \returns a reference to the value at i. 
     */
   realT & operator[](size_t i /**< [in] the index of the spectral point*/)
   {
      return _spectrum[i];
   }
   
   ///Access a single point in the spectrum, specified by its vector index.
   /**
     * \returns the value of point i. 
     */
   const realT operator[](size_t i /**< [in] the index of the spectral point*/) const
   {
      return _spectrum[i];
   } 
    
   ///Multiply two spectra together.
   /**
     * \returns a new baseSpectrum
     */
   template<typename compSpectrumT>
   baseSpectrum<realT> operator*( const compSpectrumT & spec /**< [in] the spectrum to multiply by */ )
   {
      baseSpectrum outVec;
      outVec._spectrum.resize( _spectrum.size() );
      
      for(int i=0; i< _spectrum.size(); ++i)
      {
         outVec[i] = _spectrum[i] * spec[i];
      }
      
      return outVec;
   }

   ///Calculate the mean value of the spectrum.
   realT mean()
   {
      realT sum = 0;
      
      for(int i=0; i< _spectrum.size(); ++i) sum += _spectrum[i];
      
      return sum/_spectrum.size();
   }
   
   ///Calculate the mean value of the spectrum when mutiplied by another.
   /** The result is normalized by the mean value of the input spectrum, equivalent to:
     * \f[
      \mu = \frac{\int T(\lambda)S(\lambda) d\lambda}{\int T(\lambda) d\lambda}
      \f]
     * For instance, use this to get the mean value of a spectrum in a filter.
     * 
     */
   template<typename compSpectrumT>
   realT mean( const compSpectrumT & T /**< [in] the spectrum to multiply by*/ )
   {
      realT sum1 = 0, sum2 = 0;
      
      for(int i=0; i< _spectrum.size(); ++i)
      {
         sum1 += _spectrum[i]*T[i];
         sum2 += T[i];
      }
      
      return sum1/sum2;
   }
   
};

///Class to manage an astronomical spectrum.
/** Details of the spectra, including units, file-location, and how they are read from disk, are specified by template parameter _spectrumT:
  * see \ref astrophot_spectra "the available spectrum definitions." 
  * Spectra are loaded and immediately interpolated onto a wavelength grid.  This is to facilitate manipulation of spectra, i.e. 
  * for multiplying by filters, etc.
  * 
  * Inherits functions and operators to manipulate and characterize spectra from \ref mx::astro::baseSpectrum "baseSpectrum".  
  * 
  * \tparam _spectrumT is the underlying spectrum type, which provides (through static interfaces) the specifications of how to read or calculate the spectrum.
  * \tparam freq specify whether this spectrum uses frequency instead of wavelength.  Default is false.
  * 
  * \ingroup astrophot
  */
template<typename _spectrumT, bool freq=false>
struct astroSpectrum : public baseSpectrum<typename _spectrumT::units::realT>
{
   typedef _spectrumT spectrumT;
   typedef typename spectrumT::units units;
   typedef typename units::realT realT;
   typedef typename spectrumT::paramsT paramsT;
   
   std::string _dataDir; ///< The directory containing the spectrum
      
   paramsT _params; ///< The parameters of the spectrum, e.g. its name or the numerical parameters needed to generate the name.
   
   ///Default c'tor
   astroSpectrum(){}
   
   ///Constructor specifying name, the enivronment will be queried for data directory.
   astroSpectrum( const paramsT & params /**< [in] The name of the spectrum */)
   {
      setParameters(params);
   }
   
   ///Constructor specifying name and data directory.
   astroSpectrum( const paramsT & params, ///< [in] The name of the spectrum 
                  const std::string & dataDir ///< [in] The directory containing the spectrum
                )
   {
      setParameters(params, dataDir);
   }
   
   ///Set the parameters of the spectrum, using the underlying spectrums parameter type.
   int setParameters( const paramsT & params /**< [in] The name of the spectrum */)
   {
      _params = params;
      
      if(spectrumT::dataDirEnvVar)
      {
         _dataDir = getEnv(spectrumT::dataDirEnvVar);
      }
      
      return 0;
   }
   
   ///Set the parameters of the spectrum, using the underlying spectrums parameter type.
   /** This version also sets the data directory, instead of using the enivronment variable.
     *
     * \overload
     */ 
   int setParameters( const paramsT & params,  ///< [in] The name of the spectrum
                      const std::string & dataDir ///< [in] The directory containing the spectrum
                    )
   {
      _params = params;
      _dataDir = dataDir;
      
      return 0;
   }
   
   ///Load the spectrum and interpolate it onto a wavelength scale. 
   template<typename gridT>
   int setSpectrum( gridT & lambda )
   {
      std::vector<realT> rawLambda;
      std::vector<realT> rawSpectrum;
      
      if(_dataDir == "") _dataDir = ".";
      
      std::string path = _dataDir + "/" + spectrumT::fileName(_params);
   
      if( spectrumT::readSpectrum(rawLambda, rawSpectrum, path, _params) < 0)
      {
         return -1; ///\returns -1 on an error reading the spectrum.
      }
   
      //Unit conversions
      for(int i=0; i < rawLambda.size(); ++i)
      {
         rawLambda[i] /= spectrumT::wavelengthUnits/(units::length);
         rawSpectrum[i] /= spectrumT::fluxUnits/(units::energy/(units::time * units::length * units::length * units::length));
      }
      
      mx::gsl_interpolate(gsl_interp_linear, rawLambda, rawSpectrum, lambda, this->_spectrum);
      
      for(int i=0; i < lambda.size(); ++i)
      {
         if( !std::isnormal(this->_spectrum[i])) this->_spectrum[i] = 0;
      }
      
      return 0; ///\returns 0 on success.
   }
    
   
};




} //namespace astro

} //namespace mx





#endif //mx_astro_astroSpectrum_hpp
