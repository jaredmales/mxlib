/** \file astroSpectra.hpp
  * \author Jared R. Males
  * \brief Type definitions for the astronomical spectral libraries.
  * \ingroup astrophot
  *
  */

#ifndef mx_astro_astroSpectra_hpp
#define mx_astro_astroSpectra_hpp

#include "units.hpp"
#include "../readColumns.hpp"

namespace mx
{
namespace astro 
{

/// A basic spectrum
/** 
  * \ingroup astrophot_spectra 
  */
template<typename _units>
struct basicSpectrum
{
   typedef _units units;
   typedef typename units::realT realT;
   
   static const bool freq = false;
   
   typedef std::string paramsT; ///< The parameter is a string name.
      
   ///Specify how to convert to SI wavelength units. No conversions are performed in the basic spectrum
   static constexpr realT wavelengthUnits = static_cast<realT>(1);
   
   ///Specify how to convert to SI flux units. No conversions are performed in the basic spectrum
   static constexpr realT fluxUnits = static_cast<realT>(1);
   
   ///The data directory environment variable name.
   static constexpr const char * dataDirEnvVar = 0;
   
   ///This function should calculate the file name (without the path) for the spectrum parameters.
   static std::string fileName( const paramsT & name /**< [in] The parameters of the spectrum, in this case just its file name.*/)
   {
      return name;
   }
   
   ///This function reads the spectrum, returning its raw wavelength and spectrum points.  
   /** No unit conversions or interpolations should take place.
     */
   static std::string readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                                    std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                                    const std::string & path ///< [in] the full path to the file. 
                                  )
   {
      mx::readColumns(path, rawLambda, rawSpectrum);
   }
   
};


/// A spectrum from the astroFilt filter library
/** 
  * \ingroup astrophot_spectra
  */
template<typename _units, bool _rsr=true>
struct astroFilter
{
   typedef _units units;
   typedef typename units::realT realT;
   
   static const bool freq = false;
   
   
   typedef std::string paramsT; ///< The astroFilters are parameterized by name.
   
   ///Convert from um to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e6);
   
   ///No conversion is performed on filter transmission.
   static constexpr realT fluxUnits = static_cast<realT>(1); 
   
   static constexpr const char * dataDirEnvVar = "ASTROFILT_DATADIR";
   
   static std::string fileName( const std::string name )
   {
      return name + ".dat";
   }
   
   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path ///< [in] the full path to the file.
                          )
   {
      mx::readColumns(path, rawLambda, rawSpectrum);
      realT max = 0;
      
      for(int i=0;i<rawLambda.size(); ++i)
      {
         if(_rsr) rawSpectrum[i] = rawSpectrum[i]*rawLambda[i];
         if(rawSpectrum[i] > max) max = rawSpectrum[i];
      }
      
      for(int i=0;i<rawLambda.size(); ++i)
      {
         rawSpectrum[i] /= max;
      }
   
      return 0;
      
   }
   
};


/// A spectrum from the HST calspec library
/** See http://www.stsci.edu/hst/observatory/crds/calspec.html
  * 
  * \ingroup astrophot_spectra 
  */
template<typename _units>
struct calspecSpectrum
{
   typedef _units units;
   typedef typename units::realT realT;
   
   static const bool freq = false;
   
   typedef std::string paramsT; ///< The calspec Spectra are parameterized by star name.
      
   ///Convert from nm to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e10);
   
   ///Convert from erg s-1 cm-2 A-1 to SI W m-3
   static constexpr realT fluxUnits = static_cast<realT>(1e7) / (static_cast<realT>(1e4)*static_cast<realT>(1e10)); 
   
   static constexpr const char * dataDirEnvVar = "CALSPEC_DATADIR";
   
   ///The file name is found from the star's name.
   static std::string fileName( const std::string name )
   {
      if(name == "alpha_lyr" || name == "vega") return "alpha_lyr_stis_005.asc";
      else if(name == "1740346") return "1740346_nic_002.ascii";
      else if(name == "sun" || name == "sun_reference") return "sun_reference_stis.002.asc";
   }
   
   ///Read a CALSPEC spectrum, which is a simple two column ASCII format.
   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path ///< [in] the full path to the file.
                          )
   {
      return mx::readColumns(path, rawLambda, rawSpectrum);
   }
   
};

/// A spectrum from the Pickles library
/** 
  * \ingroup astrophot_spectra 
  */
template<typename _units>
struct picklesSpectrum
{
   typedef _units units;
   typedef typename units::realT realT;
   
   static const bool freq = false;
   
   typedef std::string paramsT; ///< The Pickles spectra are parameterized by a spectral type string.
      
   ///Convert from A to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e10);
   
   ///The Pickles spectra are dimensionless.
   static constexpr realT fluxUnits = static_cast<realT>(1); 
   
   ///The Pickles spectra location is specified by the PICKLES_DATADIR environment variable. 
   static constexpr const char * dataDirEnvVar = "PICKLES_DATADIR";
   
   ///The name of the datafile is constructed from its spectral type string.
   static std::string fileName( const std::string spt )
   {
      return "uk"+ mx::toLower(spt) + ".dat";
   }
   
   ///Read a Pickles spectrum, which for these purposes is a simple two column ASCII format.
   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path ///< [in] the full path to the file.
                          )
   {
      return mx::readColumns(path, rawLambda, rawSpectrum);
   }
   
};




} //namespace astro

} //namespace mx





#endif //mx_astro_astroSpectra_hpp
