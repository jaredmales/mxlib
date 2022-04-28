/** \file astroSpectra.hpp
  * \author Jared R. Males
  * \brief Type definitions for the astronomical spectral libraries.
  * \ingroup astrophot
  *
  */

#ifndef mx_astro_astroSpectra_hpp
#define mx_astro_astroSpectra_hpp

#include "units.hpp"
#include "../ioutils/readColumns.hpp"

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
   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                                    std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                                    const std::string & path, ///< [in] the full path to the file.
                                    const paramsT & params ///< [in] the parameters are passed in case needed to construct the spectrum
                                  )
   {
      return mx::ioutils::readColumns(path, rawLambda, rawSpectrum);
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

   static std::string fileName( const std::string & name )
   {
      return name + ".dat";
   }

   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path, ///< [in] the full path to the file.
                            const paramsT & params ///< [in] the parameters are passed in case needed to construct the spectrum
                          )
   {
      mx::ioutils::readColumns(path, rawLambda, rawSpectrum);
      realT max = 0;

      for(int i=0;i<rawLambda.size(); ++i)
      {
         //if(_rsr) rawSpectrum[i] = rawSpectrum[i]*rawLambda[i];
         if(rawSpectrum[i] > max) max = rawSpectrum[i];
      }

      for(int i=0;i<rawLambda.size(); ++i)
      {
         rawSpectrum[i] /= max;
      }

      return 0;

   }

};


/// A square-wave filter spectrum
/** Parameters specify the central wavelength, width, and sampling (all in microns) of a square-wave type filter.
  * 
  * To create this filter:
  * \code
  * typedef double realT;
  * 
  * realT lam0 = 0.5; //Central wavelength in microns
  * realT fw = 0.05; //Full-width, in microns.
  * realT dlam = 0.001; //Delta-lambda for specifying the defining points.  See note.
  * 
  * astroSpectrum<sqWaveFilter<units::si<realT>>> filt({ lam0, fw, dlam});
  * 
  * filt.setSpectrum(grid_meters); // grid_meters is a vector<realT> specifying the wavelength grid in meters
  * 
  * \endcode
  * 
  * Note that dlam specifies how sharp the filter edges are when interpolated.  Larger values will make the filter more trapezoidal.
  * 
  * \ingroup astrophot_spectra
  */
template<typename _units, bool _rsr=true>
struct sqWaveFilter
{
   typedef _units units;
   typedef typename units::realT realT;

   static const bool freq = false;

   ///The square wave is parameterized by the central wavelength, width, and sampling (all in microns).
   typedef struct
   {
      realT lam0; ///< The central Wavelength in microns
      realT fw; ///< The full width of the filter in microns
      realT dlam; ///< The wavelength sampling to use in microns.
   } paramsT;

   ///Convert from um to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e6);

   ///No conversion is performed on filter transmission.
   static constexpr realT fluxUnits = static_cast<realT>(1);

   static constexpr const char * dataDirEnvVar = 0;

   static std::string fileName( const paramsT & params )
   {
      return " "; //must not be empty to avoid error
   }

   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path, ///< [in] the full path to the file.
                            const paramsT & params ///< [in] the parameters are passed in case needed to construct the spectrum
                          )
   {
      rawLambda.resize(4);
      rawSpectrum.resize(4);

      rawLambda[0] = params.lam0 - 0.5*params.fw - 0.5*params.dlam;
      rawSpectrum[0] = 0.0;

      rawLambda[1] = params.lam0 - 0.5*params.fw + 0.5*params.dlam;
      rawSpectrum[1] = 1.0;

      rawLambda[2] = params.lam0 + 0.5*params.fw - 0.5*params.dlam;
      rawSpectrum[2] = 1.0;

      rawLambda[3] = params.lam0 + 0.5*params.fw + 0.5*params.dlam;
      rawSpectrum[3] = 0.0;

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

   ///Convert from A to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e10);

   ///Convert from erg s-1 cm-2 A-1 to SI W m-3
   static constexpr realT fluxUnits = static_cast<realT>(1e7) / (static_cast<realT>(1e4)*static_cast<realT>(1e10));

   static constexpr const char * dataDirEnvVar = "CALSPEC_DATADIR";

   ///The file name is found from the star's name.
   static std::string fileName( const std::string & name )
   {
      if(name == "alpha_lyr" || name == "vega") return "alpha_lyr_stis_005.asc";
      else if(name == "1740346") return "1740346_nic_002.ascii";
      else if(name == "sun" || name == "sun_reference") return "sun_reference_stis.002.asc";
   }

   ///Read a CALSPEC spectrum, which is a simple two column ASCII format.
   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path, ///< [in] the full path to the file.
                            const paramsT & params ///< [in] the parameters are passed in case needed to construct the spectrum
                          )
   {
      return mx::ioutils::readColumns(path, rawLambda, rawSpectrum);
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
   static std::string fileName( const std::string & spt )
   {
      return "uk"+ ioutils::toLower(spt) + ".dat";
   }

   ///Read a Pickles spectrum, which for these purposes is a simple two column ASCII format.
   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path, ///< [in] the full path to the file.
                            const paramsT & params ///< [in] the parameters are passed in case needed to construct the spectrum
                          )
   {
      return mx::ioutils::readColumns(path, rawLambda, rawSpectrum);
   }

};


/// Earth Albedo Spectra
/** The spectra can be one of:
  * - "EPOXI" returns the apparent albedo spectrum from Cowan and Strait (2013) \cite cowan_2013.
  * - "RawEarthshine" returns the unormalized albedo spectrum measured using Earthshine by Turnbull et al (2006) \cite turnbull_20016.
  * - "Earthshine" returns the Earthshine spectrum normalized to match the EPOXI result of 0.27 in the 550 nm band.
  *
  * \ingroup astrophot_spectra
  */
template<typename _units>
struct earthAlbedo
{
   typedef _units units;
   typedef typename units::realT realT;

   static const bool freq = false;

   typedef std::string paramsT; ///< The name of the spectrum can be "EPOXI", "Earthshine", or "RawEarthshine".

   ///Convert from A to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e6);

   ///The Earthshine is a dimensionless albedo.
   static constexpr realT fluxUnits = static_cast<realT>(1);

   ///The location is specified by the EARTHSHINE_DATADIR environment variable.
   static constexpr const char * dataDirEnvVar = "EARTHSHINE_DATADIR";

   ///The name of the datafile is a constant.
   static std::string fileName( const std::string & name )
   {
      if(name == "EPOXI") return "cowan_2013_EPOXI_albedo.dat";
      if(name == "Earthshine") return "earthshine_epoxi_normalized.dat";
      if(name == "RawEarthshine") return "Earthshine/F7_opt_NIR_ES_data.txt";

      mxError("earthAlbeo::fileName", MXE_INVALIDARG, "name not recognized.");

      return "";
   }

   ///Read the Earthshine albedo spectrum, which is a simple two column ASCII format.
   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path, ///< [in] the full path to the file.
                            const paramsT & params ///< [in] the parameters are passed in case needed to construct the spectrum
                          )
   {
      if(mx::ioutils::readColumns(path, rawLambda, rawSpectrum) < 0) return -1;
      return 0;
   }

};

/// Venus Spectra
/** 
  *
  * \ingroup astrophot_spectra
  */
template<typename _units>
struct venusAlbedo
{
   typedef _units units;
   typedef typename units::realT realT;

   static const bool freq = false;

   typedef std::string paramsT; ///< The name of the spectrum can be "venus"

   ///Convert from A to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e6);

   ///The Earthshine is a dimensionless albedo.
   static constexpr realT fluxUnits = static_cast<realT>(1);

   ///The location is specified by the EARTHSHINE_DATADIR environment variable.
   static constexpr const char * dataDirEnvVar = "VENUS_DATADIR";

   ///The name of the datafile is a constant.
   static std::string fileName( const std::string & name )
   {
      if(name == "Venus") return "venus_combined_albedo.dat";

      mxError("earthAlbeo::fileName", MXE_INVALIDARG, "name not recognized.");

      return "";
   }

   ///Read the Earthshine albedo spectrum, which is a simple two column ASCII format.
   static int readSpectrum( std::vector<realT> & rawLambda, ///< [out] the raw wavelength vector.  This should be an empty vector on input.
                            std::vector<realT> & rawSpectrum, ///< [out] the raw spectrum.  This should be an empty vector on input.
                            const std::string & path, ///< [in] the full path to the file.
                            const paramsT & params ///< [in] the parameters are passed in case needed to construct the spectrum
                          )
   {
      if(mx::ioutils::readColumns(path, rawLambda, rawSpectrum) < 0) return -1;
      return 0;
   }

};

} //namespace astro

} //namespace mx





#endif //mx_astro_astroSpectra_hpp
