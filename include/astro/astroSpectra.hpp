/** \file astroSpectra.hpp
  * \author Jared R. Males
  * \brief Type definitions for the astronomical spectral libraries.
  * \ingroup astrophot
  *
  */

#ifndef __mx_astro_astroSpectra_hpp__
#define __mx_astro_astroSpectra_hpp__

#include "units.hpp"
#include "../readColumns.hpp"

namespace mx
{
namespace astro 
{

/// A spectrum from the HST calspec library
/** See http://www.stsci.edu/hst/observatory/crds/calspec.html
  */
template<typename _units>
struct calspecSpectrum
{
   typedef _units units;
   typedef typename units::realT realT;
   
   static const bool freq = false;
   
   ///Convert from nm to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e10);
   
   ///Convert from erg s-1 cm-2 A-1 to SI W m-3
   static constexpr realT fluxUnits = static_cast<realT>(1e7) / (static_cast<realT>(1e4)*static_cast<realT>(1e10)); 
   
   static constexpr const char * dataDirEnvVar = "CALSPEC_DATADIR";
   
   static std::string fileName( const std::string name )
   {
      if(name == "alpha_lyr" || name == "vega") return "alpha_lyr_stis_005.asc";
      else if(name == "1740346") return "1740346_nic_002.ascii";
      else if(name == "sun" || name == "sun_reference") return "sun_reference_stis.002.asc";
   }
   
   static std::string readSpectrum( std::vector<realT> & rawLambda,
                                    std::vector<realT> & rawSpectrum,
                                    const std::string & path )
   {
      mx::readColumns(path, rawGrid, rawSpectrum);
   }
   
};






} //namespace astro

} //namespace mx





#endif //__mx_astro_astroSpectra_hpp__
