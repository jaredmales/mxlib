/** \file phoenixSpectra.hpp
  * \author Jared R. Males
  * \brief Utilities for working with spectra from phoenix code.
  * \ingroup astrophot
  *
  */

#ifndef phoenixSpectrum_hpp
#define phoenixSpectrum_hpp


#include "../math/vectorUtils.hpp"
#include "../ioutils/fileUtils.hpp"

namespace mx
{
namespace astro
{

/// A spectrum from the Phoenix model, for use with ioutils::astro::astroSpectrum.
/** \ingroup astrophot_spectra
  */
template<typename _units>
struct phoenixSpectrum
{
   typedef _units units;
   typedef typename units::realT realT;
   
   static const bool freq = false;
   
   ///Convert from A to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e10);
   
   ///Convert from erg s-1 cm-2 A-1 to SI W m-3
   static constexpr realT fluxUnits = static_cast<realT>(1e7) / (static_cast<realT>(1e4)*static_cast<realT>(1e10)); 
   
   static constexpr const char * dataDirEnvVar = "PHOENIX_DATADIR";
   
   typedef std::string paramsT; ///< The parameter is a string name.
   
   static std::string fileName( const std::string & name )
   {
      return name;
   }
   
   static int readSpectrum( std::vector<realT> & rawLambda,
                            std::vector<realT> & rawSpectrum,
                            const std::string & path ,
                            const paramsT & params ///< [in] the parameters are passed in case needed to construct the spectrum
                          )
   {
      if( ioutils::readColumns(path, rawSpectrum) < 0) return -1;
      
      if( ioutils::readColumns(ioutils::parentPath(path) + "/wavelength.dat", rawLambda) < 0) return -1;
      
      return 0;
   }
   
   static void scaleSpectrum( std::vector<realT> & spectrum,
                              realT radius,
                              realT distance
                            )
   {
      for(int i=0; i < spectrum.size(); ++i)
      {
         spectrum[i] *= pow( radius / distance * ( constants::radJupiter<units>() / constants::parsec<units>()), 2);
      }
   }
   
};

///Read in, crop, scale, and re-write a Phoenix spectrum data file.
/** For working with spectra obtained from: http://perso.ens-lyon.fr/france.allard/ 
  * 
  * The Phoenix spectra contain many columns of line data which are not often used, are
  * not necessarily sorted by wavelength, and are sometimes formated so that points with leading 
  * minus signs are not in a separate column.  We also have to change the fortan D to e.
  * 
  * We also want to apply the dilution factor and take the power of 10.
  * 
  * This function deals with these issues, producing a two-column space-delimited
  * file in a specified wavelength range.
  * 
  * References:
  * - https://phoenix.ens-lyon.fr/Grids/BT-Settl/README
  * - https://phoenix.ens-lyon.fr/Grids/FORMAT
  * 
  * NOTE: This overwrites the input file!
  * 
  * 
  * \tparam floatT the floating point type in which to work
  * 
  * \ingroup astrophot_spectra
  */
template<typename floatT>
void rewritePhoenixSpectrum( const std::string & filename, ///< [in/out] complete name of the file to rewrite
                             floatT lmin,  ///< [in]  minimum wavelength to rewrite [microns]
                             floatT lmax,  ///< [in] maximum wavelemngth to rewrite [microns]
                             int sepWavelength = 0,  ///< [in] [optional] controls how wavelength is handled.  0=> wavelength is included in file as first column.  1=> wavelength is written to a separate 'wavelength.dat' file. -1=> means don't write wavelength.
                             floatT DF = -8.0  ///< [in] [optional] the dilution factor.  See references.
                           )
{
   float lambda, flambda;
        
   std::ifstream fin;
   
   fin.open(filename);
   
   if(!fin.good())
   {
      std::cerr << "Error opening file: " << filename << "\n";
      return;
   }
   
   std::vector<floatT> lambdas;
   std::vector<floatT> flambdas;
   
   int lineSize = 1024;
   char line[lineSize];
   std::string sline;
   
   fin.getline(line, lineSize);
   sline = line;
    
   if(sline.length() < 25)
   {
      std::cerr << "Error reading file: " << filename << "\n";
   }
   
   lambda = ioutils::convertFromString<floatT>(sline.substr(1,12));
   
   //First have to diagnose if this has the column problem
   int nst;
   if(sline[13] == '-') nst = 13;
   else nst = 14;
   
   //convert the D to e
   size_t dpos = sline.find('D', nst);
   sline[dpos] = 'e';
   
   
   flambda = ioutils::convertFromString<floatT>(sline.substr(nst,12));
   
   if(lambda >= lmin*1e4 && lambda <= lmax*1e4)
   {
      lambdas.push_back(lambda);
      flambdas.push_back(flambda);      
   }

   fin.getline(line, lineSize);
   sline = line;
      
   while(fin.good())
   {
      if(sline.length() < 25) continue;
      
      dpos = sline.find('D', nst);
      sline[dpos] = 'e';
   
      lambda = ioutils::convertFromString<floatT>(sline.substr(1,12));
      flambda = ioutils::convertFromString<floatT>(sline.substr(nst,12));
      
      if(lambda >= lmin*1e4 && lambda <= lmax*1e4)
      {
         lambdas.push_back(lambda);
         flambdas.push_back(flambda);
      }
      
      fin.getline(line, lineSize);
      sline = line;
   }
   
   fin.close();
   
   std::vector<size_t> idx = math::vectorSortOrder(lambdas);
   
   std::ofstream fout;
  
   std::string fname = filename;// + ".t";
   fout.open(filename);
   fout.precision(10);
   
   for(int i=0;i<lambdas.size(); ++i)
   {
      if(sepWavelength == 0)
      {
         fout << lambdas[idx[i]] << "    ";
      }
      fout << pow(10, flambdas[idx[i]] + DF) << "\n";
   }
   fout.close();

   if(sepWavelength == 1)
   {
      fname = ioutils::parentPath(filename) + "/wavelength.dat";
         
      fout.open(fname);
      fout.precision(10);
      for(int i=0;i<lambdas.size(); ++i)
      {
         fout << lambdas[idx[i]] << "\n";
      }
      fout.close();
   }
   
} //rewritePhoenixSpectrum

///Call rewritePhoenixSpectrum for all files in a directory.
/** \ingroup astropho_spectra
  */
template<typename floatT>
void rewritePhoenixSpectrumBatch( const std::string & dir,
                                  floatT lmin,
                                  floatT lmax,
                                  floatT DF = -8.0
                                )
{
   std::vector<std::string> flist = ioutils::getFileNames(dir, "lte", "", ".7");
   
   int sepWavelength = 1;
   for(int i=0; i< flist.size(); ++i)
   {
      rewritePhoenixSpectrum( flist[i], lmin, lmax, sepWavelength, DF);
      if(i == 0) sepWavelength = -1;
   }
}

} //namespace astro
} //namespace mx

#endif //phoenixSpectrum_hpp

