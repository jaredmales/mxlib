

#ifndef __astroSpectra_hpp__
#define __astroSpectra_hpp__

//#include "astroFilter.hpp"
//#include "astroconstants.h"

#include "math/vectorUtils.hpp"

namespace mx
{
#if 0
   
/** \addtogroup astrophot
  * @{
  */





/// Read the ATRAN atmospheric transmission spectrum from Gemini
/**
  *
  * \param lambda the wavelength scale
  * \param trans the transmission
  * \param site either "cp" for Cerro Pachon or "mk" for Mauna Kea
  * \param pwv precipitable water vaper, for "cp" choices are 2.3, 4.3, 7.6, 10.0.  For "mk" choices are 1.0, 1.6, 3.0, 5.0.   
  * \param airmass airmass, choices are 1.0, 1.5, 2.0.
  * \param datadir [optional] the directory containing the Gemini ATRAN files
  *
  * \tparam dataT the type of data
  */ 
template<typename dataT>
void readATRAN( std::vector<dataT> & lambda, 
                std::vector<dataT> & trans,
                const std::string & site, 
                double pwv, 
                double airmass, 
                const std::string & datadir = "" )
{
   std::string basepath;
   if(datadir == "")
   {
      basepath = getEnv("ATRAN_DATADIR");
   }
   else
   {
      basepath = datadir;
   }
   
   std::string fname;
   
   if(site == "cp") fname = basepath + "/cptrans_zm_";
   else if (site == "mk") fname = basepath + "/mktrans_zm_";
   
   char digits[32];
   
   snprintf(digits, 32, "%i", (int) (pwv*10.0));
   
   fname += digits;
   fname += "_";

   snprintf(digits, 32, "%i", (int) (airmass*10.0));
   
   fname += digits;
   fname += ".dat";
      
   readColumns(fname, lambda, trans);
   
}

/// Read a BTRAM atmospheric transmission spectrum 
/**
  *
  * \param lambda [output] the wavelength scale
  * \param trans [output] the transmission
  * \param site [input] currently only "manqui", for Cerro Manqui at LCO, is supported.
  * \param run [input] the calculation run designation
  * \param pwv [input] precipitable water vapor   
  * \param zd [input] zenith distance
  * \param datadir [input] [optional] the directory containing the Gemini ATRAN files
  *
  * \tparam dataT the type of data
  */ 
template<typename dataT>
void readBTRAM( std::vector<dataT> & lambda, 
                std::vector<dataT> & trans,
                const std::string & site,
                const std::string & run,
                double pwv, 
                double zd, 
                const std::string & datadir = "" )
{
   std::string basepath;
   if(datadir == "")
   {
      basepath = getEnv("BTRAM_DATADIR");
   }
   else
   {
      basepath = datadir;
   }
   
   std::string fname;
   
   fname = basepath + "/" + site + "/" + run + "/" + site + "_";
   
   char digits[32];
   
   snprintf(digits, 32, "zd%0.1f_", zd);
   
   fname += digits;

   snprintf(digits, 32, "pwv%0.1f.txt", pwv);
   
   fname += digits;
      

   
   readColumns<','>(fname, lambda, trans);
   
}



///Read the Earthshine spectrum.
/** Reads the Earthshine spectrum from Turnbull et al (2006) [http://adsabs.harvard.edu/abs/2006ApJ...644..551T].  
  * 
  * It is assumed that this is a two column ascii file.
  * 
  * If datadir is empty, looks in the EARTHSHINE_DATADIR environment variable for the 
  * path to the data file.
  * 
  * \param lambda [output] the wavelength grid
  * \param albedo [output] the albedo in arbitrary units
  * \param datadir [input] the directory containing the Earthshine spectrum.  If empty will be 
  *                        filled in with environment value of EARTHSHINE_DATADIR
  * 
  * \tparam dataT the data type in which the data file is read and converted from ASCII.
  */
template<typename dataT>
void readEarthShine( std::vector<dataT> & lambda, 
                     std::vector<dataT> & albedo,
                     std::string & datadir = "" )
{
   if(datadir == "")
   {
      datadir = getEnv("EARTHSHINE_DATADIR");
      
      if(datadir == "")
      {
         std::cerr << "readEarthShine: data directory not found.  Trying pwd.\n";
      }
   }
   
   std::string fname; 
   
   fname = datadir + "/F7_opt_NIR_ES_data.txt";
   
   readColumns(fname, lambda, albedo);
}

///Read the Earthshine spectrum.
/** Reads the Earthshine spectrum from Turnbull et al (2006) [http://adsabs.harvard.edu/abs/2006ApJ...644..551T].  
  * 
  * It is assumed that this is a two column ascii file.
  * 
  * Looks in the EARTHSHINE_DATADIR environment variable for the 
  * path to the data file.
  * 
  * \param lambda [output] the wavelength grid
  * \param albedo [output] the albedo in arbitrary units
  * 
  * \tparam dataT the data type in which the data file is read and converted from ASCII.
  */
template<typename dataT>
void readEarthShine( std::vector<dataT> & lambda, 
                     std::vector<dataT> & albedo )
{
   std::string datadir;
   readEarthShine(lambda, albedo, datadir);
}




/// @}
#endif

}

#endif //__astroSpectra_hpp__

