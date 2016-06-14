/** \file astroSpectrum.hpp
  * \author Jared R. Males
  * \brief Utilities for working with astronomical filters spectra
  * \ingroup astrophot
  *
  */

#ifndef __astroSpectra_hpp__
#define __astroSpectra_hpp__

#include "astroFilter.hpp"
#include "astroconstants.h"

namespace mx
{

/** \addtogroup astrophot
  * @{
  */

/// Read a Pickles atlas spectrum
/** Reads in spectrum from the Pickles atlas, usually found via the PICKLES_DATADIR
  * environment variable.  Returned spectra are in wavelength of microns, and arbitrary f_lambda flux units
  * 
  * \param lambda [output] the wavelength scale, in microns.
  * \param flambda [output] the spectrum, in arbitrary f_lambda flux units
  * \param spt  [input] the spectral type of the star.
  * \param datadir [intput] [optional] the directory containing the Pickles atlas
  * 
  * \tparam dataT the type of data to return
  */
template<typename dataT>
void readPickles( std::vector<dataT> & lambda, 
                  std::vector<dataT> & flambda,
                  const std::string & spt,
                  const std::string & datadir = "")
{
   std::string basepath;
   if(datadir == "")
   {
      basepath = getEnv("PICKLES_DATADIR");
   }
   else
   {
      basepath = datadir;
   }
      
   std::string fname = basepath + "/uk" + mx::toLower(spt) + ".dat";
   
   mx::readColumns(fname, lambda, flambda);
   
   for(int i=0; i<lambda.size(); ++i)
   {
      lambda[i] = lambda[i]/10000.;
   }
}


/// Read an HST Calspec spectrum
/** Reads in spectrum from the HST Calspec collection, usually found via the CALSPEC_DATADIR
  * environment variable.  Returned spectra are in wavelength of microns, and flux units of
  * ergs/s/um/cm^2 (converted from the native units of Calspec).
  * 
  * \param lambda [output] the profile wavelength scale, in microns.
  * \param flambda [output] the spectrum, in ergs/s/micron/cm^2
  * \param spName  [input] the name of the Calspec star, normally the prefix to the file name.
  * \param datadir [input] [optional] the directory containing the Calspec collection
  * 
  * \tparam dataT the type of data to return
  */
template<typename dataT>
void readCalspec( std::vector<dataT> & lambda, 
                  std::vector<dataT> & flambda,
                  const std::string & spName, 
                  const std::string & datadir)
{
   std::string basepath;
   if(datadir == "")
   {
      basepath = getEnv("CALSPEC_DATADIR");
   }
   else
   {
      basepath = datadir;
   }
   
   
   std::string csName;
   
   if(spName == "alpha_lyr" || spName == "vega") csName = "alpha_lyr_stis_005.asc";
   else if(spName == "1740346") csName = "1740346_nic_002.ascii";
   else if(spName == "sun" || spName == "sun_reference") csName = "sun_reference_stis.002.asc";

   if(csName == "")
   {
      std::cerr << "calspec spectrum not found for " <<   spName << "\n";
      return;
   }

      
   std::string fname = basepath + "/" + csName;
   
   mx::readColumns(fname, lambda, flambda);
 
   for(int i=0;i<lambda.size();++i)
   {
      //Convert Flam to per micron
      flambda[i] = flambda[i]*1e4;

      //Convert wavelength to microns
      lambda[i] = lambda[i]/1e4;
   }

}

template<typename dataT>
void readCalspec( std::vector<dataT> & lambda, 
                  std::vector<dataT> & flambda,
                  const std::string & spName )
{
   std::string datadir;
   
   readCalspec(lambda, flambda, spName, datadir);
}

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

//Used by rewritePhoenixSpectrum
template<typename floatT>
struct spectralPoint
{
   floatT lambda;
   floatT flambda;
   
   spectralPoint()
   {
      lambda = 0;
      flambda = 0;
   }
   
   spectralPoint(floatT l, floatT f)
   {
      lambda = l;
      flambda = f;
   }
};

//Used by rewritePhoenixSpectrum
template<typename floatT>
bool comp_spectralPoint (const spectralPoint<floatT> & sp1, const spectralPoint<floatT> & sp2) 
{ 
   return (sp1.lambda < sp2.lambda);
}


///Read in, crop, and re-write a Phoenix spectrum data file.
/** For working with spectra obtained from: http://perso.ens-lyon.fr/france.allard/ 
  * 
  * The Phoenix spectra contain many columns of line data which are not often used, are
  * not necessarily sorted by wavelength, and are sometimes formated so that points with leading 
  * minus signs are not in a separate column.  We also have to change the fortan D to e.
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
  * \param filename [input] complete name of the file to rewrite
  * \param lmin [input] minimum wavelength to rewrite [microns]
  * \param lmax [input] maximum wavelemngth to rewrite [microns]
  * 
  * \tparam floatT the floating point type in which to work
  */
template<typename floatT>
void rewritePhoenixSpectrum( const std::string & filename,
                             floatT lmin,
                             floatT lmax)
{
   float lambda, flambda;
        
   std::ifstream fin;
   
   fin.open(filename);
   
   if(!fin.good())
   {
      std::cerr << "Error opening file: " << filename << "\n";
      return;
   }
   
   std::vector<spectralPoint<floatT>> points;
   
   
   int lineSize = 1024;
   char * line = new char[lineSize];
   std::string sline;
   
   fin.getline(line, lineSize);
   sline = line;
    
   if(sline.length() < 25)
   {
      std::cerr << "Error reading file: " << filename << "\n";
   }
   lambda = mx::convertFromString<floatT>(sline.substr(1,12));
   
   //First have to diagnose if this has the column problem
   int nst;
   if(sline[13] == '-') nst = 13;
   else nst = 14;
   
   //convert the D to e
   sline[nst + 8] = 'e';
   
   flambda = mx::convertFromString<floatT>(sline.substr(nst,12));
   //if(flambda > 2) flambda = -flambda;
   
   if(lambda >= lmin*1e4 && lambda <= lmax*1e4)
   {
      points.push_back( spectralPoint<floatT>(lambda, flambda) );
   }

   fin.getline(line, lineSize);
   sline = line;
      
   while(fin.good())
   {
      if(sline.length() < 25) continue;
      sline[nst + 8] = 'e';
      
      lambda = mx::convertFromString<floatT>(sline.substr(1,12));
      flambda = mx::convertFromString<floatT>(sline.substr(nst,12));
      //if(flambda > 2) flambda = -flambda;
   
      if(lambda >= lmin*1e4 && lambda <= lmax*1e4)
      {
         points.push_back( spectralPoint<floatT>(lambda, flambda) );
      }
      
      fin.getline(line, lineSize);
      sline = line;
   }
   
   fin.close();
   delete line;
   
   std::sort(points.begin(), points.end(), comp_spectralPoint<floatT>);

   std::ofstream fout;
  
   
   std::string fname = filename;// + ".t";
   fout.open(filename);
   
   
   for(int i=0;i<points.size(); ++i)
   {
      fout.precision(10);
      fout << points[i].lambda << "    ";
      fout.precision(5);
      fout << points[i].flambda << "\n";
   }

   fout.close();

}

/// Read a Phoenix model spectrum
/** Reads in a Phoenix model spectrum specified by a filename.  Does not prepend PHOENIX_DATADIR.
  * 
  * NOTE: The spectrum must be formatted as if had been processed with rewritePhoenixSpectrum 
  * 
  * The units of the spectra are set by the scale parameter, with the following choices
  * - 0 == un-changed, returns the raw values from the data files
  * - 1 == microns, and ergs/sec/cm^2/micron, assuming the new format
  * 
  * References:
  * - https://phoenix.ens-lyon.fr/Grids/BT-Settl/README
  * - https://phoenix.ens-lyon.fr/Grids/FORMAT
  * 
  * \param lambda [output] the profile wavelength scale, in microns.
  * \param flambda [output] the spectrum, in arbitrary f_lambda flux units
  * \param filename [input] the subdirectory, specifying which model it is.
  * \param scale [input] [optional] whether and how to convert and scale the wavelength and spectrum, default is 1.
  * \param radius [input] [optional] if radius >0 and scale !=0 then flux is scaled to 10 pc.  
  *               In units of Jupiter radii, default is 1 R_jup.
  * 
  * \tparam dataT the type of data to return
  */
template<typename dataT>
void readPhoenix( std::vector<dataT> & lambda, 
                  std::vector<dataT> & flambda,
                  const std::string & filename,
                  int scale = 1,
                  dataT radius = 1 )
{
   
   mx::readColumns(filename, lambda, flambda);
   
   if(scale==1)
   {
      for(int i=0; i<lambda.size(); ++i)
      {
         flambda[i] = pow(10.0, flambda[i] - 8.0) * 1e4;
         
         if(radius > 0)
         {
            flambda[i] *= pow( radius/10.0*(RAD_JUPITER/PC_M), 2);
         }
         
         lambda[i] = lambda[i]/10000.;
      }
   }
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

}

#endif //__astroSpectra_hpp__

