

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

#endif

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
  * \param filename [input] complete name of the file to rewrite
  * \param lmin [input] minimum wavelength to rewrite [microns]
  * \param lmax [input] maximum wavelemngth to rewrite [microns]
  * 
  * \tparam floatT the floating point type in which to work
  */
template<typename floatT>
void rewritePhoenixSpectrum( const std::string & filename,
                             floatT lmin,
                             floatT lmax,
                             int sepWavelength = 0,
                             floatT DF = -8.0
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
   
   //std::vector<spectralPoint<floatT>> points;
   std::vector<floatT> lambdas;
   std::vector<floatT> flambdas;
   
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
   
   if(lambda >= lmin*1e4 && lambda <= lmax*1e4)
   {
      //points.push_back( spectralPoint<floatT>(lambda, flambda) );
      lambdas.push_back(lambda);
      flambdas.push_back(flambda);
      
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
         //points.push_back( spectralPoint<floatT>(lambda, flambda) );
         lambdas.push_back(lambda);
         flambdas.push_back(flambda);
      }
      
      fin.getline(line, lineSize);
      sline = line;
   }
   
   fin.close();
   delete line;
   
   //std::sort(points.begin(), points.end(), comp_spectralPoint<floatT>);
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
      fname = "wavelength.dat";
      fout.open(fname);
      fout.precision(10);
      for(int i=0;i<lambdas.size(); ++i)
      {
         fout << lambdas[idx[i]] << "\n";
      }
      fout.close();
   }
   
}

#if 0

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


template<typename dataT>
int readCahoy( std::vector<dataT> & lambda,
               std::vector<dataT> & albedo,
               double sepAU,
               int metal,
               double phase,
               std::string & datadir)
{
   if(datadir == "")
   {
      datadir = getEnv("CAHOYALBEDO_DATADIR");
      
      if(datadir == "")
      {
         std::cerr << "readCahoy: data directory not found.  Trying pwd.\n";
      }
   }
   
   std::string fname; 
   
   fname = datadir + "/";
   
   if(sepAU == 0.8)
   {
      fname += "0.8";
   }
   else if(sepAU == 2.0)
   {
      fname += "2.0";
   }
   else if(sepAU == 5.0)
   {
      fname += "5";
   }
   else if(sepAU = 10.0)
   {
      fname += "10";
   }
   else
   {
      mxError("readCahoy", MXE_INVALIDARG, "invalid separation (sepAU)");
      return -1;
   }

   fname += "au_";
   fname += convertToString<int>(metal);
   fname += "x_albedos/00_Spectral_Albedos   ";  
   
   char pstr[9];
   int l = snprintf(pstr, 9, "%0.5f", phase);
   
   fname+=pstr;
   
   fname += ".txt";
   
   std::vector<int> tmp1, tmp2;
   
   readColumns(fname, tmp1, lambda, tmp2, albedo);

   if(lambda.size() == 0)
   {
      std::cerr << "None: " << fname << "\n";
   }
}

struct cahoyGrid
{
   std::vector<double> _sep;
   std::vector<double> _phase;
   
   std::vector<double> _lambda;
   std::vector<std::vector<std::vector<double>>> _ag;
   
   
   void openGrid(int metal)
   {
      _sep = {0.8, 2.0, 5.0, 10.0};
      _phase = {0.0, 0.1745, 0.3491, 0.5236, 0.6981, 0.8727, 1.0472, 1.2217, 1.3963,1.5708,1.7453,1.9199, 2.0944, 2.2689, 2.4435, 2.6180, 2.7925,2.9671, 3.139};
      
      _ag.resize( _sep.size());
      for(int i=0; i< _ag.size(); ++i) _ag[i].resize(_phase.size());
      
      std::vector<double> tmpl, tmpag;
      for(int i=0; i< _sep.size(); ++i)
      {
         for(int j=0; j < _phase.size(); ++j)
         {
            tmpl.clear();
            tmpag.clear();
           
            std::string ddir = "/home/jrmales/Data/Spectra/Models/cahoy_EGP_albedos";
            readCahoy(tmpl, tmpag, _sep[i], metal, _phase[j], ddir);
            _ag[i][j] = tmpag;
         }
      }

      _lambda = tmpl;
      
   }
   
   ///Get an interpolated spectrum at arbitrary non-grid point using bilinear interpolation.
   /** If outside the endpoints of the grid, the grid is simply extended (i.e. constant albedo).
     */
   int getAgSpectrum( std::vector<double> & lambda, std::vector<double> & spect, double sep, double phase)
   {
      int i_l, i_u, j_l, j_u;
      
      phase = fmod(phase, pi<double>());
      
      if(phase < 0) phase += pi<double>();
      
      i_l = _sep.size()-1;
      while( _sep[i_l] > sep && i_l > 0) --i_l;
      
      i_u = 0;
      while( _sep[i_u] < sep && i_u < _sep.size()-1) ++i_u;
      
      j_l = _phase.size()-1;
      while( _phase[j_l] > phase && j_l > 0) --j_l;
      
      j_u = 0;
      while( _phase[j_u] < phase && j_u < _phase.size()-1) ++j_u;
      
      lambda.resize( _lambda.size() );
      spect.resize( _lambda.size() );
      
      double x = sep - _sep[i_l];
      if(x < 0) x = 0;
      
      double y = (phase - _phase[j_l]);
      if(y < 0) y = 0;
      
      double x0, x1;
            
      for(int i=0; i< spect.size(); ++i)
      {
         lambda[i] = _lambda[i];
       
         x0 = _ag[i_l][j_l][i];
         x1 = _ag[i_u][j_l][i];
         
         if(y != 0 && j_u != j_l)
         {
            x0 += (_ag[i_l][j_u][i] - _ag[i_l][j_l][i])*y/(_phase[j_u] - _phase[j_l]);
            x1 += (_ag[i_u][j_u][i] - _ag[i_u][j_l][i])*y/(_phase[j_u] - _phase[j_l]);
         }
         if( x == 0 || i_u==i_l ) spect[i] = x0;
         else spect[i] = x0 + (x1-x0)*x/( _sep[i_u] - _sep[i_l]);
      }
      
      return 0;
   }

      
};

/// @}
#endif

}

#endif //__astroSpectra_hpp__

