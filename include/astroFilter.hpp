/** \file astroFilter.hpp
  * \author Jared R. Males
  * \brief Utilities for working with astronomical filters and spectra
  * \ingroup astrophot
  *
  */

#ifndef __astroFilter_hpp__
#define __astroFilter_hpp__

#include <iostream>

#include <vector>
#include <string>
#include <cmath>

#include <boost/units/systems/si/codata/universal_constants.hpp>

#include "stringUtils.hpp"
#include "readColumns.hpp"

#include "gslInterpolation.hpp"

#include "environment.hpp"

#include "signalWindows.hpp"

namespace mx
{

   
/** \addtogroup astrophot
  * @{
  */

/// Read a Pickles atlas spectrum
/** Reads in spectrum from the Pickles atlas, usually found via the PICKLES_DATADIR
  * environment variable.  Returned spectra are in wavelength of microns, and arbitrary f_lambda flux units
  * 
  * \param spName  the name of the Calspec star, normally the prefix to the file name.
  * \param lambda the profile wavelength scale, in microns.
  * \param flambda the spectrum, in arbitrary f_lambda flux units
  * \param datadir [optional] the directory containing the Pickles atlas
  * 
  * \tparam dataT the type of data to return
  */
template<typename dataT>
void readPickles( const std::string & spt, 
                  std::vector<dataT> & lambda, 
                  std::vector<dataT> & flambda,
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
  * \param spName  the name of the Calspec star, normally the prefix to the file name.
  * \param lambda the profile wavelength scale, in microns.
  * \param flambda the spectrum, in ergs/s/micron/cm^2
  * \param datadir [optional] the directory containing the Calspec collection
  * 
  * \tparam dataT the type of data to return
  */
template<typename dataT>
void readCalspec( const std::string & spName, 
                  std::vector<dataT> & lambda, 
                  std::vector<dataT> & flambda,
                  const std::string & datadir = "" )
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

/// Read the ATRAN atmospheric transmission spectrum from Gemini
/**
  *
  * \param site either "cp" for Cerro Pachon or "mk" for Mauna Kea
  * \param pwv precipitable water vaper, for "cp" choices are 2.3, 4.3, 7.6, 10.0.  For "mk" choices are 1.0, 1.6, 3.0, 5.0.   
  * \param airmass airmass, choices are 1.0, 1.5, 2.0.
  * \param lambda the wavelength scale
  * \param trans the transmission
  * \param datadir [optional] the directory containing the Gemini ATRAN files
  *
  * \tparam dataT the type of data
  */ 
template<typename dataT>
void readATRAN( const std::string & site, 
                double pwv, 
                double airmass, 
                std::vector<dataT> & lambda, 
                std::vector<dataT> & trans,
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

template<typename dataT>
void readEarthShine( std::vector<dataT> & lambda, 
                     std::vector<dataT> & albedo,
                     const std::string & datadir = "" )
{
   std::string basepath;
   if(datadir == "")
   {
      basepath = getEnv("EARTHSHINE_DATADIR");
   }
   else
   {
      basepath = datadir;
   }
   
   std::string fname; 
   
   fname = basepath + "/F7_opt_NIR_ES_data.txt";
   
   readColumns(fname, lambda, albedo);
}

/// Normalize a transmission profile to have peak value of 1.0
/**
  * \param trans the transmission profile to be normalized. Is altered.
  *
  * \tparam dataT the type of data to return
  */ 
template<typename dataT>
void astrofiltNormalize( std::vector<dataT> & trans )
{
   dataT max = 0;
   for(int i=0; i < trans.size(); ++i) if(trans[i] > max) max = trans[i];
   
   for(int i=0; i < trans.size(); ++i) trans[i] /= max;
}

/// Erase points at beginning and end to trim unwanted data from a spectrum.
/**
  * \param lambda_0 the minimum wavelength to be included
  * \param lambda_f the maximum wavelength to be include
  * \param lambda the wavelength points of the spectrum
  * \param spectrum the spectrum 
  *  
  * \tparam dataT the type of the data 
  */
template<typename dataT>
void astrofiltTrimWavelength( dataT lambda_0,
                              dataT lambda_f,
                              std::vector<dataT> & lambda,
                              std::vector<dataT> & spectrum )
{
   typename std::vector<dataT>::iterator it_lam, it_spect;
   
   it_lam = lambda.begin();
   it_spect = spectrum.begin();
   
   int i =0;
   
   while(lambda[i] < lambda_0 && i < lambda.size())
   {
      ++i;
      ++it_lam;
      ++it_spect;
   }
   
   lambda.erase(lambda.begin(), it_lam);
   spectrum.erase(spectrum.begin(), it_spect);
   
   
   //Now find and erase from lambda_f to end
   it_lam = lambda.begin();
   it_spect = spectrum.begin();
   
   i =0;
   
   while(lambda[i] <= lambda_f && i < lambda.size())
   {
      ++i;
      ++it_lam;
      ++it_spect;
   }
   
   lambda.erase(it_lam, lambda.end());
   spectrum.erase(it_spect, spectrum.end());
   
}

/// Erase points at beginning and end to trim unwanted data from a spectrum.
/**
  * \param lambda_ref the wavelength grid to use as the reference for min and max.
  * \param lambda the wavelength points of the spectrum
  * \param spectrum the spectrum 
  *  
  * \tparam dataT the type of the data 
  */
template<typename dataT>
void astrofiltTrimWavelength( std::vector<dataT> & lambda_ref,
                              std::vector<dataT> & lambda,
                              std::vector<dataT> & spectrum )
{
   astrofiltTrimWavelength( lambda_ref[0], lambda_ref[lambda_ref.size()-1], lambda, spectrum);
}

/// Convert a transmission profile to relative spectral response.
/** The relative spectral response is given by
 * 
  * \f$ T_{RSR}(\lambda) = T(\lambda)\lambda / \mbox{max}[T(\lambda)\lambda] \f$
  * 
  * where \f$ T(\lambda) \f$ is the energy transmission profile.
  * This is appropriate for photon detectors (like CCDs).  See Bessell, PASP, 112, 961 (2000).
  * 
  * \param lambda is the wavelength scale
  * \param trans is the transmission profile
  *  
  * \tparam dataT the type of the data 
  */
template<typename dataT>
void astrofiltRSR( std::vector<dataT> & lambda,
                   std::vector<dataT> & trans )
{
   dataT max = 0;
      
   for(int i=0;i<lambda.size(); ++i)
   {
      trans[i] = trans[i]*lambda[i];
      if(trans[i] > max) max = trans[i];
   }
      
   for(int i=0;i<lambda.size(); ++i)
   {
      trans[i] /= max;
   }
      
}

/// Read a filter transmission curve from the astrofilt collection.
/**
  * Reads in a filter transmission curve from the astrofilt collection, usually found via the ASTROFILT_DATADIR
  * environment variable.  Astrofilt profiles are in wavelength of microns, and transmission normalized to peak 
  * of 1 and airmass 0.
  * 
  * \param filtName  the name of the filter, with no path or extension, e.g. 'J_NIRC2' or 'ip_VisAO'
  * \param lam the profile wavelength scale, in microns.
  * \param trans the transmission at each lambda, normalized to a peak of 1 at 0 airmass.
  * \param rsr  if true converts to relative spectral response, appropriate for photon detectors (like CCDs).  See Bessell, PASP, 112, 961 (2000).
  * \param datadir specifies the data directory.  if not set or empty, the environment is queried for then value of ASTROFILT_DATADIR.
  * 
  * \tparam dataT the type of data to return
  */
template<typename dataT>
void astrofilt( const std::string & filtName, 
                std::vector<dataT> & lambda, 
                std::vector<dataT> & trans,
                bool rsr = false,
                const std::string & datadir = "" )
{
   std::string basepath;
   if(datadir == "")
   {
      basepath = getEnv("ASTROFILT_DATADIR");
   }
   else
   {
      basepath = datadir;
   }
      
   std::string fname = basepath + "/" + filtName + ".dat";
   
   mx::readColumns(fname, lambda, trans);
   
   if(rsr)
   {
      astrofiltRSR(lambda, trans);
   }
}

/// Characterize the wavelength of a filter transmission curve.
/** For a transmission curve given by \f$ T(\lambda ) \f$  The central wavelength is defined as
  * \f$ \lambda_0 = \frac{1}{w_{eff}}\int T(\lambda ) \lambda d\lambda \f$
  * where the effective width is defined by
  * \f$ w_{eff} = \int T(\lambda ) d\lambda \f$
  * 
  * The full-width at half-maximum, FWHM, is the distance between the points at 50% of maximum \f$ T(\lambda) \f$.
  * 
  * \param lambda  the wavelength grid of the profile
  * \param trans   the transmission of the profile at each wavelength
  * \param lambda0 the central wavelength of the filter 
  * \param fwhm the full-width at half-maximum of the filter profile
  * \param weff the effective width of the filter
  *
  * \tparam dataT the data type in which calculations are performed.
  */ 
template<typename dataT>
void astrofiltCharTrans( std::vector<dataT> & lambda,
                         std::vector<dataT> & trans,
                         dataT & lambda0,
                         dataT & fwhm,
                         dataT & weff )
{
   
   weff = 0;
   
   for(int i=0; i< lambda.size()-1; ++i) weff += trans[i]*(lambda[i+1]-lambda[i]);
   weff += trans[trans.size()-1]* (lambda[lambda.size()-1] -lambda[lambda.size()-2]);
      
   lambda0 = 0;
   for(int i=0; i< lambda.size()-1; ++i) lambda0 += lambda[i]*trans[i]*(lambda[i+1]-lambda[i]);
   lambda0 += lambda[trans.size()-1]* trans[trans.size()-1] * (lambda[lambda.size()-1] -lambda[lambda.size()-2]);

   lambda0 /= weff;
   
   
   dataT max, left, right;
   
   max = 0;
   for(int i=0; i< lambda.size(); ++i)
   {
      if(trans[i] > max) max = trans[i];
   }
   
   int i=1;
   while(trans[i] < 0.5*max && i < lambda.size()) ++i;
   left = lambda[i-1] + (0.5 - trans[i-1])* ( (lambda[i] - lambda[i-1]) / (trans[i]-trans[i-1]));
   
   while(trans[i] < max && i < lambda.size()) ++i;
   
   while(trans[i] > 0.5*max && i < lambda.size()) ++i;
   right = lambda[i-1] + (0.5 - trans[i-1])* ( (lambda[i] - lambda[i-1]) / (trans[i]-trans[i-1]));
   
   //std::cout << left << " " << right << "\n";
   
   fwhm = right-left;
   
   
}

/// Characterize the fluxes of a filter transmission curve w.r.t. a spectrum
/** 
  * 
  * \param lambda  the wavelength grid of the profile
  * \param trans   the transmission of the profile at each wavelength
  * \param flambda0 the flux of the star at \f$ \lambda_0 \f$ in ergs/sec/um/cm^2
  * \param fnu0 the flux of the star at \f$ \lambda_0 \f$  in Jy
  * \param fphot0 the flux of the star at \f$ \lambda_0 \f$  in photons/sec/m^2
  * 
  * \tparam dataT the data type in which calculations are performed.
  */ 
template<typename dataT>
void astrofiltCharFlux( std::vector<dataT> & lambda_f,
                        std::vector<dataT> & trans,
                        std::vector<dataT> & lambda_s,
                        std::vector<dataT> & flambda,
                        dataT & flambda0,
                        dataT & fnu0,
                        dataT & fphot0)
{
   //std::vector<dataT> lambda_s;
   //std::vector<dataT> flambda;
   
   std::vector<dataT> trans_terp;
   
   mx::gsl_interpolate(gsl_interp_linear, lambda_f, trans, lambda_s, trans_terp);
      
   dataT min_l = 1e30;
   dataT max_l = 0;
   
   for(int i=0; i< lambda_f.size(); ++i)
   {
      if( lambda_f[i] < min_l) min_l = lambda_f[i];
      if( lambda_f[i] > max_l) max_l = lambda_f[i];
   }
   
   int min_l_i = 0;
   int max_l_i = 0;

   while( lambda_s[min_l_i] < min_l) ++min_l_i;
   max_l_i = min_l_i;
   while( lambda_s[max_l_i] <= max_l) ++max_l_i;
   
   dataT tottrans = 0;
   
   for(int i=min_l_i; i< max_l_i-1; ++i) tottrans += trans_terp[i]*(lambda_s[i+1]-lambda_s[i]);
   tottrans += trans_terp[max_l_i-1]* (lambda_s[max_l_i-1] - lambda_s[max_l_i-2]);
   
   
   // flambda0
   
   flambda0 = 0;
   
   for(int i=min_l_i; i< max_l_i-1; ++i) flambda0 += flambda[i]*trans_terp[i]*(lambda_s[i+1]-lambda_s[i]);
   flambda0 += flambda[max_l_i-1]*trans_terp[max_l_i-1]*(lambda_s[max_l_i-1]-lambda_s[max_l_i-2]);
   
   flambda0 /= tottrans;
   
   
   // fnu0
   
   fnu0 = 0;
   
   for(int i=min_l_i; i< max_l_i-1; ++i) fnu0 += flambda[i]*3.33564095E+08*pow(lambda_s[i],2)*trans_terp[i]*(lambda_s[i+1]-lambda_s[i]);
   fnu0 += flambda[max_l_i-1]*3.33564095E+08*pow(lambda_s[max_l_i-1],2)*trans_terp[max_l_i-1]*(lambda_s[max_l_i-1]-lambda_s[max_l_i-2]);
   
   fnu0 /= tottrans;
   
   //fphot0
   
   dataT h = boost::units::si::constants::codata::h / boost::units::si::joule/boost::units::si::seconds * 1e7; //ergs
   dataT c = boost::units::si::constants::codata::c / boost::units::si::meter*boost::units::si::seconds;
   
   fphot0 = 0;
 
   //Conversions to wavelength-meters and area-meters^2
   for(int i=min_l_i; i< max_l_i-1; ++i) fphot0 += flambda[i]/( (h*c)/(lambda_s[i]*1e-6))*1e4*trans_terp[i]*(lambda_s[i+1]-lambda_s[i]);
   fphot0 += flambda[max_l_i-1]/( (h*c)/(lambda_s[max_l_i-1]*1e-6))*1e4*trans_terp[max_l_i-1]*(lambda_s[max_l_i-1]-lambda_s[max_l_i-2]);
   
   fphot0 /= tottrans;
   
}

/// Characterize the fluxes of a filter transmission curve w.r.t. Vega.
/** Uses the HST Calspec spectrum of Vega.
  * 
  * 
  * \param lambda  the wavelength grid of the profile
  * \param trans   the transmission of the profile at each wavelength
  * \param flambda0 the flux of a 0 mag star at \f$ \lambda_0 \f$  in ergs/sec/um/cm^2
  * \param fnu0 the flux of a 0 mag star at \f$ \lambda_0 \f$ in Jy
  * \param fphot0 the flux of a 0 mag star at \f$ \lambda_0 \f$ in photons/sec/m^2
  * 
  * \tparam dataT the data type in which calculations are performed.
  */ 
template<typename dataT>
void astrofiltCharVega( std::vector<dataT> & lambda,
                        std::vector<dataT> & trans,
                        dataT & flambda0,
                        dataT & fnu0,
                        dataT & fphot0)
{
   std::vector<dataT> vega_lambda;
   std::vector<dataT> vega_flambda;
   
   
   readCalspec("vega", vega_lambda, vega_flambda);
   
   astrofiltCharFlux(lambda, trans, vega_lambda, vega_flambda, flambda0, fnu0, fphot0);
   

   
}

/// Multiply two spectra or filter curves.  
/** Interpolates the 2nd onto the 1st, so that the wavelength scale of the output is lambda_1.
  *
  * \param lambda_1 the first (higher sampling) wavelength scale
  * \param spectrum_1 the first (higher sampling) spectrum
  * \param lambda_2 the second wavelength scale
  * \param spectrum_2 the second spectrum, will be re-sampled onto lambda_1
  * \param spectrum_out the product of spectrum_1 and the re-sampled spectrum_2
  *
  * \tparam dataT the type of the data 
  */ 
template<typename dataT>
void astrofiltMultiply( std::vector<dataT> & lambda_1,
                        std::vector<dataT> & spectrum_1,
                        std::vector<dataT> & lambda_2,
                        std::vector<dataT> & spectrum_2,
                        std::vector<dataT> & spectrum_out )
{
   gsl_interpolate( gsl_interp_linear, lambda_2, spectrum_2, lambda_1, spectrum_out);
   
   for(int i=0; i<spectrum_out.size(); ++i) 
   {
      if( !std::isnormal(spectrum_out[i])) spectrum_out[i] = 0;
      spectrum_out[i] *= spectrum_1[i];
   }
}


// astrofiltScaleFlux

/// Scale a spectrum to a specified Vega magnitude in a filter.
/**
  * Scales the spectrum to have the specified Vega magnitude, using the HST Calspec Vega spectrum.
  * The spectrum can be scaled by \f$ f_\lambda \f$, \f$ f_\nu \f$, or \f$ f_\gamma \f$.
  *
  * \param lambda_s wavelength scale of the spectrum
  * \param spectrum the spectrum to scale, will be be altered
  * \param lambda_f the wavelength scale of the filter
  * \param trans the transmission profile of the filter
  * \param mag the desired magnitude in the filter
  * \param scaleBy specifies which flux measurement to scale by. \f$ 0 = f_\lambda \f$, \f$ 1=f_\nu \f$, or \f$ 2=f_\gamma \f$
  * 
  * \tparam dataT the type of the data
  */ 
template<typename dataT>
void astrofiltScaleFlux( std::vector<dataT> & lambda_s,
                         std::vector<dataT> & spectrum,
                         std::vector<dataT> & lambda_f,
                         std::vector<dataT> & trans,
                         dataT mag,
                         int scaleBy = 0 )
{
   dataT flambda0, fnu0, fphot0;
   dataT flambda, fnu, fphot;
   
   astrofiltCharVega(lambda_f, trans, flambda0, fnu0, fphot0);
   astrofiltCharFlux(lambda_f, trans, lambda_s, spectrum, flambda, fnu, fphot);
   
   dataT scale;
   
   switch(scaleBy)
   {
      case 0:
         scale =  (flambda0/flambda);
         break;
      case 1:
         scale =  (fnu0/fnu);
         break;
      case 2:
         scale =  (fphot0/fphot);
         break;
      default:
         scale =  (flambda0/flambda);
   }
   
   scale *= pow(10., -0.4*mag);
         
   for(int i=0; i < spectrum.size(); ++i) spectrum[i] *= scale;
   
}

template<typename dataT>
void astrofiltSqWave( dataT minLam,
                      dataT maxLam,
                      dataT dLam,
                      std::vector<dataT> & lambda,
                      std::vector<dataT> & trans,
                      dataT alpha = 0.0 )
{
   int N = (maxLam - minLam)/dLam + 1.5;
   
   lambda.resize(N);
   trans.resize(N);
   
   for(int i=0; i<N; ++i) 
   {
      lambda[i] = minLam + i*dLam;
   }
   
   tukey1d( trans.data(), trans.size(), alpha);
   
   //trans[0] = 0;
   //trans[N-1] = 0;
}

/// @}



} //namespace mx

#endif //__astroFilter_hpp__


