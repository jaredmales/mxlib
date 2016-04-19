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
      basepath = getEnv("PICKLES_ATLAS_DATADIR");
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
      
   mx::readColumns(fname, lambda, trans);
   
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
}

/// Characterize the wavelength of a filter transmission curve.
/** For a transmission curve given by \f$ T(\lambda ) \f$  The central wavelength is defined as
  * \f$ \lambda_0 = \frac{1}{w_{eff}}\int T(\lambda ) \lambda d\lambda \f$
  * where the effective width is defined by
  * \f$ w_{eff} = \int T(\lambda ) d\lambda \f$
  * 
  * The full-width at half-maximum, FWHM, is the distance between the 50% of maximum \f$ T(\lambda) \f$.
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


/// Characterize the fluxes of a filter transmission curve w.r.t. Vega.
/** Uses the HST Calspec spectrum of Vega.
  * 
  * 
  * \param lambda  the wavelength grid of the profile
  * \param trans   the transmission of the profile at each wavelength
  * \param flambda0 the flux of a 0 mag star in ergs/sec/um/cm^2
  * \param fnu0 the flux of a 0 mag star in Jy
  * \param fphot0 the flux of a 0 mag star in photons/sec/m^2
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
   
   std::vector<dataT> trans_terp;
   
   readCalspec("vega", vega_lambda, vega_flambda);
   
   mx::gsl_interpolate(gsl_interp_linear, lambda, trans, vega_lambda, trans_terp);
   
   
   dataT min_l = 1e30;
   dataT max_l = 0;
   
   for(int i=0; i< lambda.size(); ++i)
   {
      if( lambda[i] < min_l) min_l = lambda[i];
      if( lambda[i] > max_l) max_l = lambda[i];
   }
   
   int min_l_i = 0;
   int max_l_i = 0;

   while( vega_lambda[min_l_i] < min_l) ++min_l_i;
   max_l_i = min_l_i;
   while( vega_lambda[max_l_i] <= max_l) ++max_l_i;
   
   dataT tottrans = 0;
   
   for(int i=min_l_i; i< max_l_i-1; ++i) tottrans += trans_terp[i]*(vega_lambda[i+1]-vega_lambda[i]);
   tottrans += trans_terp[max_l_i-1]* (vega_lambda[max_l_i-1] - vega_lambda[max_l_i-2]);
   
   
   // flambda0
   
   flambda0 = 0;
   
   for(int i=min_l_i; i< max_l_i-1; ++i) flambda0 += vega_flambda[i]*trans_terp[i]*(vega_lambda[i+1]-vega_lambda[i]);
   flambda0 += vega_flambda[max_l_i-1]*trans_terp[max_l_i-1]*(vega_lambda[max_l_i-1]-vega_lambda[max_l_i-2]);
   
   flambda0 /= tottrans;
   
   
   // fnu0
   
   fnu0 = 0;
   
   for(int i=min_l_i; i< max_l_i-1; ++i) fnu0 += vega_flambda[i]*3.33564095E+08*pow(vega_lambda[i],2)*trans_terp[i]*(vega_lambda[i+1]-vega_lambda[i]);
   fnu0 += vega_flambda[max_l_i-1]*3.33564095E+08*pow(vega_lambda[max_l_i-1],2)*trans_terp[max_l_i-1]*(vega_lambda[max_l_i-1]-vega_lambda[max_l_i-2]);
   
   fnu0 /= tottrans;
   
   //fphot0
   
   dataT h = boost::units::si::constants::codata::h / boost::units::si::joule/boost::units::si::seconds * 1e7; //ergs
   dataT c = boost::units::si::constants::codata::c / boost::units::si::meter*boost::units::si::seconds;
   
   fphot0 = 0;
 
   //Conversions to wavelength-meters and area-meters^2
   for(int i=min_l_i; i< max_l_i-1; ++i) fphot0 += vega_flambda[i]/( (h*c)/(vega_lambda[i]*1e-6))*1e4*trans_terp[i]*(vega_lambda[i+1]-vega_lambda[i]);
   fphot0 += vega_flambda[max_l_i-1]/( (h*c)/(vega_lambda[max_l_i-1]*1e-6))*1e4*trans_terp[max_l_i-1]*(vega_lambda[max_l_i-1]-vega_lambda[max_l_i-2]);
   
   fphot0 /= tottrans;
   
}




/// @}



} //namespace mx

#endif //__astroFilter_hpp__


