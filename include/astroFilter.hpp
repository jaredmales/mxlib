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

#include "astroSpectrum.hpp"

namespace mx
{

   
//Forward Decls.
template<typename dataT>
void readCalspec( std::vector<dataT> & lambda, 
                  std::vector<dataT> & flambda,
                  const std::string & spName, 
                  const std::string & datadir);

template<typename dataT>
void readCalspec( std::vector<dataT> & lambda, 
                  std::vector<dataT> & flambda,
                  const std::string & spName );


/** \addtogroup astrophot
  * @{
  */

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
  * \param lambda [out] the wavelength points of the spectrum
  * \param spectrum [out] the spectrum 
  * \param lambda_0 [in] the minimum wavelength to be included
  * \param lambda_f [in] the maximum wavelength to be include
  *  
  * \tparam dataT the type of the data 
  */
template<typename dataT>
void astrofiltTrimWavelength( std::vector<dataT> & lambda,
                              std::vector<dataT> & spectrum,
                              dataT lambda_0,
                              dataT lambda_f )
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
  * 
  * \param lambda [out] the wavelength points of the spectrum
  * \param spectrum [out] the spectrum 
  * \param lambda_ref [in] the wavelength grid to use as the reference for min and max.
  * 
  * \tparam dataT the type of the data 
  */
template<typename dataT>
void astrofiltTrimWavelength( std::vector<dataT> & lambda,
                              std::vector<dataT> & spectrum,
                              std::vector<dataT> & lambda_ref )
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
  * \param lambda [in] The wavelength scale
  * \param trans [in] The transmission profile
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
  * \param lam [out] the profile wavelength scale, in microns.  Is cleared.
  * \param trans [out] the transmission at each lambda, normalized to a peak of 1 at 0 airmass. Is cleared.
  * \param filtName  [in] the name of the filter, with no path or extension, e.g. 'J_NIRC2' or 'ip_VisAO'
  * \param rsr [in] [optional] if true converts to relative spectral response, appropriate for photon detectors (like CCDs).  See Bessell, PASP, 112, 961 (2000).
  * \param datadir [in] [optional] specifies the data directory.  if not set or empty, the environment is queried for then value of ASTROFILT_DATADIR.
  * 
  * \tparam dataT the type of data to return
  */
template<typename dataT>
void astrofilt( std::vector<dataT> & lambda, 
                std::vector<dataT> & trans,
                const std::string & filtName, 
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
   
      
   lambda.clear();
   trans.clear();
   
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
  * \param lambda0 [out] the central wavelength of the filter 
  * \param fwhm [out] the full-width at half-maximum of the filter profile
  * \param weff [out] the effective width of the filter
  * \param lambda [in] the wavelength grid of the profile
  * \param trans  [in] the transmission of the profile at each wavelength
  *
  * \tparam dataT the data type in which calculations are performed.
  */ 
template<typename dataT>
void astrofiltCharTrans( dataT & lambda0,
                         dataT & fwhm,
                         dataT & weff,
                         std::vector<dataT> & lambda,
                         std::vector<dataT> & trans )
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
  * \param flambda0 [out] the flux of the star at \f$ \lambda_0 \f$ in ergs/sec/um/cm^2
  * \param fnu0 [out] the flux of the star at \f$ \lambda_0 \f$  in Jy
  * \param fphot0 [out] the flux of the star at \f$ \lambda_0 \f$  in photons/sec/m^2
  * \param lambda_f [in] the wavelength grid of the filter
  * \param trans  [in] the transmission of the filter at each wavelength
  * \param lambda_s [in] the wavelength scale of the spectrum
  * \param flambda [in] the spectrum
  * 
  * \tparam dataT the data type in which calculations are performed.
  */ 
template<typename dataT>
void astrofiltCharFlux( dataT & flambda0,
                        dataT & fnu0,
                        dataT & fphot0,
                        std::vector<dataT> & lambda_f,
                        std::vector<dataT> & trans,
                        std::vector<dataT> & lambda_s,
                        std::vector<dataT> & flambda)
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
  * \param flambda0 [out] the flux of a 0 mag star at \f$ \lambda_0 \f$  in ergs/sec/um/cm^2
  * \param fnu0 [out] the flux of a 0 mag star at \f$ \lambda_0 \f$ in Jy
  * \param fphot0 [out] the flux of a 0 mag star at \f$ \lambda_0 \f$ in photons/sec/m^2
  * \param lambda [in] the wavelength grid of the profile
  * \param trans [in] the transmission of the profile at each wavelength
  * 
  * 
  * \tparam dataT the data type in which calculations are performed.
  */ 
template<typename dataT>
void astrofiltCharVega( dataT & flambda0,
                        dataT & fnu0,
                        dataT & fphot0,
                        std::vector<dataT> & lambda,
                        std::vector<dataT> & trans )
{
   std::vector<dataT> vega_lambda;
   std::vector<dataT> vega_flambda;
   
   
   readCalspec(vega_lambda, vega_flambda, "vega");
   
   astrofiltCharFlux(flambda0, fnu0, fphot0, lambda, trans, vega_lambda, vega_flambda);
   
   
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

/// Multiply a spectrum or filter curve by a scalar.  
/** 
  *
  * \param spectrum the spectrum
  * \param scale the scale to multiply by
  * \param spectrum_out the product of spectrum_1 and the re-sampled spectrum_2
  *
  * \tparam dataT the type of the data 
  */ 
template<typename dataT>
void astrofiltMultiply( std::vector<dataT> & spectrum,
                        dataT scale,
                        std::vector<dataT> & spectrum_out )
{
   spectrum_out.resize(spectrum.size());
   
   for(int i=0; i<spectrum_out.size(); ++i) 
   {
      spectrum_out[i] * spectrum[i] * scale;
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
void astrofiltScaleMag( std::vector<dataT> & lambda_s,
                        std::vector<dataT> & spectrum,
                        std::vector<dataT> & lambda_f,
                        std::vector<dataT> & trans,
                        dataT mag,
                        int scaleBy = 0 )
{
   dataT flambda0, fnu0, fphot0;
   dataT flambda, fnu, fphot;
   
   astrofiltCharVega(flambda0, fnu0, fphot0, lambda_f, trans);
   astrofiltCharFlux(flambda, fnu, fphot, lambda_f, trans, lambda_s, spectrum);
   
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

/// Scale a spectrum to a specified mean value in a filter.
/**
  * Scales the spectrum to have the specified mean value.
  * The spectrum can be scaled by \f$ f_\lambda \f$, \f$ f_\nu \f$, or \f$ f_\gamma \f$.
  *
  * \param lambda_s wavelength scale of the spectrum
  * \param spectrum the spectrum to scale, will be be altered
  * \param lambda_f the wavelength scale of the filter
  * \param trans the transmission profile of the filter
  * \param mean the desired mean value in the filter
  * \param scaleBy specifies which flux measurement to scale by. \f$ 0 = f_\lambda \f$, \f$ 1=f_\nu \f$, or \f$ 2=f_\gamma \f$
  * 
  * \tparam dataT the type of the data
  */ 
template<typename dataT>
void astrofiltScaleMean( std::vector<dataT> & lambda_s,
                         std::vector<dataT> & spectrum,
                         std::vector<dataT> & lambda_f,
                         std::vector<dataT> & trans,
                         dataT mean,
                         int scaleBy = 0 )
{
   dataT flambda, fnu, fphot;
   
   astrofiltCharFlux(lambda_f, trans, lambda_s, spectrum, flambda, fnu, fphot);
   
   dataT scale;
   
   switch(scaleBy)
   {
      case 0:
         scale =  (mean/flambda);
         break;
      case 1:
         scale =  (mean/fnu);
         break;
      case 2:
         scale =  (mean/fphot);
         break;
      default:
         scale =  (mean/flambda);
   }
            
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


