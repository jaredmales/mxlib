/** \file astroSpectrum.hpp
  * \author Jared R. Males
  * \brief A class for working with astronomical spectra.
  * \ingroup astrophot
  *
  */

#ifndef mx_astro_astroSpectrum_hpp
#define mx_astro_astroSpectrum_hpp

#include <cmath>


#include "../sys/environment.hpp"
#include "../ioutils/readColumns.hpp"
#include "../math/gslInterpolation.hpp"


namespace mx
{
namespace astro
{

///Base spectrum class which provides manipulation and characterization functionality.
/**
  * \ingroup astrophot
  */
template<typename realT>
struct baseSpectrum
{
   std::vector<realT> _spectrum; ///< Contains the spectrum after it is set.

   /// Get the current size of the spectrum
   /**
     * \returns the number of points in the spectrum
     */ 
   size_t size()
   {
      return _spectrum.size();
   }
   
   ///Access a single point in the spectrum, specified by its vector index.
   /**
     * \returns a reference to the value at i.
     */
   realT & operator[](size_t i /**< [in] the index of the spectral point*/)
   {
      return _spectrum[i];
   }

   ///Access a single point in the spectrum, specified by its vector index.
   /**
     * \returns the value of point i.
     */
   const realT operator[](size_t i /**< [in] the index of the spectral point*/) const
   {
      return _spectrum[i];
   }

   ///Multiply two spectra together.
   /**
     * \returns a new baseSpectrum
     */
   template<typename compSpectrumT>
   baseSpectrum<realT> operator*( const compSpectrumT & spec /**< [in] the spectrum to multiply by */ )
   {
      baseSpectrum outVec;
      outVec._spectrum.resize( _spectrum.size() );

      for(int i=0; i< _spectrum.size(); ++i)
      {
         outVec[i] = _spectrum[i] * spec[i];
      }

      return outVec;
   }

   ///Calculate the mean value of the spectrum.
   /**
     * \returns the mean value of the spectrum
     */
   realT mean()
   {
      realT sum = 0;

      for(int i=0; i< _spectrum.size(); ++i) sum += _spectrum[i];

      return sum/_spectrum.size();
   }

   ///Calculate the mean value of the spectrum when mutiplied by another.
   /** The result is normalized by the mean value of the input spectrum, equivalent to:
     * \f[
      \mu = \frac{\int T(\lambda)S(\lambda) d\lambda}{\int T(\lambda) d\lambda}
      \f]
     * For instance, use this to get the mean value of a spectrum in a filter.
     *
     * \returns the mean value of the multiplied spectrum.
     *
     * \tparam compSpectrumT the vector-like type of the comparison spectrum.
     */
   template<typename compSpectrumT>
   realT mean( const compSpectrumT & T /**< [in] the spectrum to multiply by*/ )
   {
      realT sum1 = 0, sum2 = 0;

      for(int i=0; i< _spectrum.size(); ++i)
      {
         sum1 += _spectrum[i]*T[i];
         sum2 += T[i];
      }

      return sum1/sum2;
   }


   /// Characterize the spectrum as a filter transmission curve.
   /** For a transmission curve given by \f$ T(\lambda ) \f$  The central wavelength is defined as
     * \f[ 
       \lambda_0 = \frac{1}{w_{eff}}\int \frac{T(\lambda )}{T_{max}} \lambda d\lambda 
       \f]
     * where the effective width is defined by
       \f[ 
        w_{eff} = \int \frac{T(\lambda )} { T_{max}} d\lambda 
        \f]
     *
     * The full-width at half-maximum, FWHM, is the distance between the points at 50% of maximum \f$ T(\lambda) \f$.
     *
     */
   void charTrans( realT & lambda0, ///< [out] the central wavelength of the filter
                   realT & weff, ///< [out] the effective width of the filter
                   realT & max, ///< [out] the maximum value of the transmission curve
                   realT & fwhm, ///< [out] the full-width at half-maximum of the filter profile
                   std::vector<realT> & lambda ///< [in] the wavelength scale, should correspond to the spectrum.
                 )
   {
      size_t N = lambda.size();
      realT dl = lambda[1] - lambda[0];
      realT half = static_cast<realT>(0.5);
      
      
      max = 0;
      for(int i=0; i< N; ++i)
      {
         if(_spectrum[i] > max) max = _spectrum[i];
      }
      
      weff = half * _spectrum[0];
      for(int i=1; i< N-1; ++i) weff += _spectrum[i];
      weff += half * _spectrum[N-1];
      weff *= dl/max;
      
      lambda0 = half*lambda[0]*_spectrum[0];
      for(int i=1; i< N-1; ++i) lambda0 += lambda[i]*_spectrum[i];
      lambda0 += half*lambda[N-1]*_spectrum[N-1];
      lambda0 *= dl/max/weff;

      realT left, right;

      int i=1;
      while(i < lambda.size())
      {
         if(_spectrum[i] >= 0.5*max) break;
         ++i;
      }

      //Interpolate
      left = lambda[i-1] + (0.5*max - _spectrum[i-1])* ( dl / (_spectrum[i]-_spectrum[i-1]));

      i = N-2;
      while( _spectrum[i] < 0.5*max )
      {
         --i;
         if(i < 0 ) break;
      }

      right = lambda[i] + (_spectrum[i]-0.5*max)* ( dl / (_spectrum[i]-_spectrum[i+1]));

      fwhm = right-left;
   }

   /// Characterize the flux densities of the spectrum w.r.t. a filter transmission curve
   /** To obtain the flux (e.g. W/m^2) multiply these quantities by the effective width calculated using
     * \ref charTrans.
     * 
     * \warning this only produces correct fphot0 for a spectrum in W/m^3.  DO NOT USE FOR ANYTHING ELSE.
     *
     * \todo use unit conversions to make it work for everything.
     * \todo check on integration method, should it be trap?
     *
     */
   void charFlux( realT & flambda0, ///< [out] the flux of the star at \f$ \lambda_0 \f$ in W/m^3
                  realT & fnu0,     ///< [out] the flux of the star at \f$ \lambda_0 \f$  in Jy 
                  realT & fphot0,   ///< [out] the flux of the star at \f$ \lambda_0 \f$  in photons/sec/m^3
                  std::vector<realT> & lambda, ///< [in] the wavelength scale of this spectrum.
                  std::vector<realT> & trans  ///< [in] the filter transmission curve over which to characterize, on the same wavelength grid.
                )
   {
      size_t N = lambda.size();
      realT dl = lambda[1] - lambda[0];
      realT half = static_cast<realT>(0.5);
      
      realT tottrans = half*trans[0];
      for(int i=1; i< N-1; ++i) tottrans += trans[i];
      tottrans += half*trans[N-1];
      
      tottrans *= dl;
      
      // flambda0

      flambda0 = half*_spectrum[0]*trans[0];
      for(int i=1; i< N-1; ++i) flambda0 += _spectrum[i]*trans[i];
      flambda0 += half*_spectrum[N-1]*trans[N-1];
      flambda0  *= dl/tottrans;


      // fnu0

      fnu0 = half*_spectrum[0]*pow(lambda[0],2)*trans[0];
      for(int i=1; i< N-1; ++i) fnu0 += _spectrum[i]*pow(lambda[i],2)*trans[i];
      fnu0 += half*_spectrum[N-1]*pow(lambda[N-1],2)*trans[N-1];
      
      fnu0 *= 3.33564095E+17 * dl / tottrans;
      
      //fphot0

      constexpr realT h = constants::h<units::si<realT>>();
      constexpr realT c = constants::c<units::si<realT>>();

      fphot0 = half*_spectrum[0] * lambda[0]*trans[0];

      for(int i=1; i< N-1; ++i) fphot0 += _spectrum[i] * lambda[i] *trans[i];
      fphot0 += half*_spectrum[N-1] * lambda[N-1]*trans[N-1];

      fphot0 *= dl / (h*c*tottrans);

   }
};

///Class to manage an astronomical spectrum.
/** Details of the spectra, including units, file-location, and how they are read from disk, are specified by template parameter _spectrumT:
  * see \ref astrophot_spectra "the available spectrum definitions."
  * Spectra are loaded and immediately interpolated onto a wavelength grid.  This is to facilitate manipulation of spectra, i.e.
  * for multiplying by filters, etc.
  *
  * Inherits functions and operators to manipulate and characterize spectra from \ref mx::astro::baseSpectrum "baseSpectrum".
  *
  * \tparam _spectrumT is the underlying spectrum type, which provides (through static interfaces) the specifications of how to read or calculate the spectrum.
  * \tparam freq specify whether this spectrum uses frequency instead of wavelength.  Default is false.
  *
  * \ingroup astrophot
  */
template<typename _spectrumT, bool freq=false>
struct astroSpectrum : public baseSpectrum<typename _spectrumT::units::realT>
{
   typedef _spectrumT spectrumT;
   typedef typename spectrumT::units units;
   typedef typename units::realT realT;
   typedef typename spectrumT::paramsT paramsT;

   std::string _dataDir; ///< The directory containing the spectrum

   paramsT _params; ///< The parameters of the spectrum, e.g. its name or the numerical parameters needed to generate the name.

   ///Default c'tor
   astroSpectrum(){}

   ///Constructor specifying name, the enivronment will be queried for data directory.
   explicit astroSpectrum( const paramsT & params /**< [in] The name of the spectrum */)
   {
      setParameters(params);
   }

   ///Constructor specifying name and data directory.
   astroSpectrum( const paramsT & params, ///< [in] The name of the spectrum
                  const std::string & dataDir ///< [in] The directory containing the spectrum
                )
   {
      setParameters(params, dataDir);
   }

   ///Set the parameters of the spectrum, using the underlying spectrums parameter type.
   int setParameters( const paramsT & params /**< [in] The name of the spectrum */)
   {
      _params = params;

      if(spectrumT::dataDirEnvVar)
      {
         _dataDir = sys::getEnv(spectrumT::dataDirEnvVar);
      }

      return 0;
   }

   ///Set the parameters of the spectrum, using the underlying spectrums parameter type.
   /** This version also sets the data directory, instead of using the enivronment variable.
     *
     * \overload
     */
   int setParameters( const paramsT & params,  ///< [in] The name of the spectrum
                      const std::string & dataDir ///< [in] The directory containing the spectrum
                    )
   {
      _params = params;
      _dataDir = dataDir;

      return 0;
   }

   ///Load the spectrum and interpolate it onto a wavelength scale.
   template<typename gridT>
   int setSpectrum( gridT & lambda )
   {
      std::vector<realT> rawLambda;
      std::vector<realT> rawSpectrum;

      if(_dataDir == "") _dataDir = ".";

      std::string path = _dataDir + "/" + spectrumT::fileName(_params);

      if( spectrumT::readSpectrum(rawLambda, rawSpectrum, path, _params) < 0)
      {
         return -1; ///\returns -1 on an error reading the spectrum.
      }

      //Unit conversions
      for(int i=0; i < rawLambda.size(); ++i)
      {
         rawLambda[i] /= spectrumT::wavelengthUnits/(units::length);
         rawSpectrum[i] /= spectrumT::fluxUnits/(units::energy/(units::time * units::length * units::length * units::length));
      }

      mx::gsl_interpolate(gsl_interp_linear, rawLambda, rawSpectrum, lambda, this->_spectrum);

      for(int i=0; i < lambda.size(); ++i)
      {
         if( !std::isnormal(this->_spectrum[i])) this->_spectrum[i] = 0;
      }

      return 0; ///\returns 0 on success.
   }


};




} //namespace astro

} //namespace mx





#endif //mx_astro_astroSpectrum_hpp
