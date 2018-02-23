/** \file psdUtils.hpp
  * \brief Tools for working with PSDs
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef psdUtils_hpp
#define psdUtils_hpp

#include "../fft/fft.hpp"
#include "../math/vectorUtils.hpp"

namespace mx
{
namespace sigproc 
{
   
/** \ingroup psds
  * @{
  */

/// Calculates the frequency sampling for a grid given maximum dimension and maximum frequency.
/** The freq_sampling is
  * @f$ \Delta f = f_{max}/ (0.5*dim) @f$
  * where @f$ f_{max} = 1/(2\Delta t) @f$ is the maximum frequency and @f$ dim @f$ is the size of the grid.
  *
  * \param [in] dim is the size of the grid
  * \param [in] f_max is the maximum frequency of the grid
  * 
  * \returns the sampling interval @f$ \delta f @f$
  * 
  * \tparam realT is the real floating point type used for calculations.
  * 
  */
template<class realT> 
realT freq_sampling( size_t dim, 
                     realT f_max )
{
   return (f_max/(0.5*dim));
}

///Create a 1-D frequency grid
/**
  * \param [out] vec the pre-allocated Eigen-type 1xN or Nx1 array, on return contains the frequency grid
  * \param [in] dt the temporal sampling of the time series
  * \param [in] inverse [optional] if true
  * 
  * \tparam eigenArr the Eigen-like array type
  */ 
template<typename eigenArr>
void frequency_grid1D( eigenArr & vec,
                       typename eigenArr::Scalar dt,
                       bool inverse = false )
{
   typename eigenArr::Index dim, dim_1, dim_2;
   typename eigenArr::Scalar df;
   
   dim_1 = vec.rows();
   dim_2 = vec.cols();
   
   dim = std::max(dim_1, dim_2);
   
   df = freq_sampling(dim, 0.5/dt);
      
   if( !inverse )
   {
      for(int ii=0; ii < ceil(0.5*(dim-1) + 1); ++ii)
      {
         vec(ii) = ii*df;
      }
   
      for(int ii=ceil(0.5*(dim-1)+1); ii < dim_1; ++ii)
      {
         vec(ii) = (ii-dim)*df;
      }
   }
   else
   {
      for(int ii=0; ii < dim; ++ii)
      {
         vec(ii) = df * ii / dim;
      }
   }
}

///Create a frequency grid
template<typename eigenArr>
void frequency_grid( eigenArr & arr,
                     typename eigenArr::Scalar dt)
{
   typename eigenArr::Index dim_1, dim_2;
   typename eigenArr::Scalar f_1, f_2, df;
   
   dim_1 = arr.rows();
   dim_2 = arr.cols();
   
   df = freq_sampling(std::max(dim_1, dim_2), 0.5/dt);

   for(int ii=0; ii < 0.5*(dim_1-1) + 1; ++ii)
   {
      f_1 = ii*df;
      for(int jj=0; jj < 0.5*(dim_2-1)+1; ++jj)
      {
         f_2 = jj*df;
         
         arr(ii, jj) = sqrt(f_1*f_1 + f_2*f_2);
      }
                  
      for(int jj=0.5*(dim_2-1)+1; jj < dim_2; ++jj)
      {
         f_2 = (jj-dim_2) * df;
         
         arr(ii, jj) = sqrt(f_1*f_1 + f_2*f_2);
      }
   }
   
   for(int ii=0.5*(dim_1-1)+1; ii < dim_1; ++ii)
   {
      f_1 = (ii-dim_1)*df;
      for(int jj=0; jj < 0.5*(dim_2-1) + 1; ++jj)
      {
         f_2 = jj*df;
         
         arr(ii, jj) = sqrt(f_1*f_1 + f_2*f_2);
      }
                  
      for(int jj=0.5*(dim_2-1)+1; jj < dim_2; ++jj)
      {
         f_2 = (jj-dim_2) * df;
         
         arr(ii, jj) = sqrt(f_1*f_1 + f_2*f_2);
      }
   }
}

///Calculate the normalization for a 1-D @f$ 1/|f|^\alpha @f$ PSD.
/**
  * \param [in] fmin is the minimum non-zero absolute value of frequency
  * \param [in] fmax is the maximum absolute value of frequencey
  * \param [in] alpha is the power-law exponent, by convention @f$ \alpha > 0 @f$.
  * 
  * \returns the normalization for a 2-sided power law PSD.
  * 
  * \tparam realT is the real floating point type used for calculations.
  */
template<typename realT>
realT oneoverf_norm(realT fmin, realT fmax, realT alpha)
{
   realT integ = 2*(pow(fmax, -1.0*alpha + 1.0) - pow(fmin, -1.0*alpha + 1.0))/(-1.0*alpha + 1.0);
   
   return 1/integ;
}

///Calculate the normalization for a 2-D @f$ 1/|k|^\alpha @f$ PSD.
/**
  * \param [in] kmin is the minimum non-zero absolute value of frequency
  * \param [in] kmax is the maximum absolute value of frequencey
  * \param [in] alpha is the power-law exponent, by convention @f$ \alpha > 0 @f$.
  * 
  * \returns the normalization for a 2-D, 2-sided power law PSD.
  * 
  * \tparam realT is the real floating point type used for calculations.
  */
template<typename realT>
realT oneoverk_norm(realT kmin, realT kmax, realT alpha)
{
   realT integ = 2*(pow(kmax, -1*alpha + 2.0) - pow(kmin, -1.0*alpha + 2.0))/(-1*alpha + 2.0);
   
   return 1/integ;
}

/// Generates a @f$ 1/|f|^\alpha @f$ power spectrum
/** 
  * Populates an Eigen array  with
  * \f[
  *  P(|f| = 0) = 0 
  * \f]
  * \f[
  *  P(|f| > 0) = \frac{\beta}{|f|^{\alpha}} 
  * \f]
  * 
  *
  * \param [out] psd is the array to populate
  * \param [in] freq is a frequency grid, must be the same logical size as psd
  * \param [in] alpha is the power law exponent, by convention @f$ alpha > 0 @f$.
  * \param [in] beta [optional is a normalization constant to multiply the raw spectrum by.  If beta==-1 (default) then 
  *                           the PSD is normalized using \ref oneoverf_norm.
  *
  * \tparam eigenArrp is the Eigen-like array type of the psd 
  * \tparam eigenArrf is the Eigen-like array type of the frequency grid 
  */
template<typename eigenArrp, typename eigenArrf> 
void oneoverf_psd( eigenArrp  & psd,
                   eigenArrf & freq,
                   typename eigenArrp::Scalar alpha,
                   typename eigenArrp::Scalar beta = -1 )
{   
   typedef typename eigenArrp::Scalar Scalar;
   
   typename eigenArrp::Index dim_1, dim_2;
   Scalar f_x, f_y, p;
   
   dim_1 = psd.rows();
   dim_2 = psd.cols();
   
   
   
   if(beta==-1)
   {
      Scalar fmin;
      Scalar fmax;
      
      fmax = freq.abs().maxCoeff();
      
      //Find minimum non-zero Coeff.
      fmin = (freq.abs()>0).select(freq.abs(), freq.abs()+fmax).minCoeff();
      
      beta = oneoverf_norm(fmin, fmax, alpha);
   }
   
   for(int ii =0; ii < dim_1; ++ii)
   {
      for(int jj=0; jj < dim_2; ++jj)
      {
         if(freq(ii,jj) == 0)
         {
            p = 0;
         }
         else
         {
            p = beta / std::pow(std::abs(freq(ii,jj)), alpha);
         }
         psd(ii,jj) = p;
      }
   }   
}
/// Generate a 1-D von Karman power spectrum
/** 
  * Populates an Eigen array  with
  *
  * \f[
  *  P(f) = \frac{\beta}{ (f^2 + (1/T_0)^2)^{\alpha/2}} e^{ - f^2 t_0^2}
  * \f]
  * 
  * If you set \f$ T_0 \le 0 \f$ and \f$ t_0 = 0\f$ this reverts to a simple \f$ 1/f^\alpha \f$ law (i.e. 
  * it treats this as infinite outer scale and inner scale).  
  *
  * \tparam floatT a floating point 
  */
template<typename floatT> 
int vonKarmanPSD( std::vector<floatT> & psd, ///< [out] the PSD vector, will be resized.
                  std::vector<floatT> & f,   ///< [in] the frequency vector
                  floatT beta,               ///< [in] the scaling constant
                  floatT T0,                 ///< [in] the outer scale 
                  floatT t0,                 ///< [in] the inner scale
                  floatT alpha               ///< [in] the exponent, by convention @f$ alpha > 0 @f$.
                )
{   
   
   
   floatT T02;
   if(T0 > 0) T02 = 1.0/(T0*T0);
   else T02 = 0;
   
   floatT sqrt_alpha = 0.5*alpha;
   
   psd.resize(f.size());
   
   for(int i=0; i< f.size(); ++i)
   {
      floatT p = beta / pow( pow(f[i],2) + T02, sqrt_alpha);
      if(t0 > 0 ) p *= exp(-1*pow( f[i]*t0, 2)); 
      psd[i] = p;
   }   
   
   return 0;
}

/// Generates a von Karman power spectrum
/** 
  * Populates an Eigen array  with
  *
  * \f[
  *  P(k) = \frac{\beta}{ (k^2 + (1/L_0)^2)^{\alpha/2}} e^{ - k^2 l_0^2}
  * \f]
  * 
  * If you set \f$ L_0 \le 0 \f$ and \f$ l_0 = 0\f$ this reverts to a simple \f$ 1/f^\alpha \f$ law (i.e. 
  * it treats this as infinite outer scale and inner scale).  
  *
  * \param [out] psd is the array to populate, allocated.
  * \param [in] freq is a frequency grid, must be the same logical size as psd
  * \param [in] alpha is the power law exponent, by convention @f$ alpha > 0 @f$.
  * \param [in] L0 [optional] is the outer scale. 
  * \param [in] l0 [optional] is the inner scale. 
  * \param [in] beta [optional] is a normalization constant to multiply the raw spectrum by.  If beta==-1 (default) then 
  *                           the PSD is normalized using \ref oneoverf_norm.
  *
  * \tparam eigenArrp is the Eigen array type of the psd 
  * \tparam eigenArrf is the Eigen array type of the frequency grid 
  */
template<typename eigenArrp, typename eigenArrf> 
void vonKarman_psd( eigenArrp  & psd,
                   eigenArrf & freq,
                   typename eigenArrp::Scalar alpha,
                   typename eigenArrp::Scalar L0 = 0,
                   typename eigenArrp::Scalar l0 = 0,
                   typename eigenArrp::Scalar beta = -1 )
{   
   typedef typename eigenArrp::Scalar Scalar;
   
   typename eigenArrp::Index dim_1, dim_2;
   Scalar f_x, f_y, p;
   
   dim_1 = psd.rows();
   dim_2 = psd.cols();
   
   
   
   if(beta==-1)
   {
      Scalar fmin;
      Scalar fmax;
      
      fmax = freq.abs().maxCoeff();
      
      //Find minimum non-zero Coeff.
      fmin = (freq.abs()>0).select(freq.abs(), freq.abs()+fmax).minCoeff();
      
      beta = oneoverf_norm(fmin, fmax, alpha);
   }
   
   
   Scalar L02;
   if(L0 > 0) L02 = 1.0/(L0*L0);
   else L02 = 0;
   
   Scalar sqrt_alpha = 0.5*alpha;//std::sqrt(alpha);
   
   for(int ii =0; ii < dim_1; ++ii)
   {
      for(int jj=0; jj < dim_2; ++jj)
      {
         if(freq(ii,jj) == 0 && L02 == 0)
         {
            p = 0;
         }
         else
         {
            p = beta / pow( pow(freq(ii,jj),2) + L02, sqrt_alpha);
            if(l0 > 0 ) p *= exp(-1*pow( freq(ii,jj)*l0, 2)); 
         }
         psd(ii,jj) = p;
      }
   }   
}


   
///Calculate the average periodogram from a time series for a specified averaging interval and overlap.
/** The time series should be mean subtracted before passing to this function.
  * 
  * The frequency scale of the output periodogram is 1/(2*pgram.size()*dt), where the factor of 2 is due to the one-sided-ness of the result.
  * 
  * If a window is supplied, the PSD is normalized so that the one-sided integral is equal to the variance of the input time-series.
  *
  * 
  */
template<typename realT>
void averagePeriodogram( std::vector<realT> & pgram,  ///< [out] the resultant periodogram.
                         std::vector<realT> & ts,  ///< [in] is input the time-series.
                         realT dt,   ///< [in] is the sampling time of time-series.
                         realT avgLen,    ///< [in] is the length of the averaging interval, same units as dt.
                         realT olap,    ///< [in] is the length of the overlap region, same units as avgLen.
                         std::vector<realT> & w   ///< [in] a vector of length ( (int) avgLen/dt) containing a window.  If empty then then the square window is used.
                       ) 
{
   int Nper = avgLen/dt;
   int Nover = (avgLen-olap)/dt;

   pgram.resize(Nper/2., 0);
   
   std::vector<std::complex<realT> > fftwork, cwork;
   
   fftwork.resize(Nper);
   cwork.resize(Nper);
   
   int Navg = ts.size()/Nover;
   
   while(Navg*Nover + Nper > ts.size()) --Navg;

   if(Navg < 1) Navg = 1; //Always do at least 1!
   
   mx::fftT<std::complex<realT>, std::complex<realT>, 1, 0> fft(Nper);

   for(int i=0;i<Navg;++i)
   {
      realT v;
      
      for(int j=0;j<Nper;++j) 
      {
         v = ts[i*Nover + j];
 
         if(w.size() == Nper) v *= w[j];
         
         cwork[j] = std::complex<realT>(v, 0);
      }
         
      fft( fftwork.data(), cwork.data()); 
       
      for(int j=0;j<pgram.size();++j) pgram[j] += norm(fftwork[j]); //pow(abs(fftwork[j]),2);
   }

   realT varNorm = 1;
   if(w.size() == Nper)
   {
      realT sum = 0;
      realT df = 1.0/(2.0*pgram.size()*dt); //factor of 2 since it's a one-sided PSD 
      for(int j =0; j< pgram.size(); ++j) sum += pgram[j]*df;
      varNorm = mx::math::vectorVariance(ts) / sum;
   }
   
   for(int j=0;j<pgram.size();++j) pgram[j] *= varNorm / (Nper*Navg);
}

///Augment a 1-sided PSD to standard 2-sided FFT form.
/** Allocates psdTwoSided to hold a flipped copy of psdOneSided.
  * Default assumes that psdOneSided[0] corresponds to 0 frequency,
  * but this can be changed by setting zeroFreq to a non-zero value.
  * In this case psdTwoSided[0] is set to 0, and the augmented psd
  * is shifted by 1.
  * 
  * Example:
  * 
  * {1,2,3,4,5} --> {0,1,2,3,4,5,-4,-3,-2,-1}
  * 
  * Entries in psdOneSided are cast to the value_type of psdTwoSided,
  * for instance to allow for conversion to complex type.
  *
  */ 
template<typename vectorTout, typename vectorTin>
void augment1SidedPSD( vectorTout & psdTwoSided, ///< [out] on return contains the FFT storage order copy of psdOneSided.
                       vectorTin  & psdOneSided, ///< [in] the one-sided PSD to augment
                       bool zeroFreq = false,       ///< [in] [optional] set to true if psdOneSided does not contain a zero frequency component.
                       typename vectorTin::value_type scale = 1  ///< [in] [optional] value to scale the input by when copying to the output.  Use 0.5 if power needs to be re-normalized for 2-sided PSD.
                     )
{
   typedef typename vectorTout::value_type outT;
   
   int needZero = 1;
   
   size_t N; 
   
   if( zeroFreq == 0 ) 
   {
      needZero = 0;
      N = 2*psdOneSided.size()-2;
   }
   else
   {
      N = 2*psdOneSided.size();
   }
   
   psdTwoSided.resize(N);
   
   //First set the 0-freq point
   if( needZero)
   {
      psdTwoSided[0] = outT(0.0);
   }
   else
   {
      psdTwoSided[0] = outT(psdOneSided[0] * scale);
   }
   
   //Now set all the rest.
   int i;
   for(i=0; i < psdOneSided.size() - 1 - (1-needZero); ++i)
   {
      psdTwoSided[i + 1] = outT(psdOneSided[i + (1-needZero)] * scale);
      psdTwoSided[i + psdOneSided.size()+ needZero] = outT(psdOneSided[ psdOneSided.size() - 2 - i] * scale);
   }
   psdTwoSided[i + 1] = outT(psdOneSided[i + (1-needZero) ] * scale);

}
   
///Augment a 1-sided frequency scale to standard FFT form.
/** Allocates freqTwoSided to hold a flipped copy of freqOneSided.
  * If freqOneSided[0] is not 0, freqTwoSided[0] is set to 0, and the augmented
  * frequency scale is shifted by 1.
  * 
  * Example:
  * 
  * {1,2,3,4,5} --> {0,1,2,3,4,5,-4,-3,-2,-1}
  * 
  */ 
template<typename T>
void augment1SidedPSDFreq( std::vector<T> & freqTwoSided, ///< [out] on return contains the FFT storage order copy of freqOneSided.
                           std::vector<T> & freqOneSided  ///< [in] the one-sided frequency scale to augment
                         )
{
   int needZero = 1;
   
   size_t N; 
   
   if(freqOneSided[0] != 0) 
   {
      N = 2*freqOneSided.size();
   }
   else
   {
      needZero = 0;
      N = 2*freqOneSided.size()-2;
   }
   
   freqTwoSided.resize(N);
   
   if( needZero)
   {
      freqTwoSided[0] = 0.0;
   }
   else
   {
      freqTwoSided[0] = freqOneSided[0]; //0
   }
   
   int i;
   for(i=0; i < freqOneSided.size() - 1 - (1-needZero); ++i)
   {
      freqTwoSided[i + 1] = freqOneSided[i + (1-needZero) ];
      freqTwoSided[i + freqOneSided.size()+ needZero] = -freqOneSided[ freqOneSided.size() - 2 - i];
   }
   freqTwoSided[i + 1] = freqOneSided[i + (1-needZero) ];

}

///@}

} //namespace sigproc 
} //namespace mx

#endif //psdUtils_hpp
