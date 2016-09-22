/** \file psds.hpp
  * \brief Tools for working with PSDs
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup psds
  *
  */

#ifndef __psds_hpp__
#define __psds_hpp__

#include "fft.hpp"
#include "eigenImage.hpp"

namespace mx
{
   
/** \ingroup psds
  * @{
  */

/// Calculates the frequency sampling for a grid given maximum dimension and maximum frequency.
/** The freq_sampling is
  * @f$ \Delta f = f_{max}/ (0.5*dim) @f$
  * where @f$ f_{max} = 1/(2\Delta t) @f$ is the maximum frequency and @f$ dim @f$ is the size of the grid.
  *
  * \param dim is the size of the grid
  * \param f_max is the maximum frequency of the grid
  * \returns the sampling interval @f$ \delta f @f$
  */
template<class realT> 
realT freq_sampling( size_t dim, 
                     realT f_max
                   )
{
   return (f_max/(0.5*dim));
}

///Create a 1-D frequency grid
/**
  * \param vec [out] the pre-allocated Eigen-type 1xN or Nx1 array, on return contains the frequency grid
  * \param dt [in] the temporal sampling of the time series
  * \param inverse [in] [optional] if true
  * 
  * \tparam eigenArr the Eigen array type
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
  * \param fmin is the minimum non-zero absolute value of frequency
  * \param fmax is the maximum absolute value of frequencey
  * \param alpha is the power-law exponent, by convention @f$ \alpha > 0 @f$.
  * 
  * \returns the normalization for a 2-sided power law PSD.
  */
template<typename realT>
realT oneoverf_norm(realT fmin, realT fmax, realT alpha)
{
   realT integ = 2*(pow(fmax, -1.0*alpha + 1.0) - pow(fmin, -1.0*alpha + 1.0))/(-1.0*alpha + 1.0);
   
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
  * \param psd [output] is the array to populate
  * \param freq [input] is a frequency grid, must be the same logical size as psd
  * \param alpha [input] is the power law exponent, by convention @f$ alpha > 0 @f$.
  * \param beta [optional input] is a normalization constant to multiply the raw spectrum by.  If beta==-1 (default) then 
  *                           the PSD is normalized using \ref oneoverf_norm.
  *
  * \tparam eigenArrp is the Eigen array type of the psd 
  * \tparam eigenArrf is the Eigen array type of the frequency grid 
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

/// Generates a von Karman power spectrum
/** 
  * Populates an Eigen array  with
  *
  * \f[
  *  P(f) = \frac{\beta}{ (f^2 + (1/L_0)^2)^{\alpha/2}} 
  * \f]
  * 
  * If you set \f$ L_0 \le 0 \f$ this reverts to standard \f$ 1/f^\alpha \f$ law (i.e. 
  * it treats this as infinite outer scale).  
  *
  * \param psd [output] is the array to populate, allocated.
  * \param freq [input] is a frequency grid, must be the same logical size as psd
  * \param alpha [input] is the power law exponent, by convention @f$ alpha > 0 @f$.
  * \param L0 [input] is the outer scale.  
  * \param beta [optional input] is a normalization constant to multiply the raw spectrum by.  If beta==-1 (default) then 
  *                           the PSD is normalized using \ref oneoverf_norm.
  *
  * \tparam eigenArrp is the Eigen array type of the psd 
  * \tparam eigenArrf is the Eigen array type of the frequency grid 
  */
template<typename eigenArrp, typename eigenArrf> 
void vonKarman_psd( eigenArrp  & psd,
                   eigenArrf & freq,
                   typename eigenArrp::Scalar alpha,
                   typename eigenArrp::Scalar L0,
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
            p = beta / std::pow(   std::pow(freq(ii,jj),2) + L02, sqrt_alpha);
         }
         psd(ii,jj) = p;
      }
   }   
}

///Filter a noise vector or array by the given PSD using Fourier convolution.
/**
  * \param noise [input/output] is the noise array, an Eigen type. On output it is replaced with the filtered noise.
  * \param psd [input] is the 2-sided PSD with which to filter the noise.
  * 
  */
template<typename eigenArrn, typename eigenArrp, typename carrT>
void psd_filter( eigenArrn & noise,
                 eigenArrp & psd ,
                 carrT & ft
               )
{
   typedef std::complex<typename eigenArrn::Scalar> complexT;
   
   //Eigen::Array< complexT, Eigen::Dynamic, Eigen::Dynamic> ft(noise.rows(), noise.cols());
   Eigen::Array< complexT, Eigen::Dynamic, Eigen::Dynamic> ft2(noise.rows(), noise.cols());
   Eigen::Array< complexT, Eigen::Dynamic, Eigen::Dynamic> cnoise(noise.rows(), noise.cols());
   
   //Make noise a complex number
   for(int ii=0;ii<noise.rows();++ii)
   {
      for(int jj=0; jj<noise.cols(); ++jj)
      {
         cnoise(ii,jj) = complexT(noise(ii,jj),0);// noise(ii,jj));
      }
   }
   
   mx::fftT<complexT, complexT, 2, 0> fft(noise.rows(), noise.cols());
   mx::fftT<complexT, complexT, 2, 0> fftR(noise.rows(), noise.cols(), FFTW_BACKWARD);
   
   fft.fft( cnoise.data(),ft.data());
   
   ft *= psd.sqrt();//*exp(complexT(0, 0.5*3.14159));
        
   fftR.fft(ft.data(), ft2.data()); //in-place
   
   noise = ft2.real()/(noise.rows()*noise.cols());
   
}


///Calculate the average periodogram from a time series for a specified averaging interval and overlap.
/**
  * \param pgram [out] the resultant periodogram
  * \param ts [in] the time series
  * \param dt [in] the sampling time of ts
  * \param avgLen [in] the length of the averaging interval, same units as dt
  * \param olap [in] the length of the overlap region, same units as avgLen
  * \param w  [in] a vector of length ( (int) avgLen/dt) or empty, which is a window.
  */
template<typename floatT>
void averagePeriodogram( std::vector<floatT> & pgram, 
                         std::vector<floatT> & ts, 
                         floatT dt, 
                         floatT avgLen, 
                         floatT olap,
                         std::vector<floatT> & w)
{
   int Nper = avgLen/dt;
   int Nover = (avgLen-olap)/dt;
   
   pgram.resize(Nper/2., 0);
   
   std::vector<floatT> work;
   work.resize(Nper);
   
   std::vector<std::complex<floatT> > fftwork;
   fftwork.resize(Nper);
   
   int Navg = ts.size()/Nover;
   
   while(Navg*Nover + Nper > ts.size()) --Navg;

   typename std::vector<floatT>::iterator first, last;
   
   mx::fftT<std::complex<floatT>, std::complex<floatT>, 1, 0> fft(Nper);

   for(int i=0;i<Navg;++i)
   {
      first = ts.begin() + (i*Nover);
      last = first + Nper;
      
      work.assign(first, last);
      floatT mean = 0;
      for (int j=0;j<work.size();++j) mean += work[j];
      //mean = std::accumulate(work.begin(), work.end(), 0);// / work.size();
      mean/=work.size();
      
      
      for(int j=0;j<work.size();++j) 
      {
         work[j] = (work[j] - mean);//
         if(w.size() == Nper) work[j] *= w[j];
      }
            
      fft.fft(work.data(), fftwork.data()); 
       
      for(int j=0;j<pgram.size();++j) pgram[j] += norm(fftwork[j]); //pow(abs(fftwork[j]),2);

   }

   for(int j=0;j<pgram.size();++j) pgram[j] /= (work.size()*Navg);
   
   //fftw_destroy_plan(p);

}





///Augment a 1-sided PSD to standard FFT form.
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
  * Note: entries in psdOneSided are cast to the value_type of psdTwoSided,
  * for instance to allow for conversion to complex type.
  * 
  * \param psdTwoSided [out] on return contains the FFT storage order copy of psdOneSided.
  * \param psdOneSided [in] the one-sided PSD to augment
  * \param zeroFreq [in] [optional] set to a non-zero value if psdOneSided does not contain a zero frequency component.
  */ 
template<typename vectorTout, typename vectorTin, typename dataT>
void augment1SidedPSD( vectorTout & psdTwoSided, 
                       vectorTin  & psdOneSided,
                       dataT zeroFreq = 0 )
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
   
   if( needZero)
   {
      psdTwoSided[0] = outT(0.0);
   }
   else
   {
      psdTwoSided[0] = outT(psdOneSided[0]);
   }
   
   int i;
   for(i=0; i < psdOneSided.size() - 1 - (1-needZero); ++i)
   {
      psdTwoSided[i + 1] = outT(psdOneSided[i + (1-needZero) ]);
      psdTwoSided[i + psdOneSided.size()+ needZero] = outT(psdOneSided[ psdOneSided.size() - 2 - i]);
   }
   psdTwoSided[i + 1] = outT(psdOneSided[i + (1-needZero) ]);

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
  * \param freqTwoSided [out] on return contains the FFT storage order copy of freqOneSided.
  * \param freqOneSided [in] the one-sided frequency scale to augment
  */ 
template<typename T>
void augment1SidedPSDFreq( std::vector<T> & freqTwoSided,
                           std::vector<T> & freqOneSided )
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

} //namespace mx

#endif //__psds_hpp__
