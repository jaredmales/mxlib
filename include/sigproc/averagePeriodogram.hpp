/** \file averagePeriodogram.hpp
  * \brief A class to manage calculation of periodograms from timeseries data.
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef averagePeriodogram_hpp
#define averagePeriodogram_hpp

#include "psdUtils.hpp"


namespace mx
{
namespace sigproc
{

/** \ingroup psds
  * @{
  */


/// Calculate the average periodogram of a time-series.
/**
  * Implements the overlapped average periodogram, cf. pp 843 to 845 of \cite 10.5555.1795494. 
  * Optionally includes a window (by default no window is used, which is equivalent to the rectangle window). 
  * Can also be used to calculate the unaveraged non-overlapped periodogram of the entire time-series.
  * 
  * \todo this is currently a mess.  We don't need dt.  If avgLen is the same as sz, or 0 it is size, then just do the non-overlapped case.
  */
template<typename realT>
class averagePeriodogram
{
   
protected:
   
   /** \name Configuration Parameters
     *
     *@{
     */
   
   size_t m_avgLen {0};
   size_t m_overlap {0};
   
   realT m_dt {1}; ///< The time sampling.  Only used for normalization and calculation of the frequency scale.

   std::vector<realT> m_win;
      
   ///@}
   
   size_t m_size {0} ///< The size of the periodogram vector, calculated as m_avgLen/2 + 1;
   
   int m_nOver {0}; ///< The number of overlapping segments.  Calculated from m_avgLen and m_overlap;
   
   realT m_df {1};

   
   mx::fftT< realT, std::complex<realT>, 1, 0> m_fft;
   
   realT * m_tsWork {nullptr};
   size_t m_tsWorkSize {0};
   
   std::complex<realT> * m_fftWork {nullptr};
   size_t m_fftWorkSize {0};
   
public:
   
   /// C'tor which sets up the non overlapped periodogram of the timeseries.
   averagePeriodogram( size_t avgLen /**< [in] The length of averaging.*/ );
   
   /// C'tor setting up an arbitrary overlap.
   averagePeriodogram( size_t avgLen,
                       size_t olap
                     );
   
   ~averagePeriodogram();
   
   int resize( size_t avgLen, 
               size_t olap 
             );
      
   std::vector<realT> & win();
   
   int resizeWin();
   
   size_t pgramSize();
   
   void operator()( realT * pgram,
                    const realT * ts,
                    size_t sz
                  );

   
   void operator()( std::vector<realT> & pgram,
                    const std::vector<realT> & ts
                  );
   
   std::vector<realT> operator()(std::vector<realT> & ts);

};

template<typename realT>
averagePeriodogram<realT>::averagePeriodogram(size_t avgLen)
{
   resize(avgLen, 0);
}

template<typename realT>
averagePeriodogram<realT>::averagePeriodogram( size_t avgLen,
                                               size_t olap
                                             )
{
   resize(avgLen, olap);
}

template<typename realT>
averagePeriodogram<realT>::~averagePeriodogram()
{
   if(m_tsWork) fftw_free(m_tsWork);
   if(m_fftWork) fftw_free(m_fftWork);
}

template<typename realT>
int averagePeriodogram<realT>::resize( size_t avgLen, 
                                       size_t olap 
                                     )
{
   m_avgLen = avgLen;
   m_overlap = olap;
   
   m_nOver = (m_avgLen-m_overlap);
   
   m_fft.plan(m_avgLen, 0, MXFFT_FORWARD, false);
   
   if(m_tsWork) fftw_free(m_tsWork);
   m_tsWork = fftw_malloc<realT>( m_avgLen );
   
   if(m_fftWork) fftw_free(m_fftWork);
   m_fftWork = fftw_malloc<std::complex<realT>>( (m_avgLen/2 + 1) );
   
   return 0;
}
   
template<typename realT>
std::vector<realT> & averagePeriodogram<realT>::win()
{
   return m_win;
}

template<typename realT>
int averagePeriodogram<realT>::resizeWin()
{
   m_win.resize(m_avgLen, static_cast<realT>(1));
   
   return 0;
}

template<typename realT>
size_t averagePeriodogram<realT>::pgramSize()
{
   return  m_avgLen/2 + 1;
}

template<typename realT>
void averagePeriodogram<realT>::operator()( realT * pgram,
                                            const realT * ts,
                                            size_t sz
                                          )
{
   if( m_win.size() > 0 && m_win.size() != m_avgLen )
   {
      std::cerr << "averagePeriodogram: Window size not correct.\n";
   }

   size_t pgSize = pgramSize();
   
   int Navg = sz/m_nOver;

   while(Navg*m_nOver + m_avgLen > sz) --Navg;

   if(Navg < 1) Navg = 1; //Always do at least 1!

   //mx::fftT< realT, std::complex<realT>, 1, 0> fft(m_avgLen);

   for(int i=0;i<Navg;++i)
   {
      if(m_win.size() == m_avgLen) //no if inside the for loop
      {
         for(size_t j=0;j<m_avgLen;++j)
         {
            m_tsWork[j] = ts[i*m_nOver + j] * m_win[j];
         }
      }
      else
      {
         for(size_t j=0;j<m_avgLen;++j)
         {
            m_tsWork[j] = ts[i*m_nOver + j];
         }
      }

      m_fft( m_fftWork, m_tsWork);

      for(size_t j=0;j<pgSize;++j) pgram[j] += norm(m_fftWork[j]);
   }

   //This is what you'd do to normalize for dt=1
   //for(size_t j=0;j<pgSize;++j) pgram[j] /= (m_avgLen*Navg);

   //but we will just always normalize:

   realT df = 1.0/(2.0*pgSize*m_dt); //factor of 2 since it's a one-sided PSD
   
   realT pgramVar = psdVar(df, pgram, pgSize);
   
   realT tsVar = mx::math::vectorVariance(ts, sz);

   for(size_t j =0; j< pgSize; ++j) pgram[j] *= tsVar/pgramVar; //*df;
   
   
}

template<typename realT>
void averagePeriodogram<realT>::operator()( std::vector<realT> & pgram,
                                            const std::vector<realT> & ts
                                          )
{
   pgram.resize(pgramSize());
   operator()( pgram.data(), ts.data(), ts.size());
   
}

template<typename realT>
std::vector<realT> averagePeriodogram<realT>::operator()(std::vector<realT> & ts)
{
   std::vector<realT> pgram;
   operator()(pgram, ts);
   return pgram;
}
   
#if 0
///Calculate the average periodogram from a time series for a specified averaging interval and overlap.
/** The time series should be mean subtracted before passing to this function.
  *
  * The frequency scale of the output periodogram is 1/(2*pgram.size()*dt), where the factor of 2 is due to the one-sided-ness of the result.
  *
  * If a window is supplied, the PSD is normalized so that the one-sided integral is equal to the variance of the input time-series.
  * 
  * If you just want the FFT of a timeseries with no overlap, set dt=1, avgLen=\<length of time series\>, olap=0.  This will still window
  * and normalize if needed.
  */
template<typename realT>
void averagePeriodogram( std::vector<realT> & pgram, ///< [out] the resultant periodogram.
                         std::vector<realT> & ts,    ///< [in] is input the time-series.
                         realT dt,                   ///< [in] is the sampling time of time-series.
                         realT avgLen,               ///< [in] is the length of the averaging interval, same units as dt.
                         realT olap,                 ///< [in] is the length of the overlap region, same units as avgLen.
                         std::vector<realT> & w      ///< [in] a vector of length ( (int) avgLen/dt) containing a window.  If empty then then the square window is used.
                       )
{
   size_t Nper = avgLen/dt;
   int Nover = (avgLen-olap)/dt;

   if( w.size() > 0 && w.size() != Nper )
   {
      std::cerr << "averagePeriodogram: Window size not correct.\n";
   }


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

      for(size_t j=0;j<Nper;++j)
      {
         v = ts[i*Nover + j];

         if(w.size() == Nper) v *= w[j];

         cwork[j] = std::complex<realT>(v, 0);
      }

      fft( fftwork.data(), cwork.data());

      for(size_t j=0;j<pgram.size();++j) pgram[j] += norm(fftwork[j]); //pow(abs(fftwork[j]),2);
   }

   for(size_t j=0;j<pgram.size();++j) pgram[j] /= (Nper*Navg);

   //realT varNorm = 1;
   if(w.size() == Nper)
   {
      realT df = 1.0/(2.0*pgram.size()*dt); //factor of 2 since it's a one-sided PSD
      realT pgramVar = psdVar(df, pgram);

      realT tsVar = mx::math::vectorVariance(ts);

      for(size_t j =0; j< pgram.size(); ++j) pgram[j] *= tsVar/pgramVar; //*df;

   }


};

#endif

///@}

} //namespace sigproc
} //namespace mx

#endif //averagePeriodogram_hpp
