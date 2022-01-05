/** \file averagePeriodogram.hpp
  * \brief A class to manage calculation of periodograms from time series data.
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

/// Calculate the average periodogram of a time-series.
/**
  * Implements the overlapped average periodogram, cf. pp 843 to 845 of \cite 10.5555.1795494. This produces
  * the variance normalized 1-sided PSD estimate of the periodogram.  
  * 
  * Optionally includes a window (by default no window is used, which is equivalent to the rectangle window). 
  * 
  * Can also be used to calculate the unaveraged non-overlapped periodogram of the entire time series.
  * 
  * Example:
  \code
   typedef float realT;
   
   mx::fftwEnvironment<realT> fftwEnv; //for efficient fft planning
   
   std::vector<realT> ts = getTimeSeries(); //Some function to populate ts.
   
   math::vectorMeanSub(ts); //The time series should be overall mean subtracted.
   
   realT dt = 0.1; //The sampling of ts is 0.1 seconds.
   
   averagePeriodogram<realT> avgPgram(2.0/dt, dt); //This sets the averaging length to 2 seconds (20 samples), with a 1 second overlap (set automatically)
   
   avgPgram.win(window::hann); //Set the Hann window
   
   std::vector<realT> pgram = avgPgram(ts); //Calculate the half-overlapped periodogram estimate for ts.
   
   //Now output the periodogram estimate.  This will give two columns, the first frequency, the second power at that frequency.
   for(size_t n=0; n< avgPgram.size(); ++n)
   {
      std::cout << avgPgram[n] << " " << pgram[n] << "\n"; 
   }
  
  \endcode
  * 
  * \ingroup psds
  */
template<typename realT>
class averagePeriodogram
{
   
protected:
   
   /** \name Configuration Parameters
     *
     *@{
     */
   
   size_t m_avgLen {0};  ///< The number of samples in each periodgoram estimate.
   size_t m_overlap {0}; ///< The number of samples by which to overlap.  This should almost always be 0.5*m_avgLen.  Set 0 for the non-overlapped case.
   
   realT m_dt {1}; ///< The time sampling.  Only used for normalization and calculation of the frequency scale.

   std::vector<realT> m_win; ///< The window function.  By default this is empty, which is equivalent to setting it to the rectanbular window.
      
   ///@}
   
   size_t m_size {0}; ///< The size of the periodogram vector, calculated as m_avgLen/2 + 1;
   
   int m_nOver {0}; ///< The number of overlapping segments.  Calculated from m_avgLen and m_overlap;
   
   realT m_df {1}; ///< The frequency sampling.  This is used only for normalization and frequency scale output.

   
   math::fft::fftT< realT, std::complex<realT>, 1, 0> m_fft;
   
   realT * m_tsWork {nullptr};
   size_t m_tsWorkSize {0};
   
   std::complex<realT> * m_fftWork {nullptr};
   size_t m_fftWorkSize {0};
   
public:
   
   /// C'tor which sets up the optimum overlapped periodogram of the timeseries.
   /** Sets m_overlap = 0.5*m_avgLen.  If you desire the non-overlapped case use 
     * the alternate constructor:
     * \code
       averagePeriogram p( avgLen, 0, dt);
       \endcode
     */
   explicit averagePeriodogram( size_t avgLen /**< [in] The length of averaging in samples.*/ );
   
   /// C'tor which sets up the optimum overlapped periodogram of the timeseries and sets the sampling.
   /** Sets m_overlap = 0.5*m_avgLen.  If you desire the non-overlapped case use 
     * the alternate constructor:
     * \code
       averagePeriogram p( avgLen, 0, dt);
       \endcode
     */
   averagePeriodogram( size_t avgLen, ///< [in] The length of averaging in samples. 
                       realT dt ///< [in] the sampling interval of the time-series
                     );
   
   /// C'tor setting up an arbitrary overlap.
   /** Set olap to 0 for the unoverlapped case.  
     */
   averagePeriodogram( size_t avgLen, ///< [in] The number of samples in each periodgoram estimate.
                       size_t olap,   ///< [in] The number of samples by which to overlap.  This should almost always be 0.5*m_avgLen.  Set 0 for the non-overlapped case.
                       realT dt       ///< [in] the sampling interval of the time-series
                     );
   
   /// D'tor, frees all working memory.
   ~averagePeriodogram();
   
   /// Set the sampling interval of the time-series 
   /** This also sets the frequency scale of the output.
     */
   void dt( realT ndt);
   
   /// Get the sampling interval of the time-series
   /** \returns m_dt
     */
   realT dt();
   
   /// Get the frequency interval of the periodogram
   /** \returns m_df
     */
   realT df();
   
   /// Resize the periodogram working memory, setting up the 1/2-overlapped optimum case.
   /** This sets the overlap to 0.5*avgLen.
     *
     * Also performs fft planning.
     * 
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int resize( size_t avgLen /**< [in] The number of samples in each periodgoram estimate. */);
   
   /// Resize the periodogram working memory.
   /** Also performs fft planning.
     * 
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int resize( size_t avgLen, ///< [in] The number of samples in each periodgoram estimate.
               size_t olap    ///< [in] The number of samples by which to overlap.  This should almost always be 0.5*m_avgLen.  Set 0 for the non-overlapped case.
             );
      
   /// Get a reference to the window vector.
   /** This allows population of the window with any arbitrary function.  You should call
     * resizeWin() before using this to make sure that the vector has the correct length.
     *
     * \returns a reference to m_win.
     */ 
   std::vector<realT> & win();
   
   /// Set the window vector using a function.
   /** This will resize m_win, and then call the function with m_win as the argument.
     * For example to use the hann window defined \ref signal_windows1D:
     * \code
       avgPgram.win( window::hann );
       \endcode
     * which will set the Hann window.
     */
   void win( void(*winFunc)(std::vector<realT> &) /**<[in] pointer to a function which takes a pre-allocated vector and fills it in with a window function*/);
   
   /// Resize the window and, unless not desired, initialize to the rectangular window by setting all 1s.
   /** If `setRect=false` then the window is not initialized.
     */
   void resizeWin( bool setRect=true /**< [in] [optional] if false, then the window is not initialized to 1 */);
   
   /// Calculate the periodogram for a time-series.
   /** 
     *
     */ 
   void operator()( realT * pgram,    ///< [out] a pre-allocated to size() array which will be populated with the periodogram of the time-series 
                    const realT * ts, ///< [in] the time-series
                    size_t sz         ///< [in] the length of the time-series array
                  );

   
   /// Calculate the periodogram for a time-series.
   /** \overload
     *
     */ 
   void operator()( std::vector<realT> & pgram,   ///< [out] a vector which will be allocated and populated with the periodogram of the time-series 
                    const std::vector<realT> & ts ///< [in] the time-series
                  );
   
   /// Calculate the periodogram for a time-series.
   /** \overload
     *
     * \returns the periodogram as a vector.
     */ 
   std::vector<realT> operator()(std::vector<realT> & ts /**< [in] the time-series*/ );
   
   /// Return the size of the periodogram.
   /** 
     *
     * \returns the size of the periodogram estimate for the current setup.
     */ 
   size_t size();
   
   /// Get the frequency at a given point in the periodogram.
   /** 
     * \returns the frequency at point n in the periodogram.
     */ 
   realT operator[]( size_t n /**<[in] the point in the periodogram at which the frequency is desired*/ );

};

template<typename realT>
averagePeriodogram<realT>::averagePeriodogram( size_t avgLen )
{
   resize(avgLen, 0.5*avgLen);
   dt(1);
}

template<typename realT>
averagePeriodogram<realT>::averagePeriodogram( size_t avgLen,
                                               realT ndt
                                             )
{
   resize(avgLen, 0.5*avgLen);
   dt(ndt);
} 

template<typename realT>
averagePeriodogram<realT>::averagePeriodogram( size_t avgLen,
                                               size_t olap,
                                               realT ndt
                                             )
{
   resize(avgLen, olap);
   dt(ndt);
}

template<typename realT>
averagePeriodogram<realT>::~averagePeriodogram()
{
   if(m_tsWork) fftw_free(m_tsWork);
   if(m_fftWork) fftw_free(m_fftWork);
}

template<typename realT>
void averagePeriodogram<realT>::dt( realT ndt)
{
   m_dt = ndt;
   m_df = 1.0/(m_avgLen*m_dt);
}

template<typename realT>
realT averagePeriodogram<realT>::dt()
{
   return m_dt;
}
   
template<typename realT>
realT averagePeriodogram<realT>::df()
{
   return m_df;
}

   
template<typename realT>
int averagePeriodogram<realT>::resize( size_t avgLen )
{
   return resize(avgLen, 0.5*avgLen);
}

template<typename realT>
int averagePeriodogram<realT>::resize( size_t avgLen, 
                                       size_t olap 
                                     )
{
   m_avgLen = avgLen;
    
   m_size = m_avgLen/2 + 1;
   
   m_df = 1.0/(m_avgLen*m_dt);
    
   m_overlap = olap;
   
   m_nOver = (m_avgLen-m_overlap);
   
   m_fft.plan(m_avgLen, MXFFT_FORWARD, false);
   
   if(m_tsWork) fftw_free(m_tsWork);
   m_tsWork = math::fft::fftw_malloc<realT>( m_avgLen );
   
   if(m_fftWork) fftw_free(m_fftWork);
   
   m_fftWork = math::fft::fftw_malloc<std::complex<realT>>( (m_avgLen/2 + 1) );
   
   return 0;
}
   
template<typename realT>
std::vector<realT> & averagePeriodogram<realT>::win()
{
   return m_win;
}

template<typename realT>
void averagePeriodogram<realT>::win( void(*winFunc)(std::vector<realT> &) )
{
   resizeWin(false);
   winFunc(m_win);
}
   
template<typename realT>
void averagePeriodogram<realT>::resizeWin(bool setRect)
{
   if(setRect)
   {
      m_win.resize(m_avgLen, static_cast<realT>(1));
   }
   else
   {
      m_win.resize(m_avgLen);
   }
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

      for(size_t j=0;j<m_size;++j) pgram[j] += norm(m_fftWork[j]);
   }

   //This is what you'd do to normalize for dt=1
   //for(size_t j=0;j<m_size;++j) pgram[j] /= (m_avgLen*Navg);

   //but we will just always normalize:

   realT pgramVar = psdVar1sided(m_df, pgram, m_size);
   
   realT tsVar = mx::math::vectorVariance(ts, sz);

   for(size_t j =0; j< m_size; ++j) pgram[j] *= tsVar/pgramVar; 
   
}

template<typename realT>
void averagePeriodogram<realT>::operator()( std::vector<realT> & pgram,
                                            const std::vector<realT> & ts
                                          )
{
   pgram.resize(m_size);
   operator()( pgram.data(), ts.data(), ts.size());
   
}

template<typename realT>
std::vector<realT> averagePeriodogram<realT>::operator()(std::vector<realT> & ts)
{
   std::vector<realT> pgram;
   operator()(pgram, ts);
   return pgram;
}
  
template<typename realT>
size_t averagePeriodogram<realT>::size()
{
   return m_size;
}

template<typename realT>
realT averagePeriodogram<realT>::operator[]( size_t n )
{
   return n*m_df;
}
   

} //namespace sigproc
} //namespace mx

#endif //averagePeriodogram_hpp
