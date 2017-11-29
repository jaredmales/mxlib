/** \file autocorrelation.hpp
  * \brief Tools for working with autocorrelations
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
  *
  */

#ifndef autocorrelation_hpp
#define autocorrelation_hpp

#include "../fft/fft.hpp"

namespace mx
{
   
namespace sigproc 
{
   

///Calculate the autocorrelation of a time-series
/**
  * \tparam T is the real type of the data and autocorrelation. 
  * 
  * \ingroup signal_processing
  */
template<typename T>
void autocorrelation( T * ac, ///< [out] is the pre-allocated array of length Nac which will contain the autocorrelation on return
                      size_t Nac, ///< [in] is the length of ac 
                      T * sig, ///< [in] is the input time-series (signal)
                      size_t Nsig ///< [in] is the length of the input time-series
                    )
{

   #pragma omp parallel for
   for(int i=0; i< Nac; ++i)
   {
      ac[i] = 0;
      
      for(int j=0; j<Nsig-i; ++j)
      {
         ac[i] += sig[j]*sig[j+i];
      }
   }
   
   T norm = ac[0];
   for(int i=0; i<Nac; ++i) ac[i] /= norm;
}

///Calculate the autocorrelation of a time-series
/**
  * \tparam T is the real type of the data and autocorrelation. 
  * 
  * \ingroup signal_processing
  */
template<typename T>
void autocorrelation( std::vector<T> & ac, ///< [out] will contain the autocorrelation on return.  If ac.size()==0 then it is resized to sig.size().
                      std::vector<T> & sig  ///< [in] is the input time-series (signal)
                    )
{
   if(ac.size()==0) ac.resize(sig.size());
   autocorrelation( ac.data(), ac.size(), sig.data(), sig.size());
}

/// Functor for calculating the autocorrelation given a PSD 
/** Stores the fftT object and related working memory so that
  * repeated calls do not re-allocate or re-plan the FFT.
  *
  * \tparam T is the real type of the PSD and resultig A.C. 
  *
  * \ingroup signal_processing 
  */
template<typename T>
struct autocorrelationFromPSD
{
   std::vector<std::complex<T>> fftOut;
   std::vector<std::complex<T>> fftIn;
   
   mx::fftT<std::complex<T>, std::complex<T>, 1, 0> fft;
   
   /// Calculate the A.C. as the inverse FFT of the PSD 
   /** This calculates the circular autocorrelation from the PSD.
     * 
     */
   void operator()( T * ac, ///<  [out] pre-allocated array, on output contains the first Nac points of the autocorrelation
                    size_t Nac, ///< [in] the allocated size of ac.
                    T * psd, ///< [in] the 2-sided FFT storage order PSD 
                    size_t Npsd ///< [in] the number of points in the PSD 
                  )
   {
      fft.plan( Npsd, MXFFT_FORWARD);
      
      fftOut.resize(Npsd);
      fftIn.resize(Npsd);
      
      for(int i=0; i< Npsd; ++i) fftIn[i] = psd[i];
      
      fft( fftOut.data(),fftIn.data() );
      
      T norm = fftOut[0].real();
      for(int i=0; i < Npsd && i < Nac; ++i) ac[i] = fftOut[i].real()/norm;
   
   }

   /// Calculate the A.C. as the inverse FFT of the PSD 
   /** This calculates the circular autocorrelation from the PSD.
     * 
     */
   void operator() ( std::vector<T> & ac, ///< [out] On output contains the autocorrelation.  If ac.size()==0, it is allocated to psd.size().
                     std::vector<T> & psd  ///<  [in] the 2-sided FFT storage order PSD 
                   )
   {
      if(ac.size() == 0) ac.resize(psd.size());
      operator()( ac.data(), ac.size(), psd.data(), psd.size() );
   }
};



} //namespace sigproc 
} //namespace mx

#endif //autocorrelation_hpp
