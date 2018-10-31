/** \file speckleAmpPSD.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Calculates the PSD of speckle intensity given a modified Fourier mode amplitude PSD.
  * \ingroup mxAO_files
  *
  */

#ifndef speckleAmpPSD_hpp
#define speckleAmpPSD_hpp

#include <vector>

#include <Eigen/Dense>

#include "../../randomT.hpp"
#include "../../math/vectorUtils.hpp"

#include "../../sigproc/signalWindows.hpp"
#include "../../sigproc/psdUtils.hpp"
#include "../../sigproc/psdFilter.hpp"
#include "../../sigproc/psdVarMean.hpp"
using namespace mx::sigproc;


namespace mx
{
namespace AO
{
namespace analysis
{

/// Calculate and accumulate the means of a timeseries in bins of various sizes.
/** Used to, say, calculate the variance of the mean as a function of sample size.
  *
  * \returns 0 on success.
  */ 
template<typename vectorT, typename binVectorT>
int binMeans( std::vector<vectorT> & means, ///< [out] the means in each bin.  Not cleared, but Will be resized with new means appended.
              binVectorT & bins, ///< [in] the bin sizes in which to calculate the means
              vectorT & v ///< [in] the input vector to bin .
            )
{
   vectorT binv;

   for(size_t i=0; i< bins.size(); ++i)
   {
      mx::math::vectorRebin(binv, v, bins[i], true);

      means[i].resize(means[i].size() + binv.size());
      for(size_t j=0; j< binv.size(); ++j)  means[i][ means[i].size() - binv.size() + j] = binv[j];
   }
   
   return 0;
}

/// Calculate the variance of the mean vs bin size in a speckle time-series
/**
  * Does so by generating a time-series from the PSD using Fourier-domain convolution, then
  * calculates the binned-variances in the generated time-series.
  *
  */
template<typename realT>
int speckleAmpVarMean( std::vector<realT> & vars,  ///< [out] [optional] The binned variances of the time series generated.
                       std::vector<realT> & bins,  ///< [in]  [optional] The bin sizes to use in calculating vars.
                       std::vector<realT> & freq,  ///< [in] The Frequency grid of the input PSD
                       std::vector<realT> & fmPSD, ///< [in] The Fourier mode PSD.
                       int N                       ///< [in] The number of trials to use in calculating the amplitude PSD.
                     )
{
   std::vector<realT> vpsd2;

   bool hasZero = false;
   if(freq[0] == 0) hasZero = true;

   //First augment to two-sided DFT form
   mx::sigproc::augment1SidedPSD( vpsd2, fmPSD, !hasZero, 0.5);

   Eigen::Array<realT, -1,-1> psd2( vpsd2.size(), 1);
   for(size_t i=0; i< vpsd2.size(); ++i) psd2(i,0) = vpsd2[i];

   int Nwd = 0.5*psd2.rows();
   int NwdStart = 0.5*psd2.rows() - 0.5*Nwd;
   
   int Nsamp = 0.5*psd2.rows();
   int NsampStart = 0.5*Nwd - 0.5*Nsamp;
      
   std::vector<std::vector<realT>> means;
   means.resize(bins.size());

   #pragma omp parallel
   {
      mx::sigproc::psdFilter<realT> filt;
      filt.psd(psd2);

      mx::normDistT<realT> normVar;
      normVar.seed();

      //The two modal amplitude series
      Eigen::Array<realT, -1, -1> n( psd2.rows(), 1), nm(psd2.rows(),1);

      //The speckle amplitude
      std::vector<realT> vnl( Nwd );
      std::vector<realT> vn( Nsamp );

      #pragma omp for
      for(int k=0; k < N; ++k)
      {
         //Generate the two time-series
         //Note don't use the 2-for-1 option of psdFilter, it will produce correlated noise series.
         for(int i=0; i< psd2.rows(); ++i) n(i,0) = normVar;
         filt(n);

         for(int i=0; i< psd2.rows(); ++i) nm(i,0) = normVar;
         filt(nm);

         //Calculate the speckle amplitude
         for(int i= 0; i< Nwd; ++i)
         {
            vnl[i] = (pow(n(i+NwdStart,0),2) + pow(nm(i+NwdStart,0),2));
         }

         //Get the middle sample
         for(int i=0; i<Nsamp; ++i) vn[i] = vnl[i+NsampStart];
                  
         //Accumulate
         #pragma omp critical
         {
            binMeans( means, bins, vn);
         }
      }
   }//pragma omp parallel

   //Calculate binned variances.
   vars.resize(bins.size());
   for(size_t i=0; i< bins.size(); ++i)
   {
      vars[i] = mx::math::vectorVariance(means[i]) ;
   }

   return 0;
}

/// Calculate the PSD of the speckle intensity given the PSD of Fourier mode amplitude
/**
  * Does so by generating a time-series from the PSD using Fourier-domain convolution, and calculating
  * the average periodogram.
  *
  * A Hann window is used.
  *
  * Statistics are normalized so that the mean contrast is equal to the variance of the input Fourier mode.
  *
  * Will also calculate the binned-variances in the generated time-series, if the arguments vars and bins are not null pointers.
  *
  * \todo Figure out what's going on with psdFilter normalization!
  * \todo probably need to not do overlapped averaging in periodogram.  Could drop mean-sub then.
  */
template<typename realT>
int speckleAmpPSD( std::vector<realT> & spFreq,         ///< [out] The frequency grid of the output PSD
                   std::vector<realT> & spPSD,          ///< [out] The speckle amplitude PSD corresponding to the freq coordinates.  Will be resized.
                   std::vector<realT> & freq,           ///< [in] The Frequency grid of the input PSD
                   std::vector<realT> & fmPSD,          ///< [in] The Fourier mode PSD.
                   int N,                               ///< [in] The number of trials to use in calculating the amplitude PSD.
                   std::vector<realT> * vars = nullptr, ///< [out] [optional] The binned variances of the time series generated.
                   std::vector<realT> * bins = nullptr  ///< [in]  [optional] The bin sizes to use in calculating vars.
                 )
{

   std::vector<realT> vpsd2;

   bool hasZero = false;
   if(freq[0] == 0) hasZero = true;

   //First augment to two-sided DFT form
   mx::sigproc::augment1SidedPSD( vpsd2, fmPSD, !hasZero, 0.5);

   Eigen::Array<realT, -1,-1> psd2( vpsd2.size(), 1);
   for(size_t i=0; i< vpsd2.size(); ++i) psd2(i,0) = vpsd2[i];

   //The time sampling
   realT dt = 1./(2*freq.back());

   int Nwd = 0.5*psd2.rows();
   int NwdStart = 0.5*psd2.rows() - 0.5*Nwd;
   
   int Nsamp = 0.5*psd2.rows();
   std::cerr << "Nspeck = " << Nsamp << " (" << Nsamp*dt << ")\n";
   int NsampStart = 0.5*Nwd - 0.5*Nsamp;
      

   spPSD.resize(Nsamp/2,0);

   //Make a Hann window
   std::vector<realT> w;
   w.resize(Nsamp);
   mx::sigproc::tukey1d(w.data(), w.size(), 1.0);


   std::vector<std::vector<realT>> means;
   if( bins != nullptr && vars != nullptr) means.resize(bins->size());

   //Calculate the Fourier Mode variance for normalization
   //realT fmVar = psdVar( freq[1]-freq[0], fmPSD);


   #pragma omp parallel
   {
      mx::sigproc::psdFilter<realT> filt;
      filt.psd(psd2);

      mx::normDistT<realT> normVar;
      normVar.seed();

      //The two modal amplitude series
      Eigen::Array<realT, -1, -1> n( psd2.rows(), 1), nm(psd2.rows(),1);

      //The speckle amplitude
      std::vector<realT> vnl( Nwd );
      std::vector<realT> vn( Nsamp );

      //The temporary periodogram
      std::vector<realT> tpgram(spPSD.size());

      #pragma omp for
      for(int k=0; k < N; ++k)
      {
         //Generate the two time-series
         //Note don't use the 2-for-1 option of psdFilter, it will produce correlated noise series.
         for(int i=0; i< psd2.rows(); ++i) n(i,0) = normVar;
         filt(n);

         for(int i=0; i< psd2.rows(); ++i) nm(i,0) = normVar;
         filt(nm);

         //Calculate the speckle amplitude and mean-subtract
         realT mn = 0;
         for(int i= 0; i< Nwd; ++i)
         {
            vnl[i] = (pow(n(i+NwdStart,0),2) + pow(nm(i+NwdStart,0),2));

            mn += vnl[i];
         }
         mn /= vnl.size();

         
         //Mean subtract, and normalize by the Fourier-Mode variance, which should be the mean contrast.
         //for(size_t i=0; i< vnl.size(); ++i) vnl[i] = (vnl[i]-mn) * fmVar/mn ;

         for(int i=0; i<Nsamp; ++i) vn[i] = vnl[i+NsampStart];
         
         //Calculate PSD of the speckle amplitude
         mx::sigproc::averagePeriodogram<realT>( tpgram, vn, 1, Nsamp, 0, w);

         //std::cerr << spPSD.size() << " " << tpgram.size() << "\n";
         //Accumulate
         #pragma omp critical
         {
            if( bins != nullptr && vars != nullptr) binMeans( means, *bins, vn);
            for(size_t i=0; i< spPSD.size(); ++i) spPSD[i] += tpgram[i];
         }
      }
   }//pragma omp parallel

   //Normalize pgram by number of trials
   for(size_t i=0; i< spPSD.size(); ++i) spPSD[i] /= N;

   //-- Calculate frequency grid (which is maybe not identical to input freq)
   realT df = (1.0)/(2.0*spPSD.size()*dt);
   spFreq.resize(spPSD.size());
   for(size_t i =0; i< spFreq.size(); ++i) spFreq[i] = i*df;

   //Calculate binned variances.
   if( bins != nullptr && vars != nullptr)
   {
      vars->resize(bins->size());
      for(size_t i=0; i< bins->size(); ++i)
      {
         (*vars)[i] = mx::math::vectorVariance(means[i]) ;
      }
   }

   return 0;
}

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //speckleAmpPSD_hpp
