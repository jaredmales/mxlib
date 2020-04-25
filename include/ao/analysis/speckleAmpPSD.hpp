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
#include "../../sigproc/averagePeriodogram.hpp"

#include "jincFuncs.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{



/// Calculate the variance of the mean vs bin size in a speckle time-series
/**
  * Does so by generating a time-series from the PSD using Fourier-domain convolution, then
  * calculates the binned-variances in the generated time-series.
  *
  */
template<typename realT>
int speckleAmpVarMean( std::vector<realT> & vars,                    ///< [out] The binned variances of the time series generated.
                       std::vector<realT> & bins,                    ///< [in] The bin sizes to use to calculate the variances.
                       std::vector<realT> & freq,                    ///< [in] The Frequency grid of the input PSD
                       std::vector<realT> & fmPSD,                   ///< [in] The open-loop Fourier mode PSD.
                       std::vector<realT> & fmDeltaPhase,            ///< [in] The phase angle between the two Fourier components, as a function of freq.
                       std::vector<std::complex<realT>> & fmXferFxn, ///< [in] The complex error transfer function, as a function of freq
                       std::vector<realT> & nPSD,                    ///< [in] The open-loop noise PSD
                       std::vector<std::complex<realT>> & nXferFxn,  ///< [in] The noise transfer function, as a function of freq
                       int m,
                       int n,
                       int N                       ///< [in] The number of trials to use in calculating the amplitude time-series.
                     )
{
   std::vector<realT> vpsd2, dphase2, npsd2;
   std::vector<std::complex<realT>> xfer2, nxfer2;

   bool hasZero = false;
   if(freq[0] == 0) hasZero = true;

   //First augment to two-sided DFT form
   sigproc::augment1SidedPSD( vpsd2, fmPSD, !hasZero, 0.5); //<<<<<<<<****** This needs to be the open loop PSD.

   //And augment the h+ to h- phase shift to two sided form
   sigproc::augment1SidedPSD( dphase2, fmDeltaPhase, !hasZero, 1.0);
   
   //And augment the xfer fxn to two sided form
   sigproc::augment1SidedPSD( xfer2, fmXferFxn, !hasZero, 1.0);
   
   //Now augment the noise PSD to two-sided DFT form
   sigproc::augment1SidedPSD( npsd2, nPSD, !hasZero, 0.5);
   
   //And augment the noise xfer fx
   sigproc::augment1SidedPSD( nxfer2, nXferFxn, !hasZero, 1.0);
   
   
   
   Eigen::Array<realT, -1,-1> psd2( vpsd2.size(), 1);
   for(size_t i=0; i< vpsd2.size(); ++i) psd2(i,0) = vpsd2[i];

   int Nwd = 0.5*psd2.rows();
   int NwdStart = 0.5*psd2.rows() - 0.5*Nwd;
   
   int Nsamp = 0.5*psd2.rows();
   int NsampStart = 0.5*Nwd - 0.5*Nsamp;
      
   std::vector<std::vector<realT>> means;
   

   realT Jp, Jm;
            
   jincFuncs( Jp, Jm, m, n, m, n);
            
   //std::cerr << "jfuncs: " << Jp << " " << Jm << "\n";
   
   #pragma omp parallel
   {
      mx::sigproc::psdFilter<realT> filt;
      filt.psd(psd2); 

      fftT<realT, std::complex<realT>, 1, 0> fft(vpsd2.size());
      fftT<std::complex<realT>, realT, 1, 0> fftB(vpsd2.size(), MXFFT_BACKWARD);
      
      std::vector<std::complex<realT>> tform1(vpsd2.size());
      std::vector<std::complex<realT>> tform2(vpsd2.size());
   
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
         for(int i=0; i< psd2.rows(); ++i) 
         {
            n(i,0) = normVar;
            //nm(i,0) = normVar;
         }
         
         filt.filter(n);//, &nm);

         fft(tform1.data(), n.data());
            
         //<<<<<<<<****** Apply the phase shift to the second one.
         for(size_t m=0;m<tform1.size();++m)
         {
            // Apply the phase shift to form the 2nd time series
            tform2[m] = tform1[m]*exp( std::complex<realT>(0, dphase2[m] ));//+ 0.02*std::complex<realT>(nm(n,0), nm(n,0));

            //Xpply the augmented ETF to two time-series
            tform1[m] *= xfer2[m]/std::complex<realT>(tform1.size(),0) ;
            tform2[m] *= xfer2[m]/std::complex<realT>(tform1.size(),0) ;
         }

         //<<<<<<<<****** This needs to be the open loop PSD.
         fftB(n.data(), tform1.data());
         fftB(nm.data(), tform2.data());
         
         //Calculate the speckle amplitude
         for(int i= 0; i< Nwd; ++i)
         {
            vnl[i] = (pow(n(i+NwdStart,0),2)   + pow(nm(i+NwdStart,0),2));//*(Jp*Jp + Jm*Jm) + 8*n(i+NwdStart,0)*nm(i+NwdStart,0)*Jp*Jm;
         }

         //Get the middle sample
         for(int i=0; i<Nsamp; ++i) vn[i] = vnl[i+NsampStart];
                  
         //std::cerr << math::vectorVariance(vn) << "\n";
         //Accumulate
         #pragma omp critical
         {
            math::vectorBinMeans( means, bins, vn);
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

      sigproc::averagePeriodogram<realT> avgPgram(Nsamp, 0, 1);
      avgPgram.win(sigproc::window::hann);
      
      //The temporary periodogram
      std::vector<realT> tpgram(spPSD.size());

      #pragma omp for
      for(int k=0; k < N; ++k)
      {
         //Generate the two time-series
         //Note don't use the 2-for-1 option of psdFilter, it will produce correlated noise series.
         for(int i=0; i< psd2.rows(); ++i) 
         {
            n(i,0) = normVar;
            nm(i,0) = normVar;
         }
         
         //filt.filter(n, &nm);
         filt(n);
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
         math::vectorMeanSub(vn);
         
         //Calculate PSD of the speckle amplitude
         //mx::sigproc::averagePeriodogram<realT>( tpgram, vn, 1, Nsamp, 0, w);
         avgPgram(tpgram, vn);
         
         //std::cerr << spPSD.size() << " " << tpgram.size() << "\n";
         //Accumulate
         #pragma omp critical
         {
            if( bins != nullptr && vars != nullptr) math::vectorBinMeans( means, *bins, vn);
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
