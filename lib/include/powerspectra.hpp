/** \file powerspectra.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Tools for working with Fourier transforms and power spectra
  * \ingroup psds
  * 
*/

#ifndef __mx_powerspectra__
#define __mx_powerspectra__

#include <cmath>
#include "randomT.hpp"

namespace mx
{

/** \addtogroup psds 
  */
//@{


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
template<typename eigenVec>
void frequency_grid1D( eigenVec & vec,
                       typename eigenVec::Scalar dt)
{
   typename eigenVec::Index dim, dim_1, dim_2;
   typename eigenVec::Scalar df;
   
   dim_1 = vec.rows();
   dim_2 = vec.cols();
   
   dim = std::max(dim_1, dim_2);
   
   df = freq_sampling(dim, 0.5/dt);
      
   for(int ii=0; ii < ceil(0.5*(dim-1) + 1); ++ii)
   {
      vec(ii) = ii*df;
   }
   
   for(int ii=ceil(0.5*(dim-1)+1); ii < dim_1; ++ii)
   {
      vec(ii) = (ii-dim)*df;
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
   
template<typename eigenArrp, typename eigenArrs>
void calcPSD( eigenArrp & psd,
              eigenArrs & sample
            )
{
   typedef typename eigenArrp::Scalar realT;
   
   Eigen::Array<std::complex<realT>, Eigen::Dynamic, Eigen::Dynamic> ft;
   
   ft.resize(sample.rows(), sample.cols());
   ft.setZero();
   ft.real() = sample;
   
   fft(ft, ft);
   
   psd = (ft*conj(ft)).real();
}

template<typename eigenArra, typename eigenArrp>
void calcAC( eigenArra & ac,
             eigenArrp & psd
           )
{
   typedef typename eigenArra::Scalar realT;
   
   Eigen::Array<std::complex<realT>, Eigen::Dynamic, Eigen::Dynamic> ft;
   
   ft.resize(psd.rows(), psd.cols());
   ft.setZero();
   
   ft.real() = psd;
   
   fft(ft, ft);
   
   ac = ft.real();
}   
   
   
///Calculate the normalization for a @f$ 1/|f|^\alpha @f$ PSD.
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
/** \ingroup psds
  * Populates an Eigen array  with
  * \f[
  *  P(|f| = 0) = 0 
  * \f]
  * \f[
  *  P(|f| > 0) = c/|f|^{\alpha} 
  * \f]
  * 
  *
  * \param psd [output] is the array to populate
  * \param freq [input] is a frequency grid, must be the same logical size as psd
  * \param alpha [input] is the power law exponent, by convention @f$ alpha > 0 @f$.
  * \param c [optional input] is a normalization constant to multiply the raw spectrum by.  If c==-1 (default) then 
  *                           the PSD is normalized using \ref oneoverf_norm.
  *
  * \tparam eigenArrp is the Eigen array type of the psd 
  * \tparam eigenArrf is the Eigen array type of the frequency grid 
  */
template<typename eigenArrp, typename eigenArrf> 
void oneoverf_psd( eigenArrp  & psd,
                   eigenArrf & freq,
                   typename eigenArrp::Scalar alpha,
                   typename eigenArrp::Scalar c = -1 )
{   
   typedef typename eigenArrp::Scalar Scalar;
   
   typename eigenArrp::Index dim_1, dim_2;
   Scalar f_x, f_y, p;
   
   dim_1 = psd.rows();
   dim_2 = psd.cols();
   
   
   if(c==-1)
   {
      Scalar fmin;
      Scalar fmax;
      
      fmax = freq.abs().maxCoeff();
      
      //Find minimum non-zero Coeff.
      fmin = (freq.abs()>0).select(freq.abs(), freq.abs()+fmax).minCoeff();
      
      c = oneoverf_norm(fmin, fmax, alpha);
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
            p = c/pow(abs(freq(ii,jj)), alpha);
         }
         psd(ii,jj) = p;
      }
   }   
}


template<typename eigenArrn, typename eigenArrp>
void psd_filter( eigenArrn & noise,
                 eigenArrp & psd )
{
   typedef std::complex<typename eigenArrn::Scalar> complexT;
   
   Eigen::Array< complexT, Eigen::Dynamic, Eigen::Dynamic> ft(noise.rows(), noise.cols());
   Eigen::Array< complexT, Eigen::Dynamic, Eigen::Dynamic> cnoise(noise.rows(), noise.cols());
   
   //Make noise a complex 
   for(int ii=0;ii<noise.rows();++ii)
   {
      for(int jj=0; jj<noise.cols(); ++jj)
      {
         cnoise(ii,jj) = complexT(noise(ii,jj), 0);
      }
   }
   
   fft(ft, cnoise);
   
   ft *= psd.sqrt();
   fft(ft, ft); //in-place
   
   noise = ft.real()/(noise.rows()*noise.cols());
}

template<typename arithT>
arithT sinc(arithT x)
{
   if(x == 0) return 1.;
   
   return sin(x)/x;
}

template<typename eigenArrs, typename ranT, typename realT>
void generateSubHarmonics( eigenArrs & subh,
                           int nsub,
                           ranT & rand,
                           realT df,
                           realT alpha,
                           realT c
                         )
{
   typedef typename eigenArrs::Scalar complexT;
   
   
   subh.setZero();
   
   if(nsub == 0) return; 
  
   realT Nk, Nl, fks, fls, fs, W, amp, phase, rcos, isin; 

   int dim=2;
   
   if(subh.rows() == 1 || subh.cols() == 1) dim=1;
   
   Nk = subh.rows();
   Nl = subh.cols();
      
   for(int n=0; n < nsub; ++n)
   {
      W = sqrt(1./pow(3., dim*(n+1)));

      for(int ks = 0; ks < 2; ++ks)
      {
         if(subh.rows() == 1 && ks != 0) continue;
                     
         fks = ks/pow(3., n+1);
         
         for(int ls = 0; ls < 2; ++ls)
         {
            if(subh.cols() == 1 && ls != 0) continue;

            if(ks == 0 && ls == 0) continue;
                  
            fls = ls/pow(3., n+1);
                  
            fs = sqrt(fks*fks + fls*fls);//*df;
     
            //Calculate amp
            // -- c is the normalization constant
            // -- W is the subharmonic weight
            // -- (ks+ls) multiplies by sqrt(2) if in the (1,1) corner
            // -- The rest is power law
            amp = sqrt(c)* W * sqrt(ks + ls) * (1./pow(abs(fs*df), 0.5*alpha));
            
            //Get next random number
            phase = rand;
            
            
            for(int k=0; k < subh.rows(); ++k)
            {     
               for(int l=0; l < subh.cols(); ++l)
               {
                  rcos = amp*cos(D2PI*(fks*k/Nk + fls*l/Nl + phase));
                  isin = amp*sin(D2PI*(fks*k/Nk + fls*l/Nl + phase));
                  subh(k,l) += complexT(rcos, isin);    
               }
            }
         }
      }
   }
   
}
  
template<typename eigenArrn, typename ranT>
void oneoverf_noise( eigenArrn & noise, 
                     ranT & rand, 
                     typename eigenArrn::Scalar dt,
                     typename eigenArrn::Scalar alpha,
                     typename eigenArrn::Scalar c=-1,
                     int nsub = 0
                   )
{
   typedef typename eigenArrn::Scalar Scalar;
   typedef typename eigenArrn::Index Index;
   typedef std::complex<typename eigenArrn::Scalar> complexT;
   
   Scalar df;
   Index dim_1, dim_2;
   
   eigenArrn freq(noise.rows(), noise.cols());
   eigenArrn psd(noise.rows(), noise.cols());
   Eigen::Array<complexT, Dynamic, Dynamic> subh(noise.rows(), noise.cols());
   
   dim_1 = noise.rows();
   dim_2 = noise.cols();
   
   df = freq_sampling(std::max(dim_1, dim_2), 0.5/dt);
   
   
   frequency_grid( freq, dt);
   
   
   if(c==-1)
   {
      Scalar fmin;
      Scalar fmax;
      
      fmax = freq.abs().maxCoeff();
      
      //Find minimum non-zero Coeff.
      fmin = (freq.abs()>0).select(freq.abs(), freq.abs()+fmax).minCoeff();
      
      c = oneoverf_norm(fmin, fmax, alpha);
   }
   
   oneoverf_psd( psd, freq, alpha, c);
   
   
   Eigen::Array< complexT, Eigen::Dynamic, Eigen::Dynamic> ft(noise.rows(), noise.cols());
   Eigen::Array< complexT, Eigen::Dynamic, Eigen::Dynamic> cnoise(noise.rows(), noise.cols());
   
   //Make noise a complex 
   for(int ii=0;ii<noise.rows();++ii)
   {
      for(int jj=0; jj<noise.cols(); ++jj)
      {
         cnoise(ii,jj) = psd.sqrt()(ii,jj)*exp(complexT(0,D2PI*rand));
      }
   }
   
   //cnoise *= (psd.sqrt()*complexT(1,0) + subh);
   //cnoise += subh;
   //fft(ft, cnoise);
   
   //ft = ft*psd.sqrt() + subh;
//   ft = psd.sqrt()*cnoise;// + subh;
   fft(ft, cnoise, FFTW_BACKWARD); //in-place

   generateSubHarmonics(subh, nsub, rand, df, alpha, c);
   
   ft = ft + subh;
   
   for(int ii=0;ii<noise.rows();++ii)
   {
      for(int jj=0; jj<noise.cols(); ++jj)
      {
         noise(ii,jj) = ft(ii,jj).real();
      }
   }
}


template<typename eigenArrw>
void tukeyWindow1D( eigenArrw & w,
                  typename eigenArrw::Index N,
                  typename eigenArrw::Scalar alpha
                )
{
   typedef typename eigenArrw::Index Index;
   typedef typename eigenArrw::Scalar Scalar;
   
   w.resize(N,1);

   //Avoid divide by 0, etc.,if we don't need it 
   if(alpha == 0.)
   {
      w.setOnes();
      return;
   }

   Index t = 0;
   for(t;t <= 0.5*alpha*(N-1); ++t)
   {
      w(t) = 0.5*(1. + cos(DPI*(2*t/(alpha*N-1) - 1)));
   }
   for(t; t < (N-1)*(1-0.5*alpha); ++t)
   {
      w(t) = 1.;
   }
   for(t; t<N; ++t)
   {
      w(t) = 0.5*(1. + cos(DPI*(2*t/(alpha*N-1) - 2./alpha + 1))) ;
   }
}

/*group psds*/
//@} 

}//namespace mx

template<typename realT>
void test_psd_and_ac(realT N, int simNx, realT alpha, std::string fname)
{
   Eigen::ArrayXXd noise, sample, w, psd, ac, psdavg, acavg;
   int Ntrials = 1000;

   uni_distd rand;

   noise.resize(N*simNx, 1);

   sample.resize(N,1);
   psd.resize(N,1);
   ac.resize(N,1);

   psdavg.resize(N,1);
   psdavg.setZero();
   acavg.resize(N,1);
   acavg.setZero();

   w.resize(N,1);
   tukeyWindow1D(w, N, 1.0);

   for(int n = 0; n<Ntrials; ++n)
   {
      oneoverf_noise(noise, rand, 1./(24.*3.), alpha, -1, 0);
      //noise = noise-noise.mean();
      //pout(noise.square().mean());
      //noise = noise/sqrt(noise.square().mean());

      sample = noise.block(.5*(N*simNx-1)-.5*N, 0, N, 1);

      calcPSD(psd,  (sample-sample.mean())*w );

      psdavg += psd; 
      
      calcPSD(psd,  sample*w);
      calcAC(ac, psd);

      double mx = ac(0);//.maxCoeff();

      ac /= mx;

      acavg += ac;
      
      
      if( n == floor(0.5*Ntrials) )
      {
         std::ofstream fout;
         fout.open(fname+"_sample.txt");
         for(int i=0; i<N; ++i)
         {
            fout << sample(i) << "\n";
         }

         fout.close();
         
         fout.open(fname+"_noise.txt");
         for(int i=0; i<N*simNx; ++i)
         {
            fout << noise(i) << "\n";
         }

         fout.close();
      }
   }
   std::ofstream fout;
   fout.open(fname+"_psdac.txt");

   for(int i=0; i<N; ++i)
   {
      fout << psdavg(i)/Ntrials << " " << acavg(i)/Ntrials << "\n";
   }

   fout.close();

}  
   
   
   
   
   
#endif //#ifdef __mx_powerspectra__

