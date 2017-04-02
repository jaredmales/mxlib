
#ifndef __linearPredictor_hpp__
#define __linearPredictor_hpp__

#include <vector>
#include <complex>

#include <mx/eigenUtils.hpp>
#include <mx/timeUtils.hpp>

#include "nrToeplitz.hpp"

template<typename _realT>
struct linearPredictor
{
   typedef _realT realT;
   
   Eigen::Array<realT, -1, -1> _c;
   
   realT _setCondition;
   realT _actCondition;
   int _nRejected;
   
   void calcCoefficients( std::vector<realT> & ac,
                            size_t Nc,
                            realT condition = 0 )
   {
      
      if(condition == 0)
      {
         return calcCoefficientsToep( ac, Nc );
      }
      
      
      Eigen::Array<realT,-1,-1> Rmat, Rvec, PInv, LPcoeff;   
      
      Rmat.resize(Nc, Nc);
      Rvec.resize(1, Nc);
   
      for(int i=0; i<Nc; ++i)
      {
         for(int j=0; j<Nc; ++j)
         {
            Rmat(i,j) = ac[ fabs(i-j)];
         }
         
         Rvec(0,i) = ac[i+1];
      }  

      realT tmpCond = condition;
      
      _setCondition = condition;
      double t0 = mx::get_curr_time();
      
      mx::eigenPseudoInverse(PInv, tmpCond, _nRejected,  Rmat, condition);
      
      
      std::cerr << mx::get_curr_time() - t0 << "\n";
      _actCondition = tmpCond;
      
      _c = Rvec.matrix()*PInv.matrix();

   }

   void calcCoefficientsToep( std::vector<realT> & ac,
                              size_t Nc )
   {
      std::vector<realT> r, x, y;
      
      r.resize(2.*Nc-1);
      x.resize(Nc);
      y.resize(Nc);
      
      for(int i=0; i< Nc; ++i) r[i] = ac[Nc-i - 1];
      for(int i=Nc; i< 2*Nc-1; ++i) r[i] = ac[i-Nc+1];
      
      for(int i=0;i<Nc; ++i) y[i] = ac[i+1];
      
      nrToeplitz(r.data(), x.data(), y.data(), Nc);
      
      _c.resize(1, Nc);
      for(int i=0; i< Nc; ++i) _c(0,i) = x[i];
      
   }
   
   realT c(size_t i)
   {
      return _c(0,i);
   }
   
   realT predict( std::vector<realT> & hist,
                  size_t idx = 0 )
   {
      realT x = 0;
      
      for(int i=0; i< _c.cols(); ++i)
      {
         x += _c[i]*hist[idx];
         ++idx;
         if(idx >= hist.size()) idx = 0; 
      }
   }
   
   realT spectralResponse( realT f, realT fs)
   {
      int n = _c.cols();
      
      std::complex<realT> He = 0;
      for(int j=0; j < n; ++j)
      {
         realT s = (j+1.0)* 2*3.14159;
         He += _c(0,j) * exp(s*std::complex<realT>(0,-1.0)*f/fs);
      }

      realT one = 1.0;
      return std::norm( one/ (one - He));
   }
   
};

#endif //__linearPredictor_hpp__

