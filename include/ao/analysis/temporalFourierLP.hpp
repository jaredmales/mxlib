
#ifndef temporalFourierLP_hpp
#define temporalFourierLP_hpp


#include <vector>

#include <mx/sigproc/autocorrelation.hpp>

#include <mx/sigproc/psdUtils.hpp>

#include <mx/math/geo.hpp>

#include <mx/sigproc/linearPredictor.hpp>

namespace mx
{
namespace AO
{
namespace analysis 
{
   
template<typename realT>
struct temporalFourierLP
{
   std::vector<realT> PSDtn;
   
   std::vector<realT> psd2s;
   
   std::vector<realT> ac;
   
   sigproc::autocorrelationFromPSD<realT> acpsd;

   sigproc::linearPredictor<realT> lp;
   
   bool _windFact;
   
   temporalFourierLP()
   {
      _windFact = false;
   }
   
   void doit( std::vector<realT> & PSDt, 
              std::vector<realT> & PSDn, 
              int Nc,
              realT _k_dir,
              realT _dir_wind )
   {
      PSDtn.resize(PSDt.size());
      
      realT max= 0, tot = 0;
      
      for(int i=0; i< PSDt.size(); ++i)
      {
         if(PSDt[i] > max) max = PSDt[i];
         if(PSDt[i] > PSDn[0]) tot += PSDt[i]*(i+1);
      }
      
      //std::cerr << max << " " << tot << "\n";
      
      realT nFact = 1.0;
      
      if(_windFact)
      {
         realT dang = abs(math::angleDiff<1>( _k_dir, _dir_wind));
         //if(dang > 0.5*pi<realT>()) dang -= 0.5*pi<realT>();
         nFact = abs(1.0 - dang/ pi<realT>() * 1.8);
         if(nFact < 0.1) nFact = 0.1;
      
      }
      
      for(int i=0; i< PSDt.size(); ++i)
      {
         PSDtn[i] = PSDt[i] +  nFact*PSDn[i];// + 0;// 1e-9;
      }   

      sigproc::augment1SidedPSD( psd2s, PSDtn,1);

      ac.resize(psd2s.size());
   
      //#pragma omp critical
      acpsd(ac, psd2s);

      lp.calcCoefficients(ac, Nc);
   
      
   }
};

}//namespace analysis
}//namespace ao
}//namespace mx

#endif // temporalFourierLP_hpp
