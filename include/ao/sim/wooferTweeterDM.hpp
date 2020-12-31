/** \file 
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief 
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef _wooferTweeterDM_hpp__
#define _wooferTweeterDM_hpp__

#include "deformableMirror.hpp"

namespace mx
{

namespace AO
{

namespace sim
{

struct wooferTweeterDMSpec
{
   std::string name;
   std::string basisName;
   
   deformableMirrorSpec woofer;
   deformableMirrorSpec tweeter;
   
   int wooferModes;
   
   wooferTweeterDMSpec()
   {
      wooferModes = 0;
   }
};

// template<typename floatT>
// struct wooferTweeterCommand
// {
//    bool woofer;
//    Eigen::Array< floatT, Eigen::Dynamic, Eigen::Dynamic> wooferVect;
//    
//    bool tweeter;
//    Eigen::Array< floatT, Eigen::Dynamic, Eigen::Dynamic> tweeterVect;
// };

template<typename _floatT>
class wooferTweeterDM
{
public:   

   typedef _floatT floatT;
     
   typedef std::complex<floatT> complexT;
   
   ///The wavefront data type
   typedef wavefront<floatT>  wavefrontT;
   
   ///The pupil image type
   typedef Eigen::Array< floatT, Eigen::Dynamic, Eigen::Dynamic> imageT;
 
   typedef wooferTweeterDMSpec specT;
   
   //typedef wooferTweeterCommand<floatT> commandT;
   //typedef Eigen::Array< floatT, -1, -1> commandT;
   
   deformableMirror<floatT> woofer;
   deformableMirror<floatT> tweeter;
   
   int _wooferModes;

   std::string _name;
   std::string _basisName;

public:

   ///Default c'tor.
   wooferTweeterDM();
   
   int initialize( specT & spec, 
                   const std::string & pupil); 
      
   
   
   ///Get the calibration amplitude.
//   void calAmp(const std::vector<floatT> & ca);
   void calAmp(floatT ca); 
   
   void applyMode( wavefrontT & wf, 
                   int modeNo, 
                   floatT amp, 
                   floatT lambda );
   
   template<typename commandT>
   void setShape( commandT & commandV );
      
   void applyShape(wavefrontT & wf,  floatT lambda);

   int nModes()
   {
      return woofer.nModes() + tweeter.nModes();
   }
   
   double t0, t1, t_mm, t_sum;
   
   std::string name() {return _name;}
   std::string basisName() {return _basisName;}
   
};


template<typename _floatT>
wooferTweeterDM<_floatT>::wooferTweeterDM()
{
   ds9_interface_set_title(&woofer.ds9i_shape, "Woofer_Shape");
   ds9_interface_set_title(&woofer.ds9i_phase, "Woofer_Phase");
   ds9_interface_set_title(&woofer.ds9i_acts, "Woofer_Acts");
   ds9_interface_set_title(&tweeter.ds9i_shape, "Tweeter_Shape");
   ds9_interface_set_title(&tweeter.ds9i_phase, "Tweeter_Phase");
   ds9_interface_set_title(&tweeter.ds9i_acts, "Tweeter_Acts");
   
   t_mm = 0;
   t_sum = 0;
}

template<typename _floatT>
int wooferTweeterDM<_floatT>::initialize( specT & spec, 
                                           const std::string & pupil )
{

   _name = spec.name;
   _basisName = spec.basisName;
   
   woofer.initialize(spec.woofer, pupil);
      
   tweeter.initialize(spec.tweeter, pupil);
   
   _wooferModes = spec.wooferModes;
}
   
   
   


// template<typename _floatT>
// void wooferTweeterDM<_floatT>::calAmp(const std::vector<_floatT> & ca)
// {
//    woofer.calAmp(ca[0]);
//    tweeter.calAmp(ca[1]);
// }

template<typename _floatT>
void wooferTweeterDM<_floatT>::calAmp(_floatT ca)
{
   woofer.calAmp(ca);
   tweeter.calAmp(ca);
}


template<typename _floatT>
void wooferTweeterDM<_floatT>::applyMode( wavefrontT & wf, 
                                          int modeNo, 
                                          floatT amp, 
                                          floatT lambda )
{  

   if(modeNo < _wooferModes)
   {
      woofer.applyMode(wf, modeNo, amp, lambda);
   }
   else
   {
      tweeter.applyMode(wf, modeNo - _wooferModes, amp, lambda);
   }
}
   

template<typename _floatT>
template<typename commandT>
void wooferTweeterDM<_floatT>::setShape( commandT & commandV )
{
   //static int called = 0;
   static commandT avgWoofV;
   
   commandT woofV, tweetV;
   
   BREAD_CRUMB;
   woofV.measurement = commandV.measurement.block(0, 0, 1, _wooferModes);
   woofV.iterNo = commandV.iterNo;
   
   BREAD_CRUMB;
//    std::cerr << _wooferModes << "\n";
//    std::cerr << commandV.cols() << "\n";
   
   tweetV.measurement = commandV.measurement.block(0, _wooferModes, 1,commandV.measurement.cols()-_wooferModes);
   tweetV.iterNo = commandV.iterNo;
   
   BREAD_CRUMB;
   
   woofer.setShape(woofV);
         
//    if(called == 0)
//    {
//       avgWoofV.measurement = woofV.measurement;
//       ++called;
//    }
//    else if(called == 1)
//    {
//       avgWoofV.measurement += woofV.measurement;
//       avgWoofV.measurement /= (called + 1);
//       
//       woofer.setShape(avgWoofV.measurement);
//       called = 0;
//    }
   
   
   BREAD_CRUMB;
   tweeter.setShape(tweetV);
   BREAD_CRUMB;
}

   
   
   
//    static int called = 0;
//    
//    if(commandV.woofer)
//    {
//       woofer.setShape(commandV.wooferVect);
//    }
//    
//    //if(commandV.tweeter)
//    if(called > 50)
//    {
//       tweeter.setShape(commandV.tweeterVect);
//    }
//    
//    ++called;
//}
   
template<typename _floatT>
void wooferTweeterDM<_floatT>::applyShape(wavefrontT & wf,  floatT lambda)
{
   
   woofer.applyShape(wf, lambda);
   tweeter.applyShape(wf, lambda);
   
}

} //sim
} //AO
} //namespace mx

#endif //__deformableMirror_hpp__

