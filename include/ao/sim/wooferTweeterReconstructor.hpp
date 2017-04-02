#ifndef __wooferTweeterReconstructor_hpp__
#define __wooferTweeterReconstructor_hpp__

#pragma GCC system_header
#include <Eigen/Dense>

#include "wooferTweeterDM.hpp"

namespace mx
{
namespace AO
{
namespace sim
{

template< typename wooferReconT, typename tweeterReconT >
struct wooferTweeterReconstructorSpec
{
   typename wooferReconT::specT woofer;
   typename tweeterReconT::specT tweeter;
};



template< typename wooferReconT, typename tweeterReconT > 
class wooferTweeterReconstructor
{
public:
   
   typedef typename wooferReconT::floatT floatT;
   
   ///The type of the measurement (i.e. the slope vector)
   //typedef wooferTweeterCommand<floatT> measurementT;
   
   ///The type of the WFS image
   typedef Eigen::Array<floatT, -1, -1> imageT;
    
   typedef wooferTweeterReconstructorSpec<wooferReconT, tweeterReconT> specT;
   

   
   wooferReconT _woofer;
   tweeterReconT _tweeter;
   
public:   
   ///Default c'tor
   wooferTweeterReconstructor();
   
   template<typename AOSysT>
   void initialize(AOSysT & AOSys, specT & spec)
   {
      _woofer.initialize(AOSys, spec.woofer);
      _tweeter.initialize(AOSys, spec.tweeter);
      
   }

   ///Get the calibration amplitude used in response matrix acquisition (_calAmp)
   std::vector<floatT> calAmp();
   
#if 0   
   ///Calculate the slope measurement
   /**
     * \param slopes [out] a (_measurementSize X 2)  array of slopes
     * \param wfsImage [in] the WFS image from which to measure the slopes
     */   
   void calcMeasurement(measurementT & slopes, imageT & wfsImage);
#endif

   ///Reconstruct the wavefront from the input image, producing the modal amplitude vector 
   template<typename measurementT, typename wfsImageT>
   void reconstruct( measurementT & commandVect, 
                     wfsImageT & wfsImage);

#if 0  
   ///Initialize the response matrix for acquisition
   /** 
     * \param nmodes the number of modes 
     * \param calamp the calibration amplitude
     * \param detrows the number of detector rows
     * \param detcols the number of detector columns
     */ 
   void initializeRMat(int nmodes, floatT calamp, int detrows,int detcols);
   
   ///Accumalte the next measurement in the response matrix
   /** 
     * \param i the measurement index
     * \param measureVec is the i-th measurement vector
     */ 
   void accumulateRMat(int i, measurementT &measureVec);
   void accumulateRMat(int i, measurementT &measureVec, imageT & wfsImage);
   
   ///Write the accumulated response matrix to disk
   /**
     * \param fname the name, including path, of the response matrix
     */ 
   void saveRMat(std::string fname);
      
   void saveRImages(std::string fname);
      
#endif   
   
};


template< typename wooferReconT, typename tweeterReconT >
wooferTweeterReconstructor<wooferReconT, tweeterReconT>::wooferTweeterReconstructor()
{
}


template< typename wooferReconT, typename tweeterReconT >
std::vector<typename wooferReconT::floatT> wooferTweeterReconstructor<wooferReconT, tweeterReconT>::calAmp()
{
   std::vector<typename wooferReconT::floatT> ca(2);
   
   ca[0] = _woofer.calAmp();
   ca[1] = _tweeter.calAmp();
   
   return ca;
   
}


#if 0
template< typename wooferReconT, typename tweeterReconT >
void wooferTweeterReconstructor<wooferReconT, tweeterReconT>::calcMeasurement(measurementT & slopes, imageT & wfsImage)
{
}
#endif

template< typename wooferReconT, typename tweeterReconT >
template< typename measurementT, typename wfsImageT >
void wooferTweeterReconstructor<wooferReconT, tweeterReconT>::reconstruct(measurementT & commandVect, wfsImageT & wfsImage)
{
   measurementT woofV;
   _woofer.reconstruct( woofV, wfsImage);
   
   measurementT tweetV;
   _tweeter.reconstruct( tweetV, wfsImage);
   
   commandVect.measurement.resize(1, woofV.measurement.cols() + tweetV.measurement.cols());
   
   
   for(int i=0; i< woofV.measurement.cols(); ++i)
   {
      commandVect.measurement(0,i) = woofV.measurement(0,i);
   }
   
   for(int i= 0; i< tweetV.measurement.cols(); ++i)
   {
      commandVect.measurement(0, woofV.measurement.cols()+i) = tweetV.measurement(0,i);
   }
   
   
}

#if 0
template< typename wooferReconT, typename tweeterReconT >
void wooferTweeterReconstructor<wooferReconT, tweeterReconT>::initializeRMat(int nModes, floatT calamp, int detRows, int detCols)
{
}

template< typename wooferReconT, typename tweeterReconT >
void wooferTweeterReconstructor<wooferReconT, tweeterReconT>::accumulateRMat(int i, measurementT &measureVec)
{

}

template< typename wooferReconT, typename tweeterReconT >
void wooferTweeterReconstructor<wooferReconT, tweeterReconT>::accumulateRMat(int i, measurementT &measureVec, imageT & wfsImage)
{

}

template< typename wooferReconT, typename tweeterReconT >
void wooferTweeterReconstructor<wooferReconT, tweeterReconT>::saveRMat(std::string fname)
{

}

template< typename wooferReconT, typename tweeterReconT >
void wooferTweeterReconstructor<wooferReconT, tweeterReconT>::saveRImages(std::string fname)
{
}
#endif

} //namespace sim 
} //namespace AO
} //namespace mx

#endif //__wooferTweeterReconstructor_hpp__

