/** \file pyramidSensor.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declaration and definition of a standard 4 quadrant pyramid WFS.
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef __pyramidSensor_hpp__
#define __pyramidSensor_hpp__


#include <mx/wfp/imagingUtils.hpp>
#include <mx/wfp/fraunhoferPropagator.hpp>
#include <mx/timeUtils.hpp>
#include <mx/improc/fitsFile.hpp>
#include <mx/improc/ds9Interface.hpp>

#include "wavefront.hpp"

#include <mx/improc/imageTransforms.hpp>

#ifdef DEBUG
#define BREAD_CRUMB std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n"; 
#else
#define BREAD_CRUMB 
#endif


namespace mx
{
namespace AO 
{
namespace sim 
{
 
template<typename _realT>
struct wfsImageT
{
   typedef _realT realT;
   
   unsigned iterNo;
   
   ///The wavefront sensor detector image type
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;
   
   imageT image;
   
   imageT tipImage;
};

template<typename _realT, typename _detectorT>
class pyramidSensor
{
public:   
   
   typedef _realT realT;
     
   typedef std::complex<realT> complexT;
   
   ///The wavefront data type
   typedef wavefront<realT> wavefrontT;

   ///The wavefront complex field type
   typedef mx::imagingArray<std::complex<realT>,fftwAllocator<std::complex<realT> >, 0> complexFieldT;
   
   ///The wavefront sensor detector image type
   //typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> wfsImageT;
   
   typedef _detectorT detectorT;
   
protected:   
   
   /* Standard WFS Interface: */
   int _wfSz; ///< Size of the wavefront in pixels

   int _detRows; ///<The number of rows of the WFS detector.  After forming the image the WFS detector plane is re-binned to this.
   int _detCols; ///<The number of columns of the WFS detector.  After forming the image the WFS detector plane is re-binned to this.

   realT _lambda; ///< Central wavelength, in meters

public:
   std::vector<realT> _wavelengths; ///< Vector of wavelengths in the WFS bandpass
   std::vector<realT> _wavelengthWeights; ///< The relative weights of the wavelengths
   
protected:
   
   int _iTime; ///<Integration time in loop steps
         
   int _roTime; ///<Readout time in loop steps

   realT _simStep; ///<The simulation stepsize in seconds.

   
   /* PyWFS Specific: */
   
   realT _wfPS; ///< Wavefront pixel scale, in meters/pixel
   
   realT _D; ///< Telescope diameter, in meters
   
   int _modSteps; ///<Number of steps in the modulation simulation

   realT _modRadius; ///<Radius of the modulation in pixels
   
   int _quadSz; ///<The size of the PyWFS quadrant
   
   
   wfp::fraunhoferPropagator<complexFieldT> fi;
   
   bool _opdMaskMade;
   complexFieldT _opdMask;
   
   bool tiltsMade;
   std::vector<complexFieldT> tilts;

   int _iTime_counter; 
   
   int _reading;
   
   int _roTime_counter;
   
   std::vector<wavefrontT> _wavefronts;
   
   int _lastWavefront;
   
   improc::ds9Interface ds9i;
   
public:   
   ///Default c'tor
   pyramidSensor();

   /* The standard WFS interface: */
   
   detectorT detector;
   
   ///The image on the detector, resized from wfsImage
   wfsImageT<realT> detectorImage;

   ///Get the wavefront size in pixels
   /**
     * \returns the wavefront size in pixels
     */ 
   int wfSz();
   
   ///Set the wavefront size in pixels.
   /**
     * \param sz is the new size
     */ 
   void wfSz( int sz /**< */ );
   
   ///Get the detector rows  in pixels
   /**
     * \returns _detRows
     */ 
   int detRows();
      
   ///Get the detector columns  in pixels
   /**
     * \returns _detCols
     */ 
   int detCols();
   
   ///Set the detector columns in pixels.
   /**
     * \param sz is the new size
     */ 
   void detSize( int nrows, ///< 
                 int ncols ///<
               );
   
   ///Get the PyWFS central wavelength
   /**
     * \returns the central wavelength in meters
     */
   realT lambda();
   
   ///Set the PyWFS central wavelength
   /**
     * \param d is the new central wavelength in meters
     */
   void lambda(realT l /**< */ );
   
   ///Get the PyWFS integration time, in time steps
   int iTime(); 
   
   ///Set the PyWFS integration time, in time steps
   void iTime(int it /**< */ );
   
   ///Get the PyWFS detector readout time, in time steps
   int roTime(); 
   
   ///Set the PyWFS detector readout time, in time steps
   void roTime(int rt /**< */ );
   
   ///Get the simulation step-size, in seconds.
   realT simStep(); 
   
   ///Set the simulation step-size, in seconds.
   void simStep(realT st /**< */ );
   
   template<typename AOSysT>
   void linkSystem(AOSysT & AOSys /**< */ );
   
   
   ///Sense the wavefront aberrations
   /** Returns true if a new wavefront measurement is ready.
     * Retruns false if still integrating.
     */
   bool senseWavefront(wavefrontT & pupilPlane /**< */ );
   
   ///Sense the wavefront aberrations in a calibration mode
   /** Allows for faster calibrations.
     */
   bool senseWavefrontCal(wavefrontT & pupilPlane /**< */ );
   
   
   
   
   
   /* The PyWFS  Specific Interface */
   
   ///Get the wavefront pixel scale in meters per pixel
   /**
     * \returns the wavefront pixel scale in meters/pixel
     */ 
   int wfPS();
   
   ///Set the wavefront pixel scale in meters per pixel
   /**
     * \param ps is the pixel scale
     */ 
   void wfPS(realT ps /**< */ );
   
   ///Get the telescope diameter
   /**
     * \returns the telescope diameter in meters
     */
   realT D();
   
   ///Set the telescope diameter
   /**
     * \param d is the new size in meters
     */
   void D(realT d /**< */ );
   
   
   ///Get the number of modulation steps
   /**
     * \returns _modSteps;
     */ 
   int modSteps();
   
   ///Set the number of modulation steps
   /**
     * \param mSt is the new number of modulation steps
     */ 
   void modSteps(int mSt /**< */ );

   ///Get the radius of modulation
   /**
     * \returns _modRadius;
     */
   realT modRadius();

   ///Set the modulation radius
   /**
     * \param mR is the new modulation radius in lambda/D
     */    
   void modRadius(realT mR /**< */ );
   
   ///Get the quadrant size in pixels
   /** This is the size of the quadrant in un-binned wavefront space
     * 
     * 
     * \returns _quadSz
     */ 
   int quadSz();
   
   ///Set the quadrant size in pixels.
   /** This is the size of the quadrant in un-binned wavefront space
     * It should be at least the size of the Pupil.  Make larger than the pupil
     * to have smaller pupil images on the PyWFS detector.
     * 
     * \param sz is the new size
     */ 
   void quadSz(int sz /**< */ );


   template<typename pupilT>
   void makeRefTipImage( pupilT & pupil);
   
   int ref;
   
protected:
   
   ///The image formed by the WFS
   wfsImageT<realT> wfsImage;
   
   wfsImageT<realT> wfsTipImage;
   
   wfsImageT<realT> refTipImage;
   
   void makeOpdMask();
   
   void makeTilts();
   
   void doSenseWavefront();
   void doSenseWavefrontNoMod(wavefrontT &  /**< */ );

   bool firstRun;
   
};

template<typename _realT,  typename _detectorT>
pyramidSensor<_realT, _detectorT>::pyramidSensor()
{
   _wfSz = 0;
   _detRows = 0;
   _detCols = 0;
   _lambda = 0;
   
   
   _modSteps = 16;
   _modRadius = 16;
   
   
   iTime(1);
   _iTime_counter = 0;
   
   _reading =0;
   _roTime = 1;
   _roTime_counter = 0;
   
   _simStep = 0.001;
   
   _opdMaskMade = false;
   tiltsMade = false;
   
   firstRun = true;
   
   ds9i.title("PyWFS");
   
   ref = 1;
}

template<typename _realT,  typename _detectorT>
int pyramidSensor<_realT, _detectorT>::wfSz()
{
   return _wfSz;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::wfSz(int sz)
{
   if( _wfSz == sz) return;
   
   _wfSz = sz;
               
   fi.setWavefrontSizePixels(_wfSz);

   tiltsMade = false;
   _opdMaskMade = false;
}

template<typename _realT,  typename _detectorT>
int pyramidSensor<_realT, _detectorT>::detRows()
{
   return _detRows;
}


template<typename _realT,  typename _detectorT>
int pyramidSensor<_realT, _detectorT>::detCols()
{
   return _detCols;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::detSize(int nrows, int ncols)
{
   if( _detRows == nrows && _detCols == ncols) return;
   
   _detRows = nrows;
   _detCols = ncols;
   
   detector.setSize(_detRows, _detCols);
   detectorImage.image.resize(_detRows, _detCols);
   
   
}

template<typename _realT,  typename _detectorT>
_realT pyramidSensor<_realT, _detectorT>::lambda()
{
   return _lambda;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::lambda(_realT l)
{
   _lambda = l;
}


template<typename _realT,  typename _detectorT>
template<typename AOSysT>
void pyramidSensor<_realT, _detectorT>::linkSystem(AOSysT & AOSys)
{
   AOSys.wfs.wfPS(AOSys._wfPS);
   AOSys.wfs.D(AOSys._D);

}

template<typename _realT,  typename _detectorT>
int pyramidSensor<_realT, _detectorT>::wfPS()
{
   return _wfPS;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::wfPS(_realT ps)
{
   _wfPS = ps;
}

template<typename _realT,  typename _detectorT>
_realT pyramidSensor<_realT, _detectorT>::D()
{
   return _D;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::D(_realT d)
{
   _D = d;
}



template<typename _realT,  typename _detectorT>
int pyramidSensor<_realT, _detectorT>::modSteps()
{
   return _modSteps;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::modSteps(int mSt)
{
   _modSteps = mSt;
   
   tiltsMade = false;
}

template<typename _realT,  typename _detectorT>
_realT pyramidSensor<_realT, _detectorT>::modRadius()
{
   return _modRadius;// * mx::fftPlateScale(_wfSz, _wfPS, _lambda)/(_lambda/D);
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::modRadius(_realT mR)
{
   _modRadius = mR ;// (_lambda/D)/mx::fftPlateScale(_wfSz, _wfPS, _lambda);
   
   tiltsMade = false;
}

template<typename _realT,  typename _detectorT>
int pyramidSensor<_realT, _detectorT>::quadSz()
{
   return _quadSz;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::quadSz(int sz)
{
   if( _quadSz == sz) return;
      
   _quadSz = sz;
   
   wfsImage.image.resize(2*_quadSz, 2*_quadSz);
   
   _opdMaskMade = false;
}

template<typename _realT,  typename _detectorT>
template<typename pupilT>
void pyramidSensor<_realT, _detectorT>::makeRefTipImage(pupilT & pupil)
{
   wavefrontT currWF;
   currWF.setAmplitude( pupil );
   currWF.setPhase( pupil*0 );

   refTipImage.image.resize(0,0);
   
   senseWavefrontCal(currWF);
   
   refTipImage = wfsTipImage;
   
   //ds9_interface_display_raw( &ds9i, 1, refTipImage.image.data(), refTipImage.image.rows(), refTipImage.image.cols(),1, mx::getFitsBITPIX<realT>());
   
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::makeOpdMask()
{
   complexFieldT opdMaskQ;

   _opdMask.resize(_wfSz, _wfSz);
   opdMaskQ.resize(  _wfSz, _wfSz);
   
   opdMaskQ.set(std::complex<_realT>(0,1));
   wfp::tiltWavefront(opdMaskQ, 0.5*_quadSz, 0.5*_quadSz);
   wfp::extractBlock(_opdMask, 0, 0.5*_wfSz, 0, 0.5*_wfSz, opdMaskQ, 0 , 0);

   opdMaskQ.set(std::complex<_realT>(0,1));
   wfp::tiltWavefront( opdMaskQ, -0.5*_quadSz, -0.5*_quadSz); 
   wfp::extractBlock(_opdMask, 0.5*_wfSz, 0.5*_wfSz, 0.5*_wfSz, 0.5*_wfSz, opdMaskQ, 0 , 0);
   
   opdMaskQ.set(std::complex<_realT>(0,1));
   wfp::tiltWavefront( opdMaskQ, 0.5*_quadSz, -0.5*_quadSz); 
   wfp::extractBlock(_opdMask, 0, 0.5*_wfSz, 0.5*_wfSz, 0.5*_wfSz, opdMaskQ, 0 , 0);
   
   opdMaskQ.set(std::complex<_realT>(0,1));
   wfp::tiltWavefront( opdMaskQ, -0.5*_quadSz, 0.5*_quadSz);
   wfp::extractBlock(_opdMask, 0.5*_wfSz, 0.5*_wfSz, 0, 0.5*_wfSz, opdMaskQ, 0 , 0);
   
   _opdMaskMade = true;

   
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::makeTilts()
{
   _realT pi = boost::math::constants::pi<_realT>();
   
   _realT dang = 2*pi/(_modSteps);
   _realT dx, dy;

   tilts.resize(_modSteps);
   
   std::cout << "WF Size: " << _wfSz << "\n";
   std::cout << "WF PS:   " << _wfPS << "\n";
   std::cout << "Lambda:  " << _lambda << "\n";
   
   std::cout << "Pyr. PS: " << wfp::fftPlateScale<realT>(_wfSz, _wfPS, _lambda)*206265. << " (mas/pix)\n";
   std::cout << "Mod rad: " << _modRadius * (_lambda/_D)/wfp::fftPlateScale<realT>(_wfSz, _wfPS, _lambda) << " (pixels)\n";
   
   for(int i=0; i < _modSteps; ++i)
   { 
      dx = _modRadius * (_lambda/_D) / wfp::fftPlateScale<realT>(_wfSz, _wfPS, _lambda) * cos(0.5*dang+dang * i);
      dy = _modRadius * (_lambda/_D) /  wfp::fftPlateScale<realT>(_wfSz, _wfPS, _lambda) * sin(0.5*dang+dang * i);

      tilts[i].resize(_wfSz, _wfSz);
      tilts[i].set(std::complex<_realT>(0,1));
    
      wfp::tiltWavefront(tilts[i], dx, dy);
   }
   
   tiltsMade = true;
}

template<typename _realT,  typename _detectorT>
int pyramidSensor<_realT, _detectorT>::iTime()
{
   return _iTime;
}
   
template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::iTime(int it)
{
   if(it < 1)
   {
      std::cerr << "iTime must be >= 1.  Correcting\n";
      it = 1;
   }
   
   _iTime = it;
   
   _wavefronts.resize(_iTime+2);
   _lastWavefront = -1;
   
   detector.expTime(_simStep*_iTime);
   
}

template<typename _realT,  typename _detectorT>
int pyramidSensor<_realT, _detectorT>::roTime()
{
   return roTime;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::roTime(int rt)
{
   if(rt < 1)
   {
      std::cerr << "roTime must be >= 1.  Correcting\n";
      rt = 1;
   }
   
   _roTime = rt;
   

}

template<typename _realT,  typename _detectorT>
_realT pyramidSensor<_realT, _detectorT>::simStep()
{
   return simStep;
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::simStep(_realT st)
{
   
   _simStep = st;
   
   detector.expTime(_simStep*_iTime);
   

}


template<typename _realT,  typename _detectorT>
bool pyramidSensor<_realT, _detectorT>::senseWavefront(wavefrontT & pupilPlane)
{
   
   ++_lastWavefront;
   if(_lastWavefront >= _wavefronts.size()) _lastWavefront = 0;
   _wavefronts[_lastWavefront].amplitude = pupilPlane.amplitude;
   _wavefronts[_lastWavefront].phase = pupilPlane.phase;
   _wavefronts[_lastWavefront].iterNo = pupilPlane.iterNo;
   
   //Always skip the first one for averaging to center of iTime.
   if(firstRun)
   {
      firstRun = false;
      return false;
   }
   
   ++_iTime_counter;

   
   bool rv = false;
   
   if(_reading)
   {
      ++_roTime_counter;
      
      if(_roTime_counter >= _roTime)
      {
         detector.exposeImage(detectorImage.image, wfsImage.image);
         
         detectorImage.tipImage = wfsTipImage.image;
         detectorImage.iterNo = wfsImage.iterNo;
         
         //ds9_interface_display_raw( &ds9i, 1, detectorImage.image.data(), detectorImage.image.rows(), detectorImage.image.cols(),1, mx::getFitsBITPIX<realT>());
         
         _roTime_counter = 0;
         _reading=0;
         rv = true;
      }
   }
   
   if( _iTime_counter >= _iTime)
   {
      doSenseWavefront();
      _iTime_counter = 0;
      
      _reading = 1;
      _roTime_counter = 0;
   }



   return rv;
   
}

template<typename _realT,  typename _detectorT>
bool pyramidSensor<_realT, _detectorT>::senseWavefrontCal(wavefrontT & pupilPlane)
{
   
   _lastWavefront =1;

   _wavefronts[0].amplitude = pupilPlane.amplitude;
   _wavefronts[0].phase = pupilPlane.phase;
   
   _wavefronts[1].amplitude = pupilPlane.amplitude;
   _wavefronts[1].phase = pupilPlane.phase;
   

   BREAD_CRUMB;
   
   doSenseWavefront();
      
   BREAD_CRUMB;
   
   detector.exposeImage(detectorImage.image, wfsImage.image);

   BREAD_CRUMB;
   detectorImage.tipImage = wfsTipImage.image;
   
   BREAD_CRUMB;
    //ds9_interface_display_raw( &ds9i, 1, wfsImage.data(), wfsImage.rows(), wfsImage.cols(),1, mx::getFitsBITPIX<realT>());
   //ds9_interface_display_raw( &ds9i, 1, detectorImage.data(), detectorImage.rows(), detectorImage.cols(),1, mx::getFitsBITPIX<realT>());
   
   return true;
   
}

template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::doSenseWavefront()
{ 
   if(tiltsMade == false) makeTilts();

   BREAD_CRUMB;
   
   wavefrontT pupilPlane;
   
   /* Here make average wavefront for now */
   int _firstWavefront = _lastWavefront - _iTime;
   if(_firstWavefront < 0) _firstWavefront += _wavefronts.size();
   
   pupilPlane.amplitude = _wavefronts[_firstWavefront].amplitude;
   pupilPlane.phase = _wavefronts[_firstWavefront].phase;
   
   realT avgIt = _wavefronts[_firstWavefront].iterNo;
   
   //std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
   BREAD_CRUMB;
   
   std::cerr << pupilPlane.amplitude.rows() << " " << _wavefronts[_firstWavefront].amplitude.rows() << std::endl;
   for(int i=0; i<_iTime; ++i)
   {
      ++_firstWavefront;
      if( (size_t) _firstWavefront >= _wavefronts.size()) _firstWavefront = 0;
      
      pupilPlane.amplitude += _wavefronts[_firstWavefront].amplitude;
      pupilPlane.phase += _wavefronts[_firstWavefront].phase;
      avgIt += _wavefronts[_firstWavefront].iterNo;
   }
   
   BREAD_CRUMB;
   
   pupilPlane.amplitude /= (_iTime+1);
   pupilPlane.phase /= (_iTime+1);
   
   avgIt /= (_iTime + 1.0);
   
   /*=====================================*/
   
   if(_modRadius == 0) return doSenseWavefrontNoMod(pupilPlane);
   
   BREAD_CRUMB;
   
   wfsImage.image.resize(2*_quadSz, 2*_quadSz);
   
   wfsImage.image.setZero();
         
   wfsTipImage.image.resize(64, 64);
   wfsTipImage.image.setZero();
   
   //---------------------------------------
   //  Check if wavelength vector is filled
   //---------------------------------------
   if( _wavelengths.size() == 0 )
   {
      _wavelengths.resize(1, _lambda);
      _wavelengthWeights.resize(1, 1.0);
   }
   
   
   if(!_opdMaskMade) makeOpdMask();
   
   
   
   for(size_t l = 0; l<_wavelengths.size(); ++l)
   {
      
   complexFieldT pupilPlaneCF;
   
   pupilPlane.lambda = _lambda;
   pupilPlane.getWavefront(pupilPlaneCF, _wavelengths[l], _wfSz); //_wavelengths[i]
   


   #pragma omp parallel 
   {
      complexFieldT tiltedPlane;
      complexFieldT tiltedPlaneDF;
      
      complexFieldT focalPlane;
      complexFieldT focalPlaneDefocus;
      
      complexFieldT sensorPlane;
      typename wfsImageT<realT>::imageT pyramidImage;
      
      typename wfsImageT<realT>::imageT sensorImage;
         
      tiltedPlane.resize(_wfSz, _wfSz);
      tiltedPlaneDF.resize(_wfSz, _wfSz);
      
      focalPlane.resize(_wfSz, _wfSz);
      focalPlaneDefocus.resize(_wfSz, _wfSz);
      
      sensorPlane.resize(_wfSz, _wfSz);
      
      pyramidImage.resize( 64, 64);
      pyramidImage.setZero();
      
      sensorImage.resize(_quadSz*2, _quadSz*2);
      sensorImage.setZero();

      #pragma omp for 
      for(int i=0; i < _modSteps; ++i)
      { 
         int nelem = _wfSz*_wfSz;

         complexT * tp_data = tiltedPlane.data();  
         complexT * tpDF_data = tiltedPlaneDF.data();
         
         complexT * pp_data = pupilPlaneCF.data();
         complexT * ti_data = tilts[i].data();
         complexT * fp_data = focalPlane.data();
         complexT * opd_data = _opdMask.data();
         
         //---------------------------------------------
         //Apply the modulationg tip 
         //---------------------------------------------
         for(int ii=0; ii< nelem; ++ii)
         {
            tp_data[ii] = pp_data[ii]*ti_data[ii];
            tpDF_data[ii] = tp_data[ii];
         }
         
         //---------------------------------------------
         //Propagate to Pyramid tip 
         //---------------------------------------------
         fi.propagatePupilToFocal(focalPlane, tiltedPlane);            
         
   /*      
         for(int ii =0; ii < _wfSz; ++ii)
         {
            for(int jj=0; jj < _wfSz; ++jj)
            {
               tiltedPlaneDF(ii,jj) = tiltedPlaneDF(ii,jj)*defocus(ii,jj);
            }
         }
         
         fi.propagatePupilToFocal(focalPlaneDefocus, tiltedPlaneDF);            
         
         extractIntensityImageAccum(pyramidImage, 0, 64, 0, 64, focalPlaneDefocus, 0.5*(_wfSz-1) - 0.5*(64-1), 0.5*(_wfSz-1) - 0.5*(64-1));
     */    
         
         
         
         
         //---------------------------------------------
         //Now apply the pyramid OPD 
         //---------------------------------------------
         for(int ii=0; ii< nelem; ++ii)
         {
            fp_data[ii] = fp_data[ii]*opd_data[ii];
         }
   
         //---------------------------------------------
         //Propagate to sensor plane
         //---------------------------------------------
         fi.propagateFocalToPupil(sensorPlane, focalPlane);
         
         //---------------------------------------------
         //Extract the image.
         //---------------------------------------------
         wfp::extractIntensityImageAccum(sensorImage, 0, 2*_quadSz, 0, 2*_quadSz, sensorPlane, 0.5*_wfSz-_quadSz, 0.5*_wfSz-_quadSz);
         
   
      }//for
      
      BREAD_CRUMB;
                  
      #pragma omp critical
      {
         wfsImage.image += sensorImage * _wavelengthWeights[l];
        // wfsTipImage.image += pyramidImage * _wavelengthWeights[l];
         
      }
      
      
   }//#pragma omp parallel
   
   } //l for wavelength

         
   
   
   BREAD_CRUMB;
   
   wfsImage.image /= _modSteps;
   wfsImage.iterNo = avgIt;
   
   
  /* wfsTipImage.image /= _modSteps;
   
   wfsTipImage.image /= wfsTipImage.image.sum();
   
   
   if (refTipImage.image.rows() == wfsTipImage.image.rows() && refTipImage.image.cols() == wfsTipImage.image.cols())
   {
      wfsTipImage.image -= refTipImage.image;
      
      if(!ref)
      {
         wfsTipImage.image *= 0;
      }
      else
      {
         wfsImage.image *= 0;
      }
      //wfsTipImage.image/=refTipImage.image;
   }
*/
   //ds9_interface_display_raw( &ds9i, 1, wfsTipImage.image.data(), wfsTipImage.image.rows(), wfsTipImage.image.cols(),1, mx::getFitsBITPIX<realT>());
   
}


template<typename _realT,  typename _detectorT>
void pyramidSensor<_realT, _detectorT>::doSenseWavefrontNoMod(wavefrontT & pupilPlane)
{
   BREAD_CRUMB;
   
   wfsImage.image.resize(2*_quadSz, 2*_quadSz);
   wfsImage.image.setZero();
         
   complexFieldT pupilPlaneCF;
   
   pupilPlane.getWavefront(pupilPlaneCF, _wfSz);
   
   if(!_opdMaskMade) makeOpdMask();
   
   complexFieldT tiltedPlane;
   complexFieldT focalPlane;
   complexFieldT sensorPlane;
     
   tiltedPlane.resize(_wfSz, _wfSz);
   focalPlane.resize(_wfSz, _wfSz);
   sensorPlane.resize(_wfSz, _wfSz);
      

   int nelem = _wfSz*_wfSz;
         
   complexT * tp_data = tiltedPlane.data();
   complexT * pp_data = pupilPlaneCF.data();
   complexT * opd_data = _opdMask.data();
   complexT * fp_data = focalPlane.data();
   
   BREAD_CRUMB;
   
   for(int ii=0; ii< nelem; ++ii)
   {
      tp_data[ii] = pp_data[ii];
   }
         
   BREAD_CRUMB;
         
   //---------------------------------------------
   //Propagate to Pyramid tip 
   //---------------------------------------------
   fi.propagatePupilToFocal(focalPlane, tiltedPlane);          

   
   BREAD_CRUMB;
   
   //---------------------------------------------
   //Now apply the pyramid OPD 
   //---------------------------------------------
   for(int ii=0; ii< nelem; ++ii)
   {
      fp_data[ii] = fp_data[ii]*opd_data[ii];
   }
   
   
   BREAD_CRUMB;
   
   //---------------------------------------------
   //Propagate to sensor plane
   //---------------------------------------------
   fi.propagateFocalToPupil(sensorPlane, focalPlane);
         
   BREAD_CRUMB;

   //---------------------------------------------
   //Extract the image.
   //---------------------------------------------
   wfp::extractIntensityImageAccum(wfsImage.image, 0, 2*_quadSz, 0, 2*_quadSz, sensorPlane, 0.5*_wfSz-_quadSz, 0.5*_wfSz-_quadSz);
   
   BREAD_CRUMB;
   
}

} //namespace sim 
} //namespace AO
} //namespace mx

#endif //__pyramidSensor_hpp__

