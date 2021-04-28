/** \file 
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief 
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef __pyramidSensor_hpp__
#define __pyramidSensor_hpp__


#include "mx/imagingUtils.hpp"
#include "mx/fraunhoferImager.hpp"
#include "mx/timeUtils.hpp"
#include "mx/fitsFile.hpp"
//#include "mx/ds9_interface.h"
#include "../../math/constants.hpp"

#include "wavefront.hpp"

#include "mx/imageTransforms.hpp"

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
   
template<typename _floatT, typename _detectorT>
class pyramidSensor
{
public:   
   
   typedef _floatT floatT;
     
   typedef std::complex<floatT> complexT;
   
   ///The wavefront data type
   typedef wavefront<floatT> wavefrontT;

   ///The wavefront complex field type
   typedef mx::imagingArray<std::complex<floatT>,fftwAllocator<std::complex<floatT> >, 0> complexFieldT;
   
   ///The wavefront sensor detector image type
   typedef Eigen::Array< floatT, Eigen::Dynamic, Eigen::Dynamic> wfsImageT;
   
   typedef _detectorT detectorT;
   
protected:   
   
   /* Standard WFS Interface: */
   int _wfSz; ///< Size of the wavefront in pixels

   int _detRows; ///<The number of rows of the WFS detector.  After forming the image the WFS detector plane is re-binned to this.
   int _detCols; ///<The number of columns of the WFS detector.  After forming the image the WFS detector plane is re-binned to this.

   floatT _lambda; ///< Central wavelength, in meters

   int _iTime; ///<Integration time in loop steps
         
   int _roTime; ///<Readout time in loop steps

   floatT _simStep; ///<The simulation stepsize in seconds.

   
   /* PyWFS Specific: */
   
   floatT _wfPS; ///< Wavefront pixel scale, in meters/pixel
   
   floatT _D; ///< Telescope diameter, in meters
   
   int _modSteps; ///<Number of steps in the modulation simulation

   floatT _modRadius; ///<Radius of the modulation in pixels
   
   int _quadSz; ///<The size of the PyWFS quadrant
   
   
   mx::fraunhoferImager<complexFieldT> fi;
   
   bool tiltsMade;
   std::vector<complexFieldT> tilts;

   int _iTime_counter; 
   
   int _reading;
   
   int _roTime_counter;
   
   std::vector<wavefrontT> _wavefronts;
   
   int _lastWavefront;
   
   ds9_interface ds9i;
   
public:   
   ///Default c'tor
   pyramidSensor();

   /* The standard WFS interface: */
   
   detectorT detector;
   
   ///The image on the detector, resized from wfsImage
   wfsImageT detectorImage;

   ///Get the wavefront size in pixels
   /**
     * \returns the wavefront size in pixels
     */ 
   int wfSz();
   
   ///Set the wavefront size in pixels.
   /**
     * \param sz is the new size
     */ 
   void wfSz(int sz);
   
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
   void detSize(int nrows, int ncols);
   
   ///Get the PyWFS central wavelength
   /**
     * \returns the central wavelength in meters
     */
   floatT lambda();
   
   ///Set the PyWFS central wavelength
   /**
     * \param d is the new central wavelength in meters
     */
   void lambda(floatT l);
   
   ///Get the PyWFS integration time, in time steps
   int iTime(); 
   
   ///Set the PyWFS integration time, in time steps
   void iTime(int it);
   
   ///Get the PyWFS detector readout time, in time steps
   int roTime(); 
   
   ///Set the PyWFS detector readout time, in time steps
   void roTime(int rt);
   
   ///Get the simulation step-size, in seconds.
   floatT simStep(); 
   
   ///Set the simulation step-size, in seconds.
   void simStep(floatT st);
   
   template<typename AOSysT>
   void linkSystem(AOSysT & AOSys);
   
   
   ///Sense the wavefront aberrations
   /** Returns true if a new wavefront measurement is ready.
     * Retruns false if still integrating.
     */
   bool senseWavefront(wavefrontT & pupilPlane);
   
   ///Sense the wavefront aberrations in a calibration mode
   /** Allows for faster calibrations.
     */
   bool senseWavefrontCal(wavefrontT & pupilPlane);
   
   
   
   
   
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
   void wfPS(floatT ps);
   
   ///Get the telescope diameter
   /**
     * \returns the telescope diameter in meters
     */
   floatT D();
   
   ///Set the telescope diameter
   /**
     * \param d is the new size in meters
     */
   void D(floatT d);
   
   
   ///Get the number of modulation steps
   /**
     * \returns _modSteps;
     */ 
   int modSteps();
   
   ///Set the number of modulation steps
   /**
     * \param mSt is the new number of modulation steps
     */ 
   void modSteps(int mSt);

   ///Get the radius of modulation
   /**
     * \returns _modRadius;
     */
   floatT modRadius();

   ///Set the modulation radius
   /**
     * \param mR is the new modulation radius in lambda/D
     */    
   void modRadius(floatT mR);
   
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
   void quadSz(int sz);


protected:
   
   ///The image formed by the WFS
   wfsImageT wfsImage;
   
   void makeTilts();
   
   void doSenseWavefront();
   void doSenseWavefrontNoMod(wavefrontT &);

   bool firstRun;
   
};

template<typename _floatT,  typename _detectorT>
pyramidSensor<_floatT, _detectorT>::pyramidSensor()
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
   
   tiltsMade = false;
   
   firstRun = true;
   
   ds9_interface_init(&ds9i);
   ds9_interface_set_title(&ds9i, "PyWFS");
}

template<typename _floatT,  typename _detectorT>
int pyramidSensor<_floatT, _detectorT>::wfSz()
{
   return _wfSz;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::wfSz(int sz)
{
   if( _wfSz == sz) return;
   
   _wfSz = sz;
               
   fi.setWavefrontSizePixels(_wfSz);

   tiltsMade = false;
}

template<typename _floatT,  typename _detectorT>
int pyramidSensor<_floatT, _detectorT>::detRows()
{
   return _detRows;
}


template<typename _floatT,  typename _detectorT>
int pyramidSensor<_floatT, _detectorT>::detCols()
{
   return _detCols;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::detSize(int nrows, int ncols)
{
   if( _detRows == nrows && _detCols == ncols) return;
   
   _detRows = nrows;
   _detCols = ncols;
   
   detector.setSize(_detRows, _detCols);
   detectorImage.resize(_detRows, _detCols);
   
   
}

template<typename _floatT,  typename _detectorT>
_floatT pyramidSensor<_floatT, _detectorT>::lambda()
{
   return _lambda;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::lambda(_floatT l)
{
   _lambda = l;
}


template<typename _floatT,  typename _detectorT>
template<typename AOSysT>
void pyramidSensor<_floatT, _detectorT>::linkSystem(AOSysT & AOSys)
{
   AOSys.wfs.wfPS(AOSys._wfPS);
   AOSys.wfs.D(AOSys._D);

}

template<typename _floatT,  typename _detectorT>
int pyramidSensor<_floatT, _detectorT>::wfPS()
{
   return _wfPS;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::wfPS(_floatT ps)
{
   _wfPS = ps;
}

template<typename _floatT,  typename _detectorT>
_floatT pyramidSensor<_floatT, _detectorT>::D()
{
   return _D;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::D(_floatT d)
{
   _D = d;
}



template<typename _floatT,  typename _detectorT>
int pyramidSensor<_floatT, _detectorT>::modSteps()
{
   return _modSteps;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::modSteps(int mSt)
{
   _modSteps = mSt;
   
   tiltsMade = false;
}

template<typename _floatT,  typename _detectorT>
_floatT pyramidSensor<_floatT, _detectorT>::modRadius()
{
   return _modRadius;// * mx::fftPlateScale(_wfSz, _wfPS, _lambda)/(_lambda/D);
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::modRadius(_floatT mR)
{
   _modRadius = mR ;// (_lambda/D)/mx::fftPlateScale(_wfSz, _wfPS, _lambda);
   
   tiltsMade = false;
}

template<typename _floatT,  typename _detectorT>
int pyramidSensor<_floatT, _detectorT>::quadSz()
{
   return _quadSz;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::quadSz(int sz)
{
   if( _quadSz == sz) return;
      
   _quadSz = sz;
   
   wfsImage.resize(2*_quadSz, 2*_quadSz);
}


template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::makeTilts()
{
   constexpr _floatT pi = math::pi<_floatT>();
   
   _floatT dang = 2*pi/_modSteps;
   _floatT dx, dy;

   tilts.resize(_modSteps);
   
   std::cout << "WF Size: " << _wfSz << "\n";
   std::cout << "WF PS:   " << _wfPS << "\n";
   std::cout << "Lambda:  " << _lambda << "\n";
   
   std::cout << "Pyr. PS: " << mx::fftPlateScale<floatT>(_wfSz, _wfPS, _lambda)*206265. << " (mas/pix)\n";
   std::cout << "Mod rad: " << _modRadius * (_lambda/_D)/mx::fftPlateScale<floatT>(_wfSz, _wfPS, _lambda) << " (pixels)\n";
   
   for(int i=0; i < _modSteps; ++i)
   { 
      dx = _modRadius * (_lambda/_D) / mx::fftPlateScale<floatT>(_wfSz, _wfPS, _lambda) * cos(0.5*dang+dang * i);
      dy = _modRadius * (_lambda/_D) /  mx::fftPlateScale<floatT>(_wfSz, _wfPS, _lambda) * sin(0.5*dang+dang * i);

      tilts[i].resize(_wfSz, _wfSz);
      tilts[i].set(std::complex<_floatT>(0,1));
    
      tiltWavefront(tilts[i], dx, dy);
   }
   
   tiltsMade = true;
}

template<typename _floatT,  typename _detectorT>
int pyramidSensor<_floatT, _detectorT>::iTime()
{
   return _iTime;
}
   
template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::iTime(int it)
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

template<typename _floatT,  typename _detectorT>
int pyramidSensor<_floatT, _detectorT>::roTime()
{
   return roTime;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::roTime(int rt)
{
   if(rt < 1)
   {
      std::cerr << "roTime must be >= 1.  Correcting\n";
      rt = 1;
   }
   
   _roTime = rt;
   

}

template<typename _floatT,  typename _detectorT>
_floatT pyramidSensor<_floatT, _detectorT>::simStep()
{
   return simStep;
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::simStep(_floatT st)
{
   
   _simStep = st;
   
   detector.expTime(_simStep*_iTime);
   

}


template<typename _floatT,  typename _detectorT>
bool pyramidSensor<_floatT, _detectorT>::senseWavefront(wavefrontT & pupilPlane)
{
   
   ++_lastWavefront;
   if(_lastWavefront >= _wavefronts.size()) _lastWavefront = 0;
   _wavefronts[_lastWavefront].amplitude = pupilPlane.amplitude;
   _wavefronts[_lastWavefront].phase = pupilPlane.phase;
   
   
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
         std::cerr << "PyWFS: reading\n";
         detector.exposeImage(detectorImage, wfsImage);
         ds9_interface_display_raw( &ds9i, 1, detectorImage.data(), detectorImage.rows(), detectorImage.cols(),1, mx::getFitsBITPIX<floatT>());
         
         _roTime_counter = 0;
         _reading=0;
         rv = true;
      }
   }
   
   if( _iTime_counter >= _iTime)
   {
      std::cerr << "PyWFS: sensing\n";
      doSenseWavefront();
      _iTime_counter = 0;
      
      _reading = 1;
      _roTime_counter = 0;
   }



   return rv;
   
}

template<typename _floatT,  typename _detectorT>
bool pyramidSensor<_floatT, _detectorT>::senseWavefrontCal(wavefrontT & pupilPlane)
{
   
   _lastWavefront =1;

   _wavefronts[0].amplitude = pupilPlane.amplitude;
   _wavefronts[0].phase = pupilPlane.phase;
   
   _wavefronts[1].amplitude = pupilPlane.amplitude;
   _wavefronts[1].phase = pupilPlane.phase;
   
   doSenseWavefront();
      
   detector.exposeImage(detectorImage, wfsImage);

   //ds9_interface_display_raw( &ds9i, 1, detectorImage.data(), detectorImage.rows(), detectorImage.cols(),1, mx::getFitsBITPIX<floatT>());
   
   return true;
   
}

template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::doSenseWavefront()
{
   constexpr _floatT pi = math::pi<_floatT>();
 
   if(tiltsMade == false) makeTilts();


   //std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
   BREAD_CRUMB;
   
   wavefrontT pupilPlane;
   
   /* Here make average wavefront for now */
   int _firstWavefront = _lastWavefront - _iTime;
   if(_firstWavefront < 0) _firstWavefront += _wavefronts.size();
   
   pupilPlane.amplitude = _wavefronts[_firstWavefront].amplitude;
   pupilPlane.phase = _wavefronts[_firstWavefront].phase;
   
   //std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
   BREAD_CRUMB;
   
   for(int i=0; i<_iTime; ++i)
   {
      ++_firstWavefront;
      if(_firstWavefront >= _wavefronts.size()) _firstWavefront = 0;
      
      pupilPlane.amplitude += _wavefronts[_firstWavefront].amplitude;
      pupilPlane.phase += _wavefronts[_firstWavefront].phase;      
   }
   
   pupilPlane.amplitude /= (_iTime+1);
   pupilPlane.phase /= (_iTime+1);
   
   /*=====================================*/

   
   if(_modRadius == 0) return doSenseWavefrontNoMod(pupilPlane);
   
   //std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
   BREAD_CRUMB;
   
   wfsImage.resize(2*_quadSz, 2*_quadSz);
   wfsImage.setZero();
         
   complexFieldT pupilPlaneCF;
   pupilPlane.getWavefront(pupilPlaneCF, _wfSz);
   
   complexFieldT focalPlaneCF;
   focalPlaneCF.resize(_wfSz, _wfSz);
   
   #pragma omp parallel 
   {
      complexFieldT tiltedPlane;
      complexFieldT focalPlane;
      complexFieldT sensorPlane;
      complexFieldT temp_ul;
      complexFieldT temp_ur;
      complexFieldT temp_ll;
      complexFieldT temp_lr;
      
      wfsImageT sensorImage;
         
      tiltedPlane.resize(_wfSz, _wfSz);
      focalPlane.resize(_wfSz, _wfSz);
      sensorPlane.resize(_wfSz, _wfSz);
      
      temp_ul.resize(_wfSz, _wfSz);
      temp_ul.set((complexT)0);
      
      temp_ur.resize(_wfSz, _wfSz);
      temp_ur.set((complexT)0);
      
      temp_ll.resize(_wfSz, _wfSz);
      temp_ll.set((complexT)0);
      
      temp_lr.resize(_wfSz, _wfSz);
      temp_lr.set((complexT)0);
      
      
      //sensorImage.resize(_wfSz, _wfSz);
      sensorImage.resize(_quadSz*2, _quadSz*2);
      sensorImage.setZero();

      int dSz = (0.5*_wfSz-0.5*_quadSz + 2*_quadSz) - _wfSz;
      if(dSz  < 0) dSz = 0;
      
      
      #pragma omp for 
      for(int i=0; i < _modSteps; ++i)
      { 
         int nelem = _wfSz*_wfSz;
         
         complexT * tp_data = tiltedPlane.data();
         complexT * pp_data = pupilPlaneCF.data();
         complexT * ti_data = tilts[i].data();
                  
         for(int ii=0; ii< nelem; ++ii)
         {
            tp_data[ii] = pp_data[ii]*ti_data[ii];
         }
         
         fi.propagatePupilToFocal(focalPlane, tiltedPlane);            
         temp_ul.set(complexT(0,0));
         extractBlock(temp_ul, 0,0.5*_wfSz, 0, 0.5*_wfSz, focalPlane, 0, 0);
         
         fi.propagateFocalToPupil(sensorPlane, temp_ul);
         extractIntensityImageAccum(sensorImage, 0,2*_quadSz-dSz, 0,2*_quadSz-dSz, sensorPlane, 0.5*_wfSz-0.5*_quadSz, 0.5*_wfSz-0.5*_quadSz);

         temp_ur.set( complexT(0,0));
         extractBlock(temp_ur, 0.5*_wfSz,0.5*_wfSz, 0, 0.5*_wfSz, focalPlane, 0.5*_wfSz, 0);
         fi.propagateFocalToPupil(sensorPlane, temp_ur);
         extractIntensityImageAccum(sensorImage, dSz,2*_quadSz-dSz, 0, 2*_quadSz-dSz, sensorPlane, 0.5*_wfSz-1.5*_quadSz+dSz, 0.5*_wfSz-0.5*_quadSz);

         temp_ll.set( complexT(0,0));
         extractBlock(temp_ll, 0, 0.5*_wfSz, 0.5*_wfSz, 0.5*_wfSz, focalPlane, 0, 0.5*_wfSz);
         fi.propagateFocalToPupil(sensorPlane, temp_ll);
         extractIntensityImageAccum(sensorImage, 0,2*_quadSz-dSz, dSz, 2*_quadSz-dSz, sensorPlane, 0.5*_wfSz-0.5*_quadSz, 0.5*_wfSz-1.5*_quadSz+dSz);

         temp_lr.set( complexT(0,0));
         extractBlock(temp_lr, 0.5*_wfSz,0.5*_wfSz, 0.5*_wfSz, 0.5*_wfSz, focalPlane, 0.5*_wfSz, 0.5*_wfSz);
         fi.propagateFocalToPupil(sensorPlane, temp_lr);
         extractIntensityImageAccum(sensorImage, dSz,2*_quadSz-dSz, dSz, 2*_quadSz-dSz, sensorPlane, 0.5*_wfSz-1.5*_quadSz+dSz, 0.5*_wfSz-1.5*_quadSz+dSz);
         
         
         
      }//for
      
      BREAD_CRUMB;
                  
      //std::cerr << wfsImage.rows() << " " << wfsImage.cols() << "\n";
      //std::cerr << sensorImage.rows() << " " << sensorImage.cols() << "\n";

      #pragma omp critical
      {
         wfsImage += sensorImage;
      }
      
      
   }//#pragma omp parallel
   
   BREAD_CRUMB;
   
   wfsImage /= _modSteps;
}


template<typename _floatT,  typename _detectorT>
void pyramidSensor<_floatT, _detectorT>::doSenseWavefrontNoMod(wavefrontT & pupilPlane)
{
   constexpr _floatT pi = math::pi<_floatT>();
 

   BREAD_CRUMB;
   
   wfsImage.resize(2*_quadSz, 2*_quadSz);
   wfsImage.setZero();
         
   complexFieldT pupilPlaneCF;
   pupilPlane.getWavefront(pupilPlaneCF, _wfSz);
   
   complexFieldT focalPlaneCF;
   focalPlaneCF.resize(_wfSz, _wfSz);
   
      complexFieldT tiltedPlane;
      complexFieldT focalPlane;
      complexFieldT sensorPlane;
      complexFieldT temp_ul;
      complexFieldT temp_ur;
      complexFieldT temp_ll;
      complexFieldT temp_lr;
      
      wfsImageT sensorImage;
         
      tiltedPlane.resize(_wfSz, _wfSz);
      focalPlane.resize(_wfSz, _wfSz);
      sensorPlane.resize(_wfSz, _wfSz);
      
      temp_ul.resize(_wfSz, _wfSz);
      temp_ul.set((complexT)0);
      
      temp_ur.resize(_wfSz, _wfSz);
      temp_ur.set((complexT)0);
      
      temp_ll.resize(_wfSz, _wfSz);
      temp_ll.set((complexT)0);
      
      temp_lr.resize(_wfSz, _wfSz);
      temp_lr.set((complexT)0);
      
      sensorImage.resize(_quadSz*2, _quadSz*2);
      sensorImage.setZero();

      int dSz = (0.5*_wfSz-0.5*_quadSz + 2*_quadSz) - _wfSz;
      if(dSz  < 0) dSz = 0;
         int nelem = _wfSz*_wfSz;
         
         complexT * tp_data = tiltedPlane.data();
         complexT * pp_data = pupilPlaneCF.data();
         
         BREAD_CRUMB;
         for(int ii=0; ii< nelem; ++ii)
         {
            tp_data[ii] = pp_data[ii];//*ti_data[ii];
         }
         
         BREAD_CRUMB;
         
         fi.propagatePupilToFocal(focalPlane, tiltedPlane);            
         temp_ul.set(complexT(0,0));
         extractBlock(temp_ul, 0,0.5*_wfSz, 0, 0.5*_wfSz, focalPlane, 0, 0);
         
         fi.propagateFocalToPupil(sensorPlane, temp_ul);
         extractIntensityImageAccum(sensorImage, 0,2*_quadSz-dSz, 0,2*_quadSz-dSz, sensorPlane, 0.5*_wfSz-0.5*_quadSz, 0.5*_wfSz-0.5*_quadSz);

         temp_ur.set( complexT(0,0));
         extractBlock(temp_ur, 0.5*_wfSz,0.5*_wfSz, 0, 0.5*_wfSz, focalPlane, 0.5*_wfSz, 0);
         fi.propagateFocalToPupil(sensorPlane, temp_ur);
         extractIntensityImageAccum(sensorImage, dSz,2*_quadSz-dSz, 0, 2*_quadSz-dSz, sensorPlane, 0.5*_wfSz-1.5*_quadSz+dSz, 0.5*_wfSz-0.5*_quadSz);

         temp_ll.set( complexT(0,0));
         extractBlock(temp_ll, 0, 0.5*_wfSz, 0.5*_wfSz, 0.5*_wfSz, focalPlane, 0, 0.5*_wfSz);
         fi.propagateFocalToPupil(sensorPlane, temp_ll);
         extractIntensityImageAccum(sensorImage, 0,2*_quadSz-dSz, dSz, 2*_quadSz-dSz, sensorPlane, 0.5*_wfSz-0.5*_quadSz, 0.5*_wfSz-1.5*_quadSz+dSz);

         temp_lr.set( complexT(0,0));
         extractBlock(temp_lr, 0.5*_wfSz,0.5*_wfSz, 0.5*_wfSz, 0.5*_wfSz, focalPlane, 0.5*_wfSz, 0.5*_wfSz);
         fi.propagateFocalToPupil(sensorPlane, temp_lr);
         extractIntensityImageAccum(sensorImage, dSz,2*_quadSz-dSz, dSz, 2*_quadSz-dSz, sensorPlane, 0.5*_wfSz-1.5*_quadSz+dSz, 0.5*_wfSz-1.5*_quadSz+dSz);
         
         
         
      BREAD_CRUMB;
      
      
      wfsImage += sensorImage;
   
      
      
   
   BREAD_CRUMB;
   
}

} //namespace sim 
} //namespace AO
} //namespace mx

#endif //__pyramidSensor_hpp__

