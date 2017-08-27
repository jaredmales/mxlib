/** \file directPhaseSensor.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declaration and definition of a direct phase sensor.
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef directPhaseSensor_hpp
#define directPhaseSensor_hpp

#include "../../improc/fitsFile.hpp"
#include "../../improc/eigenCube.hpp"
#include "../../improc/ds9Interface.hpp"

//#include <mx/fraunhoferImager.hpp>
#include "../../sigproc/psdFilter.hpp"

#include "wavefront.hpp"


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
};
   
using namespace boost::math::constants;
   
template<typename _realT, typename _detectorT>
class directPhaseSensor
{
public:   
   typedef _realT realT;
     
   typedef std::complex<_realT> complexT;
   
   ///The wavefront data type
   typedef wavefront<_realT> wavefrontT;
   
   ///The wavefront complex field type
   typedef mx::imagingArray<std::complex<_realT>,fftwAllocator<std::complex<_realT> >, 0> complexFieldT;
   
   ///The wavefront sensor detector image type
   //typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> wfsImageT;
   
   typedef _detectorT detectorT;
   
protected:   
   
   /* Standard WFS Interface: */
   int _wfSz; ///< Size of the wavefront in pixels

   int _detRows; ///<The number of rows of the WFS detector.  After forming the image the WFS detector plane is re-binned to this.
   int _detCols; ///<The number of columns of the WFS detector.  After forming the image the WFS detector plane is re-binned to this.

   realT _lambda; ///< Central wavelength, in meters
   
   int _iTime; ///<Integration time in loop steps
         
   int _roTime; ///<Readout time in loop steps

   realT _simStep; ///<The simulation stepsize in seconds.

   
   /* Direct Phase Sensor specific: */
   
   int _iTime_counter; 
   
   int _reading;
   
   int _roTime_counter;

   bool firstRun;

   std::vector<wavefrontT> _wavefronts;
   
   int _lastWavefront;
   
   ///The image formed by the WFS
   wfsImageT<realT> wfsImage;

   improc::ds9Interface ds9;
   
public:   
   ///Default c'tor
   directPhaseSensor();

   /* The standard WFS interface: */
   
   detectorT detector;
   
   ///The image on the detector, resized from wfsImage
   wfsImageT<realT> detectorImage;
   //detectorImageT detectorImage;
   
   ///Get the wavefront size in pixels
   /**
     * \returns the wavefront size in pixels
     */ 
   int wfSz();
   
   ///Set the wavefront size in pixels.
   void wfSz( int sz /**< [in] The new size */);
   
   ///Get the detector rows  in pixels
   /**
     * \returns _detRows
     */ 
   int detRows();
   
   ///Set the detector rows in pixels.
   /**
     * \param sz is the new size
     */ 
   void detRows(int sz);
   
   ///Get the detector Cols  in pixels
   /**
     * \returns _detCols
     */ 
   int detCols();
   
   ///Set the detector columns in pixels.
   /**
     * \param sz is the new size
     */ 
   void detCols(int sz);
   
   ///Set the detector columns in pixels.
   /**
     * \param sz is the new size
     */ 
   void detSize(int nrows, int ncols);
   
   ///Get the PyWFS central wavelength
   /**
     * \returns the central wavelength in meters
     */
   realT lambda();
   
   ///Set the PyWFS central wavelength
   /**
     * \param d is the new central wavelength in meters
     */
   void lambda(realT l);
   
   
   ///Get the PyWFS integration time, in time steps
   int iTime(); 
   
   ///Set the PyWFS integration time, in time steps
   void iTime(int it);
   
   ///Get the PyWFS detector readout time, in time steps
   int roTime(); 
   
   ///Set the PyWFS detector readout time, in time steps
   void roTime(int rt);
   
   ///Get the simulation step-size, in seconds.
   realT simStep(); 
   
   ///Set the simulation step-size, in seconds.
   void simStep(realT st);
   
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
   
   void doSenseWavefront();
   
   
   
   sigproc::psdFilter<realT> _filter;
   
   bool applyFilter;
   
   void setFilter( int width );
   
   
};

template<typename _realT, typename _detectorT>
directPhaseSensor<_realT, _detectorT>::directPhaseSensor()
{
   _wfSz = 0;
   _detRows = 0;
   _detCols = 0;
   _lambda = 0;
   
   iTime(1);
   _iTime_counter = 0;
   
   _reading = 0;
   _roTime = 1;   
   _roTime_counter = 0;
   
   _simStep = 0.001;
   
   firstRun = true;
   
   ds9.title("DPWFS");
   
   applyFilter = false;

}


template<typename _realT, typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::wfSz()
{
   return _wfSz;
}

template<typename _realT, typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::wfSz(int sz)
{
   if( _wfSz == sz) return;
   
   _wfSz = sz;
}

template<typename _realT, typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::detRows()
{
   return _detRows;
}

template<typename _realT, typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::detCols()
{
   return _detCols;
}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::detSize(int nrows, int ncols)
{
   if( _detRows == nrows && _detCols == ncols) return;
   
   _detRows = nrows;
   _detCols = ncols;
   
   detector.setSize(_detRows, _detCols);
   detectorImage.image.resize(_detRows, _detCols);
   
   
}

template<typename _realT, typename _detectorT>
_realT directPhaseSensor<_realT, _detectorT>::lambda()
{
   return _lambda;
}

template<typename _realT, typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::lambda(_realT l)
{
   _lambda = l;
}

template<typename _realT,  typename _detectorT>
template<typename AOSysT>
void directPhaseSensor<_realT, _detectorT>::linkSystem(AOSysT & AOSys)
{
 
}

template<typename _realT,  typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::iTime()
{
   return _iTime;
}
   
template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::iTime(int it)
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
int directPhaseSensor<_realT, _detectorT>::roTime()
{
   return roTime;
}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::roTime(int rt)
{
   if(rt < 1)
   {
      std::cerr << "roTime must be >= 1.  Correcting\n";
      rt = 1;
   }
   
   _roTime = rt;
   

}

template<typename _realT,  typename _detectorT>
_realT directPhaseSensor<_realT, _detectorT>::simStep()
{
   return simStep;
}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::simStep(_realT st)
{
   
   _simStep = st;
   
   detector.expTime(_simStep*_iTime);
   

}


template<typename _realT, typename _detectorT>
bool directPhaseSensor<_realT, _detectorT>::senseWavefront(wavefrontT & pupilPlane)
{
   ++_lastWavefront;
   if(_lastWavefront >= _wavefronts.size()) _lastWavefront = 0;
   
   
   wavefrontT pPlane = pupilPlane;
   
   _wavefronts[_lastWavefront].amplitude = pPlane.amplitude;
   _wavefronts[_lastWavefront].phase = pPlane.phase;
   _wavefronts[_lastWavefront].iterNo = pPlane.iterNo;
   
   //std::cerr << _lastWavefront << " " << pPlane.iterNo << "\n";
   //Always skip the first one for averaging to center of iTime.
   if(firstRun)
   {
      firstRun = false;
      return false;
   }
   
   ++_iTime_counter;

   
   bool rv = false;
   
   if(0)//_reading)
   {
      ++_roTime_counter;
      
      if(_roTime_counter >= _roTime)
      {
         detectorImage.image = wfsImage.image.block( 0.5*(wfsImage.image.rows()-1) - 0.5*(detectorImage.image.rows()-1), 0.5*(wfsImage.image.cols()-1) - 0.5*(detectorImage.image.cols()-1), detectorImage.image.rows(), detectorImage.image.cols());
         
         detectorImage.iterNo = wfsImage.iterNo;
         
         //ds9_interface_display_raw( &ds9i, 1, detectorImage.data(), detectorImage.rows(), detectorImage.cols(),1, mx::getFitsBITPIX<realT>());
         
         _roTime_counter = 0;
         _reading=0;
         rv = true;
      }
   }
   
   if( _iTime_counter >= _iTime)
   {
      //std::cerr << "DPWFS: sensing\n";
      doSenseWavefront();
      
      _iTime_counter = 0;
      
      _reading = 1;
      _roTime_counter = 0;
      
      //Just do the read
      detectorImage.image = wfsImage.image.block( 0.5*(wfsImage.image.rows()-1) - 0.5*(detectorImage.image.rows()-1), 0.5*(wfsImage.image.cols()-1) - 0.5*(detectorImage.image.cols()-1), detectorImage.image.rows(), detectorImage.image.cols());
         
      
      if(applyFilter)
      {
         _filter.filter(detectorImage.image);
      }
      
      detectorImage.iterNo = wfsImage.iterNo;
       
      //ds9_interface_display_raw( &ds9i, 1, detectorImage.image.data(), detectorImage.image.rows(), detectorImage.image.cols(),1, mx::getFitsBITPIX<realT>());
      
        _roTime_counter = 0;
      _reading=0;
      rv = true;
      
   }

   //std::cerr << "DPWFS: " << _iTime_counter << " " << _reading << " " << _roTime_counter << " " << rv << "\n";

   return rv;
   
}
   
template<typename _realT,  typename _detectorT>
bool directPhaseSensor<_realT, _detectorT>::senseWavefrontCal(wavefrontT & pupilPlane)
{
   
   _lastWavefront = 1;

   _wavefronts[0].amplitude = pupilPlane.amplitude;
   _wavefronts[0].phase = pupilPlane.phase;
   
   _wavefronts[1].amplitude = pupilPlane.amplitude;
   _wavefronts[1].phase = pupilPlane.phase;
   
   doSenseWavefront();
      
   
   detectorImage.image = wfsImage.image.block( 0.5*(wfsImage.image.rows()-1) - 0.5*(detectorImage.image.rows()-1), 0.5*(wfsImage.image.cols()-1) - 0.5*(detectorImage.image.cols()-1),             detectorImage.image.rows(), detectorImage.image.cols());

   
   return true;
   
}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::doSenseWavefront()
{
 
   wavefrontT pupilPlane;
#if 0
   BREAD_CRUMB;
   
   /* Here make average wavefront for now */
   int _firstWavefront = _lastWavefront - _iTime+1;
   if(_firstWavefront < 0) _firstWavefront += _wavefronts.size();
   
   pupilPlane.amplitude = _wavefronts[_firstWavefront].amplitude;
   pupilPlane.phase = _wavefronts[_firstWavefront].phase;
   
   pupilPlane.iterNo = _wavefronts[_firstWavefront].iterNo;
   
   //std::cerr << "DPS Averaging: " << _wavefronts[_firstWavefront].iterNo << " " ;
   BREAD_CRUMB;
   
   for(int i=0; i < _iTime-1; ++i)
   {
      ++_firstWavefront;
      if(_firstWavefront >= _wavefronts.size()) _firstWavefront = 0;
      
      pupilPlane.amplitude += _wavefronts[_firstWavefront].amplitude;
      pupilPlane.phase += _wavefronts[_firstWavefront].phase;      
      pupilPlane.iterNo += _wavefronts[_firstWavefront].iterNo;
      
      //std::cerr << _wavefronts[_firstWavefront].iterNo << " ";
   }
   //std::cerr  << " = " << _iTime << "\n";
   
   pupilPlane.amplitude /= (_iTime);
   pupilPlane.phase /= (_iTime);
   pupilPlane.iterNo /= (_iTime);
   
   /*=====================================*/

   wfsImage.image = pupilPlane.phase;
   wfsImage.iterNo  = pupilPlane.iterNo;
#else

   double _firstWavefront = _lastWavefront;// - 0.5*_iTime;
   if(_firstWavefront < 0) 
   {
      _firstWavefront += _wavefronts.size();
   }
   
   wfsImage.image = _wavefronts[(int)floor(_firstWavefront)].phase;
   wfsImage.iterNo = _wavefronts[(int)floor(_firstWavefront)].iterNo;
   
   //std::cerr << "DPS Averagxing: " << _lastWavefront << " " << _firstWavefront << " " <<  wfsImage.iterNo << " 1\n"; 
   
#endif
      
   
}
  
   
template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::setFilter( int width )
{
   int nr = detectorImage.image.rows();
   int nc = detectorImage.image.cols();
   
   typename wfsImageT<_realT>::imageT filterMask;
   
   filterMask.resize( nr, nc );
   filterMask.setZero();
   
   filterMask.block(0,0, width+1, width+1).setConstant(1.0);
   filterMask.block(0, nc - width, width+1, width).setConstant(1.0);
   filterMask.block(nr - width, 0, width, width+1).setConstant(1.0);
   filterMask.block(nr - width, nc - width, width, width).setConstant(1.0);
   
   realT v;
   
//    for(int i=0; i< 0.5*nr; ++i)
//    {
//       for(int j=0; j< 0.5*nc; ++j)
//       {
//          if( i <= width && j <=width)
//          {
//             v = 1;
//          }
//          else
//          {
//             v = 0;
//             if( i <= width )
//             {
//                v = (10 - (j-width))/10.;
//             }
//             else if(j <= width)
//             {
//                v = (10 - (i-width))/10.;
//             }
//             else
//             {
//                v = (10 - sqrt( pow(i-width,2) + pow(j-width,2)))/10.0;
//             }
//             
//             if(v < 0) v = 0;
//          }
//          
//          filterMask(i,j) = v;
//          filterMask(i, nc-1-j) = v;
//          filterMask(nr-1-i, j) = v;
//          filterMask(nr-1-i, nc-1-j) = v;
//             
//       }
//    }
   
   
   //ds9_interface_display_raw( &ds9i, 1, filterMask.data(), nr, nc,1, mx::getFitsBITPIX<realT>());
   
   //exit(0);
   _filter.psdSqrt(filterMask);
   

   
}


} //namespace sim 
   
} //namespace AO   
   
} //namespace mx

#endif //directPhaseSensor_hpp

