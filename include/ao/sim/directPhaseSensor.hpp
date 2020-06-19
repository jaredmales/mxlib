/** \file directPhaseSensor.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declaration and definition of a direct phase sensor.
  * \ingroup mxAO_sim_files
  *
  */

#ifndef directPhaseSensor_hpp
#define directPhaseSensor_hpp

#include "../../randomT.hpp"

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

   realT iterNo;

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

   typedef _detectorT detectorT;

   typedef mx::randomT<realT, std::mt19937_64, std::normal_distribution<realT> > norm_distT;

protected:

   /* Standard WFS Interface: */
   int m_wfSz {0}; ///< Size of the wavefront in pixels

   int m_detRows {0}; ///<The number of rows of the WFS m_detector.  After forming the image the WFS detector plane is re-binned to this.
   int m_detCols {0}; ///<The number of columns of the WFS m_detector.  After forming the image the WFS detector plane is re-binned to this.

   realT m_lambda {0}; ///< Central wavelength, in meters

   int m_iTime {0}; ///<Integration time in loop steps

   int m_roTime {1}; ///<Readout time in loop steps

   realT m_simStep {0.001}; ///<The simulation stepsize in seconds.


   /* Direct Phase Sensor specific: */

   int m_iTime_counter {0};

   int m_reading {0};

   int m_roTime_counter {0};

   bool m_firstRun {true};

   std::vector<wavefrontT> m_wavefronts;

   int m_lastWavefront;

   ///The image formed by the WFS
   wfsImageT<realT> m_wfsImage;


   
   improc::ds9Interface m_ds9;
   improc::ds9Interface m_ds9f;

public:
   ///Default c'tor
   directPhaseSensor();

   /* The standard WFS interface: */

   detectorT m_detector;

   ///The image on the detector, resized from wfsImage
   wfsImageT<realT> m_detectorImage;
  
   ///Get the wavefront size in pixels
   /**
     * \returns the wavefront size in pixels
     */
   int wfSz();

   ///Set the wavefront size in pixels.
   void wfSz( int sz /**< [in] The new size */);

   ///Get the detector rows  in pixels
   /**
     * \returns m_detRows
     */
   int detRows();

   ///Set the detector rows in pixels.
   /**
     * \param sz is the new size
     */
   void detRows(int sz);

   ///Get the detector Cols  in pixels
   /**
     * \returns m_detCols
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

   /** \name Spatial Filtering
     * @{
     */
protected:
   sigproc::psdFilter<realT> m_filter; ///< The spatial filter class object.

   bool m_applyFilter {false}; ///< Flag to control whether or not the spatial filter is applied.
   
   int m_filterWidth {0}; ///< The half-width of the filter, in Fourier domain pixels corresponding to 1/D spatial frequency units.

public:
   
   /// Set the flag controlling whether or not the spatial filter is applied
   void applyFilter(bool af /**< [in] true to apply the filter, false to not apply*/); 
   
   /// Get the current value of the flag controlling whether or not the spatial filter is applied
   /**
     * \returns the current value of m_applyFilter
     */ 
   bool applyFilter();
   
   /// Set the spatial filter half-width
   /** This sets up the filter and loads it into m_filter.  The filter is a hard-edged square of width 2*width+1 pixels in the Fourier plane.
     * 
     * m_filterWidth must still be set to true after calling this function.
     */ 
   void filterWidth( int width /**< [in] the desired half-width of the filter*/);
   
   /// Get the current value of the filter half-width
   /**
     * \returns the current value of m_filterWidth
     */
   int filterWidth();
   
   ///@}
   
   /** \name Photon Noise
     *
     * @{
     */ 
protected:
   
   norm_distT m_normVar; ///< Gets normal-distributed variates
   
   /// The photon noise senstivity parameter.  
   /** If 0 (default) no noise is applied.  A value of 1 represents the ideal interferomter.
     * Note that this is constant vs. spatial frequency.  See Guyon (2005) \cite guyon_2005 for more information about \f$ \beta_p \f$.
     */ 
   realT m_beta_p {0};
   
public:
   
   /// Set the new value of the photon noise sensitivity parameter, beta_p.
   void beta_p(realT bp /**< [in] the new value of beta_p */);
   
   /// Get the current value of the photon noise sensitivity parameter.
   /**
     * \returns the current value of m_beta_p.
     */  
   realT beta_p();
   
   ///@}
};

template<typename _realT, typename _detectorT>
directPhaseSensor<_realT, _detectorT>::directPhaseSensor()
{
   iTime(1);



   m_ds9.title("DPWFS");
   m_ds9f.title("DPWFS_Filtered");


   m_normVar.seed();
}


template<typename _realT, typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::wfSz()
{
   return m_wfSz;
}

template<typename _realT, typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::wfSz(int sz)
{
   if( m_wfSz == sz) return;

   m_wfSz = sz;
}

template<typename _realT, typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::detRows()
{
   return m_detRows;
}

template<typename _realT, typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::detCols()
{
   return m_detCols;
}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::detSize(int nrows, int ncols)
{
   if( m_detRows == nrows && m_detCols == ncols) return;

   m_detRows = nrows;
   m_detCols = ncols;

   m_detector.setSize(m_detRows, m_detCols);
   m_detectorImage.image.resize(m_detRows, m_detCols);


}

template<typename _realT, typename _detectorT>
_realT directPhaseSensor<_realT, _detectorT>::lambda()
{
   return m_lambda;
}

template<typename _realT, typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::lambda(_realT l)
{
   m_lambda = l;
}

template<typename _realT,  typename _detectorT>
template<typename AOSysT>
void directPhaseSensor<_realT, _detectorT>::linkSystem(AOSysT & AOSys)
{
   static_cast<void>(AOSys);
}

template<typename _realT,  typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::iTime()
{
   return m_iTime;
}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::iTime(int it)
{
   if(it < 1)
   {
      std::cerr << "iTime must be >= 1.  Correcting\n";
      it = 1;
   }

   m_iTime = it;

   m_wavefronts.resize(m_iTime+2+100);
   m_lastWavefront = -1;

   m_detector.expTime(m_simStep*m_iTime);

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

   m_roTime = rt;


}

template<typename _realT,  typename _detectorT>
_realT directPhaseSensor<_realT, _detectorT>::simStep()
{
   return simStep;
}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::simStep(_realT st)
{

   m_simStep = st;

   m_detector.expTime(m_simStep*m_iTime);


}


template<typename _realT, typename _detectorT>
bool directPhaseSensor<_realT, _detectorT>::senseWavefront(wavefrontT & pupilPlane)
{
   ++m_lastWavefront;
   if((size_t) m_lastWavefront >= m_wavefronts.size()) m_lastWavefront = 0;


   wavefrontT pPlane = pupilPlane;

   m_wavefronts[m_lastWavefront].amplitude = pPlane.amplitude;
   m_wavefronts[m_lastWavefront].phase = pPlane.phase;
   m_wavefronts[m_lastWavefront].iterNo = pPlane.iterNo;

   //std::cerr << m_lastWavefront << " " << pPlane.iterNo << "\n";
   //Always skip the first one for averaging to center of iTime.
   if(m_firstRun)
   {
      m_firstRun = false;
      return false;
   }

   ++m_iTime_counter;


   bool rv = false;

   if(0)//m_reading)
   {
      ++m_roTime_counter;

      if(m_roTime_counter >= m_roTime)
      {
         m_detectorImage.image = m_wfsImage.image.block( 0.5*(m_wfsImage.image.rows()-1) - 0.5*(m_detectorImage.image.rows()-1), 0.5*(m_wfsImage.image.cols()-1) - 0.5*(m_detectorImage.image.cols()-1), m_detectorImage.image.rows(), m_detectorImage.image.cols());

         m_detectorImage.iterNo = m_wfsImage.iterNo;

         //ds9_interface_display_raw( &ds9i, 1, m_detectorImage.data(), m_detectorImage.rows(), m_detectorImage.cols(),1, mx::getFitsBITPIX<realT>());

         m_roTime_counter = 0;
         m_reading=0;
         rv = true;
      }
   }

   if( m_iTime_counter >= m_iTime)
   {
      //std::cerr << "DPWFS: sensing\n";
      doSenseWavefront();

      m_iTime_counter = 0;

      m_reading = 1;
      m_roTime_counter = 0;

      //Just do the read
      m_detectorImage.image = m_wfsImage.image.block( 0.5*(m_wfsImage.image.rows()-1) - 0.5*(m_detectorImage.image.rows()-1), 0.5*(m_wfsImage.image.cols()-1) - 0.5*(m_detectorImage.image.cols()-1), m_detectorImage.image.rows(), m_detectorImage.image.cols());

      //*** Adding Noise:
      if(m_beta_p > 0)
      {
         //1) Subtract the min
         realT phaseMin = m_detectorImage.image.minCoeff();
         typename wfsImageT<realT>::imageT noiseIm = m_detectorImage.image - phaseMin;
         
         //2) Set to 0 outside the pupil
         for(int r=0;r<noiseIm.rows();++r)
         {
            for(int c=0;c<noiseIm.cols();++c)
            {
               if(pupilPlane.amplitude(r,c) == 0) noiseIm(r,c)=0;
            }
         }
         
         //3) Normalize to total photons for an interferometer
         realT phaseSum = noiseIm.sum();
         realT totalPhots = 2*pupilPlane.amplitude.square().sum()*m_detector.expTime();
         
         
         //4) Add Poisson noise
         
         for(int r=0;r<noiseIm.rows();++r)
         {
            for(int c=0;c<noiseIm.cols();++c)
            {
               if(noiseIm(r,c) == 0)
               {
                  continue;
               }
               else
               {
                  m_detectorImage.image(r,c) += m_beta_p*sqrt(noiseIm(r,c)*phaseSum/totalPhots)*m_normVar;
               }
            }
         }
         
      } //if(m_beta_p > 0)

      //ds9(m_detectorImage.image);

      if(m_applyFilter)
      {
         m_filter.filter(m_detectorImage.image);
      }

      m_detectorImage.iterNo = m_wfsImage.iterNo;

      m_roTime_counter = 0;
      m_reading=0;
      rv = true;

   }

   //std::cerr << "DPWFS: " << m_iTime_counter << " " << m_reading << " " << m_roTime_counter << " " << rv << "\n";

   return rv;

}

template<typename _realT,  typename _detectorT>
bool directPhaseSensor<_realT, _detectorT>::senseWavefrontCal(wavefrontT & pupilPlane)
{

   m_lastWavefront = 1;

   m_wavefronts[0].amplitude = pupilPlane.amplitude;
   m_wavefronts[0].phase = pupilPlane.phase;

   m_wavefronts[1].amplitude = pupilPlane.amplitude;
   m_wavefronts[1].phase = pupilPlane.phase;

   doSenseWavefront();

   m_detectorImage.image = m_wfsImage.image.block( 0.5*(m_wfsImage.image.rows()-1) - 0.5*(m_detectorImage.image.rows()-1), 0.5*(m_wfsImage.image.cols()-1) - 0.5*(m_detectorImage.image.cols()-1), m_detectorImage.image.rows(), m_detectorImage.image.cols());

   //ds9(m_detectorImage.image);

   if(m_applyFilter)
   {
      m_filter.filter(m_detectorImage.image);

      //ds9f(m_detectorImage.image);
   }



   return true;

}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::doSenseWavefront()
{
   wavefrontT pupilPlane;

   BREAD_CRUMB;
   
   /* Here make average wavefront for now */
   int _firstWavefront = m_lastWavefront - m_iTime;
   if(_firstWavefront < 0) _firstWavefront += m_wavefronts.size();

   //std::cerr << m_lastWavefront << " " << m_iTime_counter << " " << _firstWavefront << "\n";
   
   pupilPlane.amplitude = 0.5*m_wavefronts[_firstWavefront].amplitude;
   pupilPlane.phase = 0.5*m_wavefronts[_firstWavefront].phase;

   pupilPlane.iterNo = 0.5*m_wavefronts[_firstWavefront].iterNo;


   //std::cerr << "DPS Averaging: " << m_wavefronts[_firstWavefront].iterNo << " " ;
   BREAD_CRUMB;

   if(m_wavefronts[_firstWavefront].iterNo < m_iTime)
   {
      m_wfsImage.image = pupilPlane.phase;
      m_wfsImage.iterNo  = pupilPlane.iterNo;
      return;
   }

   for(int i=0; i < m_iTime - 1; ++i)
   {
      ++_firstWavefront;
      if( (size_t) _firstWavefront >= m_wavefronts.size()) _firstWavefront = 0;

      //std::cerr << m_lastWavefront << " " << m_iTime_counter << " " << _firstWavefront << "\n";

      pupilPlane.amplitude += m_wavefronts[_firstWavefront].amplitude;
      pupilPlane.phase += m_wavefronts[_firstWavefront].phase;
      pupilPlane.iterNo += m_wavefronts[_firstWavefront].iterNo;

   }

   ++_firstWavefront;
   if( (size_t) _firstWavefront >= m_wavefronts.size()) _firstWavefront = 0;

   //std::cerr << m_lastWavefront << " " << m_iTime_counter << " " << _firstWavefront << "\n";

   pupilPlane.amplitude += 0.5*m_wavefronts[_firstWavefront].amplitude;
   pupilPlane.phase += 0.5*m_wavefronts[_firstWavefront].phase;
   pupilPlane.iterNo += 0.5*m_wavefronts[_firstWavefront].iterNo;
      
   pupilPlane.amplitude /= (m_iTime);
   pupilPlane.phase /= (m_iTime);
   pupilPlane.iterNo /= (m_iTime);

   m_wfsImage.image = pupilPlane.phase;
   m_wfsImage.iterNo  = pupilPlane.iterNo;

}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::applyFilter(bool af)
{
   m_applyFilter = af;
}

template<typename _realT,  typename _detectorT>
bool directPhaseSensor<_realT, _detectorT>::applyFilter()
{
   return m_applyFilter;
}
   

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::filterWidth( int width )
{
   int nr = m_detectorImage.image.rows();
   int nc = m_detectorImage.image.cols();

   typename wfsImageT<_realT>::imageT filterMask;

   filterMask.resize( nr, nc );
   filterMask.setZero();

   filterMask.block(0,0, width+1, width+1).setConstant(1.0);
   filterMask.block(0, nc - width, width+1, width).setConstant(1.0);
   filterMask.block(nr - width, 0, width, width+1).setConstant(1.0);
   filterMask.block(nr - width, nc - width, width, width).setConstant(1.0);

   m_filter.psdSqrt(filterMask);

   m_filterWidth = width;

}

template<typename _realT,  typename _detectorT>
int directPhaseSensor<_realT, _detectorT>::filterWidth()
{
   return m_filterWidth;
}

template<typename _realT,  typename _detectorT>
void directPhaseSensor<_realT, _detectorT>::beta_p(realT bp)
{
   m_beta_p = bp;
}

template<typename realT,  typename detectorT>
realT directPhaseSensor<realT, detectorT>::beta_p()
{
   return m_beta_p;
}

} //namespace sim

} //namespace AO

} //namespace mx

#endif //directPhaseSensor_hpp
