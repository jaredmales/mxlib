/** \file pyramidSensor.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declaration and definition of a standard 4 quadrant pyramid WFS.
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef __pyramidSensor_hpp__
#define __pyramidSensor_hpp__


#include "../../wfp/imagingUtils.hpp"
#include "../../wfp/fraunhoferPropagator.hpp"
#include "../../sys/timeUtils.hpp"
#include "../../ioutils/fits/fitsFile.hpp"
#include "../../improc/imageTransforms.hpp"
#include "../../math/constants.hpp"

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
   typedef wfp::imagingArray<std::complex<realT>, wfp::fftwAllocator<std::complex<realT> >, 0> complexFieldT;
   
   ///The wavefront sensor detector image type
   //typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> wfsImageT;
   
   typedef _detectorT detectorT;
   
protected:   
   
   /* Standard WFS Interface: */
   int m_wfSz; ///< Size of the wavefront in pixels

   int _detRows; ///<The number of rows of the WFS detector.  After forming the image the WFS detector plane is re-binned to this.
   int _detCols; ///<The number of columns of the WFS detector.  After forming the image the WFS detector plane is re-binned to this.

   realT m_lambda; ///< Central wavelength, in meters

public:
   
   /** \todo when the filter should be set with astrospectrum, and should re-calculate the central wavelength.
     * \todo need to verify that the wavefront propagation is appropriately chromatic
     */ 
   std::vector<realT> m_wavelengths; ///< Vector of wavelengths in the WFS bandpass
   std::vector<realT> _wavelengthWeights; ///< The relative weights of the wavelengths
   
protected:
   
   int _iTime; ///<Integration time in loop steps
         
   int _roTime; ///<Readout time in loop steps

   realT _simStep; ///<The simulation stepsize in seconds.

   
   /* PyWFS Specific: */
   
   realT m_wfPS; ///< Wavefront pixel scale, in meters/pixel
   
   realT _D; ///< Telescope diameter, in meters
   
   /** \name Modulation
     * Pyramid sensor modulation setup
     *
     * @{
     */ 
protected:
   
   int m_modSteps {20}; ///<Number of steps in the modulation simulation

   realT m_perStep {1}; ///< The minimum number of lamba/D per step to take.
   
   realT m_modRadius {3.0}; ///<Radius of the modulation in pixels
   
   bool m_wpMod {false}; ///< If true, then the whole pixel modulation path is used, with shifts instead of F.T.s
   
   realT m_effRad, m_rmsRad, m_bestRad, m_bestAng;
   
   std::vector<realT> m_modShift_x; ///< the x-coords of the continuous modulation path
   std::vector<realT> m_modShift_y; ///< the y-coords of the continuous modulation path
   
   std::vector<int> m_modShiftWP_x; ///< the x-coords of the whole-pixel modulation path
   std::vector<int> m_modShiftWP_y; ///< the y-coords of the whole-pixel modulation path
   
public:   
   
   ///Get the minimum number of modulation steps
   /**
     * \returns m_perStep;
     */ 
   realT perStep();
   
   ///Set the minimum number of modulation steps
   /**
     * \param mSt is the new number of modulation steps
     */ 
   void perStep(realT prStp /**< The minimum number of lamba/D per step to take*/ );

   ///Get the number of modulation steps
   /**
     * \returns m_modSteps, which is defined by perStep.
     */ 
   int modSteps();
   
   ///Get the radius of modulation
   /**
     * \returns m_modRadius;
     */
   realT modRadius();

   ///Set the modulation radius
   /**
     * \param mR is the new modulation radius in lambda/D
     */    
   void modRadius(realT mR /**< [in] the new value of modulation radius */ );
   
   bool wholePixelModulation();
   
   void wholePixelModulation( bool wpMod /**< [in] set whether whole pixel modulation is used (true) or not (false)*/);
   
   void wholePixelModulation( bool wpMod,   ///< [in] set whether whole pixel modulation is used (true) or not (false)
                              realT mR,     ///< [in] the new value of modulation radius
                              realT perStep ///< [in] The minimum number of lamba/D per step to take
                            );
   
   /// Find the optimum whole-pixel modulation path
   /** Finds optimum whole-pixel modulation path with the minimum
     * least-squares difference from the continous path. 
     * 
     * \returns 0 on success
     */ 
   static int optWPMod( realT & effRad,                  ///< [out] the average or effective radius after pixelation
                        realT & rmsRad,                  ///< [out] the rms error from the desired radius (in lambda/D)
                        realT & bestRad,                 ///< [out] the modulation radius which gives the best results
                        realT & bestAng,                 ///< [out] the best fit angle to rotate the modulation steps by (radians)
                        std::vector<realT> & modShift_x, ///< [out] x-coords of the continuous modulation path 
                        std::vector<realT> & modShift_y, ///< [out] y-coords of the continuous modulation path
                        std::vector<int> & modShiftWP_x, ///< [out] x-coords of the pixelated modulation path
                        std::vector<int> & modShiftWP_y, ///< [out] y-coords of the pixelated modulation path
                        realT desRad,                    ///< [in] the desired modulation radius, in lam/D units
                        realT lD,                        ///< [in] pixels per lam/D
                        realT perStep                    ///< [in] the minimum lam/D separation per modulation step
                      );
   
   
   ///@}
   
   int m_quadSz; ///<The size of the PyWFS quadrant
   
   
   wfp::fraunhoferPropagator<complexFieldT> fi;
   
   bool m_opdMaskMade {false};
   complexFieldT m_opdMask;
   
   bool m_tiltsMade {false};
   std::vector<complexFieldT> m_tilts;

   bool m_preAllocated {false};
   complexFieldT m_pupilPlaneCF;
   
   //Pre-allocated working memory:
   
   std::vector<complexFieldT> m_th_tiltedPlane; ///< Thread-local modulated wavefront
   
   complexFieldT m_focalPlane; ///< Global tip wavefront, used for whole-pixel shifting
   
   std::vector<complexFieldT> m_th_focalPlane; ///< Thread-local tip wavefront, used for FFT tilting
   
   std::vector<complexFieldT> m_th_sensorPlane; ///< Thread-local sensor-pupil-plane wavefront 
   
   std::vector<typename wfsImageT<realT>::imageT> m_th_sensorImage; ///< Thread-local sensor-pupil-plane intensity image  
   
   
   int _iTime_counter; 
   
   int _reading;
   
   int _roTime_counter;
   
   std::vector<wavefrontT> _wavefronts;
   
   int _lastWavefront;
   
public:   
   ///Default c'tor
   pyramidSensor();

   /* The standard WFS interface: */
   
   detectorT detector;
   
   ///The image on the detector, resized from m_wfsImage
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
   
   
   
   ///Get the quadrant size in pixels
   /** This is the size of the quadrant in un-binned wavefront space
     * 
     * 
     * \returns m_quadSz
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
   
//protected:
   
   ///The image formed by the WFS
   wfsImageT<realT> m_wfsImage;
   
   wfsImageT<realT> wfsTipImage;
   
   wfsImageT<realT> refTipImage;
   
   void makeOpdMask();
   
   void makeTilts();
   
   void allocThreadMem();
   
   void preAllocate();
   
public:
   void doSenseWavefront();
   void doSenseWavefront(wavefrontT &  /**< */);
   void doSenseWavefrontNoMod(wavefrontT &  /**< */ );

   bool firstRun;
   
};

template<typename realT, typename detectorT>
pyramidSensor<realT, detectorT>::pyramidSensor()
{
   m_wfSz = 0;
   _detRows = 0;
   _detCols = 0;
   m_lambda = 0;
   
   
   m_modSteps = 16;
   m_modRadius = 16;
   
   
   iTime(1);
   _iTime_counter = 0;
   
   _reading =0;
   _roTime = 1;
   _roTime_counter = 0;
   
   _simStep = 0.001;
   
   m_opdMaskMade = false;
   
   firstRun = true;
   
   ref = 1;
}

template<typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::wfSz()
{
   return m_wfSz;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::wfSz(int sz)
{
   if( m_wfSz == sz) return;
   
   m_wfSz = sz;
               
   fi.setWavefrontSizePixels(m_wfSz);

   m_tiltsMade = false;
   m_opdMaskMade = false;
   m_preAllocated = false;
}

template<typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::detRows()
{
   return _detRows;
}


template<typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::detCols()
{
   return _detCols;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::detSize(int nrows, int ncols)
{
   if( _detRows == nrows && _detCols == ncols) return;
   
   _detRows = nrows;
   _detCols = ncols;
   
   detector.setSize(_detRows, _detCols);
   detectorImage.image.resize(_detRows, _detCols);
   
   
}

template<typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::lambda()
{
   return m_lambda;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::lambda(realT l)
{
   m_lambda = l;
   
   //---------------------------------------
   //  Check if wavelength vector is filled
   //---------------------------------------
   if( m_wavelengths.size() == 0 )
   {
      m_wavelengths.resize(1, m_lambda);
      _wavelengthWeights.resize(1, 1.0);
   }
}


template<typename realT, typename detectorT>
template<typename AOSysT>
void pyramidSensor<realT, detectorT>::linkSystem(AOSysT & AOSys)
{
   AOSys.wfs.wfPS(AOSys.m_wfPS);
   AOSys.wfs.D(AOSys._D);

}

template<typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::wfPS()
{
   return m_wfPS;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::wfPS(realT ps)
{
   m_wfPS = ps;
   m_tiltsMade = false;
}

template<typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::D()
{
   return _D;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::D(realT d)
{
   _D = d;
}

template<typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::perStep()
{
   return m_perStep;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::perStep(realT prStp)
{
   m_perStep = prStp;
   
   if(m_modRadius <= 0)
   {
      m_modSteps = 0;
      return;
   }
   
   realT radPerStep = m_perStep/m_modRadius;
   
   //Get the minimum number of steps to meet m_perStep while having symmetry for the quadrants
   m_modSteps = 1;
   while( math::half_pi<realT>()/m_modSteps > radPerStep ) ++m_modSteps;
   m_modSteps *= 4;
   
   m_tiltsMade = false;
}

template<typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::modSteps()
{
   return m_modSteps;
}

template<typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::modRadius()
{
   return m_modRadius;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::modRadius(realT mR)
{
   m_modRadius = mR ;
   perStep(m_perStep); //to calculate m_modSteps;
   
   m_tiltsMade = false;
}

template<typename realT,  typename detectorT>
void pyramidSensor<realT, detectorT>::wholePixelModulation( bool wpMod )
{
   m_wpMod = wpMod;
   
   if(m_wpMod == false)
   {
      m_modShift_x.clear();
      m_modShift_y.clear();
      m_modShiftWP_x.clear();
      m_modShiftWP_y.clear();
      return;
   }
   
   m_tiltsMade = false;
   
   
   
   
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::wholePixelModulation( bool wpMod,
                                                            realT mRad,
                                                            realT prStp
                                                          )
{
   modRadius(mRad);
   perStep(prStp);
   
   wholePixelModulation(wpMod);
}
   
template<typename realT,  typename detectorT>
int pyramidSensor<realT, detectorT>::optWPMod( realT & effRad,                  
                                               realT & rmsRad,                  
                                               realT & bestRad,                 
                                               realT & bestAng,                 
                                               std::vector<realT> & modShift_x, 
                                               std::vector<realT> & modShift_y, 
                                               std::vector<int> & modShiftWP_x, 
                                               std::vector<int> & modShiftWP_y, 
                                               realT desRad,                    
                                               realT lD,                        
                                               realT perStep                    
                                             )
{
   int bestSteps;
   
   realT modRad = desRad;
   
   realT radPerStep = perStep/desRad;
   
   realT minerr = 1e9;
   
   int nT = 200;
   for(int t = 0; t < nT+1; ++t)
   {
      modRad = desRad * (1.1 - 0.2/nT*t);
      
      //Ensure we end up on the edges
      int nSteps = 1;
      while( math::half_pi<realT>()/nSteps > radPerStep ) 
      {
         ++nSteps;
      }
      nSteps *= 4;
      
      modShift_x.resize(nSteps);
      modShift_y.resize(nSteps);
      modShiftWP_x.resize(nSteps);
      modShiftWP_y.resize(nSteps);
   
      for(int a = 0; a < 100; ++a)
      {
         realT dang = math::half_pi<realT>() * a/100.0;
         
         realT avErr = 0;
         for(int n=0;n<nSteps; ++n)
         {
      
            realT ang = math::two_pi<realT>()/nSteps * n+dang;
            modShift_x[n] = modRad*lD*cos(ang);
            modShift_y[n] = modRad*lD*sin(ang);
      
            modShiftWP_x[n] = std::round(modShift_x[n]);
            modShiftWP_y[n] = std::round(modShift_y[n]);
         
            modShift_x[n] = desRad*lD*cos(ang);
            modShift_y[n] = desRad*lD*sin(ang);
         
            avErr += ( pow(modShiftWP_x[n] - modShift_x[n],2) + pow(modShiftWP_y[n]-modShift_y[n],2));
         }
      
         avErr /= nSteps;
      
         if(avErr <= minerr)
         {
            minerr = avErr;
            bestRad = modRad;
            bestAng = dang;
            bestSteps = nSteps;
         }
         
      }
   }
   
   //Now do it at best value:
   
   modRad = bestRad;
      
   int nSteps = bestSteps;
   
   modShift_x.resize(nSteps);
   modShift_y.resize(nSteps);
   modShiftWP_x.resize(nSteps);
   modShiftWP_y.resize(nSteps);

   realT avErr = 0;
   effRad = 0;
   for(int n=0;n<nSteps; ++n)
   {
      realT ang = math::two_pi<realT>()/nSteps * n+ bestAng;
      modShift_x[n] = modRad*lD*cos(ang);
      modShift_y[n] = modRad*lD*sin(ang);
   
      modShiftWP_x[n] = std::round(modShift_x[n]);
      modShiftWP_y[n] = std::round(modShift_y[n]);
      
      modShift_x[n] = desRad*lD*cos(ang);
      modShift_y[n] = desRad*lD*sin(ang);
      
      avErr += ( pow(modShiftWP_x[n] - modShift_x[n],2) + pow(modShiftWP_y[n]-modShift_y[n],2));
      effRad += sqrt( pow(modShiftWP_x[n],2) + pow(modShiftWP_y[n],2));
   }

   effRad = effRad/nSteps/lD;
   rmsRad = sqrt(avErr/nSteps)/lD;

   
   
   return 0;
   
}



template<typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::quadSz()
{
   return m_quadSz;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::quadSz(int sz)
{
   if( m_quadSz == sz) return;
      
   m_quadSz = sz;
   
   m_wfsImage.image.resize(2*m_quadSz, 2*m_quadSz);
   
   m_opdMaskMade = false;
   m_preAllocated = false;
}

template<typename realT, typename detectorT>
template<typename pupilT>
void pyramidSensor<realT, detectorT>::makeRefTipImage(pupilT & pupil)
{
   wavefrontT currWF;
   currWF.setAmplitude( pupil );
   currWF.setPhase( pupil*0 );

   refTipImage.image.resize(0,0);
   
   senseWavefrontCal(currWF);
   
   refTipImage = wfsTipImage;
   
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::makeOpdMask()
{
   complexFieldT opdMaskQ;

   m_opdMask.resize(m_wfSz, m_wfSz);
   opdMaskQ.resize(  m_wfSz, m_wfSz);
   
   opdMaskQ.set(std::complex<realT>(0,1));
   wfp::tiltWavefront(opdMaskQ, 0.5*m_quadSz, 0.5*m_quadSz);
   wfp::extractBlock(m_opdMask, 0, 0.5*m_wfSz, 0, 0.5*m_wfSz, opdMaskQ, 0 , 0);

   opdMaskQ.set(std::complex<realT>(0,1));
   wfp::tiltWavefront( opdMaskQ, -0.5*m_quadSz, -0.5*m_quadSz); 
   wfp::extractBlock(m_opdMask, 0.5*m_wfSz, 0.5*m_wfSz, 0.5*m_wfSz, 0.5*m_wfSz, opdMaskQ, 0 , 0);
   
   opdMaskQ.set(std::complex<realT>(0,1));
   wfp::tiltWavefront( opdMaskQ, 0.5*m_quadSz, -0.5*m_quadSz); 
   wfp::extractBlock(m_opdMask, 0, 0.5*m_wfSz, 0.5*m_wfSz, 0.5*m_wfSz, opdMaskQ, 0 , 0);
   
   opdMaskQ.set(std::complex<realT>(0,1));
   wfp::tiltWavefront( opdMaskQ, -0.5*m_quadSz, 0.5*m_quadSz);
   wfp::extractBlock(m_opdMask, 0.5*m_wfSz, 0.5*m_wfSz, 0, 0.5*m_wfSz, opdMaskQ, 0 , 0);
   
   m_opdMaskMade = true;

   
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::makeTilts()
{
   if(!m_wpMod)
   {
      constexpr realT pi = math::pi<realT>();
      
      realT dang = 2*pi/(m_modSteps);
      realT dx, dy;
      
      m_tilts.resize(m_modSteps);
      
      std::cout << "WF Size: " << m_wfSz << "\n";
      std::cout << "WF PS:   " << m_wfPS << "\n";
      std::cout << "Lambda:  " << m_lambda << "\n";
      std::cout << "Pyr. PS: " << wfp::fftPlateScale<realT>(m_wfSz, m_wfPS, m_lambda)*206265. << " (mas/pix)\n";
      std::cout << "Mod. steps: " << m_modSteps << "\n";
      std::cout << "Mod rad: " << m_modRadius * (m_lambda/_D)/wfp::fftPlateScale<realT>(m_wfSz, m_wfPS, m_lambda) << " (pixels)\n";
      
      for(int i=0; i < m_modSteps; ++i)
      { 
         dx = m_modRadius * (m_lambda/_D) / wfp::fftPlateScale<realT>(m_wfSz, m_wfPS, m_lambda) * cos(0.5*dang+dang * i);
         dy = m_modRadius * (m_lambda/_D) /  wfp::fftPlateScale<realT>(m_wfSz, m_wfPS, m_lambda) * sin(0.5*dang+dang * i);
      
         m_tilts[i].resize(m_wfSz, m_wfSz);
         m_tilts[i].set(std::complex<realT>(0,1));
       
         wfp::tiltWavefront(m_tilts[i], dx, dy);
      }
   }
   else
   {
      realT lD = m_wfSz / ( _D / m_wfPS );
   
      optWPMod( m_effRad,  m_rmsRad, m_bestRad, m_bestAng, m_modShift_x, m_modShift_y, m_modShiftWP_x, m_modShiftWP_y, m_modRadius, lD, m_perStep);
      m_modSteps = m_modShift_x.size();
      
      std::cout << "WF Size:       " << m_wfSz << "\n";
      std::cout << "WF PS:         " << m_wfPS << "\n";
      std::cout << "Lambda:        " << m_lambda << "\n";
      std::cout << "Mod. steps:    " << m_modSteps << "\n";
      std::cout << "Eff. Mod. Rad: " << m_effRad << " +/- " << m_rmsRad << " lambda/D\n";
      std::cout << "Angle shift:   " << m_bestAng * 180./math::pi<realT>() << " deg\n";
   }
   
   m_tiltsMade = true;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::allocThreadMem()
{
   m_pupilPlaneCF.resize(m_wfSz, m_wfSz);

   m_focalPlane.resize(m_wfSz, m_wfSz);
   
   int maxTh = omp_get_max_threads();
   m_th_tiltedPlane.resize(maxTh);
   
   m_th_focalPlane.resize(maxTh);
   
   m_th_sensorPlane.resize(maxTh);
   
   m_th_sensorImage.resize(maxTh);  
   
   for(int nTh=0;nTh<maxTh; ++nTh)
   {
      m_th_tiltedPlane[nTh].resize(m_wfSz, m_wfSz);
      
      m_th_focalPlane[nTh].resize(m_wfSz, m_wfSz);
      
      m_th_sensorPlane[nTh].resize(m_wfSz, m_wfSz);
      
      m_th_sensorImage[nTh].resize(m_quadSz*2, m_quadSz*2);
   }
   
   m_preAllocated = true;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::preAllocate()
{
   makeTilts();
   makeOpdMask();
   allocThreadMem();
}

template<typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::iTime()
{
   return _iTime;
}
   
template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::iTime(int it)
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

template<typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::roTime()
{
   return roTime;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::roTime(int rt)
{
   if(rt < 1)
   {
      std::cerr << "roTime must be >= 1.  Correcting\n";
      rt = 1;
   }
   
   _roTime = rt;
   

}

template<typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::simStep()
{
   return simStep;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::simStep(realT st)
{
   
   _simStep = st;
   
   detector.expTime(_simStep*_iTime);
   

}


template<typename realT, typename detectorT>
bool pyramidSensor<realT, detectorT>::senseWavefront(wavefrontT & pupilPlane)
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
         detector.exposeImage(detectorImage.image, m_wfsImage.image);
         
         detectorImage.tipImage = wfsTipImage.image;
         detectorImage.iterNo = m_wfsImage.iterNo;
         
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

template<typename realT, typename detectorT>
bool pyramidSensor<realT, detectorT>::senseWavefrontCal(wavefrontT & pupilPlane)
{
   
   _lastWavefront =1;

   _wavefronts[0].amplitude = pupilPlane.amplitude;
   _wavefronts[0].phase = pupilPlane.phase;
   
   _wavefronts[1].amplitude = pupilPlane.amplitude;
   _wavefronts[1].phase = pupilPlane.phase;
   

   BREAD_CRUMB;
   
   doSenseWavefront();
      
   BREAD_CRUMB;
   
   detector.exposeImage(detectorImage.image, m_wfsImage.image);

   BREAD_CRUMB;
   detectorImage.tipImage = wfsTipImage.image;
   
   BREAD_CRUMB;
   return true;
   
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::doSenseWavefront()
{ 
   

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
   doSenseWavefront(pupilPlane);
   
   m_wfsImage.iterNo = avgIt;
}

template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::doSenseWavefront(wavefrontT & pupilPlane)
{    
   if(m_modRadius == 0) return doSenseWavefrontNoMod(pupilPlane);
   
   if(!m_tiltsMade) makeTilts();
   if(!m_opdMaskMade) makeOpdMask();
   if(!m_preAllocated) preAllocate();
   
   BREAD_CRUMB;
   
   m_wfsImage.image.resize(2*m_quadSz, 2*m_quadSz);
   m_wfsImage.image.setZero();
            
   for(size_t l = 0; l<m_wavelengths.size(); ++l)
   {
      pupilPlane.lambda = m_lambda;
      pupilPlane.getWavefront(m_pupilPlaneCF, m_wavelengths[l], m_wfSz);
      
      if(m_wpMod)
      {
         //Do propagation to tip here
         fi.propagatePupilToFocal(m_focalPlane, m_pupilPlaneCF, false);
      }
      
      #pragma omp parallel 
      {
         int nTh = omp_get_thread_num();
                     
         m_th_sensorImage[nTh].setZero();
      
         complexT * tp_data;
         complexT * pp_data; 
         complexT * ti_data; 
         complexT * fp_data; 
         complexT * opd_data;
            
         pp_data = m_pupilPlaneCF.data();
         tp_data = m_th_tiltedPlane[nTh].data();  
         fp_data = m_th_focalPlane[nTh].data();
         opd_data = m_opdMask.data();
         
         int nelem = m_wfSz*m_wfSz;
         
         #pragma omp for 
         for(int i=0; i < m_modSteps; ++i)
         { 
            if(!m_wpMod)
            {
               
               ti_data = m_tilts[i].data();
               
               //---------------------------------------------
               //Apply the modulating tip 
               //---------------------------------------------
               for(int ii=0; ii< nelem; ++ii)
               {
                  tp_data[ii] = pp_data[ii]*ti_data[ii];
               }
               
               //---------------------------------------------
               //Propagate to Pyramid tip 
               //---------------------------------------------
               
               fi.propagatePupilToFocal(m_th_focalPlane[nTh], m_th_tiltedPlane[nTh], false);   
               
               //---------------------------------------------
               //Now apply the pyramid OPD 
               //---------------------------------------------
               for(int ii=0; ii< nelem; ++ii)
               {
                  fp_data[ii] = fp_data[ii]*opd_data[ii];
               }
            }
            else
            {
               //---------------------------------------------
               //Whole-pixel shift of tip image
               //---------------------------------------------
               improc::imageShiftWP(m_th_focalPlane[nTh], m_focalPlane, m_opdMask, m_modShiftWP_x[i], m_modShiftWP_y[i]);
            }
         
            
      
            //---------------------------------------------
            //Propagate to sensor plane
            //---------------------------------------------
            fi.propagateFocalToPupil(m_th_sensorPlane[nTh], m_th_focalPlane[nTh],false);
            
            //---------------------------------------------
            //Extract the image.
            //---------------------------------------------
            wfp::extractIntensityImageAccum(m_th_sensorImage[nTh], 0, 2*m_quadSz, 0, 2*m_quadSz, m_th_sensorPlane[nTh], 0.5*m_wfSz-m_quadSz, 0.5*m_wfSz-m_quadSz);
            
      
         }//for
         
         BREAD_CRUMB;
                     
         #pragma omp critical
         {
            m_wfsImage.image += m_th_sensorImage[nTh] * _wavelengthWeights[l];
           // wfsTipImage.image += pyramidImage * _wavelengthWeights[l];
            
         }
      }//#pragma omp parallel
      
   } //l for wavelength

         
   
   
   BREAD_CRUMB;
   
   m_wfsImage.image /= m_modSteps;

   
   
}


template<typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::doSenseWavefrontNoMod(wavefrontT & pupilPlane)
{
   BREAD_CRUMB;
   
   m_wfsImage.image.resize(2*m_quadSz, 2*m_quadSz);
   m_wfsImage.image.setZero();
         
   complexFieldT m_pupilPlaneCF;
   
   pupilPlane.getWavefront(m_pupilPlaneCF, m_wfSz);
   
   if(!m_opdMaskMade) makeOpdMask();
   
   complexFieldT tiltedPlane;
   complexFieldT focalPlane;
   complexFieldT sensorPlane;
     
   tiltedPlane.resize(m_wfSz, m_wfSz);
   focalPlane.resize(m_wfSz, m_wfSz);
   sensorPlane.resize(m_wfSz, m_wfSz);
      

   int nelem = m_wfSz*m_wfSz;
         
   complexT * tp_data = tiltedPlane.data();
   complexT * pp_data = m_pupilPlaneCF.data();
   complexT * opd_data = m_opdMask.data();
   complexT * fp_data = focalPlane.data();
   
   BREAD_CRUMB;
   
   //---------------------------------------------
   //Apply modulator tilt 
   //---------------------------------------------
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
   wfp::extractIntensityImageAccum(m_wfsImage.image, 0, 2*m_quadSz, 0, 2*m_quadSz, sensorPlane, 0.5*m_wfSz-m_quadSz, 0.5*m_wfSz-m_quadSz);
   
   BREAD_CRUMB;
   
}

} //namespace sim 
} //namespace AO
} //namespace mx

#endif //__pyramidSensor_hpp__

