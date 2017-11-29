/** \file 
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief 
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef simulatedAOSystem_hpp
#define simulatedAOSystem_hpp

#include <iostream>
#include <fstream>

#include "../../improc/imagePads.hpp"
#include "../../wfp/imagingUtils.hpp"
#include "../../wfp/fraunhoferPropagator.hpp"
#include "../../improc/ds9Interface.hpp"

#include "../../improc/fitsFile.hpp"
#include "../../improc/fitsUtils.hpp"

#include "../../improc/eigenImage.hpp"
#include "../../improc/eigenCube.hpp"

#include "../../timeUtils.hpp"
#include "../../sigproc/signalWindows.hpp"

#include "wavefront.hpp"
#include "../aoPaths.hpp"

#include <mx/improc/ds9Interface.hpp>

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

template<typename complexWavefrontT, typename realPupilT>
void idealCoronagraph(complexWavefrontT & wf, realPupilT & pupil)
{
   typedef typename complexWavefrontT::Scalar complexT;
   typedef typename realPupilT::Scalar        realT;
      
   Eigen::Map<Eigen::Array<complexT,-1,-1> > eigWf(wf.data(), wf.cols(), wf.rows());
   Eigen::Map<Eigen::Array<realT,-1,-1> > eigPup(pupil.data(), pupil.cols(), pupil.rows());
      
   
   eigWf = eigWf - ((eigWf*eigPup).sum()/(eigPup*eigPup).sum())*eigPup;
   
   
}

   

      
      
///A simulated AO system.
/**
  *
  * Minimum requirement for _turbSeqT:
  \code
  {
     int wfPS(realT ps); //Sets the wavefront platescale (only needed for cascaded systems).
     int frames(); //Returns the number of frames, sets the number of iterations to run.
     int nextWF(wavefront<realT> & wf); //Fills in the wavefront with phase and amplitude
  };
  \endcode 
  * 
  */
template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
class simulatedAOSystem
{
public:
   
   bool ds9off;
   
   typedef _realT realT; ///< The floating point type used for all calculations
   
   typedef mx::AO::sim::wavefront<realT> wavefrontT; ///< The wavefront type
   
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imageT; ///<The real image type
   
   typedef mx::imagingArray<std::complex<realT>, mx::fftwAllocator<std::complex<realT> >, 0> complexImageT; ///<The complex image type
   
   typedef _wfsT wfsT;
   typedef _reconT reconT;
   typedef _filterT filterT;
   typedef _dmT dmT;
   typedef _turbSeqT turbSeqT;
   
   std::string _sysName; ///< The system name for use in mx::AO::path 
   std::string _wfsName; ///< The WFS name for use in the mx::AO::path
   std::string _pupilName; ///< The pupil name for use in the mx::AO::path
   std::string _pupilMaskName; ///< The pupil mask name for use in the mx::AO::path
   
   
   wfsT wfs;
   reconT recon;
   filterT filter;
   dmT dm;
   turbSeqT turbSeq;
   
   realT _simStep;
   
   int _wfSz;
   realT _wfPS;
   
   void wfPS(realT wps)
   {
   }
   
   realT _D;
   
   realT _wfsLambda;
   realT _sciLambda;
   
   long _frameCounter;

   bool _loopClosed;
   int _loopClosedDelay;
   int _lowOrdersDelay;
   
   std::string _rmsFile;
   std::ofstream _rmsOut;
   bool _rmsUnwrap;
   
   //imageT pupilTop, pupilBot, pupilLeft, pupilRight;
   
   std::string _ampFile;
   std::ofstream _ampOut;
   
   
   /** Image outputs
     * @{ 
     */

   bool _writeIndFrames;
   
   long _saveFrameStart;
   
   std::string _psfFileBase;
   mx::improc::eigenCube<realT> _psfs;
   imageT _psfOut;
   
   bool _doCoron;
   mx::improc::eigenCube<realT> _corons;
   imageT _coronOut;
   
   int _nPerFile;
   int _currImage;
   
   int _currFile;
   
   complexImageT _complexPupil;
   complexImageT _complexFocal;
   imageT _realFocal;

   complexImageT _complexPupilCoron;
   complexImageT _complexFocalCoron;
   imageT _realPupil;
   imageT _realPhase;
   imageT _realAmp;
   imageT _realFocalCoron;
   
   wfp::fraunhoferPropagator<complexImageT> _fi;
         
   
   std::vector<typename filterT::commandT> _delayedCommands;
   int _commandDelay;
   std::vector<int> _goodCommands;
   
   /** @}
     */
   
   ///Default c'tor
   simulatedAOSystem();
   
   ///Destructor
   ~simulatedAOSystem();
   
   
   ///Initialize the system
   /**
     * Initializes the basic parts of the system.
     * 
     * \returns 0 on success
     * \returns -1 on an error (simulation should not continue if this happens).
     * 
     */ 
   int initSystem( const std::string & sysName,       ///< [in] System name.
                   typename dmT::specT & dmSpec,      ///< [in] DM Specification.
                   const std::string & wfsName,       ///< [in] WFS Name.
                   const std::string & pupilName,     ///< [in] Name of the system pupil.
                   const std::string & pupilMaskName, ///< [in] Name of the pupil mask, possibly apodized.  If empty string "", then pupilName is used.
                   const int & wfSz                   ///< [in] Size of the wavefront used for propagation.
                 );
   
   void initSim( typename reconT::specT & reconSpec,
                 realT simStep,
                 int commandDelay );

   imageT _pupil; ///< The system pupil.  This is generally a binary mask and will be applied at various points in the propagation.
   imageT _pupilMask; ///< A pupil mask which is applied once at the beginning of propagation.  Could be apodized and/or different from _pupil.
   
   imageT _postMask;
   imageT _coronPhase;
   
   
   int _npix;
   
   ///Measure the system response matrix
   /** System should be initialized with initSystemCal.
     *
     * \param amp 
     * \param rmatName
     * \param nmodes
     */ 
   void takeResponseMatrix(realT amp, std::string rmatName, int nmodes=0);

   int frames();
   
   void calcOpenLoopAmps(wavefrontT & wf);
   
   void nextWF(wavefrontT & wf);

   void runTurbulence();
   
   
};

template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::simulatedAOSystem()
{
   _frameCounter = 0;
   _loopClosed = false;      
   _loopClosedDelay = 0;
   _lowOrdersDelay = 0;
   
   _rmsUnwrap = false;

   _writeIndFrames = false;
   _saveFrameStart = 0;
   _doCoron = false;
   _nPerFile = 500;
   _currImage = 0;
   _currFile = 0;
   
   _npix = 0;
   
}
   
template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::~simulatedAOSystem()   
{
   if(_ampOut.is_open()) _ampOut.close();
   
   if(_rmsOut.is_open()) _rmsOut.close();
#if 1
   //Write out the psf and coron cubes one last time
   if(_psfFileBase !="" && _currImage > 0 && _writeIndFrames)
   {
      mx::improc::fitsFile<realT> ff;
         
      std::string fn = _psfFileBase + "_psf_" + mx::convertToString(_currFile) + ".fits";
   
      BREAD_CRUMB;
         
      ff.write(fn, _psfs);//_psfs.data(), _psfs.rows(), _psfs.cols(), _currImage);
      
      if(_doCoron)
      {
         std::string fn = _psfFileBase + "_coron_" + mx::convertToString(_currFile) + ".fits";
            
         BREAD_CRUMB;
         
         ff.write(fn, _corons);// _corons.data(), _corons.rows(), _corons.cols(), _currImage);
      }
   }
#endif

}

template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
int simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::initSystem( const std::string & sysName,
                                                                                        typename _dmT::specT & dmSpec,
                                                                                        const std::string & wfsName,
                                                                                        const std::string & pupilName,
                                                                                        const std::string & pupilMaskName,
                                                                                        const int & wfSz )
{
   _sysName = sysName;
   _wfsName = wfsName;
   _pupilName = pupilName;
   
   
   if(pupilMaskName == "")
   {
      _pupilMaskName = _pupilName;
   }
   else
   {
      _pupilMaskName = pupilMaskName;
   }
   
   
   improc::fitsFile<realT> ff;
   improc::fitsHeader head;

   //Initialize the pupil.
         
   std::string pupilFile = mx::AO::path::pupil::pupilFile(_pupilName);
         
   ff.read(_pupil, head, pupilFile);
   
   _D = head["PUPILD"].Value<realT>(); //pupilD;
   _wfPS = head["SCALE"].Value<realT>();
   
   //Pupil Mask Initialization
   pupilFile = mx::AO::path::pupil::pupilFile(_pupilMaskName);
   
   std::cerr << pupilMaskName << " " << _pupilMaskName << " " << pupilFile << "\n";
   
   ff.read(_pupilMask, pupilFile);

   if( _pupilMask.rows() != _pupil.rows() || _pupilMask.cols() != _pupil.cols())
   {
      mxError("simulatedAOSystem::initSystem", MXE_SIZEERR, "pupil and pupilMask must be the same size.");
      return -1;
   }
   
   
   //Set the wavefront size.
   _wfSz = wfSz;
   wfs.wfSz(_wfSz);
   
   

   
   
   
   
   
   //Turbulence sequence
   turbSeq._pupil = &_pupil;
   turbSeq.wfPS(_wfPS);

   //DM Initialization.
   dm.initialize( dmSpec, _pupilName);
   
   wfs.linkSystem(*this);
   
   _loopClosed = false;
   
   
   return 0;
   
}

template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
void simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::initSim( typename _reconT::specT & reconSpec,
                                                                                     realT simStep,
                                                                                     int commandDelay )
{
   
   improc::fitsFile<realT> ff;
   improc::fitsHeader head;

   
   
   _simStep = simStep;
   wfs.simStep(_simStep);

   
   filter.initialize(dm.nModes());
   
   recon.initialize(*this, reconSpec);
   dm.calAmp(recon.calAmp());
   
   _loopClosed = false;
   
   _commandDelay = commandDelay;
   
   _delayedCommands.resize((_commandDelay+1)*5);
   _goodCommands.resize((_commandDelay+1)*5);
   
}


// template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
// void simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::initSystemCal( const std::string & sysName,
//                                                                                            const std::string & dmName,
//                                                                                            const std::string & wfsName,
//                                                                                            const std::string & pupilName,
//                                                                                            const std::string & basisName,
//                                                                                            const bool & basisOrtho,
//                                                                                            const int & wfSz )
// {
//    _sysName = sysName;
//    _dmName = dmName;
//    _wfsName = wfsName;
//    _pupilName = pupilName;
//    _basisName = basisName;
//    _basisOrtho = basisOrtho;
//    
//    mx::fitsFile<realT> ff;
//    mx::fitsHeader head;
// 
//          
//    std::string pupilFile = mx::AO::path::pupil::pupilFile(_pupilName);
//          
//    ff.read(pupilFile, _pupil, head);
//    
//    //_D = head["SCALE"].Value<realT>(); //pupilD;
//    _D = head["PUPILD"].Value<realT>(); //pupilD;
//    _wfPS = head["SCALE"].Value<realT>();
//    
//    _wfSz = wfSz;
//    wfs.wfSz(_wfSz);
//    
//    _wfPS = _D/std::max(_pupil.rows(), _pupil.cols());
//    turbSeq.wfPS(_wfPS);
// 
//    //DM Initialization.
//    dm.initialize( _dmName, _basisName, _basisOrtho, _pupilName);
//    
//    //dm.loadModes(basisSet, pupil);
// 
//    wfs.linkSystem(*this);
//    
//    _loopClosed = false;
//    
//    
// }

template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
void simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::takeResponseMatrix( realT amp, 
                                                                                                std::string rmatID, 
                                                                                                int nmodes )
{
   complexImageT cpup;
   cpup.resize(_wfSz, _wfSz);
   
   recon.initializeRMat(dm.nModes(), amp, wfs.detRows(), wfs.detCols());
      
   typename filterT::commandT measureVec;
   
   wavefrontT currWF;
   currWF.setAmplitude(_pupil);
    
   wfs.iTime(1);
   wfs.roTime(1);
   wfs.detector.noNoise(true);

   double tO,tF, t0, t1, t_applyMode = 0, t_senseWF = 0, t_calcMeas = 0, t_accum = 0;

   std::cerr << dm.nModes() << "\n";
   
   tO = get_curr_time();
   
   improc::ds9Interface ds9;
      
   for(int i=0;i< dm.nModes(); ++i)
   {
      BREAD_CRUMB;
      
      currWF.setPhase(_pupil*0);
      
      realT s_amp = amp;
      if(nmodes>0 && i>= nmodes) 
      {
         s_amp =0;
         wfs.detectorImage.image.setZero();
      }
      else
      {
         t0 = get_curr_time();
         dm.applyMode(currWF, i, s_amp, 0.8e-6);
         t1 = get_curr_time();
         t_applyMode += t1-t0;
      
         BREAD_CRUMB;
         
         t0 = get_curr_time();
         wfs.senseWavefrontCal(currWF);
         t1 = get_curr_time();
         
         t_senseWF += t1-t0;
      }
      
      BREAD_CRUMB;
      
      t0 = get_curr_time();
      
      recon.calcMeasurement(measureVec, wfs.detectorImage);
      
      //ds9(wfs.detectorImage.image);
      
      t1 = get_curr_time();
 
      t_calcMeas += t1-t0;
      
      BREAD_CRUMB;

      t0 = get_curr_time();
      recon.accumulateRMat(i, measureVec, wfs.detectorImage);
      t1 = get_curr_time();
      t_accum += t1-t0;
      
      BREAD_CRUMB;
      
      //ds9_display(1, wfs.detectorImage.data(), wfs.detectorImage.cols(), wfs.detectorImage.rows(),1, mx::getFitsBITPIX<realT>());
   }
   tF = get_curr_time();
   std::cout << ( (realT) dm.nModes())/(tF - tO) << " Hz\n";
   
   std::cout << t_applyMode/dm.nModes() << " " << t_senseWF/dm.nModes() << " " << t_calcMeas/dm.nModes() << " " << t_accum/dm.nModes();
   std::cout << " " << dm.t_mm/dm.nModes() << " " << dm.t_sum/dm.nModes() << " " << "\n";
   
   std::string fname;
   fname = mx::AO::path::sys::cal::rMat(_sysName, dm.name(), _wfsName, _pupilName, dm.basisName(), rmatID, true);
   
   std::cout << fname << "\n";
   recon.saveRMat(fname);
   
   fname = mx::AO::path::sys::cal::rImages(_sysName, dm.name(), _wfsName, _pupilName, dm.basisName(), rmatID, true);
   recon.saveRImages(fname);
   
}
   
template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
int simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::frames()
{
   return turbSeq.frames();
}

/*
template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
void simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::calcOpenLoopAmps(wavefrontT & wf)
{
    BREAD_CRUMB;
   
   int npix = _pupil->sum();
   
   imageT olAmps.resize(1, dm.nModes());
   
   #pragma omp parallel for
   for(int j=0; j< dm.nModes; ++j)
   {
      realT amp;
      
      
      amp = (wf.phase*dm._infF->image(j)).sum()/ npix;
      
      olAmps(0,j) = amp;
   }

   
}
*/

template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
void simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::nextWF(wavefrontT & wf)
{
   
   BREAD_CRUMB;
      
   if(_npix == 0)
   {
      _npix = _pupil.sum();
   }
   
   BREAD_CRUMB;
   
   turbSeq.nextWF(wf);
   wf.iterNo = _frameCounter;
   
   //std::cout << "Input Photons: " << wf.amplitude.square().sum() << "\n";
      
   BREAD_CRUMB;
   
   //Mean subtraction on the system pupil.  
   realT mn = (wf.phase * _pupil).sum()/_npix;
   
   
   //Apply the pupil mask just once.
   wf.phase = (wf.phase-mn)*_pupilMask;


   BREAD_CRUMB;
  /* 
   if(_openLoopAmps)
   {
      calcOpenLoopAmps(wf);
   }*/
   
   
   BREAD_CRUMB;
   
   realT rms_ol = sqrt(  wf.phase.square().sum()/ _npix );
   
   BREAD_CRUMB;
   
   typename filterT::commandT measuredAmps, commandAmps;
   
   BREAD_CRUMB;
   
   filter.initMeasurements(measuredAmps, commandAmps);
   
   if(_frameCounter == 0)
   {
      for(int i=0; i<_delayedCommands.size();++i) 
      {
         _goodCommands[i] = 0;
      }
   }
   
   
   BREAD_CRUMB;

//   _wfsLambda = 0.78e-6;;
   
   BREAD_CRUMB;
   if(_loopClosed) dm.applyShape(wf, _wfsLambda);

   
   
   
   
   BREAD_CRUMB;
   if(_loopClosed)
   {
      bool newCV;
      
      BREAD_CRUMB;
   
      newCV = wfs.senseWavefront(wf);
         
      if(newCV)
      {         
         BREAD_CRUMB;
         
         recon.reconstruct(measuredAmps, wfs.detectorImage);
         
         BREAD_CRUMB;
         
         if(_ampFile != "")
         {
            if(!_ampOut.is_open())
            {
               _ampOut.open(_ampFile);
               _ampOut << std::scientific;
            }
            
            _ampOut << measuredAmps.iterNo << "> ";
            for(int i=0;i<measuredAmps.measurement.cols(); ++i)
            {
               _ampOut << measuredAmps.measurement(0,i) << " ";
            }
            _ampOut << std::endl;
         }
            
         
         BREAD_CRUMB;
         
  
         int nAmps = ( _frameCounter % _delayedCommands.size() );
                  
         _delayedCommands[nAmps].measurement = measuredAmps.measurement;
         _delayedCommands[nAmps].iterNo = measuredAmps.iterNo;
         
         _goodCommands[nAmps] = 1;
         
         
         //std::cerr << "nAmps: " << nAmps <<  " " << 1 << "\n";
      }
      else
      {
         int nAmps = ( _frameCounter % _delayedCommands.size() );
         _goodCommands[nAmps] = 0;
         //std::cerr << "nAmps: " << nAmps <<  " " << 0 << "\n";
      }
      
      
      int _currCommand = ( _frameCounter % _delayedCommands.size() ) - _commandDelay;
      if(_currCommand < 0) _currCommand += _delayedCommands.size();
      
      
      //std::cerr << "\t\t Current Command: " << ( _frameCounter % _delayedCommands.size() ) << " " << _currCommand << " " << _goodCommands[_currCommand] << "\n";
      if(_goodCommands[_currCommand])
      {
         BREAD_CRUMB;
         
         filter.filterCommands(commandAmps, _delayedCommands[_currCommand], _frameCounter);
            
         BREAD_CRUMB;
         
         dm.setShape(commandAmps);
      }
      
         
   }
   else
   {
      int nAmps = ( _frameCounter % _delayedCommands.size() );
      _goodCommands[nAmps] = 0;
      //std::cerr << "nAmps: " << nAmps <<  " " << 0 << "\n";
   }
   
   
   //**** If the _postMask isn't set, set it to _pupil ****//
   if(_postMask.rows() != _pupil.rows())
   {
      _postMask = _pupil;
   }
   
   
    //Mean subtraction on the system pupil.  
   mn = (wf.phase * _pupil).sum()/_npix;
   //Apply the pupil mask just once.
   wf.phase = (wf.phase-mn)*_pupil;//_pupilMask;
   
   //wf.phase *= _postMask;
   //wf.amplitude *= _postMask;
   
   //**** Calculate RMS phase ****//
   realT rms_cl;
   
   imageT uwphase = wf.phase;
   mn = uwphase.sum()/_postMask.sum();
   rms_cl = sqrt( uwphase.square().sum()/ _postMask.sum() );
   
   std::cout << _frameCounter << " WFE: " << rms_ol << " " << rms_cl << " [rad rms phase]\n";

   wfs._filter.filter(wf.phase);
   
   if(_rmsFile != "")
   {
      if(! _rmsOut.is_open() )
      {
         _rmsOut.open(_rmsFile);
         _rmsOut << "#open-loop-wfe    closed-loop-wfe  [rad rms phase]\n";
      }
      _rmsOut << rms_ol << " " << rms_cl << std::endl;
   }
      
   if(_psfFileBase != "" && _frameCounter > _saveFrameStart)
   {
      //if(_psfs.planes() == 0)
      if(_psfOut.rows() == 0 || (_psfs.planes() == 0 && _writeIndFrames))
      {
         if(_writeIndFrames) _psfs.resize(_wfSz, _wfSz, _nPerFile);
         _psfOut.resize(_wfSz, _wfSz);
         _psfOut.setZero();
         
         _complexPupil.resize( _wfSz, _wfSz);
         _complexFocal.resize( _wfSz, _wfSz);
         _realFocal.resize(_wfSz, _wfSz);
         
         if(_doCoron) 
         {
            if(_writeIndFrames) _corons.resize(_wfSz, _wfSz, _nPerFile);
            _coronOut.resize(_wfSz, _wfSz);
            _coronOut.setZero();
            
            _complexPupilCoron.resize( _wfSz, _wfSz);
            _complexFocalCoron.resize( _wfSz, _wfSz);
            _realFocalCoron.resize(_wfSz, _wfSz);
            
            _realPhase.resize(_wfSz, _wfSz);
            _realAmp.resize(_wfSz, _wfSz);
            
            
            //Create Coronagraph pupil.
            //improc::padImage(_realPupil, _postMask, 0.5*(_wfSz-_pupil.rows()),0);
            improc::padImage(_realPupil, _pupil, 0.5*(_wfSz-_pupil.rows()),0);
            
//             improc::ds9Interface ds9;
//             ds9(_pupil,1);
//             ds9(_realPupil,2);
         }            
      }

      
      if(_doCoron)
      {
         
         wavefrontT cwf;
         cwf.phase = wf.phase;
         cwf.amplitude = wf.amplitude;
         
         //improc::ds9Interface ds9;
         //ds9(cwf.phase,1);
         //ds9(cwf.amplitude,2);
         
         bool ideal = true;
         if(_coronPhase.rows() > 0)
         {
            cwf.phase += _coronPhase;
            ideal = false;
         }
         
         cwf.getWavefront(_complexPupil, _wfSz);
         
         //Haven't implemented assignment in imagingArray
         Eigen::Map<Eigen::Array<std::complex<realT>,-1,-1> > mt(_complexPupilCoron.data(), _complexPupilCoron.rows(), _complexPupilCoron.cols());
         mt = Eigen::Map<Eigen::Array<std::complex<realT>,-1,-1> >(_complexPupil.data(), _complexPupil.rows(), _complexPupil.cols());

         
         if(ideal) idealCoronagraph(_complexPupilCoron, _realPupil);
         
         _fi.propagatePupilToFocal(_complexFocalCoron, _complexPupilCoron);
      
         wfp::extractIntensityImage(_realFocalCoron,0,_complexFocalCoron.rows(),0,_complexFocalCoron.cols(),_complexFocalCoron,0,0);
         
         if(_writeIndFrames) _corons.image(_currImage) = _realFocalCoron;
         
         _coronOut += _realFocalCoron;
      }

      wf.getWavefront(_complexPupil, _wfSz);
#if 0
      imageT nm;
      nm.resize(_complexPupil.rows(), _complexPupil.cols());
      for(int i=0;i<nm.rows();++i)
      {
         for(int j=0; j<nm.cols();++j)
         {
            nm(i,j) = real(_complexPupil(i,j));
         }
      }
      improc::ds9Interface ds9;
      ds9(nm);
#endif

//       for(int ni=0; ni< _complexPupilCoron.rows(); ++ni)
//       {
//          for(int nj=0; nj< _complexPupilCoron.cols(); ++nj)
//          {
//             _complexPupilCoron(ni,nj) = _complexPupilCoron(ni,nj)*_coronMask(ni,nj);
//          }
//       }
            
      _fi.propagatePupilToFocal(_complexFocal, _complexPupil);
      wfp::extractIntensityImage(_realFocal,0,_complexFocal.rows(),0,_complexFocal.cols(),_complexFocal,0,0);
      if(_writeIndFrames) _psfs.image(_currImage) = _realFocal;
      _psfOut += _realFocal;
      
      #if 0
      ds9_display((ds9off*3) + 4, _realFocal.data(), _realFocal.rows(),_realFocal.cols(),1, mx::getFitsBITPIX<realT>());
      #endif
      
      ++_currImage;
#if 1
      if(_currImage >= _nPerFile && _writeIndFrames)
      {
         mx::improc::fitsFile<realT> ff;
         
         std::string fn = _psfFileBase + "_psf_" + mx::convertToString(_currFile) + ".fits";
   
         BREAD_CRUMB;
         
         ff.write(fn, _psfs);// _psfs.data(), _psfs.rows(), _psfs.cols(), _psfs.planes());
      
         if(_doCoron)
         {
            std::string fn = _psfFileBase + "_coron_" + mx::convertToString(_currFile) + ".fits";
            
            BREAD_CRUMB;
            
            ff.write(fn, _corons);//_corons.data(), _corons.rows(), _corons.cols(), _corons.planes());
         }

         ++_currFile;
         _currImage = 0;
      }
#endif    
   }//if(_psfFileBase != "" && _frameCounter > _saveFrameStart)
   
   ++_frameCounter;
   
}//void simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::nextWF(wavefrontT & wf)

template<typename _realT, typename _wfsT, typename _reconT, typename _filterT, typename _dmT, typename _turbSeqT>
void simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::runTurbulence()
{   
   wavefrontT currWF;


   _wfsLambda = wfs.lambda();
   
   for(int i=0;i<turbSeq.frames();++i)
   { 
      //std::cout << i << "/" << turbSeq.frames() << "\n" ;
   
      
      if(i == 0) 
      {
         turbSeq._loopClosed = true;
         //wfs.applyFilter = false;
      }
      
      if(i == _loopClosedDelay) _loopClosed = true;
      
      
   
      BREAD_CRUMB;
      
      nextWF(currWF);
      
      BREAD_CRUMB;
      
   }
   
   if(_psfFileBase != "")
   {
      improc::fitsFile<realT> ff;
      std::string fn = _psfFileBase + "_psf.fits";
      ff.write(fn, _psfOut);
      
      if(_doCoron)
      {
         fn = _psfFileBase + "_coron.fits";
         ff.write(fn, _coronOut);
      }
   }
   
      
}//void simulatedAOSystem<_realT, _wfsT, _reconT, _filterT, _dmT, _turbSeqT>::runTurbulence()
   
} //namespace sim 
} //namespace AO
} //namespace mx

#endif //simulatedAOSystem_hpp
