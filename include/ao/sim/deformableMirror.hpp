/** \file 
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief 
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef __deformableMirror_hpp__
#define __deformableMirror_hpp__

#include <cmath>

#include <mx/imagingUtils.hpp>
#include <mx/fraunhoferImager.hpp>
#include <mx/timeUtils.hpp>
#include <mx/fitsFile.hpp>
#include <mx/ds9_interface.h>

#include "../aoPaths.hpp"
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
   
using namespace boost::math::constants;

struct deformableMirrorSpec
{
   std::string name;
   std::string basisName;
};

template<typename _realT>
class deformableMirror
{
public:   

   typedef _realT realT;
     
   typedef std::complex<realT> complexT;
   
   ///The wavefront data type
   typedef wavefront<realT>  wavefrontT;
   
   ///The pupil image type
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;
 
   typedef deformableMirrorSpec specT;
   
protected:
 
   std::string _name;
   
   std::string _basisName;
   
   std::string _pupilName;
   
   ///Time to move from old shape to new shape, in loop time steps
   int _settleTime;

   int _settling;
   
   int _settleTime_counter;

   ///The amplitude used when measuring the response matrix of the current basis set.
   float _calAmp;
   
   ///The system pupil
   imageT _pupil;
      
public:   
   ds9_interface ds9i;
   
   //The modes-2-command matrix for the basis
   Eigen::Array<realT, -1, -1> _m2c;
   
   //The mirror influence functions
   mx::eigenCube<realT> _infF;
   
   
protected:
   //The current shape of the mirror
   imageT _shape;   

   //The shape of the mirror when movement begins
   imageT _oldShape;
   
   //The next shape of the mirror after movement ends.
   imageT _nextShape;

   
public:

   ///Default c'tor.
   deformableMirror();
   
   ~deformableMirror()
   {
      if(_commandFileOpen)
      {
         _commandFout.close();
      }
   }
         
   int initialize( specT & spec, 
                   const std::string & pupil); 
      
   std::string name()
   {
      return _name;
   }
   
   std::string basisName()
   {
      return _basisName;
   }
   
   
   ///Get the settling time of the DM
   int settleTime();
   
   ///Set the settling time of the DM
   void settleTime(int st);

   ///Get the calibration amplitude.
   /** The modal commands are relative to this value.
     */
   realT calAmp();
   
   ///Set the calibration amplitude.
   void calAmp(realT ca);
      
   ///Get the number of modes in the  M2C.
   int nModes();
   
   ///Apply a single mode.
   void applyMode(wavefrontT & wf, int modeNo, realT amp, realT lambda);
   
   unsigned _settlingIter;
   unsigned _settledIter;
   
   template<typename commandT>
   void setShape(commandT commandV);
      
   void applyShape(wavefrontT & wf,  realT lambda);

   double t0, t1, t_mm, t_sum;
   
   bool _writeCommands;
   bool _commandFileOpen;
   std::string _commandFile;
   std::ofstream _commandFout;
   
};


template<typename _realT>
deformableMirror<_realT>::deformableMirror()
{
   _settleTime = 1;
   _settling = 0;
   _settleTime_counter = 0;
   _calAmp = 1e-6;
   
   ds9_interface_init(&ds9i);
   ds9_interface_set_title(&ds9i, "DM");
   
   t_mm = 0;
   t_sum = 0;
   
   _writeCommands = false;
   _commandFileOpen = false;
}

// template<typename _realT>
// int deformableMirror<_realT>::initialize( const std::string & name,
//                                            const std::string & basis,
//                                            bool ortho,
//                                            const std::string & pupil )
template<typename _realT>
int deformableMirror<_realT>::initialize( specT & spec, 
                                           const std::string & pupil )
{
   _name = spec.name;
   _basisName = spec.basisName;
   _pupilName = pupil;
   
   mx::fitsFile<_realT> ff;

   std::string pName;
   pName = mx::AO::path::pupil::pupilFile(_pupilName);
   ff.read(pName, _pupil);
   
   
   std::string ifName;
   ifName = mx::AO::path::dm::influenceFunctions(_name);
   ff.read(ifName, _infF);


   
   std::string m2cName;

   m2cName = mx::AO::path::dm::M2c( _name, _basisName );
 
   ff.read(m2cName, _m2c);
   
   //std::cerr << "_infF: " << _infF.rows() << " " <<  _infF.cols() << "\n";
   _shape.resize(_infF.rows(), _infF.cols());
   
   _shape.setZero();   
   _nextShape = _shape;
   _oldShape = _shape;
   

}
   
   
   
template<typename _realT>
int deformableMirror<_realT>::settleTime()
{
   return settleTime;
}
   
template<typename _realT>
void deformableMirror<_realT>::settleTime(int st)
{
   if(st < 1)
   {
      std::cerr << "DM: settleTime must be > 1.  Correcting.\n";
      st = 1;
   }
   
   _settleTime = st;
}

template<typename _realT>
_realT deformableMirror<_realT>::calAmp()
{
   return _calAmp;
}

template<typename _realT>
void deformableMirror<_realT>::calAmp(realT ca)
{
   _calAmp = ca;
}
   
   
template<typename _realT>
int deformableMirror<_realT>::nModes()
{
   return _m2c.cols();
}
 

template<typename _realT>
void deformableMirror<_realT>::applyMode(wavefrontT & wf, int modeNo, realT amp, realT lambda)
{  
   
   Eigen::Array<_realT,-1,-1> commandV(1, nModes());
   commandV.setZero();

   commandV(0, modeNo) = 1;

   Eigen::Array<_realT, -1, -1> c;

   t0 = get_curr_time();
   c = _m2c.matrix() * commandV.matrix().transpose();
   t1 = get_curr_time();
   t_mm += t1-t0;
   
   imageT shape( _infF.rows(), _infF.cols());
   
   shape = c(0,0)*_infF.image(0);

   t0 = get_curr_time();
   #pragma omp parallel
   {
      Eigen::Array<_realT, -1, -1> tmp;
      //tmp.resize(_infF.rows(), _infF.cols());
      //tmp.setZero();
      tmp.setZero(_infF.rows(), _infF.cols());
      
      #pragma omp for schedule(static)
      for(int i=1;i < _infF.planes(); ++i)
      {
         tmp +=  c(i,0) * _infF.image(i);
      }
      #pragma omp critical
      shape += tmp;
   }
   
   t1 = get_curr_time();
   t_sum += t1-t0;
  // ds9_interface_display_raw( &ds9i, 1, shape.data(), shape.rows(), shape.cols(),1, mx::getFitsBITPIX<realT>());
   
   wf.phase += 2*amp*shape*_pupil*two_pi<realT>()/lambda;

}

template<typename _realT>
template<typename commandT>
void deformableMirror<_realT>::setShape(commandT commandV)
{
   
   if(_settling)
   {
      std::cerr << "DM: new command received while still settling.\n";
      return;
   }
   
   
   Eigen::Array<_realT, -1, -1> c;
    
   //c = -1*_calAmp*_m2c.matrix() * commandV.measurement.matrix().transpose();
   c = -1*_m2c.matrix() * commandV.measurement.matrix().transpose();
   
   if( _writeCommands )
   {
      if(!_commandFileOpen)
      {
         _commandFout.open( _commandFile );
         _commandFout << std::scientific;
         _commandFileOpen = true;
      }
      _commandFout << commandV.iterNo << "> ";
      for(int i=0; i<c.rows(); ++i)
      {
         _commandFout << commandV.measurement(0,i) << " ";
      }
      _commandFout << std::endl;
      
      _commandFout << commandV.iterNo << "> ";
      for(int i=0; i<c.rows(); ++i)
      {
         _commandFout << c(i,0) << " ";
      }
      _commandFout << std::endl;
      
//      _commandFout << c.transpose() << "\n\n";
   }
   c*=_calAmp;
   
   
   bool skipFrame = 0;
   
   ///\todo Should check for command limits here.
   for(int i=0; i< _infF.planes(); ++i)
   {
      if( std::isnan( c(i,0) ) || !std::isfinite(c(i,0))) 
      {
         skipFrame = true;
         break;
      }
   }
   
   if(skipFrame)
   {
      std::cerr << "SKIP FRAME\n";
      return;
   }
   
   //_nextShape.setZero();
   _nextShape = c(0,0)*_infF.image(0);
   #pragma omp parallel
   {
      Eigen::Array<_realT, -1, -1> tmp ;
      tmp.resize(_infF.rows(), _infF.cols());
      tmp.setZero();
      
      #pragma omp for
      for(int i=1;i < _infF.planes(); ++i)
      {
         tmp +=  c(i,0) * _infF.image(i);
      }
      #pragma omp critical
      _nextShape += tmp;
   }
      
   _oldShape = _shape;
#if 0
   _settling = 1;
   _settleTime_counter = 0;
   _settlingIter = commandV.iterNo;
#endif
//Ignoring settling time.   
   _shape = _nextShape;
   
}
   
template<typename _realT>
void deformableMirror<_realT>::applyShape(wavefrontT & wf,  realT lambda)
{
   
#if 0
   BREAD_CRUMB;
   
   if(_settling)
   {
      BREAD_CRUMB;
       _shape = _oldShape + (_nextShape-_oldShape)*( (realT) _settleTime_counter+1.0)/( (realT) _settleTime);
      
       ++_settleTime_counter;
      
       if(_settleTime_counter >= _settleTime) 
       {
          _settling = 0;
          _settledIter = _settlingIter;
       }
   }
#endif
   //imageT dshape = _shape*(2.*3.14159/900e-6)*2.;
   
   //ds9_display(4, dshape.data(), _shape.rows(), _shape.cols(), 1, mx::getFitsBITPIX<realT>());
   
   
   BREAD_CRUMB;
   
   wf.phase += 2*_shape*_pupil*two_pi<realT>()/lambda;
   
   BREAD_CRUMB;
   
   std::cerr << "DM Shape applied: " << wf.iterNo << " " << _settledIter << "\n";
}

} //sim
} //AO
} //namespace mx

#endif //__deformableMirror_hpp__

