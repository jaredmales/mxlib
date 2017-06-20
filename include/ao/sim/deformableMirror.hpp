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
   
   ///The system pupil, possibly apodized, etc.
   imageT _pupil;
   
      
public:   
   ds9_interface ds9i_shape, ds9i_phase, ds9i_acts;
   
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
   
   int display_phase;
   int display_phase_counter;
   int display_shape;
   int display_shape_counter;
   
   int display_acts;
   int display_acts_counter;
   
   realT _commandLimit;
   
   Eigen::Array<double,-1,-1> _pos, _map;
   
};


template<typename _realT>
deformableMirror<_realT>::deformableMirror()
{
   _settleTime = 1;
   _settling = 0;
   _settleTime_counter = 0;
   _calAmp = 1e-6;
   
   ds9_interface_init(&ds9i_shape);
   ds9_interface_set_title(&ds9i_shape, "DM_Shape");
   ds9_interface_init(&ds9i_phase);
   ds9_interface_set_title(&ds9i_phase, "DM_Phase");
   
   ds9_interface_init(&ds9i_acts);
   ds9_interface_set_title(&ds9i_acts, "DM_Actuators");
   
   t_mm = 0;
   t_sum = 0;
   
   _writeCommands = false;
   _commandFileOpen = false;
   
   display_phase = 0;
   display_phase_counter = 0;
   display_shape = 0;
   display_shape_counter = 0;
   
   display_acts = 0;
   display_acts_counter = 0;
   
   _commandLimit = 0;
}


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
   
   
   if(_name == "modalDM")
   {
      std::cerr << "modalDM\n";
      
      std::string ifName;
      ifName = mx::AO::path::basis::modes(_basisName);
      ff.read(ifName, _infF);
      
      _m2c.resize( _infF.planes(), _infF.planes());
      _m2c.setZero();
      
      for(int i=0;i<_m2c.rows();++i) _m2c(i,i) = 1.0;
      
   }
   else
   {
      std::string ifName;
      ifName = mx::AO::path::dm::influenceFunctions(_name);
      
      
      mx::eigenCube<realT> infFLoad;
      ff.read(ifName, infFLoad);
   
      realT c = 0.5*(infFLoad.rows()-1);
      realT w = 0.5*(_pupil.rows()-1);
      
      _infF.resize( _pupil.rows(), _pupil.cols(), infFLoad.planes());
      
      for(int i=0;i<infFLoad.planes(); ++i)
      {
         _infF.image(i) = infFLoad.image(i).block( c-w, c-w, _pupil.rows(), _pupil.rows());
      }
      
      std::string m2cName;

      m2cName = mx::AO::path::dm::M2c( _name, _basisName );
 
      ff.read(m2cName, _m2c);
      
      
      std::string posName =  mx::AO::path::dm::actuatorPositions(_name, true);
   
      ff.read(posName, _pos);

   }
   
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

template<typename realT>
void makeMap( Eigen::Array<realT, -1, -1> & map,  Eigen::Array<realT, -1, -1> & pos, Eigen::Array<realT, -1, -1> & act)
{
   
   realT minx = pos.row(0).minCoeff();
   realT maxx = pos.row(0).maxCoeff();
   
   int i=0;
   realT dx = 0;
   while(dx == 0)
   {
      dx = fabs(pos(0,i)- pos(0,0));
      ++i;
   }
   
   realT miny = pos.row(1).minCoeff();
   realT maxy = pos.row(1).maxCoeff();
   
   i = 0;
   realT dy = 0;
   
   while(dy == 0)
   {
      dy = fabs(pos(1,i)- pos(1,0));
      ++i;
   }
   
   int nx = (maxx-minx)/dx + 1;
   int ny = (maxy-miny)/dy + 1;
   
   
   map.resize(nx, ny);
   map.setZero();
   
   realT x, y;
   
   for(int j=0;j < pos.cols(); ++j)
   {
      x = (pos(0,j) - minx)/dx;
      
      y = ( pos(1,j) - miny ) / dy;
      
      map(x,y) = act(j,0);
   }
   
   
   
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
   c = -1*_calAmp*_m2c.matrix() * commandV.measurement.matrix().transpose();
   
   if(_commandLimit > 0)
   {
      realT r1 = sqrt( pow(_pos(0,1) - _pos(0,0),2) + pow(_pos(1,1) - _pos(1,0),2));
   

      realT r;
   
      int nLimited = 0;
      
      for(int i=0; i< _pos.cols(); ++i)
      {
         for(int j=i+1; j< _pos.cols(); ++j)
         {
            r = sqrt( pow(_pos(0,j) - _pos(0,i),2) + pow(_pos(1,j) - _pos(1,i),2));
         
            if( fabs(r1 - r) < .01 )
            {
               realT iact = fabs( c(i,0) - c(j,0) ); 
               if( iact > _commandLimit )
               {
                  std::cerr << "Limited Actuators " << i << " " << j << "\n";
                  ++nLimited;
                  c(i,0) *= _commandLimit/iact;
                  c(j,0) *= _commandLimit/iact;
               }
            }
         }
      }
      if(nLimited > 0) std::cerr << nLimited << " strokes limited\n";
      
   }
   
   
   
   
   if( display_acts > 0)
   {
      ++display_acts_counter;
      
      if(display_acts_counter >= display_acts)
      {
         makeMap( _map, _pos, c);

         ds9_interface_display_raw( &ds9i_acts, 1, _map.data(), _map.rows(), _map.cols(),1, mx::getFitsBITPIX<realT>());
         display_acts_counter = 0;
      }
   }

   
#if 0   
   if(_commandLimit > 0 )
   {
      for(int i=0; i < c.rows(); ++i)
      {
         if(c(i,0) > _commandLimit ) c(i,0) = _commandLimit;
         if(c(i,0) < -1*_commandLimit ) c(i,0) = -1*_commandLimit;
      }
   }
#endif


   
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
         _commandFout << c(i,0) << " ";
      }
      _commandFout << std::endl;
   }
   
   
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
    
    
    _nextShape *= _pupil;
   
    _nextShape = (_nextShape - _nextShape.sum()/_pupil.sum())*_pupil;

   
   _oldShape = _shape;
#if 1
   _settling = 1;
   _settleTime_counter = 0;
   _settlingIter = commandV.iterNo;
#else
//Ignoring settling time.   
   _shape = _nextShape;
#endif
      
}
   
template<typename _realT>
void deformableMirror<_realT>::applyShape(wavefrontT & wf,  realT lambda)
{
   
#if 1
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

   if( display_shape > 0)
   {
      ++display_shape_counter;
      
      if(display_shape_counter >= display_shape)
      {
         ds9_interface_display_raw( &ds9i_shape, 1, _shape.data(), _shape.rows(), _shape.cols(),1, mx::getFitsBITPIX<realT>());
         display_shape_counter = 0;
      }
   }
   
   BREAD_CRUMB;
   
   wf.phase += 2*_shape*_pupil*two_pi<realT>()/lambda;
   
   if( display_phase > 0)
   {
      ++display_phase_counter;
      
      if(display_phase_counter >= display_phase)
      {
         ds9_interface_display_raw( &ds9i_phase, 1, wf.phase.data(), wf.phase.rows(), wf.phase.cols(),1, mx::getFitsBITPIX<realT>());

         display_phase_counter = 0;
      }
   }
   
   
   BREAD_CRUMB;
   
   //std::cerr << "DM " << _infF.planes() << " Shape applied: " << wf.iterNo << " " << _settledIter << "\n";
}

} //sim
} //AO
} //namespace mx

#endif //__deformableMirror_hpp__

