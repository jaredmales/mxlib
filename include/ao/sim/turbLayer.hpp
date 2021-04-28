/** \file turbLayer.hpp
  * \brief Declaration and definition of a turbulence layer.
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup mxAO_sim_files
  *
  */

#ifndef turbLayer_hpp
#define turbLayer_hpp

#include <vector>

#include <Eigen/Dense>

#include "../../math/randomT.hpp"

#include "../../improc/imageTransforms.hpp"


namespace mx
{
namespace AO
{
namespace sim
{
   
///Simulation of a single turbulent layer
/** \todo document this
  * \todo add facility for changing interpolator
  */
template<typename _realT>
struct turbLayer
{
   typedef _realT realT;
   typedef Eigen::Array<realT, -1, -1> arrayT;
   
   size_t _wfSz {0};
   size_t _buffSz {0};
   size_t _scrnSz {0};
   
   realT _r0 {0};
   realT _L0 {0};
   realT _l0 {0};
   realT _pupD {0};
   realT _Cn2 {0};
   realT _z {0};
   realT _windV {0};
   realT _windD {0};
   
   realT _dx {0};
   realT _dy {0};

   int _nCombo {0};
   
   std::vector<_realT> _x0;
   std::vector<_realT> _y0;
   
   std::vector<int> _last_wdx;
   std::vector<int> _last_wdy;

   arrayT phase;
   std::vector<arrayT> shiftPhaseWP;
   arrayT shiftPhase;
   arrayT shiftPhaseWork;

   mx::math::uniDistT<realT> uniVar; ///< Uniform deviate, used in shiftRandom.
   
   //turbLayer();
   
   void setLayer( int wfSz,
                  int buffSz,
                  int scrnSz,
                  realT r0,
                  realT L0,
                  realT l0, 
                  realT pupD,
                  realT Cn2,
                  realT z,
                  realT windV, 
                  realT windD );
   
   void alloc();

   ///Shift to a timestep.
   /**
     * \param [in] dt is the new timestep.
     */ 
   void shift( realT dt );
   
   
   ///Seed the uniform deviation.  Call this if you intend to use shiftRandom.
   /** This only needs to be called once.
     */
   void initRandom();
   
   ///Shift by a random amount using the uniform distribution
   /** Call initRandom() once before calling this method.
     * 
     * \param [in] nofract if true then the fractional part is ignored, and only a whole-pixel shift is executed.  Default=false
     */ 
   void shiftRandom( bool nofract=false );
   
   
};

// template<typename realT>
// turbLayer<realT>::turbLayer()
// {
// }

template<typename realT>
void turbLayer<realT>::setLayer( int wfSz,
               int buffSz,
               int scrnSz,
               realT r0,
               realT L0,
               realT l0, 
               realT pupD,
               realT Cn2,
               realT z,
               realT windV, 
               realT windD )
{
   
   _wfSz = wfSz;
   _buffSz = buffSz;
   _scrnSz = scrnSz;
   
   _r0 = r0;
   _L0 = L0;
   _l0 = l0;
   
   _pupD = pupD;
   
   _Cn2 = Cn2;
   _z = z;
   _windV = windV;
   _windD = windD;
   
   _dx = _windV*cos(_windD)/(_pupD/_wfSz);
   _dy = _windV*sin(_windD)/(_pupD/_wfSz);
   
   

   alloc();
}

template<typename realT>
void turbLayer<realT>::alloc()
{      
   phase.resize(_scrnSz, _scrnSz);
   
   shiftPhaseWP.resize(_nCombo);
   for(int i=0;i<_nCombo; ++i)
   {
      shiftPhaseWP[i].resize( _wfSz+2*_buffSz, _wfSz+2*_buffSz);
   }
   
   shiftPhase.resize( _wfSz+2*_buffSz, _wfSz+2*_buffSz);
   
   shiftPhaseWork.resize( _wfSz+2*_buffSz, _wfSz+2*_buffSz);
   
   _x0.resize(_nCombo);
   _y0.resize(_nCombo);
   
   _last_wdx.resize(_nCombo);
   _last_wdy.resize(_nCombo);
   
   initRandom();
   
   for(int i=0;i<_nCombo; ++i)
   {
      _x0[i] = 0; //floor(uniVar * (_scrnSz));
      _y0[i] = 0; //floor(uniVar * (_scrnSz));
      
      _last_wdx[i] = _scrnSz + 1;
      _last_wdy[i] = _scrnSz + 1;
   
   }
   
   //fft.plan( _wfSz+2*_buffSz, _wfSz+2*_buffSz, MXFFT_FORWARD, true);
   //fftR.plan( _wfSz+2*_buffSz, _wfSz+2*_buffSz, MXFFT_BACKWARD, true);
}

template<typename realT>
void turbLayer<realT>::shift( realT dt )
{
   shiftPhase.setZero();
   
   for(int i=0; i < _nCombo; ++i)
   {
      int wdx, wdy;
      realT ddx, ddy;
   
      ddx = _x0[i] + _dx*dt;
      ddy = _y0[i] + _dy*dt;
   
      wdx = (int) trunc(ddx);
      ddx -= wdx;
      wdy = (int) trunc(ddy);
      ddy -= wdy;
   
      wdx %= _scrnSz;
      wdy %= _scrnSz;

      //Check for a new whole-pixel shift
      if(wdx != _last_wdx[i] || wdy != _last_wdy[i])
      {
         //Need a whole pixel shift
         improc::imageShiftWP(shiftPhaseWP[i], phase, wdx, wdy);
      }
   
      //Do the sub-pixel shift      
      improc::imageShift( shiftPhaseWork, shiftPhaseWP[i], ddx, ddy, improc::cubicConvolTransform<realT>(-0.5));
      shiftPhase += shiftPhaseWork;
      
      _last_wdx[i] = wdx;
      _last_wdy[i] = wdy;
   }
   if(_nCombo>1) shiftPhase /= sqrt(_nCombo);
}

template<typename realT>
void turbLayer<realT>::initRandom()
{
   uniVar.seed();
}

template<typename realT>
void turbLayer<realT>::shiftRandom( bool nofract )
{
   int wdx, wdy;
   realT ddx, ddy;


   ddx = uniVar * (_scrnSz);
   ddy = uniVar * (_scrnSz);
   
   wdx = (int) trunc(ddx);
   ddx -= wdx;
   wdy = (int) trunc(ddy);
   ddy -= wdy;
   
   if(nofract)
   {  
      ddx=0;
      ddy=0;
   }
   
   improc::imageShiftWP(shiftPhaseWP, phase, wdx, wdy);
   
   if(ddx !=0 || ddy != 0)
   {
      improc::imageShift( shiftPhase, shiftPhaseWP, ddx, ddy, improc::cubicConvolTransform<realT>(-0.5));
   }
   else
   {
      shiftPhase = shiftPhaseWP;
   }
}

} //namespace sim
} //namespace AO
} //namespace mx
#endif //turbLayer_hpp
