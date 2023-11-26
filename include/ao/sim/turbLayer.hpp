/** \file turbLayer.hpp
  * \brief Declaration and definition of a turbulence layer.
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup mxAO_sim_files
  *
  */

//***********************************************************************//
// Copyright 2023 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef mx_AO_sim_turbLayer_hpp
#define mx_AO_sim_turbLayer_hpp

#include <vector>

#include <Eigen/Dense>

#include "../../improc/imageTransforms.hpp"

namespace mx
{
namespace AO
{
namespace sim
{
   
template<typename _aoSystemT>
struct turbAtmosphere;

///Simulation of a single turbulent layer
/** \todo document this
  * \todo add facility for changing interpolator
  * 
  * \ingroup mxAOSim
  */
template<typename _aoSystemT>
struct turbLayer
{
    typedef _aoSystemT aoSystemT;
    typedef typename aoSystemT::realT realT;
    typedef Eigen::Array<realT, -1, -1> imageT;
   
    turbAtmosphere<aoSystemT> * m_parent {nullptr};
   
    uint32_t m_layerNo;

    uint32_t m_scrnSz;

    imageT m_freq;
    imageT m_psd;

    realT m_dx {0};
    realT m_dy {0};

    realT m_x0;
    realT m_y0;
   
    int m_last_wdx;
    int m_last_wdy;

    imageT m_phase;
    imageT m_shiftPhaseWP;
    imageT m_shiftPhase;
    imageT m_shiftPhaseWork;

    mx::math::uniDistT<realT> uniVar; ///< Uniform deviate, used in shiftRandom.
   
    void setLayer( turbAtmosphere<aoSystemT> * parent,
                   int layerNo,
                   uint32_t scrnSz
                 );
    
    uint32_t layerNo();

    uint32_t scrnSz();

    void alloc();

    /// Deallocate memory necessary for phase screen generation.
    void genDealloc();

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

// template<typename aoSystemT>
// turbLayer<aoSystemT>::turbLayer()
// {
// }

template<typename aoSystemT>
void turbLayer<aoSystemT>::setLayer( turbAtmosphere<aoSystemT> * parent,
                                     int layerNo,
                                     uint32_t scrnSz
                                   )
{
    m_parent = parent;
    m_layerNo = layerNo;
   
    m_scrnSz = scrnSz;

    if(m_parent == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbLayer::setLayer", "parent is not set (m_parent is nullptr)"); 
    }

    if(m_parent->aosys() == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbLayer::setLayer", "parent ao system is not set (m_parent->aosys() is nullptr)"); 
    }

    realT vwind = m_parent->aosys()->atm.layer_v_wind(m_layerNo);
    realT dwind = m_parent->aosys()->atm.layer_dir(m_layerNo);
    realT pupD = m_parent->aosys()->D();
    uint32_t wfSz = m_parent->wfSz();

    m_dx = vwind*cos(dwind)/(pupD/wfSz);
    m_dy = vwind*sin(dwind)/(pupD/wfSz);
   
    alloc();
}

template<typename aoSystemT>
uint32_t turbLayer<aoSystemT>::layerNo()
{
    return m_layerNo;
}

template<typename aoSystemT>
uint32_t turbLayer<aoSystemT>::scrnSz()
{
    return m_scrnSz;
}

template<typename aoSystemT>
void turbLayer<aoSystemT>::alloc()
{
    if(m_parent == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbLayer::alloc", "parent is not set (m_parent is nullptr)"); 
    }

    if(m_parent->aosys() == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbLayer::alloc", "parent ao system is not set (m_parent->aosys() is nullptr)"); 
    }

    m_phase.resize(m_scrnSz, m_scrnSz);
   
    uint32_t wfSz = m_parent->wfSz();

    uint32_t buffSz = m_parent->buffSz();

    m_shiftPhaseWP.resize( wfSz+2*buffSz, wfSz+2*buffSz);
   
    m_shiftPhase.resize( wfSz+2*buffSz, wfSz+2*buffSz);
   
    m_shiftPhaseWork.resize( wfSz+2*buffSz, wfSz+2*buffSz);
   
    initRandom();
   
    m_x0 = 0;
    m_y0 = 0;
    m_last_wdx = m_scrnSz + 1;
    m_last_wdy = m_scrnSz + 1;

    //Now set up for generation
    realT r0 = m_parent->aosys()->atm.r_0(m_parent->aosys()->lam_sci());
    realT L0 = m_parent->aosys()->atm.L_0(m_layerNo);
    realT l0 = m_parent->aosys()->atm.l_0(m_layerNo);

    m_psd.resize(m_scrnSz, m_scrnSz);

    m_freq.resize(m_scrnSz, m_scrnSz);
    sigproc::frequencyGrid<imageT>(m_freq, m_parent->aosys()->D()/wfSz);

    realT beta = 0.0218/pow( r0, 5./3.);
    realT sqrt_alpha = 0.5*11./3.;

    realT L02;
    if(L0 > 0) L02 = 1.0/(L0*L0);
    else L02 = 0;

    #pragma omp parallel for
    for(size_t jj =0; jj < m_scrnSz; ++jj)
    {
        for(size_t ii=0; ii < m_scrnSz; ++ii)
        {
            realT p;
            if(m_freq(ii,jj) == 0 && L02 == 0)
            {
                p = 0;
            }
            else
            {
                p = beta / pow( pow(m_freq(ii,jj),2) + L02, sqrt_alpha);
                if(l0 > 0 ) p *= exp(-1*pow( m_freq(ii,jj)*l0, 2));
            }

            realT Ppiston = 0;
            realT Ptiptilt = 0;
            if(m_parent->aosys()->psd.subPiston())
            {
                Ppiston = pow(2*math::func::jinc(math::pi<realT>() * m_freq(ii,jj) * m_parent->aosys()->D()), 2);
            }
            if(m_parent->aosys()->psd.subTipTilt())
            {
                Ptiptilt = pow(4*math::func::jincN(2,math::pi<realT>() * m_freq(ii,jj) * m_parent->aosys()->D()), 2);
            }

            m_psd(ii,jj) = sqrt(p*(1.0-Ppiston-Ptiptilt));

        }
    }

    if(m_parent->shLevel() > 0)
    {
        m_psd(0,0) = 0;
        m_psd(0,1) *= 0.5*0.5;
        m_psd(1,0) *= 0.5*0.5;
        m_psd(m_psd.rows()-1, 0) *= 0.5*0.5;
        m_psd(0, m_psd.cols()-1) *= 0.5*0.5;

        m_psd(1,1) *= 0.75*0.75;
        m_psd(m_psd.rows()-1, m_psd.cols()-1) *= 0.75*0.75;
        m_psd(m_psd.rows()-1, 1) *= 0.75*0.75;
        m_psd(1, m_psd.cols()-1) *= 0.75*0.75;
    }
}

template<typename aoSystemT>
void turbLayer<aoSystemT>::genDealloc()
{
    m_freq.resize(0,0);
    m_psd.resize(0,0);    
}

template<typename aoSystemT>
void turbLayer<aoSystemT>::shift( realT dt )
{
    int wdx, wdy;
    realT ddx, ddy;
   
    ddx = m_x0 + m_dx*dt;
    ddy = m_y0 + m_dy*dt;
   
    wdx = (int) trunc(ddx);
    ddx -= wdx;

    wdy = (int) trunc(ddy);
    ddy -= wdy;
   
    wdx %= m_scrnSz;
    wdy %= m_scrnSz;

    //Check for a new whole-pixel shift
    if(wdx != m_last_wdx || wdy != m_last_wdy)
    {
       //Need a whole pixel shift (this also extracts the m_wfSz + m_buffSz subarray)
       improc::imageShiftWP(m_shiftPhaseWP, m_phase, wdx, wdy);
    }

    //Do the sub-pixel shift      
    improc::imageShift( m_shiftPhase, m_shiftPhaseWP, ddx, ddy, improc::cubicConvolTransform<realT>(-0.5));
    
    m_last_wdx = wdx;
    m_last_wdy = wdy;
}

template<typename aoSystemT>
void turbLayer<aoSystemT>::initRandom()
{
    uniVar.seed();
}

template<typename aoSystemT>
void turbLayer<aoSystemT>::shiftRandom( bool nofract )
{
   int wdx, wdy;
   realT ddx, ddy;

   ddx = uniVar * (m_scrnSz);
   ddy = uniVar * (m_scrnSz);
   
   wdx = (int) trunc(ddx);
   ddx -= wdx;
   wdy = (int) trunc(ddy);
   ddy -= wdy;
   
   if(nofract)
   {  
      ddx=0;
      ddy=0;
   }
   
   improc::imageShiftWP(m_shiftPhaseWP, m_phase, wdx, wdy);
   
   if(ddx !=0 || ddy != 0)
   {
      improc::imageShift( m_shiftPhase, m_shiftPhaseWP, ddx, ddy, improc::cubicConvolTransform<realT>(-0.5));
   }
   else
   {
      m_shiftPhase = m_shiftPhaseWP;
   }
}

} //namespace sim
} //namespace AO
} //namespace mx

#endif //mx_AO_sim_turbLayer_hpp
