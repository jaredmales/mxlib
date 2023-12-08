/** \file turbSubHarm.hpp
  * \brief A class to manage low-frequency sub-harmonic phase screen generation in atmospheric turbulence.
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


#ifndef mx_AO_sim_turbSubHarm_hpp
#define mx_AO_sim_turbSubHarm_hpp

#include <vector>

#include "turbAtmosphere.hpp"

#include "../../base/changeable.hpp"

#include "../../improc/eigenImage.hpp"
#include "../../improc/eigenCube.hpp"

#include "../../math/constants.hpp"
#include "../../math/func/jinc.hpp"

namespace mx
{
namespace AO
{
namespace sim
{

/// A class to manage low-frequency sub-harmonic phase screen generation in atmospheric turbulence.
/** Implements the method of Johansson and Gavel (1994)\cite johansson_gavel_1994. This is 
  * needed to generate adequate tip/tilt variance with large outer scales.
  * 
  * \ingroup mxAOSim
  */
template<typename _turbAtmosphereT>
class turbSubHarmonic : public base::changeable<turbSubHarmonic<_turbAtmosphereT>>
{

public:

    typedef _turbAtmosphereT turbAtmosphereT;
    typedef typename turbAtmosphereT::realT realT;

protected:

    /** \name Configuration Parameters 
      * @{ 
      */

    turbAtmosphereT * m_turbAtmo {nullptr}; ///< Pointer to the parent atmosphere object.
    
    unsigned m_level {1}; ///< The subharmonic level to apply.

    bool m_outerSubHarmonics {true}; ///< Whether or not to include the outer subharmonics

    bool m_preCalc {false}; ///< whether or not the modes are pre-calculated.

    ///@}

    /** \name Internal State
      * @{
      */

    uint32_t m_scrnSz; ///< The wavefront screen size from the layer being simulated.

    std::vector<realT> m_noise; ///< Vector of Gaussian deviates prepared for each screen generation.

    std::vector<realT> m_m; ///< m-coordinate fractional spatial frequency indices of the subharmonics 
    std::vector<realT> m_n; ///< n-coordinate fractional spatial frequency indices of the subharmonics

    std::vector<realT> m_sqrtPSD; //the square-root of the PSD at each point
    
    improc::eigenCube<realT> m_modes; ///< the pre-calculated modes

    ///@}

public:

    /** \name Construction
      * @{ 
      */
    turbSubHarmonic();

    ///@}

    /** \name Configuration
      * @{ 
      */

    /// Set the pointer to the turbulent atmosphere
    void turbAtmo( turbAtmosphereT * atm /**< [in] the new pointer to an AO atmosphere*/ );

    /// Get the pointer to the AO atmosphere
    /**
      * \returns the the current pointer to the AO atmosphere
      */
    turbAtmosphereT * turbAtmo();

    /// Set the subharmonic level to apply
    void level( uint32_t ml /**< [in] the new level*/ );

    /// Get the subharmonic level
    /**
      * \returns the current subharmonic level
      */
    uint32_t level();

    /// Set whether or not the outer subharmonics are included
    void outerSubHarmonics( bool osh /**< [in] the new value of the \ref m_outerSubHarmonics flag */);

    /// Get whether or not the outer subharmonics are included
    /**
      * \returns the current value of the \ref m_outerSubHarmonics flag 
      */
    bool outerSubHarmonics();


    /// Set whether or not to pre-calculate the modes
    void preCalc( bool pc /**< [in] the new value of the \ref m_preCalc flag*/ );

    /// Get whether or not the modes are pre-calculated
    /**
      * \returns the current value of the preCalc flag
      */
    bool preCalc();

    ///@}

    /** \name Screen Generation
      * @{ 
      */

    /// Allocate needed memory and initialize the subharmonic transform 
    void initGrid( uint32_t layerNo );

    /// Generate a realization of the subharmonic phase screen and add it to the input screen
    void screen( improc::eigenImage<realT> & scrn /**< [in] the input phase screen to which to add the low-frequency screen.  Must be m_scrnSz X m_scrnsz*/ );

    /// Deallocate memory
    void deinit();

    ///@}
};

template<typename turbAtmosphereT>
turbSubHarmonic<turbAtmosphereT>::turbSubHarmonic()
{
}

template<typename turbAtmosphereT>
void turbSubHarmonic<turbAtmosphereT>::turbAtmo( turbAtmosphereT * turbatm )
{
    if(turbatm != m_turbAtmo)
    {
        m_turbAtmo = turbatm;
        this->changed();
    }
}

template<typename turbAtmosphereT>
turbAtmosphereT * turbSubHarmonic<turbAtmosphereT>::turbAtmo()
{
    return m_turbAtmo;
}

template<typename turbAtmosphereT>
void turbSubHarmonic<turbAtmosphereT>::level( uint32_t ml )
{
    if(ml != m_level)
    {
        m_level = ml;
        this->changed();
    }
}

template<typename turbAtmosphereT>
uint32_t turbSubHarmonic<turbAtmosphereT>::level()
{
    return m_level;
}

template<typename turbAtmosphereT>
void turbSubHarmonic<turbAtmosphereT>::outerSubHarmonics( bool osh )
{
    if(osh != m_outerSubHarmonics)
    {
        m_outerSubHarmonics = osh;
        this->changed();
    }
}

template<typename turbAtmosphereT>
bool turbSubHarmonic<turbAtmosphereT>::outerSubHarmonics()
{
    return m_outerSubHarmonics;
}

template<typename turbAtmosphereT>
void turbSubHarmonic<turbAtmosphereT>::preCalc( bool pc )
{
    if(pc != m_preCalc)
    {
        m_preCalc = pc;
        this->changed();
    }
}

template<typename turbAtmosphereT>
bool turbSubHarmonic<turbAtmosphereT>::preCalc()
{
    return m_preCalc;
}

template<typename turbAtmosphereT>
void turbSubHarmonic<turbAtmosphereT>::initGrid( uint32_t layerNo )
{
    int N;

    if(m_turbAtmo == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbSubHarmonic::initGrid", "atmosphere is not set (m_turbAtmo is nullptr)"); 
    }

    if(m_turbAtmo->aosys() == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbSubHarmonic::initGrid", "ao system is not set (m_turbAtmo->m_aosys is nullptr)"); 
    }

    if(m_turbAtmo->nLayers() <= layerNo)
    {
        mxThrowException(err::invalidconfig , "mx::AO::sim::turbSubHarmonic::initGrid", "atmosphere is not setup (m_turbAtmo->m_layers size is <= layerNo)"); 
    }

    if(m_level == 0)
    {
        N = 0;
    }
    else
    {
        N = 36.0 + 32.0*(m_level-1);

        if(m_outerSubHarmonics)
        {
            N += 20;
        }
    }

    m_m.resize(0); //resized dynamically below
    m_n.resize(0);
    m_sqrtPSD.resize(N);
    m_noise.resize(N);

    m_scrnSz = m_turbAtmo->layer(layerNo).scrnSz();

    
    if(N == 0) 
    {
        
        m_modes.clear();
        this->setChangePoint();
        return;
    }

    realT r0 = m_turbAtmo->aosys()->atm.r_0(m_turbAtmo->aosys()->lam_sci());
    realT D = m_turbAtmo->aosys()->D();
    uint32_t wfSz = m_turbAtmo->wfSz();
    
    realT beta = 0.0218/pow(r0, 5./3.)/pow((D/wfSz)*m_scrnSz,2) ;

    realT sqrt_alpha = 0.5*11./3.;

    realT L0 = m_turbAtmo->aosys()->atm.L_0(layerNo);

    realT L02;
    if(L0 > 0) L02 = 1.0/pow( L0, 2);
    else L02 = 0;

    if(m_preCalc)
    {
        m_modes.resize(m_scrnSz, m_scrnSz, N);
    }

    
    std::vector<realT> scs;

    if(m_outerSubHarmonics)
    {
        realT sc = 0.25;

        m_m = std::vector<realT>({ -1.25, -0.75, -0.25,  0.25,  0.75,  1.25, -1.25,  1.25, -1.25,  1.25, -1.25, 1.25, -1.25, 1.25, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25});
        m_n = std::vector<realT>({ -1.25, -1.25, -1.25, -1.25, -1.25, -1.25, -0.75, -0.75, -0.25, -0.25,  0.25, 0.25,  0.75, 0.75,  1.25,  1.25,  1.25, 1.25, 1.25, 1.25});
        scs.resize(m_m.size(), sc);
    }
    
    for(int nl = 1; nl <= m_level; ++nl)
    {
        realT sc = pow(3.0, -nl);

        for(int mp = -3; mp < 3; ++mp)
        {
            for(int np = -3; np < 3; ++np)
            {
                if(nl < m_level)
                {
                    if(mp == -1)
                    {
                        if(np == -1 || np == 0) continue;
                    }
                    else if(mp == 0)
                    {
                        if(np == -1 || np == 0) continue;
                    }
                }

                m_m.push_back(sc*(mp+0.5));
                m_n.push_back(sc*(np+0.5));
                scs.push_back(sc);
            }
        }   
    }
    
    int n = 0;
    for(int n = 0; n < m_m.size(); ++n)
    {
        realT k = sqrt((pow(m_m[n],2) + pow(m_n[n],2)))/((D/wfSz)*m_scrnSz);
                
        realT tpsd = beta / pow( k*k + L02, sqrt_alpha);

        realT Ppiston = 0;
        realT Ptiptilt = 0;
        if(m_turbAtmo->aosys()->psd.subPiston())
        {
            Ppiston = pow(2*math::func::jinc(math::pi<realT>() * k * D), 2);
        }

        if(m_turbAtmo->aosys()->psd.subTipTilt())
        {
            Ptiptilt = pow(4*math::func::jincN(2,math::pi<realT>() * k * D), 2);
        }

        m_sqrtPSD[n] = scs[n]*sqrt(tpsd*(1-Ppiston-Ptiptilt));
                
        if(m_preCalc)
        {
            for(int cc=0; cc < m_scrnSz; ++cc)
            {
                realT np = cc - 0.5*m_scrnSz;
                for(int rr=0; rr < m_scrnSz; ++rr)
                {
                    realT mp = rr - 0.5*m_scrnSz;                    
                    m_modes.image(n)(rr,cc) = m_sqrtPSD[n]*cos(math::two_pi<realT>()*(m_m[n]*mp + m_n[n]*np)/m_scrnSz);
                }
            }
        }
    }

    this->setChangePoint();
}

template<typename turbAtmosphereT>
void turbSubHarmonic<turbAtmosphereT>::screen(improc::eigenImage<realT> & scrn)
{
    if(m_level == 0)
    {
        return;
    }

    if(m_turbAtmo == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbSubHarmonic::screen", "atmosphere is not set (m_turbAtmo is nullptr)"); 
    }

    if(this->isChanged())
    {
        mxThrowException(err::invalidconfig, "mx::AO::sim::turbSubHarmonic::screen", "configuration has changed but not re-initialized"); 
    }

    if(scrn.rows() != m_scrnSz || scrn.cols() != m_scrnSz)
    {
        mxThrowException(err::sizeerr, "mx::AO::sim::turbSubHarmonic::screen", "input screen is not the right size"); 
    }

    //Check that we're allocated
    if(m_preCalc)
    {
        if(m_modes.rows() != scrn.rows() || m_modes.cols() != scrn.cols() || m_modes.planes() != m_noise.size())
        {
            mxThrowException(err::sizeerr, "mx::AO::sim::turbSubHarmonic::screen", "modes cube wrong size, call initGrid().");
        }
    }
    else
    {
        if( m_noise.size() != m_m.size() || m_noise.size() != m_n.size() || m_sqrtPSD.size() != m_noise.size())
        {
            mxThrowException(err::sizeerr, "mx::AO::sim::turbSubHarmonic::screen", "vectors not allocated, call initGrid().");
        }
    }

    //Now fill in the noise
    for(size_t n=0; n < m_noise.size(); ++n) 
    {
        m_noise[n] = 2*m_turbAtmo->normVar();
    }

    #pragma omp parallel for
    for(int cc = 0; cc < m_scrnSz; ++cc)
    {
        realT np = cc - 0.5*m_scrnSz;
        for(int rr=0; rr < m_scrnSz; ++rr)
        {
            realT mp = rr - 0.5*m_scrnSz;

            if(m_preCalc)
            {
                for(unsigned n =0; n < m_noise.size(); ++n)
                {
                    scrn(rr,cc) += m_noise[n]*m_modes.image(n)(rr,cc);
                }
            }
            else
            {
                for(unsigned n =0; n < m_m.size(); ++n)
                {
                    scrn(rr,cc) += m_noise[n]*m_sqrtPSD[n]*cos(math::two_pi<realT>()*(m_m[n]*mp + m_n[n]*np)/m_scrnSz);
                }
            }
        }
    }
}

template<typename turbAtmosphereT>
void turbSubHarmonic<turbAtmosphereT>::deinit()
{
    m_m.clear();
    m_n.clear();
    m_sqrtPSD.clear();
    m_noise.clear();
    m_modes.resize(0, 0, 0);

    this->changed();
}

} //namespace sim
} //namespace AO
} //namespace mx

#endif //mx_AO_sim_turbSubHarm_hpp

