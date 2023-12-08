/** \file turbAtmosphere.hpp
  * \brief Declaration and definition of a turbulent atmosphere.
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

#ifndef mx_AO_sim_turbAtmosphere_hpp
#define mx_AO_sim_turbAtmosphere_hpp

#include <vector>
#include <iostream>

#include "../../sigproc/psdFilter.hpp"
#include "../../sigproc/psdUtils.hpp"

#include "../../base/changeable.hpp"

#include "../../improc/milkImage.hpp"

#include "../../math/constants.hpp"
#include "../../math/func/jinc.hpp"
#include "../../math/randomT.hpp"

#include "../../ioutils/fits/fitsFile.hpp"
#include "../../ioutils/stringUtils.hpp"
#include "../../sys/timeUtils.hpp"

#include "turbLayer.hpp"
#include "turbSubHarmonic.hpp"
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

/// A turbulent atmosphere simulator
/** Generate multi-layer phase screens using PSD filtering with the FFT. Manages
  * the shifting and combination of the phase screens (without propagation).
  * 
  * The \ref turbSubHarmonic class provides low-frequency sub-harmonics if desired
  * 
  * The \ref turbLayer class manages the actual interpolation.
  * 
  * \todo add facility for switching interpolators
  * \todo need to include wrap detection when using subharmonics (e.g. frames needs to be right)
  * 
  * \ingroup mxAOSim
  */
template<typename _aoSystemT>
struct turbAtmosphere : public base::changeable<turbAtmosphere<_aoSystemT>>
{

public:    
    typedef _aoSystemT aoSystemT;
    typedef typename aoSystemT::realT realT;
    typedef Eigen::Array<realT, -1, -1> imageT;

protected:

    /** \name Configuration Data
      * @{ 
      */
    uint32_t m_wfSz {0}; ///< Size of the wavefront in pixels.

    uint32_t m_buffSz {0}; ///< Buffer to apply around wavefront for interpolation. 

    aoSystemT * m_aosys {nullptr}; ///< The AO system object containing all system parameters.  See \ref aoSystem.

    uint32_t m_shLevel {0}; /**< Number of subharmonic levels.  0 means no subharmonics.  Generally only 
                              *  1 level is needed to obtain good Tip/Tilt results.  See \ref turbSubHarmonic.
                              */

    bool m_outerSubHarmonics {true}; /**< Whether or not the outer subharmonics are included.
                                       *  See \ref m_turbSubHarmonic.
                                       */
    bool m_shPreCalc {true}; /**< Whether or not to pre-calculate the subharmonic modes.  This is a
                               *  trade of speed vs. memory use. See \ref turbSubHarmonic.
                               */

    bool m_retain {false}; /**< Whether or not to retain working memory after screen generation.  One 
                             *  would set to true if many screens will be calculated in monte carlo fashion.
                             *  For a standard case of one set of screens being shifted by wind velocity, 
                             *  this should be false to minimize memory use.
                             */

    bool m_forceGen {true}; /**< Force generation of new screens if true.  Note that this class does not presently
                              *  do any checks of saved screens to verify they match the current configuration.
                              *  Set to true with caution.
                              * 
                              * \todo need check for changes in parameters to decide if saved frames are valid per layer
                              */
    
    std::string m_dataDir; /**< Specifies a path where phase screens are stored.  Leaving this unset or "" is equivalent to
                             *  \ref m_forceGen` == true`.
                             */

    imageT * m_pupil {0}; ///< A pointer to the pupil mask.

    realT m_timeStep {0}; ///< Length of each iteration, in seconds.
    
    ///@}

    /** \name Internal State
      * @{ 
      */
    std::vector<turbLayer<aoSystemT>> m_layers; ///< Vector of turbulent layers.

    std::vector<turbSubHarmonic<turbAtmosphere>> m_subHarms; ///< Sub-harmonic layers

    size_t m_frames {0}; ///< Length of the turbulence sequence.

    int m_nWf {0}; ///< Number of iterations which have occurred.

    math::normDistT<realT> m_normVar; ///< Normal random deviate generator.  This seeded in the constructor.

public:

    /** \name Construction and Configuration
      * @{
      */

    ///Default c'tor
    turbAtmosphere();

    ///Setup the overall atmosphere.
    /**
      */
    void setup( uint32_t wfSz,     ///< [in] The size of the wavefront in pixels.
                uint32_t buffSz,   ///< [in] The size of the interpolation buffer to use.
                aoSystemT * aosys, ///< [in] Pointer to an AO System.  See \ref m_aosys.
                uint32_t shLevel   ///< [in] number of subharmonic levels to use.  0 turns off subharmonics.  See \ref m_shLevel.
              );

    /// Set the wavefront size
    /**
      * see \ref m_wfSz
      */
    void wfSz( uint32_t ws /**< [in] the new wavefront size */);

    /// Get the wavefront size
    /**
      * \returns the current value of \ref m_wfSz
      */
    uint32_t wfSz();

    /// Set the interpolation buffer size
    /**
      * see \ref m_buffSz
      */
    void buffSz( uint32_t bs /**< [in] the new buffer size*/);

    /// Get the interpolation buffer size
    /**
      * \returns the current value of \ref m_buffSz
      */
    uint32_t buffSz();

    /// Set the pointer to an AO system object
    /**
      * see \ref m_aosys
      */
    void aosys( aoSystemT * aos /**< [in] the new pointer to an AO system object*/);

    /// Get the pointer to an AO system object
    /**
      * \returns the current value of \ref m_aosys
      */
    aoSystemT * aosys();

    /// Set the subharmonic level
    /**
      * see \ref m_shLevel
      */
    void shLevel( uint32_t shl /**< [in] the new subharmonic level*/);

    /// Get the subharmonic level
    /**
      * \returns the current value of \ref m_shLevel
      */
    uint32_t shLevel();

    /// Set whether outer subharmonics are included
    /**
      * see \ref m_outerSubHarmonics
      */
    void outerSubHarmonics( bool osh /**< [in] the new value of flag controlling outer subharmonics */);

    /// Get whether outer subharmonics are included
    /**
      * \returns the current value of \ref m_outerSubHarmonics
      */
    bool outerSubHarmonics();

    /// Set whether subharmonic modes are pre-calculated
    /**
      * see \ref m_shPreCalc
      */
    void shPreCalc( bool shp /**< [in] the new value of flag controlling subharmonic mode precalculation */);

    /// Get whether subharmonic modes are pre-calculated 
    /**
      * \returns the current value of \ref m_shPreCalc
      */
    bool shPreCalc();

    /// Set whether memory for screen generation is retained
    /**
      * see \ref m_retain
      */
    void retain( bool rtn /**< [in] the new value of the flag controlling whether memory is retained*/);

    /// Get whether memory for screen generation is retained
    /**
      * \returns the current value of \ref m_retain
      */
    bool retain();

    /// Set whether new screen generation is forced
    /**
      * see \ref m_forceGen
      */
    void forceGen( bool fg /**< [in] the new value of the flag controlling whether screen generation is forced*/);

    /// Get whether new screen generation is forced
    /**
      * \returns the current value of m_forceGen
      */
    bool forceGen();

    /// Set the data directory for saving phase screens
    /**
      * see \ref m_dataDir
      */
    void dataDir( const std::string & dd /**< [in] the new data directory for phase screen saving*/);

    /// Get the data directory for saving phase screens
    /**
      * \returns the current value of m_dataDir
      */
    std::string dataDir();

    /// Setup the layers and prepare them for phase screen generation.
    /** The number of entries in \p scrnSz must mach the number of layers in the atmosphere
      * of \ref m_aosys. 
      */
    void setLayers( const std::vector<size_t> & scrnSz /**< [in] the screen sizes for each layer.*/);

    /// Setup the layers and prepare them for phase screen generation.
    /** This sets all layer phase screens to be the same size.
      * 
      */
    void setLayers( const size_t scrnSz /**< [in] the screen sizes for all layers.*/);

    /// Get the number of layers
    /**
      * \returns the size of the layers vector \ref m_layers 
      */
    uint32_t nLayers();

    /// Get a layer
    /**
      * \returns a reference to on entry in \ref m_layers 
      */
    turbLayer<aoSystemT> & layer(uint32_t n);

    /// Get the random number generator
    /**
      * \returns a reference to \ref m_normVar 
      */
    math::normDistT<realT> & normVar();

    ///@} - Construction and Configuration


    uint32_t scrnLengthPixels();

    uint32_t scrnLengthPixels(uint32_t n);

    realT scrnLength();

    realT scrnLength(uint32_t n);

    uint32_t maxShift( realT dt );

    uint32_t maxShift( uint32_t n,
                       realT dt 
                     );


    /** \name Screen Generation
      * @{
      */

    /// Generate all phase screens
    /** Loads them from disk if possible and \ref m_forceGen == false.
      * 
      * Deallocates working memory when done, unless \ref m_retain == true.
      */ 
    void genLayers();

    /// @}

    int shift( improc::milkImage<realT> & phase,
               realT dt );

    int frames(int f);
    size_t frames();

    int wfPS(realT ps); ///< dummy function for simulatedAOSystem.  Does nothing.

    bool _loopClosed;

    void nextWF(wavefront<realT> & wf);
};

template<typename aoSystemT>
turbAtmosphere<aoSystemT>::turbAtmosphere()
{
    m_normVar.seed();
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::setup( uint32_t ws,
                                       uint32_t bs,
                                       aoSystemT * aos,
                                       uint32_t shl
                                     ) 
{
    wfSz(ws);
    buffSz(bs);
    aosys(aos);
    shLevel(shl);
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::wfSz( uint32_t ws )
{
    if(ws != m_wfSz) 
    {
        m_wfSz = ws;
        this->changed();
    }
}

template<typename aoSystemT>
uint32_t turbAtmosphere<aoSystemT>::wfSz()
{
    return m_wfSz;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::buffSz( uint32_t bs )
{
    if(bs != m_buffSz)
    {
        m_buffSz = bs;
        this->changed();
    }
}

template<typename aoSystemT>
uint32_t turbAtmosphere<aoSystemT>::buffSz()
{
    return m_buffSz;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::aosys( aoSystemT * aos)
{
    if(aos != m_aosys)
    {
        m_aosys = aos;
        this->changed();
    }
}

template<typename aoSystemT>
aoSystemT * turbAtmosphere<aoSystemT>::aosys()
{
    return m_aosys;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::shLevel( uint32_t shl )
{
    if(shl != m_shLevel)
    {
        m_shLevel = shl;
        this->changed();
    }
}

template<typename aoSystemT>
uint32_t turbAtmosphere<aoSystemT>::shLevel()
{
    return m_shLevel;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::outerSubHarmonics( bool osh )
{
    if(osh != m_outerSubHarmonics)
    {
        m_outerSubHarmonics = osh;
        this->changed();
    }
}

template<typename aoSystemT>
bool turbAtmosphere<aoSystemT>::outerSubHarmonics()
{
    return m_outerSubHarmonics;
}


template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::shPreCalc( bool shp )
{
    if(shp != m_shPreCalc)
    {
        m_shPreCalc = shp;
        this->changed();
    }
}

template<typename aoSystemT>
bool turbAtmosphere<aoSystemT>::shPreCalc()
{
    return m_shPreCalc;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::retain( bool rtn )
{
    if(rtn != m_retain)
    {
        m_retain = rtn;
        this->changed();
    }
}

template<typename aoSystemT>
bool turbAtmosphere<aoSystemT>::retain()
{
    return m_retain;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::forceGen( bool fg )
{
    if(fg != m_forceGen)
    {
        this->changed();
    }

    m_forceGen = fg;
}

template<typename aoSystemT>
bool turbAtmosphere<aoSystemT>::forceGen()
{
    return m_forceGen;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::dataDir( const std::string & dd )
{
    if(dd != m_dataDir)
    {
        this->changed();
    }

    m_dataDir = dd;
}

template<typename aoSystemT>
std::string turbAtmosphere<aoSystemT>::dataDir()
{
    return m_dataDir;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::setLayers( const std::vector<size_t> & scrnSz)
{
    if(m_aosys == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbAtmosphere::setLayers", "ao system is not set (m_aosys is nullptr)"); 
    }

    size_t nLayers = scrnSz.size();
    
    if(nLayers != m_aosys->atm.n_layers())
    {
        mxThrowException(err::invalidarg, "mx::AO::sim::turbAtmosphere::setLayers", "Size of scrnSz vector does not match atmosphere.");
    }

    if(nLayers != m_layers.size())
    {
        this->changed();
    }

    m_layers.resize(nLayers);

    for(size_t i=0; i< nLayers; ++i)
    {
        m_layers[i].setLayer(this, i, scrnSz[i]);
    }

    if(m_shLevel > 0)
    {
        if(nLayers != m_subHarms.size())
        {
            this->changed();
        }

        m_subHarms.resize(nLayers);

        for(size_t i=0; i< nLayers; ++i)
        {
            m_subHarms[i].turbAtmo(this);
            m_subHarms[i].preCalc(m_shPreCalc);
            m_subHarms[i].level(m_shLevel);
            m_subHarms[i].outerSubHarmonics(m_outerSubHarmonics);
            m_subHarms[i].initGrid(i);
        }
    }
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::setLayers( const size_t scrnSz )
{
    if(m_aosys == nullptr)
    {
        mxThrowException(err::paramnotset, "mx::AO::sim::turbAtmosphere::setLayers", "atmosphere is not set (m_atm is nullptr)"); 
    }

    size_t n = m_aosys->atm.n_layers();

    setLayers( std::vector<size_t>(n, scrnSz));
}

template<typename aoSystemT>
uint32_t turbAtmosphere<aoSystemT>::nLayers()
{
    return m_layers.size();
}

template<typename aoSystemT>
turbLayer<aoSystemT> & turbAtmosphere<aoSystemT>::layer(uint32_t n)
{
    if(n >= m_layers.size())
    {
        mxThrowException(err::invalidarg, "mx::AO::sim::turbAtmosphere::layer", "n too large for number of layers.");
    }

    return m_layers[n];
}

template<typename aoSystemT>
math::normDistT<typename aoSystemT::realT> & turbAtmosphere<aoSystemT>::normVar()
{
    return m_normVar;
}

template<typename aoSystemT>
uint32_t turbAtmosphere<aoSystemT>::scrnLengthPixels()
{
    return 0;
}

template<typename aoSystemT>
uint32_t turbAtmosphere<aoSystemT>::scrnLengthPixels(uint32_t n)
{

    return 0;
/*    realT dx = (scrnSz-wfSz-turb.buffSz()) / cos(m_aosys->atm.layer_v_wind(n));

    pow(2*pow(scrnSz-wfSz-turb.buffSz(),2), 0.5) - 1;
    */
}

template<typename aoSystemT>
realT turbAtmosphere<aoSystemT>::scrnLength()
{
    return 0;
}

template<typename aoSystemT>
realT turbAtmosphere<aoSystemT>::scrnLength(uint32_t n)
{
    return 0;
}

template<typename aoSystemT>
uint32_t turbAtmosphere<aoSystemT>::maxShift( realT dt )
{
    return 0;
}

template<typename aoSystemT>
uint32_t turbAtmosphere<aoSystemT>::maxShift( uint32_t n,
                                              realT dt 
                                            )
{
    return 0;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::genLayers()
{
    if(m_dataDir != "" && !m_forceGen)
    {
        std::string fbase = m_dataDir;
        fbase += "/";
        fbase += "layer_";

        fits::fitsFile<realT> ff;
        std::string fname;
        for(size_t i=0; i< m_layers.size(); ++i)
        {
            fname = fbase + ioutils::convertToString<int>(i) + ".fits";
            ff.read(m_layers[i].m_phase, fname);
        }
        
        this->setChangePoint();

        return;
    }

    sigproc::psdFilter<realT,2> filt;

    for(size_t i=0; i< m_layers.size(); ++i)
    {
        if(m_layers[i].m_phase.rows() != m_layers[i].scrnSz() || m_layers[i].m_phase.cols() != m_layers[i].scrnSz())
        {
            mxThrowException(err::sizeerr, "mx::AO::sim::turbAtmosphere::genLayers", "layer phase not allocated.");
        }

        if(m_layers[i].m_freq.rows() != m_layers[i].scrnSz() || m_layers[i].m_freq.cols() != m_layers[i].scrnSz())
        {
            mxThrowException(err::sizeerr, "mx::AO::sim::turbAtmosphere::genLayers", "layer freq not allocated.");
        }

        if(m_layers[i].m_psd.rows() != m_layers[i].scrnSz() || m_layers[i].m_psd.cols() != m_layers[i].scrnSz())
        {
            mxThrowException(err::sizeerr, "mx::AO::sim::turbAtmosphere::genLayers", "layer psd not allocated.");
        }

        for(size_t jj = 0; jj < m_layers[i].m_scrnSz; ++jj)
        {
            for(size_t ii = 0; ii < m_layers[i].m_scrnSz; ++ii)
            {
                m_layers[i].m_phase(ii,jj) = m_normVar;
            }
        }

        realT dkx = m_layers[i].m_freq(1,0)-m_layers[i].m_freq(0,0);
        realT dky = m_layers[i].m_freq(0,1)-m_layers[i].m_freq(0,0);
        filt.psdSqrt(m_layers[i].m_psd, dkx, dky);

        filt(m_layers[i].m_phase);
    
        if(m_shLevel > 0)
        {
            m_subHarms[i].screen(m_layers[i].m_phase);
        }
    }
   
    if(!m_retain)
    {
        for(size_t i=0; i< m_layers.size(); ++i)
        {
            m_layers[i].genDealloc();
        }

        for(size_t i=0; i< m_subHarms.size(); ++i)
        {
            m_subHarms[i].deinit();
        }
    }

    if(m_dataDir != "")
    {
        std::string fbase = m_dataDir;
        fbase += "/";
        fbase += "layer_";

        fits::fitsFile<realT> ff;
        std::string fname;
        for(size_t i=0; i< m_layers.size(); ++i)
        {
            fname = fbase + ioutils::convertToString<int>(i) + ".fits";
            ff.write(fname, m_layers[i].m_phase);
        }
    }

    this->setChangePoint();
}

template<typename aoSystemT>
int turbAtmosphere<aoSystemT>::shift( improc::milkImage<realT> & milkPhase,
                                      realT dt 
                                    )
{
    if(this->isChanged())
    {
        mxThrowException(err::invalidconfig, "mx::AO::sim::turbAtmosphere::shift", "configuration has changed but genLayers has not been run.");
    }

    improc::eigenMap<realT> phase(milkPhase);

    if(phase.rows() != m_wfSz || phase.cols() != m_wfSz)
    {
    }

    //Don't use OMP if no multiple layers b/c it seems to make it use multiple threads in some awful way
    if(m_layers.size() > 1)
    {
        #pragma omp parallel for
        for(size_t j=0; j< m_layers.size(); ++j)
        {
            m_layers[j].shift( dt );
        }
    }
    else
    {
        for(size_t j=0; j< m_layers.size(); ++j)
        {
            m_layers[j].shift( dt );
        }
    }

    milkPhase.setWrite();
    phase.setZero();
    for(size_t j=0; j< m_layers.size(); ++j)
    {
        phase += sqrt( m_aosys->atm.layer_Cn2(j)) * m_layers[j].m_shiftPhase.block(m_buffSz, m_buffSz, m_wfSz, m_wfSz);
    }
    milkPhase.post();

    return 0;
}

template<typename aoSystemT>
int turbAtmosphere<aoSystemT>::frames(int f)
{
   m_frames = f;
   
   return 0;
}

template<typename aoSystemT>
size_t turbAtmosphere<aoSystemT>::frames()
{
   return m_frames;
}

template<typename aoSystemT>
int turbAtmosphere<aoSystemT>::wfPS(realT ps)
{
   static_cast<void>(ps);
   return 0;
}

template<typename aoSystemT>
void turbAtmosphere<aoSystemT>::nextWF(wavefront<realT> & wf)
{
    if(!m_aosys)
    {    
        mxThrowException(err::paramnotset, "mx::AO::sim::turbAtmosphere<aoSystemT>::nextWF", "the AO system pointer is not set"); 
    }

    if(!m_pupil)
    {    
        mxThrowException(err::paramnotset, "mx::AO::sim::turbAtmosphere<aoSystemT>::nextWF", "the m_pupil pointer is not set"); 
    }

    shift( wf.phase, m_nWf * m_timeStep);
    ++m_nWf;

    wf.phase *= *m_pupil;

    realT pixVal = sqrt(m_aosys->Fg())*(m_aosys->D()/m_wfSz);

    wf.amplitude = pixVal*(*m_pupil);
}

} //namespace sim
} //namespace AO
} //namespace mx

#endif //mx_AO_sim_turbAtmosphere_hpp
