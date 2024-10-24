/** \file pyramidSensor.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Declaration and definition of a standard 4 quadrant pyramid WFS.
 * \ingroup mxAO_sim_files
 *
 */

#ifndef mx_AO_sim_pyramidSensor_hpp
#define mx_AO_sim_pyramidSensor_hpp

#include "../../mxException.hpp"

#include "../../wfp/imagingUtils.hpp"
#include "../../wfp/fraunhoferPropagator.hpp"
#include "../../sys/timeUtils.hpp"

#include "../../improc/eigenImage.hpp"
#include "../../improc/imageMasks.hpp"

#include "../../math/constants.hpp"
#include "../../math/geo.hpp"

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

template <typename _realT>
struct wfsImageT
{
    typedef _realT realT;

    unsigned iterNo;

    /// The wavefront sensor detector image type
    typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

    imageT image;

    imageT tipImage;
};

/// A Pyramid Sensor Simulation
/**
 *
 * \tparam _realT is the real floating point type used for calculations
 * \tparam _detectorT is the detector used to record the PyWFS image.  Must conform to the mx::AO::sim detector
 * interface specifications.
 */
template <typename _realT, typename _detectorT>
class pyramidSensor
{
  public:
    /// The real floating point type used for calculations
    typedef _realT realT;

    /// The complex floating point type used for calculations
    typedef std::complex<realT> complexT;

    /// The wavefront data type
    typedef wavefront<realT> wavefrontT;

    /// The wavefront complex field type
    typedef wfp::imagingArray<std::complex<realT>, wfp::fftwAllocator<std::complex<realT>>, 0> complexFieldT;

    /// The wavefront sensor detector image type
    typedef _detectorT detectorT;

  public:
    /// Default c'tor
    pyramidSensor();

    /** \name Standard WFS Interface
     *
     * @{
     */

  protected:
    /* Standard WFS Interface: */
    uint32_t m_wfSz{ 0 }; ///< Size of the wavefront in pixels

    uint32_t m_detRows{ 0 }; ///< The number of rows of the WFS detector.  After forming the image the WFS detector
                             ///< plane is re-binned to this.

    uint32_t m_detCols{ 0 }; ///< The number of columns of the WFS detector.  After forming the image the WFS detector
                             ///< plane is re-binned to this.

    realT m_lambda{ 0 }; ///< Central wavelength, in meters

    /// \todo when the filter should be set with astrospectrum, and should re-calculate the central wavelength.
    /// \todo need to verify that the wavefront propagation is appropriately chromatic
    std::vector<realT> m_wavelengths;      ///< Vector of wavelengths in the WFS bandpass
    std::vector<realT> _wavelengthWeights; ///< The relative weights of the wavelengths

    int m_iTime{ 1 }; ///< Integration time in loop steps

    int m_roTime{ 1 }; ///< Readout time in loop steps

    realT m_simStep{ 0.001 }; ///< The simulation stepsize in seconds.

  public:
    /// Get the wavefront size in pixels
    /**
     * \returns the wavefront size in pixels
     */
    int wfSz();

    /// Set the wavefront size in pixels.
    /**
      */
    void wfSz(const uint32_t & sz /**< the new size*/);

    /// Get the detector rows  in pixels
    /**
     * \returns m_detRows
     */
    uint32_t detRows();

    /// Get the detector columns  in pixels
    /**
     * \returns m_detCols
     */
    uint32_t detCols();

    /// Set the detector columns in pixels.
    /**
     */
    void detSize( const uint32_t &nrows, ///< The number of rows
                  const uint32_t &ncols  ///< The number of columns
    );

    /// Get the PyWFS central wavelength
    /**
     * \returns the central wavelength in meters
     */
    realT lambda();

    /// Set the PyWFS central wavelength
    /**
     */
    void lambda( const realT &l /**< The central wavelength, in meters*/ );

    /// Get the PyWFS integration time, in time steps
    int iTime();

    /// Set the PyWFS integration time, in time steps
    void iTime( const uint32_t &it /**<  the new integration time*/ );

    /// Get the PyWFS detector readout time, in time steps
    int roTime();

    /// Set the PyWFS detector readout time, in time steps
    void roTime( const uint32_t &rt /**< the new readout time*/ );

    /// Get the simulation step-size, in seconds.
    realT simStep();

    /// Set the simulation step-size, in seconds.
    void simStep( const realT &st /**< the new simulation step size*/ );

    /// Link this WFS to an AO system simulation
    template <typename AOSysT>
    void linkSystem( AOSysT &AOSys /**< The AO system simulation to link to*/ );

    /// Sense the wavefront aberrations
    /**  \returns true if a new wavefront measurement is ready.
     *  \returns false if still integrating.
     */
    bool senseWavefront( wavefrontT &pupilPlane /**< The input wavefront to be sensed*/ );

    /// Sense the wavefront aberrations in calibration mode
    /** Allows for faster calibrations.
     */
    bool senseWavefrontCal( wavefrontT &pupilPlane /**< The input wavefront to be sensed*/ );

  public:
    /// The WFS detector
    detectorT detector;

    /// The image on the detector, resized from m_wfsImage
    wfsImageT<realT> detectorImage;

    /// @}

    /** \name Pyramid Sensor Interface
      *
      * @{
      */
protected:
    uint32_t m_nSides {4}; ///< Number of sides in the pyramid

    /// The size of the pupil in wavefront pixels.
    /** This is the maximum diameter of the pupil in wavefront pixels.
     *
     */
    uint32_t m_pupilSz{ 0 };

    /// The separation of the pupil images in fraction of a pupil.  0 <= m_pupilSep, default 1. 
    /** This sets the center-to-center separation of the pupils images in the focal plane wavefront.
      * Note that the separation in detector pixels then depends on the scaling between wavefront pixels
      * (m_wfSz) and detector pixels (m_detRows and m_detCols).
      * 
      *   
      * This sets the size of the region in the pre-detection image that each pupil image
      * takes up, and therefore the size of the pre-detection image.  
      * If the pupil (as defined in the input wavefront) is 60 pixels across
      * and m_pupilSep is set to 1.06667, then there will be a 2 pixel pad around each pupil image,
      * resulting in 4 pixels between each geometric pupil image.
      * 
      * For a standard 4-sided pyramid, the pre-detection image will be
      * 2*m_pupilSep*m_pupilSz across.  For other n-sided pyramids, m_pupilSep still specifies the size of the pupil
      * image region, but the total image size will be a function of the resultant pupil positions.
      * 
      * If m_pupilSep is less than 1, this will produce the "flattened pyramid", with overlap between
      * the pupil images.  In this case, image size will also be set by pupilSz to ensure that there are enough
      * pixels included to show all pupils.
      */
    realT m_pupilSep {1}; 

    /// The angle by which to offset the pupils, in degrees. Default is 0.
    /** If this is 0, then a 4-sided pyramid makes a square as usual.  If this is set
      * to 45 degrees, then a 4-sided pyramid makes a diamond. 
      * 
      */
    realT m_angleOffset {0};

    /// The size of the resulting PyWFS image in wavefront pixels.  
    /** If \ref m_imageSzAuto is true, this is determined by number of sides (\ref m_nSides), the pupil size (\ref m_pupilSz), and the
      * pupil separation (\ref m_pupilSep).  For a 4 sided pyramid this will be the larger of
      * 2*m_pupilSep*m_pupilSz and 2*m_pupilSz.
      * 
      * If , then this is used regardless of the optimum size.
      */
    uint32_t m_imageSz {0}; 

    bool m_imageSzAuto{ true }; ///< Flag to track if \ref m_imageSz should be set to 0.

    realT m_wfPS {0}; ///< Wavefront pixel scale, in meters/pixel

    realT m_D {0}; ///< Telescope diameter, in meters

    uint32_t m_modSteps{ 20 }; ///< Number of modulation steps in one integration.  Can be set explicitly, but will be
                               ///< calculated if \ref m_perStep is set.

    realT m_perStep{
        1 }; ///< The minimum number of lamba/D per step in the modulation.  Smaller will result in more steps.

    realT m_modRadius{ 3.0 }; ///< Radius of the modulation in pixels

  public:
    /// Get the number of pyramid sides
    /**
     * \returns the number of sides on the pyramid
     */
    int nSides();

    /// Set the number of sides on the pyramid
    /**
      */
    void nSides(const uint32_t & ns /**< The new number of sides on the pyramid*/);
        
    /// Get the minimum number of modulation steps
    /**
     * \returns m_perStep;
     */
    realT perStep();

    /// Set the minimum number of modulation steps
    /**
     * \param mSt is the new number of modulation steps
     */
    void perStep( const realT &prStp /**< The minimum number of lamba/D per step to take*/ );

    /// Get the number of modulation steps
    /**
     * \returns m_modSteps, which is defined by perStep.
     */
    int modSteps();

    /// Get the radius of modulation
    /**
     * \returns m_modRadius;
     */
    realT modRadius();

    /// Set the modulation radius
    /**
     * \param mR is the new modulation radius in lambda/D
     */
    void modRadius( const realT &mR /**< [in] the new value of modulation radius */ );

    /// Get the wavefront pixel scale in meters per pixel
    /**
     * \returns the wavefront pixel scale in meters/pixel
     */
    realT wfPS();

    /// Get the telescope diameter
    /**
     * \returns the telescope diameter in meters
     */
    realT D();

    /// Set the telescope diameter
    /**
     * \param d is the new size in meters
     */
    void D( const realT &d /**< */ );

    /// Get the pupil size in pixels
    /** This is the pupil size in un-binned wavefront space
      *
      * \returns m_pupilSz
      */
    uint32_t pupilSz();

    /// Set the pupil size in pixels.
    /** This is the size of the pupils in un-binned wavefront space.
     * See \ref m_pupilSz.
     *
     */
    void pupilSz( const uint32_t &sz /**< the new pupil size.*/ );

    /// Get the pupil separation as a fraction of pupil size
    /** This is the pupil separation in un-binned wavefront space
      *
      * \returns m_pupilSep
      */
    realT pupilSep();

    /// Set the pupil separation as a fraction of pupil size
    /** This is the separation of the pupils in un-binned wavefront space.
      * See \ref m_pupilSep.
      * 
      */
    void pupilSep(const realT & sz /**< the new pupil separation.*/);

    /// Get the angle offset
    /** See \ref m_angleOffset
      *
      * \returns m_angleOffset
      */
    realT angleOffset();

    /// Set the angle offset
    /** See \ref m_angleOffset.
      * 
      */
    void angleOffset(const realT & ao /**< the new angle offset.*/);

    /// Get the image size in wavefront pixels
    /** This is the size of the image in un-binned wavefront space
     *
     * \returns m_imageSz
     */
    uint32_t imageSz();

    /// Set the image size in wavefront pixels
    /** This is the size of the image in un-binned wavefront space
     * Setting a non-zero value also sets m_imageSizeAuto to false.
     * Setting 0 also sets m_imageSizeAuto to true.
     */
    void imageSz( const uint32_t &is );

    /// Get the value of the image size auto flag
    /** This controls whether image size is set automatically
     *
     * \returns m_imageSz
     */
    bool imageSzAuto();

    /// Set the value of the image size auto flag
    /** This controls whether image size is set automatically
     *
     */
    void imageSzAuto( const bool &ia );
    ///@}

    wfp::fraunhoferPropagator<complexFieldT> m_frProp;

    bool m_opdMaskMade{ false };
    complexFieldT m_opdMask;

    bool m_tiltsMade{ false };
    std::vector<complexFieldT> m_tilts;

    bool m_preAllocated{ false };
    complexFieldT m_pupilPlaneCF;

    // Pre-allocated working memory:

    std::vector<complexFieldT> m_th_tiltedPlane; ///< Thread-local modulated wavefront

    std::vector<complexFieldT> m_th_focalPlane; ///< Thread-local tip wavefront, used for FFT tilting

    std::vector<typename wfsImageT<realT>::imageT> m_th_focalImage; ///< Thread-local tip image

    std::vector<complexFieldT> m_th_sensorPlane; ///< Thread-local sensor-pupil-plane wavefront

    std::vector<typename wfsImageT<realT>::imageT>
        m_th_sensorImage; ///< Thread-local sensor-pupil-plane intensity image

    int m_iTime_counter{ 0 };

    int m_reading{ 0 };

    int m_roTime_counter{ 0 };

    std::vector<wavefrontT> _wavefronts;

    int m_lastWavefront{ 0 };

    /// The image formed by the WFS
    wfsImageT<realT> m_wfsImage;

  public:
    wfsImageT<realT> wfsTipImage;

  protected:
    void makeOpdMask();

    void makeTilts();

    void allocThreadMem();

    void preAllocate();

    void doSenseWavefront();
    void doSenseWavefront( wavefrontT & /**< */ );
    void doSenseWavefrontNoMod( wavefrontT & /**< */ );

    bool m_firstRun{ true };
};

template <typename realT, typename detectorT>
pyramidSensor<realT, detectorT>::pyramidSensor()
{
    iTime( m_iTime );

    m_frProp.wholePixel( 0 );
}

template <typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::wfSz()
{
    return m_wfSz;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::wfSz(const uint32_t & sz)
{
    if( m_wfSz == sz )
    {
        return;
    }

    m_wfSz = sz;

    m_tiltsMade = false;
    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT>
uint32_t pyramidSensor<realT, detectorT>::detRows()
{
    return m_detRows;
}

template <typename realT, typename detectorT>
uint32_t pyramidSensor<realT, detectorT>::detCols()
{
    return m_detCols;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::detSize( const uint32_t &nrows, const uint32_t &ncols )
{
    if( m_detRows == nrows && m_detCols == ncols )
        return;

    m_detRows = nrows;
    m_detCols = ncols;

    detector.setSize(m_detRows, m_detCols);
    detectorImage.image.resize(m_detRows, m_detCols);

    m_opdMaskMade = false; //make sure size check is run on current settings
}

template <typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::lambda()
{
    return m_lambda;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::lambda( const realT &l )
{
    m_lambda = l;

    //---------------------------------------
    //  Check if wavelength vector is filled
    //---------------------------------------
    if( m_wavelengths.size() == 0 )
    {
        m_wavelengths.resize( 1, m_lambda );
        _wavelengthWeights.resize( 1, 1.0 );
    }
}

template <typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::iTime()
{
    return m_iTime;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::iTime( const uint32_t &it )
{
    if( it < 1 )
    {
        mxThrowException( mx::err::invalidconfig, "pyramidSensor::iTime", "iTime must be >= 1" );
    }

    m_iTime = it;

    _wavefronts.resize( m_iTime + 2 );
    m_lastWavefront = -1;

    detector.expTime( m_simStep * m_iTime );
}

template <typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::roTime()
{
    return roTime;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::roTime( const uint32_t &rt )
{
    if( rt < 1 )
    {
        mxThrowException( mx::err::invalidconfig, "pyramidSensor::roTime", "roTime must be >= 1" );
    }

    m_roTime = rt;
}

template <typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::simStep()
{
    return m_simStep;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::simStep( const realT &st )
{

    m_simStep = st;

    detector.expTime( m_simStep * m_iTime );
}

template <typename realT, typename detectorT>
template <typename AOSysT>
void pyramidSensor<realT, detectorT>::linkSystem( AOSysT &AOSys )
{
    AOSys.wfs.wfPS( AOSys.m_wfPS );
    AOSys.wfs.D( AOSys.m_D );
}

template <typename realT, typename detectorT>
bool pyramidSensor<realT, detectorT>::senseWavefront( wavefrontT &pupilPlane )
{

    ++m_lastWavefront;
    if( m_lastWavefront >= _wavefronts.size() )
        m_lastWavefront = 0;
    _wavefronts[m_lastWavefront].amplitude = pupilPlane.amplitude;
    _wavefronts[m_lastWavefront].phase = pupilPlane.phase;
    _wavefronts[m_lastWavefront].iterNo = pupilPlane.iterNo;

    // Always skip the first one for averaging to center of iTime.
    if( m_firstRun )
    {
        m_firstRun = false;
        return false;
    }

    ++m_iTime_counter;

    bool rv = false;

    if( m_reading )
    {
        ++m_roTime_counter;

        if( m_roTime_counter >= m_roTime )
        {
            detector.exposeImage( detectorImage.image, m_wfsImage.image );

            detectorImage.tipImage = wfsTipImage.image;
            detectorImage.iterNo = m_wfsImage.iterNo;

            m_roTime_counter = 0;
            m_reading = 0;
            rv = true;
        }
    }

    if( m_iTime_counter >= m_iTime )
    {
        doSenseWavefront();
        m_iTime_counter = 0;

        m_reading = 1;
        m_roTime_counter = 0;
    }

    return rv;
}

template <typename realT, typename detectorT>
bool pyramidSensor<realT, detectorT>::senseWavefrontCal( wavefrontT &pupilPlane )
{

    BREAD_CRUMB;

    doSenseWavefront(pupilPlane);
    
    BREAD_CRUMB;

    detector.exposeImage( detectorImage.image, m_wfsImage.image );

    BREAD_CRUMB;

    detectorImage.tipImage = wfsTipImage.image;

    BREAD_CRUMB;

    return true;
}

/* Pyramid Specifics */

template <typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::nSides()
{
    return m_nSides;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::nSides(const uint32_t & ns)
{
    m_nSides = ns;
    m_opdMaskMade = false;
}

template <typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::wfPS()
{
    return m_wfPS;
}

template <typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::D()
{
    return m_D;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::D( const realT &d )
{
    m_D = d;
    
    if(m_pupilSz > 0) //Avoid making inf or nan so m_wfPS remains unset.  note that fast-math make detecting inf and nan hard.
    {
        m_wfPS = m_D/m_pupilSz;
    }
    else
    {
        m_wfPS = 0;
    }

    m_tiltsMade = false;
    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::perStep()
{
    return m_perStep;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::perStep( const realT &prStp )
{
    m_perStep = prStp;

    if( m_modRadius <= 0 )
    {
        m_modSteps = 0;
        return;
    }

    realT radPerStep = m_perStep / m_modRadius;

    // Get the minimum number of steps to meet m_perStep while having symmetry for the quadrants
    m_modSteps = 1;
    while( math::half_pi<realT>() / m_modSteps > radPerStep )
    {
        ++m_modSteps;
    }

    m_modSteps *= 4;

    m_tiltsMade = false;
}

template <typename realT, typename detectorT>
int pyramidSensor<realT, detectorT>::modSteps()
{
    return m_modSteps;
}

template <typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::modRadius()
{
    return m_modRadius;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::modRadius( const realT &mR )
{
    m_modRadius = mR;
    perStep( m_perStep ); // to calculate m_modSteps;

    m_tiltsMade = false;
}

template <typename realT, typename detectorT>
uint32_t pyramidSensor<realT, detectorT>::pupilSz()
{
    return m_pupilSz;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::pupilSz( const uint32_t &sz )
{
    if( m_pupilSz == sz )
    {
        return;
    }

    m_pupilSz = sz;

    if(m_pupilSz > 0) //Avoid making inf or nan so m_wfPS remains unset.  note that fast-math make detecting inf and nan hard.
    {
        m_wfPS = m_D/m_pupilSz;
    }
    else
    {
        m_wfPS = 0;
    }

    m_tiltsMade = false;
    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::pupilSep()
{
    return m_pupilSep;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::pupilSep(const realT & sz)
{
    if( m_pupilSep == sz )
    {
        return;
    }

    m_pupilSep = sz;

    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT>
realT pyramidSensor<realT, detectorT>::angleOffset()
{
    return m_angleOffset;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::angleOffset(const realT & ao)
{
    if (m_angleOffset == ao)
    {
        return;
    }

    m_angleOffset = ao;

    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT>
uint32_t pyramidSensor<realT, detectorT>::imageSz()
{
    if( !m_opdMaskMade )
    {
        makeOpdMask();
    }

    return m_imageSz;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::imageSz( const uint32_t &sz )
{
    if( m_imageSz == sz )
    {
        return;
    }

    m_imageSz = sz;

    if( m_imageSz == 0 )
    {
        m_imageSzAuto = true;
    }
    else
    {
        m_imageSzAuto = false;
    }

    m_opdMaskMade = false;
    m_preAllocated = false;
}

template <typename realT, typename detectorT>
bool pyramidSensor<realT, detectorT>::imageSzAuto()
{
    return m_imageSzAuto;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::imageSzAuto( const bool &ia )
{
    m_imageSzAuto = ia;

    m_opdMaskMade = false;

    m_preAllocated = false;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::makeOpdMask()
{
    complexFieldT opdMaskQ;

    if(m_wfPS == 0)
    {
        mxThrowException(mx::err::invalidconfig, "pyramidSensor::makeOpdMask()", "wavefront platescale (m_wfPS) is 0. Must set pupilSz and D first.");
    }

    if(!std::isfinite(m_wfPS) || !std::isnormal(m_wfPS))
    {
        mxThrowException(mx::err::invalidconfig, "pyramidSensor::makeOpdMask()", "wavefront platescale (m_wfPS) is infinite. Must set pupilSz and D first.");
    }

    std::cerr << m_wfPS << " " << m_D << "\n";

    if(m_D == 0)
    {
        mxThrowException(mx::err::invalidconfig, "pyramidSensor::makeOpdMask()", "pupil diameter is 0. Must set D > 0 first.");
    }

    //Setup the Fraunhoffer Propagator
    m_frProp.setWavefrontSizePixels(m_wfSz);

    m_opdMask.resize(m_wfSz, m_wfSz);
    opdMaskQ.resize(m_wfSz, m_wfSz);

    mx::improc::eigenImage<realT> mask;
    
    mask.resize(m_opdMask.rows(), m_opdMask.cols());
    realT dang = mx::math::two_pi<realT>()/m_nSides;

    realT minx = 0;
    realT maxx = 0;
    realT miny = 0;
    realT maxy = 0;

    realT pupilRad = m_pupilSep*m_pupilSz /(2*sin(dang/2.0));

    for(int n = 0; n < m_nSides; ++n)
    {
        realT ang = m_angleOffset*math::degreesT<realT>::radians + 0.5*dang + n * dang;

        realT dx = pupilRad * cos(ang);

        if(dx < minx) minx = dx;
        if(dx > maxx) maxx = dx;

        realT dy = pupilRad * sin(ang);

        if(dy < miny) miny = dy;
        if(dy > maxy) maxy = dy;

        opdMaskQ.set(std::complex<realT>(0, 1));
        wfp::tiltWavefront(opdMaskQ, dx, dy);
        mask.setZero();
        improc::maskWedge(mask, 0.5*(mask.rows()-1), 0.5*(mask.cols()-1), math::rtod(ang), 0.5*math::rtod(dang), 1);
        wfp::extractMaskedPixels(m_opdMask, opdMaskQ, mask);
    }

    int xsz = 2*std::max( {fabs(maxx), fabs(minx)} ) + 2*std::max({(pupilRad/2),((realT)m_pupilSz/2)});
    int ysz = 2*std::max( {fabs(maxy), fabs(miny)} ) + 2*std::max({(pupilRad/2), ((realT)m_pupilSz/2)});

    if(m_imageSzAuto)
    {
        m_imageSz = std::max(xsz, ysz);
    }


    if(m_imageSz > m_wfSz)
    {
        std::string msg = "image size (m_imageSz = " + std::to_string(m_imageSz) + ") ";
        msg += "> wavefront size (m_wfSz = " + std::to_string(m_wfSz) + "). ";
        msg += "Decrease number of sides (m_nSides = " + std::to_string(m_nSides) + ") or increase wavefront size. ";
        mxThrowException(mx::err::invalidconfig, "pyramidSensor::makeOpdMask", msg);
    }

    m_wfsImage.image.resize(m_imageSz, m_imageSz);

    if(m_detRows == 0 || m_detCols == 0)
    {
        mxThrowException(mx::err::invalidconfig, "pyramidSensor::makeOpdMask", "must set detector size");
    }

    if(m_detRows > m_imageSz || m_detCols > m_imageSz)
    {
        mxThrowException(mx::err::invalidconfig, "pyramidSensor::makeOpdMask", "detector is larger than image size");
    }

    m_opdMaskMade = true;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::makeTilts()
{
    constexpr realT pi = math::pi<realT>();

    if( m_modSteps == 0 )
    {
        mxThrowException( mx::err::invalidconfig,
                          "pyramidSensor::makeTilts()",
                          "number of modulation steps (m_modSteps) has not been set." );
    }

    if(m_wfPS == 0)
    {
        mxThrowException( mx::err::invalidconfig,
                          "pyramidSensor::makeTilts()",
                          "wavefront platescale (m_wfPS) is 0. Must set pupilSz and D first." );
    }

    if( !std::isfinite( m_wfPS ) )
    {
        mxThrowException( mx::err::invalidconfig,
                          "pyramidSensor::makeTilts()",
                          "wavefront platescale (m_wfPS) is infinite. Must set pupilSz and D first." );
    }

    if(m_D == 0)
    {
        mxThrowException(mx::err::invalidconfig, "pyramidSensor::makeTilts()", "pupil diameter is 0. Must set D > 0 first.");
    }


    if(m_D == 0)
    {
        mxThrowException(mx::err::invalidconfig, "pyramidSensor::makeTilts()", "pupil diameter is 0. Must set D > 0 first.");
    }


    realT dang = 2 * pi / ( m_modSteps );
    realT dx, dy;

    m_tilts.resize( m_modSteps );

    std::cout << "WF Size: " << m_wfSz << "\n";
    std::cout << "WF PS:   " << m_wfPS << "\n";
    std::cout << "Lambda:  " << m_lambda << "\n";
    std::cout << "Pyr. PS: " << wfp::fftPlateScale<realT>( m_wfSz, m_wfPS, m_lambda ) * 206265. << " (mas/pix)\n";
    std::cout << "Mod. steps: " << m_modSteps << "\n";
    std::cout << "Mod rad: " << m_modRadius * ( m_lambda / m_D ) / wfp::fftPlateScale<realT>( m_wfSz, m_wfPS, m_lambda )
              << " (pixels)\n";

    for( int i = 0; i < m_modSteps; ++i )
    {
        dx = m_modRadius * ( m_lambda / m_D ) / wfp::fftPlateScale<realT>( m_wfSz, m_wfPS, m_lambda ) *
             cos( 0.0 * dang + dang * i );
        dy = m_modRadius * ( m_lambda / m_D ) / wfp::fftPlateScale<realT>( m_wfSz, m_wfPS, m_lambda ) *
             sin( 0.0 * dang + dang * i );

        m_tilts[i].resize( m_wfSz, m_wfSz );
        m_tilts[i].set( std::complex<realT>( 0, 1 ) );

        wfp::tiltWavefront( m_tilts[i], dx, dy );
    }

    m_tiltsMade = true;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::allocThreadMem()
{
    if( !m_opdMaskMade )
    {
        makeOpdMask(); // Needed for m_imageSz
    }

    m_pupilPlaneCF.resize( m_wfSz, m_wfSz );

    int maxTh = omp_get_max_threads();
    m_th_tiltedPlane.resize( maxTh );

    m_th_focalPlane.resize( maxTh );

    m_th_focalImage.resize( maxTh );

    m_th_sensorPlane.resize( maxTh );

    m_th_sensorImage.resize( maxTh );

    for( int nTh = 0; nTh < maxTh; ++nTh )
    {
        m_th_tiltedPlane[nTh].resize( m_wfSz, m_wfSz );

        m_th_focalPlane[nTh].resize( m_wfSz, m_wfSz );

        m_th_focalImage[nTh].resize( m_wfSz, m_wfSz );

        m_th_sensorPlane[nTh].resize( m_wfSz, m_wfSz );

        m_th_sensorImage[nTh].resize( m_imageSz, m_imageSz );
    }

    m_preAllocated = true;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::preAllocate()
{
    if( !m_tiltsMade )
    {
        makeTilts();
    }

    if( !m_opdMaskMade )
    {
        makeOpdMask();
    }

    if( !m_preAllocated )
    {
        allocThreadMem();
    }
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::doSenseWavefront()
{
    BREAD_CRUMB;

    wavefrontT pupilPlane;

    /* Here make average wavefront for now */
    int _firstWavefront = m_lastWavefront - m_iTime;
    if( _firstWavefront < 0 )
        _firstWavefront += _wavefronts.size();

    pupilPlane.amplitude = _wavefronts[_firstWavefront].amplitude;
    pupilPlane.phase = _wavefronts[_firstWavefront].phase;

    realT avgIt = _wavefronts[_firstWavefront].iterNo;

    BREAD_CRUMB;

    for( int i = 0; i < m_iTime; ++i )
    {
        ++_firstWavefront;
        if( (size_t)_firstWavefront >= _wavefronts.size() )
            _firstWavefront = 0;

        pupilPlane.amplitude += _wavefronts[_firstWavefront].amplitude;
        pupilPlane.phase += _wavefronts[_firstWavefront].phase;
        avgIt += _wavefronts[_firstWavefront].iterNo;
    }

    BREAD_CRUMB;

    pupilPlane.amplitude /= ( m_iTime + 1 );
    pupilPlane.phase /= ( m_iTime + 1 );

    avgIt /= ( m_iTime + 1.0 );

    /*=====================================*/
    doSenseWavefront( pupilPlane );

    m_wfsImage.iterNo = avgIt;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::doSenseWavefront( wavefrontT &pupilPlane )
{
    if( m_modRadius == 0 )
    {
        return doSenseWavefrontNoMod( pupilPlane );
    }

    if( !m_preAllocated )
    {
        preAllocate();
    }

    BREAD_CRUMB;

    m_wfsImage.image.resize( m_imageSz, m_imageSz );
    m_wfsImage.image.setZero();

    wfsTipImage.image.resize( m_wfSz, m_wfSz );
    wfsTipImage.image.setZero();

    for( size_t l = 0; l < m_wavelengths.size(); ++l )
    {
        pupilPlane.lambda = m_lambda;
        pupilPlane.getWavefront( m_pupilPlaneCF, m_wavelengths[l], m_wfSz );

#pragma omp parallel
        {
            int nTh = omp_get_thread_num();

            m_th_sensorImage[nTh].setZero();
            m_th_focalImage[nTh].setZero();

            complexT *tpm_Data;
            complexT *ppm_Data;
            complexT *tim_Data;
            complexT *fpm_Data;
            complexT *opdm_Data;

            ppm_Data = m_pupilPlaneCF.data();
            tpm_Data = m_th_tiltedPlane[nTh].data();
            fpm_Data = m_th_focalPlane[nTh].data();
            opdm_Data = m_opdMask.data();

            int nelem = m_wfSz * m_wfSz;

#pragma omp for
            for( int i = 0; i < m_modSteps; ++i )
            {

                tim_Data = m_tilts[i].data();

                //---------------------------------------------
                // Apply the modulating tip
                //---------------------------------------------
                for( int ii = 0; ii < nelem; ++ii )
                {
                    tpm_Data[ii] = ppm_Data[ii] * tim_Data[ii];
                }

                //---------------------------------------------
                // Propagate to Pyramid tip
                //---------------------------------------------
                m_frProp.propagatePupilToFocal( m_th_focalPlane[nTh], m_th_tiltedPlane[nTh], true );

                //---------------------------------------------
                // Extract the tip image.
                //---------------------------------------------
                wfp::extractIntensityImageAccum(
                    m_th_focalImage[nTh], 0, m_wfSz, 0, m_wfSz, m_th_focalPlane[nTh], 0, 0 );

                //---------------------------------------------
                // Now apply the pyramid OPD
                //---------------------------------------------
                for( int ii = 0; ii < nelem; ++ii )
                {
                    fpm_Data[ii] = fpm_Data[ii] * opdm_Data[ii];
                }

                //---------------------------------------------
                // Propagate to sensor plane
                //---------------------------------------------
                m_frProp.propagateFocalToPupil( m_th_sensorPlane[nTh], m_th_focalPlane[nTh], true );

                //---------------------------------------------
                // Extract the image.
                //---------------------------------------------
                wfp::extractIntensityImageAccum( m_th_sensorImage[nTh],
                                                 0,
                                                 m_imageSz,
                                                 0,
                                                 m_imageSz,
                                                 m_th_sensorPlane[nTh],
                                                 0.5 * m_wfSz - m_imageSz / 2,
                                                 0.5 * m_wfSz - m_imageSz / 2 );

            } // for

            BREAD_CRUMB;

#pragma omp critical
            {
                m_wfsImage.image += m_th_sensorImage[nTh] * _wavelengthWeights[l];
                wfsTipImage.image += m_th_focalImage[nTh] * _wavelengthWeights[l];
            }
        } // #pragma omp parallel

    } // l for wavelength

    BREAD_CRUMB;

    m_wfsImage.image /= m_modSteps;
    wfsTipImage.image /= m_modSteps;
}

template <typename realT, typename detectorT>
void pyramidSensor<realT, detectorT>::doSenseWavefrontNoMod( wavefrontT &pupilPlane )
{
    BREAD_CRUMB;

    if( !m_opdMaskMade )
    {
        makeOpdMask();
    }

    m_wfsImage.image.resize( m_imageSz, m_imageSz );
    m_wfsImage.image.setZero();

    wfsTipImage.image.resize( m_wfSz, m_wfSz );
    wfsTipImage.image.setZero();

    complexFieldT m_pupilPlaneCF;

    pupilPlane.getWavefront( m_pupilPlaneCF, m_wfSz );

    complexFieldT tiltedPlane;
    complexFieldT focalPlane;
    complexFieldT sensorPlane;

    tiltedPlane.resize( m_wfSz, m_wfSz );
    focalPlane.resize( m_wfSz, m_wfSz );
    sensorPlane.resize( m_wfSz, m_wfSz );

    int nelem = m_wfSz * m_wfSz;

    complexT *tpm_Data = tiltedPlane.data();
    complexT *ppm_Data = m_pupilPlaneCF.data();
    complexT *opdm_Data = m_opdMask.data();
    complexT *fpm_Data = focalPlane.data();

    BREAD_CRUMB;

    //---------------------------------------------
    // Apply NO modulator tilt
    //---------------------------------------------
    for( int ii = 0; ii < nelem; ++ii )
    {
        tpm_Data[ii] = ppm_Data[ii];
    }

    BREAD_CRUMB;

    //---------------------------------------------
    // Propagate to Pyramid tip
    //---------------------------------------------
    m_frProp.propagatePupilToFocal( focalPlane, tiltedPlane, true );

    BREAD_CRUMB;

    //---------------------------------------------
    // Extract the tip image.
    //---------------------------------------------
    wfp::extractIntensityImageAccum( wfsTipImage.image, 0, m_wfSz, 0, m_wfSz, focalPlane, 0, 0 );

    BREAD_CRUMB;

    //---------------------------------------------
    // Now apply the pyramid OPD
    //---------------------------------------------
    for( int ii = 0; ii < nelem; ++ii )
    {
        fpm_Data[ii] = fpm_Data[ii] * opdm_Data[ii];
    }

    BREAD_CRUMB;

    //---------------------------------------------
    // Propagate to sensor plane
    //---------------------------------------------
    m_frProp.propagateFocalToPupil( sensorPlane, focalPlane, true );

    BREAD_CRUMB;

    //---------------------------------------------
    // Extract the image.
    //---------------------------------------------
    wfp::extractIntensityImageAccum( m_wfsImage.image,
                                     0,
                                     m_imageSz,
                                     0,
                                     m_imageSz,
                                     sensorPlane,
                                     0.5 * m_wfSz - m_imageSz / 2,
                                     0.5 * m_wfSz - m_imageSz / 2 );

    BREAD_CRUMB;
}

} // namespace sim
} // namespace AO
} // namespace mx

#endif // mx_AO_sim_pyramidSensor_hpp
